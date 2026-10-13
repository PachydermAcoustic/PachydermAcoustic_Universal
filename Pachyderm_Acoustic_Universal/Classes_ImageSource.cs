//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2025, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
//'Pachyderm-Acoustic is free software; you can redistribute it and/or modify 
//'it under the terms of the GNU General Public License as published 
//'by the Free Software Foundation; either version 3 of the License, or 
//'(at your option) any later version. 
//'Pachyderm-Acoustic is distributed in the hope that it will be useful, 
//'but WITHOUT ANY WARRANTY; without even the implied warranty of 
//'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
//'GNU General Public License for more details. 
//' 
//'You should have received a copy of the GNU General Public 
//'License along with Pachyderm-Acoustic; if not, write to the Free Software 
//'Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 

using Hare.Geometry;
using MathNet.Numerics.LinearAlgebra.Double;
using Pachyderm_Acoustic.Environment;
using Pachyderm_Acoustic.Pach_Graphics;
using Pachyderm_Acoustic.Utilities;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Pachyderm_Acoustic
{
    public class ImageSourceData : Simulation_Type
    {
        private List<Deterministic_Reflection>[] ValidPaths;
        private List<Deterministic_Reflection>[,] ThreadPaths;
        private double Speed_of_Sound;
        private int MaxOrder;
        private int BandCount = 8;
        public bool ThirdOctave { get { return BandCount == 24; } }
        private Source Src;
        private Source[] Sources;
        private Source[] ThreadSources;
        private Hare.Geometry.Point[] Rec;
        private Polygon_Scene Room;
        private int[] CurrentSrf;
        private int[] CurrentEdge;
        private DateTime ST;
        private Random[] Rnd;
        private double[] Direct_Time;
        System.Threading.Thread[] T_List;
        int processorCT;
        int elementCt;
        private TimeSpan TS;
        private int SampleCT;
        private int SampleRate;
        public int SrcNo;
        public int[] Oct_choice;
        public bool Diffraction = false;
        public bool Screencalc = false;
        public bool IncludeEdges = false;

        private ImageSourceData()
        { }

        public ImageSourceData(Source Source, Receiver_Bank Receiver, Direct_Sound Direct, Polygon_Scene Rm, int MaxOrder_in, int SourceID_in)
        : this(Source, Receiver, Direct, Rm, new int[2] { 0, Receiver.Third_Oct ? 23 : 7 }, MaxOrder_in, false, SourceID_in)
        { }

        public ImageSourceData(Source Source, Receiver_Bank Receiver, Direct_Sound Direct, Polygon_Scene Rm, int MaxOrder_in, bool ED, int SourceID_in)
            : this(Source, Receiver, Direct, Rm, new int[2] { 0, Receiver.Third_Oct ? 23 : 7 }, MaxOrder_in, ED, SourceID_in)
        { }

        /// <summary>
        /// Constructor prepares image source calculation to run.
        /// </summary>
        /// <param name="Source"></param>
        /// <param name="Receiver"></param>
        /// <param name="Direct"></param>
        /// <param name="Rm"></param>
        /// <param name="Octaves">Inclusive first/last spectral band indices in the receiver-bank resolution.</param>
        /// <param name="MaxOrder_in">The maximum order to be calculated for.</param>
        public ImageSourceData(Source Source, Receiver_Bank Receiver, Direct_Sound Direct, Polygon_Scene Rm, int[] Octaves, int MaxOrder_in, bool Edge_Diffraction, int SourceID_in)
        {
            IncludeEdges = Edge_Diffraction;
            Diffraction = Edge_Diffraction;
            int bandCount = Receiver.Third_Oct ? 24 : 8;
            if (Octaves == null || Octaves.Length != 2 || Octaves[0] < 0 || Octaves[1] >= bandCount || Octaves[0] > Octaves[1])
                throw new ArgumentOutOfRangeException(nameof(Octaves), "Provide an inclusive spectral band range from 0 through " + (bandCount - 1) + ".");
            Oct_choice = new int[Octaves[1] - Octaves[0] + 1];
            for (int i = 0; i < Oct_choice.Length; i++) Oct_choice[i] = i + Octaves[0];
            SrcNo = SourceID_in;
            ValidPaths = new List<Deterministic_Reflection>[Receiver.Count];
            Speed_of_Sound = Rm.Sound_speed(0);
            MaxOrder = MaxOrder_in;
            Src = Source.ForBandResolution(Receiver.Third_Oct);
            Sources = Src is SourceCluster cluster ? cluster.Sources.ToArray() : new Source[] { Src };
            BandCount = Src.BandCount;
            Rec = new Hare.Geometry.Point[Receiver.Count];
            SampleCT = Receiver.SampleCT;
            SampleRate = Receiver.SampleRate;
            for (int i = 0; i < Receiver.Count; i++)
            {
                Rec[i] = Receiver.Origin(i);
            }
            Room = Rm;
            Direct_Time = new double[Receiver.Count];
            for (int q = 0; q < Receiver.Count; q++)
            {
                Direct_Time[q] = Direct.Min_Time(q);
            }
        }

        /// <summary>
        /// Called by Pach_RunSim_Command. Indicates whether or not the simulation has completed.
        /// </summary>
        /// <returns>Returns running if any threads in this simulation are still running. Returns stopped if all have stopped.</returns>
        public override System.Threading.ThreadState ThreadState()
        {
            foreach (System.Threading.Thread T in T_List)
            {
                if (T != null && T.IsAlive) return System.Threading.ThreadState.Running;
            }
            return System.Threading.ThreadState.Stopped;
        }

        /// <summary>
        /// Inherited member. Divides the simulation into threads, and begins.
        /// </summary>
        public override void Begin()
        {
            processorCT = System.Environment.ProcessorCount;
            ThreadPaths = new List<Deterministic_Reflection>[Rec.Length, processorCT];
            ThreadSources = new Source[processorCT];
            int SrfCT = Room.ObjectCount;
            //
            if (Diffraction) SrfCT += Room.EdgeCount;
            //
            CurrentSrf = new int[processorCT];
            CurrentEdge = new int[processorCT];
            Rnd = new Random[processorCT];
            T_List = new System.Threading.Thread[processorCT];
            elementCt = Room.ObjectCount + ((this.Diffraction) ? Room.EdgeCount : 0);
            SrfCT *= Sources.Length;

            for (int P_I = 0; P_I < processorCT; P_I++)
            {
                for (int i = 0; i < Rec.Length; i++)
                {
                    ThreadPaths[i, P_I] = new List<Deterministic_Reflection>();
                    ValidPaths[i] = new List<Deterministic_Reflection>();
                }
                int start = (int)Math.Floor((double)P_I * SrfCT / processorCT);
                int end;
                if (P_I == processorCT - 1) end = SrfCT;
                else end = (P_I + 1) * SrfCT / processorCT;

                Calc_Params T = new Calc_Params(start, end, P_I, Room.R_Seed.Next());
                System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { Calculate(T); });
                T_List[P_I] = new System.Threading.Thread(TS);
                T_List[P_I].Start();
            }
        }

        /// <summary>
        /// Called by each thread from the begin method.
        /// </summary>
        /// <param name="i">the object is type "Calc_Params" which holds all the necessary information to run a portion of the simulation.</param>
        public void Calculate(object i)
        {
            Calc_Params Params = (Calc_Params)i;
            Rnd[Params.ThreadID] = new Random(Params.RandomSeed);
            ST = DateTime.Now;
            int[] Sequence = new int[1];
            for (CurrentSrf[Params.ThreadID] = 0; CurrentSrf[Params.ThreadID] < Params.EndIndex - Params.StartIndex; CurrentSrf[Params.ThreadID]++)
            {
                List<Hare.Geometry.Point[]> Images = new List<Hare.Geometry.Point[]>();
                List<int[]> sub = new List<int[]>();
                int element = CurrentSrf[Params.ThreadID] + Params.StartIndex;
                ThreadSources[Params.ThreadID] = Sources[element / elementCt];
                Sequence[0] = element % elementCt;
                TraverseOrders(0, Sequence, Params.ThreadID, Images, sub);
            }
        }

        /// <summary>
        /// Used by specular raytracer.
        /// </summary>
        /// <param name="Sequences">List of surface index sequences to try.</param>
        public void Lookup_Sequences(List<int[]>[] Sequences)
        {
            for (int rec = 0; rec < Rec.Length; rec++) if (ValidPaths[rec] == null) ValidPaths[rec] = new List<Deterministic_Reflection>();
            processorCT = System.Environment.ProcessorCount;
            ThreadPaths = new List<Deterministic_Reflection>[Rec.Length, processorCT];
            ThreadSources = new Source[processorCT];
            int SrfCT = Room.ObjectCount;
            CurrentSrf = new int[processorCT];
            Rnd = new Random[processorCT];
            Random R = new Random();
            for (int thread = 0; thread < processorCT; thread++) Rnd[thread] = new Random(R.Next());
            T_List = new System.Threading.Thread[processorCT];
            for (int p = 0; p < Rec.Length; p++)
            {
                for (int P_I = 0; P_I < processorCT; P_I++)
                {
                    ThreadPaths[p, P_I] = new List<Deterministic_Reflection>();
                    Calc_Params T = new Calc_Params(P_I * Sequences[p].Count / processorCT, (P_I + 1) * Sequences[p].Count / processorCT, P_I, R.Next());
                    System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate
                    {
                        foreach (Source source in Sources)
                        {
                            ThreadSources[T.ThreadID] = source;
                            for (int i = T.StartIndex; i < T.EndIndex; i++)
                                ProcessPath(Sequences[p][i], T.ThreadID, p);
                        }
                    });
                    T_List[P_I] = new System.Threading.Thread(TS);
                    T_List[P_I].Start();
                }
                do
                {
                    System.Threading.Thread.Sleep(1000);
                    if (ThreadState() != System.Threading.ThreadState.Running) break;
                } while (true);

                for (int t = 0; t < processorCT; t++)
                {
                    ValidPaths[p].AddRange(ThreadPaths[p, t]);
                }
            }
        }

        /// <summary>
        /// A string to identify the type of simulation being run.
        /// </summary>
        /// <returns></returns>
        public override string Sim_Type()
        {
            return "Image Source";
        }

        /// <summary>
        /// Called by Pach_RunSim_Command. Get a string describing the status of the simulation for display.
        /// </summary>
        /// <returns></returns>
        public override string ProgressMsg()
        {
            TS = DateTime.Now - ST;
            int Srf_CT = 0;
            int surfaceCount = Room.ObjectCount * Sources.Length;
            int edgeCount = Room.EdgeCount * Sources.Length;

            for (int i = 0; i < processorCT; i++) Srf_CT += CurrentSrf[i];

            if (Diffraction)
            {
                int Edge_CT = 0;
                for (int i = 0; i < processorCT; i++) Edge_CT += CurrentEdge[i];

                TS = new TimeSpan((long)(TS.Ticks * (((double)(surfaceCount - Srf_CT) + (edgeCount - Edge_CT)) / (Srf_CT + Edge_CT + 2))));
                return string.Format("Calculating Image/Edge Set {0} of {1}. ({2} hours,{3} minutes,{4} seconds Remaining.) Press 'Esc' to Cancel...", Srf_CT + Edge_CT, surfaceCount + edgeCount, TS.Hours, TS.Minutes, TS.Seconds);
            }
            else
            {
                TS = new TimeSpan((long)(TS.Ticks * (((double)(surfaceCount - Srf_CT)) / (Srf_CT + 1))));
                return string.Format("Calculating Image Set {0} of {1}. ({2} hours,{3} minutes,{4} seconds Remaining.) Press 'Esc' to Cancel...", Srf_CT, surfaceCount, TS.Hours, TS.Minutes, TS.Seconds);
            }
        }

        /// <summary>
        /// Consolidates output from all threads into a single set of output.
        /// </summary>
        public override void Combine_ThreadLocal_Results()
        {
            for (int i = 0; i < Rec.Length; i++)
            {
                for (int p = 0; p < processorCT; p++)
                {
                    ValidPaths[i].AddRange(ThreadPaths[i, p]);
                }
            }
        }

        /// <summary>
        /// Recursive function which carries out image source/edge source method for orders greater than 1.
        /// </summary>
        /// <param name="LastOrder">The order of the calling iteration.</param>
        /// <param name="Sequence">The input sequence, to which the next surface will be appended.</param>
        /// <param name="ThreadId">The id of the thread using the function.</param>
        /// <param name="Images">the list of source images, to which the next image will be appended.</param>
        private void TraverseOrders(int Order, int[] Sequence, int ThreadId, List<Hare.Geometry.Point[]> Images0, List<int[]> sub0)
        {
            Array.Resize(ref Sequence, Order + 1);
            int ct = Order == 0 ? 1 : elementCt;
            for (int q = 0; q < ct; q++)
            {
                List<Hare.Geometry.Point[]> Images = new List<Hare.Geometry.Point[]>();
                List<int[]> sub = new List<int[]>();
                if (Order < MaxOrder)
                {
                    if (Order > 0) Sequence[Order] = q;

                    if (Images0.Count == 0)
                    {
                        Images0.Add(new Point[0]);
                        sub0.Add(new int[0]);
                    }

                    for (int r = 0; r < Images0.Count; r++)
                    {
                        Hare.Geometry.Point[] imnext = new Point[Images0[r].Length + 1];
                        int[] subnext = new int[sub0[r].Length + 1];
                        for (int i = 0; i < Images0[r].Length; i++)
                        { 
                            imnext[i] = Images0[r][i];
                            if (i == 0 || i == sub0[r].Length) subnext[i] = sub0[r][i];
                        }
                        if (Sequence[Order] > Room.ObjectCount - 1)
                        {
                            //Its an edge...
                            if (Order > 0 && Sequence[Order] != Sequence[Sequence.Length - 1]) continue;
                            int edge_id = Sequence[Order] - Room.ObjectCount;
                            if (Room.Edge_Nodes[edge_id].EdgeSources.Count < 5) continue;
                            for (int i = 0; i < Room.Edge_Nodes[edge_id].EdgeSources.Count; i++)
                            {
                                Hare.Geometry.Point[] im = imnext.Clone() as Hare.Geometry.Point[];
                                int[] s = subnext.Clone() as int[];
                                im[Order] = Room.Edge_Nodes[edge_id].EdgeSources[i].Z_mid;
                                s[Order] = i;
                                Images.Add(im.Clone() as Hare.Geometry.Point[]);
                                sub.Add(s); //sub is the edge source id.
                            }
                        }
                        else if (Room.IsPlanar(Sequence[Order]))
                        {
                            //It's Planar...
                            if (Order > 0 && Sequence[Order] != Sequence[Sequence.Length - 1]) continue;
                            Hare.Geometry.Point[] im = imnext.Clone() as Hare.Geometry.Point[];
                            int[] s = subnext.Clone() as int[];
                            im[Order] = Room.Image(Order > 0 ? Images0[r][Images0[r].Length - 1] : ThreadSources[ThreadId].Origin, 0, Room.ObjectMembers[Sequence[Order]][0]);
                            s[Order] = -1;//no sub needed...
                            Images.Add(im);
                            sub.Add(s);
                        }
                        else
                        {
                            //It's curved
                            //Repeats allowed for this case...
                            for (int i = 0; i < Room.ObjectMeshEdges[Sequence[Order]].Length; i++)
                            {
                                Hare.Geometry.Point[] im = imnext.Clone() as Hare.Geometry.Point[];
                                int[] s = subnext.Clone() as int[];
                                im[Order] = Room.ObjectMeshEdges[Sequence[Order]][i].mid;
                                s[Order] = i;//sub is the mesh edge id
                                Images.Add(im);
                                sub.Add(s); 
                            }
                        }
                        ProcessImages(Images.ToArray(), Sequence, ThreadId, sub);
                        if (Order + 1 < MaxOrder)
                        {
                            TraverseOrders(Order + 1, Sequence, ThreadId, Images, sub);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// This function calculates the actual path of the specular/diffracted/focussed reflection.
        /// </summary>
        /// <param name="Images">The list of images.</param>
        /// <param name="Sequence">The list of surface indices for reflection.</param>
        /// <param name="Threadid">The id of the calling thread.</param>
        private void ProcessImages(Hare.Geometry.Point[][] Images, int[] Sequence, int Threadid, List<int[]> sub, int one_rec_id = -1)
        {
            Source source = ThreadSources[Threadid];
            int rec_end = one_rec_id < 0 ? rec_end = Rec.Length : rec_end = one_rec_id + 1;

            for (int rec_id = one_rec_id < 0 ? 0 : one_rec_id; rec_id < rec_end; rec_id++)
            {
                double c_sound = Room.Sound_speed(Rec[rec_id]);

                double[][] Trans_Mod = new double[Images.Length][];
                int[] Seq_Polys = new int[Sequence.Length];
                List<Hare.Geometry.Point[]> PathVertices = new List<Hare.Geometry.Point[]>();
                Hare.Geometry.Point S = source.Origin;
                Hare.Geometry.Point E = Rec[rec_id];
                double df = SampleRate * .5 / 16384;

                //Find all Path Legs from Receiver to Source
                for (int r = 0; r < Images.Length; r++)
                {
                    Trans_Mod[r] = new double[BandCount];
                    foreach (int band in Oct_choice) Trans_Mod[r][band] = 1;

                    Hare.Geometry.Point[] path = new Hare.Geometry.Point[Sequence.Length + 2];
                    path[0] = S;
                    path[path.Length - 1] = E;

                    for (int q = Sequence.Length - 1; q >= 0; q--)
                    {
                        if (Sequence[q] > Room.ObjectCount - 1)
                        {
                            //It's an edge
                            if (!OcclusionIntersectED(path[q + 2], Images[r][q], Sequence[q], sub[r][q], ref Trans_Mod[r], Threadid))
                            {
                                path = null; break; 
                            }
                            else path[q + 1] = Images[r][q];
                        }
                        else
                        {
                            //It's a plane or curve!
                            //Check sequence to determine if it is curved. These should be bundled so that it is clear which (for example) edge node is associated with which (for example) other edge or curved surface reflection.
                            if (Room.IsPlanar(Sequence[q]))
                            {
                                if (!OcclusionIntersect(path[q + 2], Images[r][q], Sequence[q], ref Trans_Mod[r], ref path[q + 1], ref Seq_Polys[q], Threadid))
                                { path = null; break; }
                            }
                            else
                            {
                                /// it's curved. Check that these are meeting the ballpark criteria for a curved reflection. Then it's business as usual.
                                if (Room.ObjectMeshEdges[Sequence[q]][r].Polys.Count != 2) { path = null; break; }
                                if (!this.Check_Edge_For_Specular(path[q + 2], path[q] == null ? Images[r][q] : path[q], Room.ObjectMeshEdges[Sequence[q]][r].Polys[0].Normal * -1, Room.ObjectMeshEdges[Sequence[q]][r].Polys[1].Normal * -1, Room.ObjectMeshEdges[Sequence[q]][r].a, Room.ObjectMeshEdges[Sequence[q]][r].b, ref path[q + 1])) { path = null; break; }
                                //if (!OcclusionIntersect(path[q + 2], path[q + 1], Sequence[q], ref Trans_Mod[r], ref path[q + 1], ref Seq_Polys[q], Threadid)) { path = null; break; }
                                if (!OcclusionIntersectED(path[q + 2], Images[r][q], Sequence[q], sub[r][q], ref Trans_Mod[r], Threadid)) { path = null; break; }
                            }
                        }
                    }
                    PathVertices.Add(path);
                }

                //Check that any path was unoccluded... if so, then record this entry. If not, move on...
                if (PathVertices.Count(item => item != null) == 0) continue;
                
                //Final Occlusion Check:
                for (int r = 0; r < PathVertices.Count; r++)
                {
                    if (PathVertices[r] == null) continue;
                    if (Sequence[0] < Room.ObjectCount)
                    {
                        for (int i = 0; i < PathVertices[r].Length - 2; i++)
                        {
                            Hare.Geometry.Vector d1 = PathVertices[r][i] - PathVertices[r][i + 1], d2 = PathVertices[r][i + 2] - PathVertices[r][i + 1], n = Room.Normal(Room.ObjectMembers[Sequence[i]][0]);
                            d1.Normalize(); d2.Normalize();
                            if (Math.Sign(Hare_math.Dot(n, d1)) != Math.Sign(Hare_math.Dot(n, d2)))
                            {
                                PathVertices[r] = null;break;
                            }
                            else if ((i == PathVertices[r].Length - 3) && FinalOcclusion(PathVertices[r][0], PathVertices[r][1], Sequence[0], ref Trans_Mod[r], Threadid))
                            { PathVertices[r] = null; break; }
                        }

                    }
                    else
                    {
                        if (!OcclusionIntersectED(PathVertices[r][0], PathVertices[r][1], Sequence[0], sub[r][0], ref Trans_Mod[r], Threadid))
                            PathVertices[r] = null;
                    }
                }

                //Check again for null(occluded) paths...
                if (PathVertices.Count(item => item != null) == 0) continue; //goto Next;

                ///Process all paths for pulse entry...
                if (PathVertices.Count == 0) continue;

                //////////////////////////////
                //Process Compound Path before storing it.
                double[][] times;
                double mintime;
                int[] elementIds;
                //For 
                double[][] B_Functions_Octave = B_Function_NoInterp(PathVertices, Sequence, Trans_Mod, c_sound, Threadid, out times, out elementIds);
                double[][][] H = H_Function(B_Functions_Octave, PathVertices, Sequence, elementIds, times, out mintime, Threadid);

                ///Enter the reflection
                PathVertices.RemoveAll(item => item == null);
                double diff_length = (PathVertices[0][1] - PathVertices[PathVertices.Count - 1][1]).Length();
                if (Sequence.Max() > Room.ObjectCount && ( diff_length < 0.025 || (diff_length < 0.05 && PathVertices.Count > 5))) continue;

                ///Process all paths for pulse entry...
                if (PathVertices.Count > 1)
                {
                    double[][] Absorption = new double[Sequence.Length][];
                    for (int q = 0; q < Sequence.Length; q++)
                        Absorption[q] = Sequence[q] >= Room.ObjectCount ? new double[BandCount] : Room.AbsorptionValue[Room.ObjectMembers[Sequence[q]][0]].Coefficient_A_Broad(ThirdOctave);
                    ThreadPaths[rec_id, Threadid].Add(new Compound_Path(PathVertices.ToArray(), source.Source_ID(), Sequence, H, mintime, Room.Rho_C(0), ref Direct_Time[rec_id], Threadid, Absorption)
                    { CoherentSpeaker = Src is SourceCluster cluster && cluster.CoherentSpeakers && source.IsSpeaker, SteeringDelaySeconds = source.SpeakerDelaySeconds });
                }
                else
                {
                    //Assumed either a convex, or planar reflection
                    double[] H_in = new double[BandCount]; for (int oct = 0; oct < BandCount; oct++) H_in[oct] = H[0][oct][0];
                    ThreadPaths[rec_id, Threadid].Add(new Specular_Path(PathVertices[0], Sequence, Seq_Polys, Room, source, c_sound, H_in, Trans_Mod[0], ref Direct_Time[rec_id], Threadid, Rnd[Threadid].Next(), Oct_choice)
                    { CoherentSpeaker = Src is SourceCluster cluster && cluster.CoherentSpeakers && source.IsSpeaker, SteeringDelaySeconds = source.SpeakerDelaySeconds });
                }
            }
        }

        public bool Check_Edge_For_Specular(Hare.Geometry.Point S, Hare.Geometry.Point R, Hare.Geometry.Vector T1_Norm, Hare.Geometry.Vector T2_Norm, Hare.Geometry.Point Z1, Hare.Geometry.Point Z2, ref Hare.Geometry.Point Spec_PT)
        {
            //if (Math.Abs(Hare.Geometry.Hare_math.Dot(T1_Norm, T2_Norm)) > 0.9999)
            //    return false;
            Hare.Geometry.Vector n = Z2 - Z1;
            double Z_range2 = n.Length();
            n /= (Z_range2);
            Z_range2 /= 2;
            return Check_Edge_For_Specular(S, R, Z1, Z2, T1_Norm, T2_Norm, (Z2 + Z1) / 2, Z_range2, ref Spec_PT);
        }

        public bool Check_Edge_For_Specular(Hare.Geometry.Point S, Hare.Geometry.Point R, Hare.Geometry.Point Z1, Hare.Geometry.Point Z2, Hare.Geometry.Vector T1_Norm, Hare.Geometry.Vector T2_Norm, Hare.Geometry.Point Z_mid, double Z_range2, ref Hare.Geometry.Point Spec_Pt)
        {
            //Determine closest point to a midpoint between reflection points before and after...
            Hare.Geometry.Vector ZDir = Z2 - Z1;
            double denom = Hare_math.Dot(ZDir, ZDir);
            Hare.Geometry.Vector SZ1 = S - Z1;
            Hare.Geometry.Vector RZ1 = R - Z1;
            double s_t = Hare_math.Dot(SZ1, ZDir);
            double r_t = Hare_math.Dot(RZ1, ZDir);
            double s_d = Math.Sqrt(Hare_math.Dot(SZ1, SZ1) - (s_t * s_t) / denom);
            double r_d = Math.Sqrt(Hare_math.Dot(RZ1, SZ1) - (r_t * r_t) / denom);
            double d = s_d + r_d;
            Hare.Geometry.Point Cpt = (S * r_d/d + R * s_d/d);

            Hare.Geometry.Vector tDir = Cpt - Z1;
            double t = Hare_math.Dot(tDir, ZDir);
            if (t+0.01 < 0)
                return false;
            else
            {
                if (t-0.01 > denom)
                    return false;
            }

            Hare.Geometry.Point Z_pt = Z1 + ZDir * t / denom;
            Hare.Geometry.Vector v1 = Z_pt - S, v2 = Z_pt - R;
            v1.Normalize(); v2.Normalize();
            Hare.Geometry.Vector TheoreticalNormal = v1 + v2;
            TheoreticalNormal -= Hare_math.Dot(TheoreticalNormal, ZDir) / Hare_math.Dot(ZDir, ZDir) * ZDir;
            TheoreticalNormal.Normalize();

            Hare.Geometry.Vector bsect = T1_Norm + T2_Norm;
            bsect.Normalize();

            double dot_Norm1 = Math.Abs(Hare_math.Dot(T1_Norm, bsect));
            double dot_Norm2 = Math.Abs(Hare_math.Dot(T2_Norm, bsect));
            double testdot = Math.Abs(Hare_math.Dot(bsect, TheoreticalNormal));
            //If it is between the two normals adjacent to these edges, then reflection is valid. Return true.
            if (testdot*1.000001 > dot_Norm1 && testdot*1.00000 > dot_Norm2)
            {
                Spec_Pt = Z_pt;
                return true;
            }
            return false;
        }

        public double Power_Recursion(int pathid, int order, List<Point[]> PathVertices, int[] Sequence, double c_sound, double[] t_limits, int Threadid, ref MathNet.Numerics.LinearAlgebra.Double.DenseMatrix W_Kurvature, ref Hare.Geometry.Vector[] W_Frame)
        {
            double h = order > 0 ? Power_Recursion(pathid, order - 1, PathVertices, Sequence, c_sound, t_limits, Threadid, ref W_Kurvature, ref W_Frame) : 1;

            bool hasedge = false;

            if (order == 0)
            {
                Hare.Geometry.Vector dir = PathVertices[pathid][order + 1] - PathVertices[pathid][order];
                W_Frame = new Hare.Geometry.Vector[2];
                W_Frame[0] = Hare.Geometry.Hare_math.Cross(dir, new Hare.Geometry.Vector(0, 0, 1));
                W_Frame[1] = Hare.Geometry.Hare_math.Cross(dir, W_Frame[0]);
                double length = dir.Length();
                W_Kurvature = new DenseMatrix(2);
                W_Kurvature[0, 0] = 1.0 / length;
                W_Kurvature[1, 1] = 1.0 / length;
            }

            if (Sequence[order] > Room.ObjectCount - 1)
            {
                //an edge
                hasedge = true;
                double m = 0, l = 0;
                double[] dm = new double[0], dl = new double[0];

                if (order == 0)
                {
                    Hare.Geometry.Vector d3 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1] - PathVertices[pathid][0]; // Direction Hare.Geometry.Vector of segment S1
                    Hare.Geometry.Vector d4 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0] - PathVertices[pathid][0]; // Direction Hare.Geometry.Vector of segment S2
                    double l3 = d3.Length();
                    double l4 = d4.Length();
                    t_limits[0] += l3 / c_sound;
                    t_limits[1] += l4 / c_sound;
                }

                h *= Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Flex_Solve(PathVertices[pathid][order], PathVertices[pathid][order + 2], ref m, ref l, ref dm, ref dl) * Math.Sqrt(W_Kurvature.Determinant());
                W_Kurvature = new MathNet.Numerics.LinearAlgebra.Double.DenseMatrix(2);
                W_Kurvature[0, 0] = 1.0 / l;
                W_Kurvature[1, 1] = 1.0 / l;
                t_limits[0] += l / c_sound; t_limits[1] += l / c_sound;

                if (order == Sequence.Length - 1)
                { 
                    Hare.Geometry.Vector d1 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0] - PathVertices[pathid][Sequence.Length + 1];//Room.ObjectCount].EdgeSources[pathid].Z_limits[0]; // Direction Hare.Geometry.Vector of segment S1
                    Hare.Geometry.Vector d2 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1] - PathVertices[pathid][Sequence.Length + 1];//Room.ObjectCount].EdgeSources[pathid].Z_limits[1]; // Direction Hare.Geometry.Vector of segment S2
                    Hare.Geometry.Vector r = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0] - Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1];
                    double a = Hare.Geometry.Hare_math.Dot(d1, d1); // Squared length of segment S1, always nonnegative
                    double e = Hare.Geometry.Hare_math.Dot(d2, d2); // Squared length of segment S2, always nonnegative
                    double c = Hare.Geometry.Hare_math.Dot(d1, r);
                    double f = Hare.Geometry.Hare_math.Dot(d2, r);

                    double b = Hare.Geometry.Hare_math.Dot(d1, d2);
                    double denom = a * e - b * b;
                    double s = denom == 0 ? 0 : Math.Min(Math.Max((b * f - c * e) / denom, 0), 1);

                    // Compute point on L2 closest to S1(s) using
                    double t = (b * s + f) / e;
                    // If t in [0,1] done.
                    if (t < 1 && t > 0)
                    {
                        //Choose the reverse...
                        Hare.Geometry.Vector d3 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1] - Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0]; // Direction Hare.Geometry.Vector of segment S1
                        Hare.Geometry.Vector d4 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0] - Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1]; // Direction Hare.Geometry.Vector of segment S2
                        t_limits[0] += d3.Length() / c_sound;
                        t_limits[1] += d4.Length() / c_sound;
                    }
                    else
                    {
                        //Use the lengths represented by this test...
                        t_limits[0] += Math.Sqrt(a) / c_sound;
                        t_limits[1] += Math.Sqrt(e) / c_sound;
                    }
                }
                else
                {
                    //Frame is arbitrary (for now)...
                    Hare.Geometry.Vector dir = PathVertices[pathid][order + 1] - PathVertices[pathid][order];
                    W_Frame = new Hare.Geometry.Vector[2];
                    W_Frame[0] = Hare.Geometry.Hare_math.Cross(dir, new Hare.Geometry.Vector(0, 0, 1));
                    W_Frame[1] = Hare.Geometry.Hare_math.Cross(dir, W_Frame[0]);
                }
            }
            else if (Room.IsPlanar(Sequence[order]))
            {
                hasedge = false;

                double l = 0;
                //a plane
                if (order == 0)
                {
                    l += (PathVertices[pathid][order] - PathVertices[pathid][order + 1]).Length();
                    t_limits[0] += l / c_sound;
                    t_limits[1] += l / c_sound;
                    l = 0;
                }

                Hare.Geometry.Vector dir1 = (PathVertices[pathid][order + 1] - PathVertices[pathid][order + 2]);
                l += dir1.Length();
                t_limits[0] += l / c_sound;
                t_limits[1] += l / c_sound;
                W_Kurvature = project_wavefront(l, W_Kurvature);

                Hare.Geometry.Vector dir2 = (PathVertices[pathid][order + 1] - PathVertices[pathid][order]);
                Hare.Geometry.Vector local_N = dir1 + dir2;
                local_N.Normalize();
                W_Frame[0] -= local_N * Hare_math.Dot(W_Frame[0], local_N) * 2;
                W_Frame[1] -= local_N * Hare_math.Dot(W_Frame[1], local_N) * 2;

                if (order < Sequence.Length - 1) return h;
            }
            else
            {
                hasedge = false;
                Hare.Geometry.Point a;
                Hare.Geometry.Point b;
                double K0 = 0 + Room.Edge_Kurvatures[Sequence[order]][pathid % Room.ObjectMeshEdges.Length][0], K1 = 0 + Room.Edge_Kurvatures[Sequence[order]][pathid % Room.ObjectMeshEdges.Length][1];
                Hare.Geometry.Vector Normal = new Hare.Geometry.Vector(Room.Edge_Normals[Sequence[order]][pathid].dx, Room.Edge_Normals[Sequence[order]][pathid].dy, Room.Edge_Normals[Sequence[order]][pathid].dz);
                Hare.Geometry.Vector dir0 = PathVertices[pathid][order + 1] - PathVertices[pathid][order];
                dir0.Normalize();
                double dot = Hare.Geometry.Hare_math.Dot(dir0, Normal);
                if (dot > 0) { K0 *= -1; K1 *= -1; Normal *= -1; }
                //K0 *= -1;
                //K1 *= -1;

                //a curve
                double[] dist = null;
                if (order == 0)
                {
                    h *= Math.Sqrt(W_Kurvature.Determinant());

                    if (K0 < 0 && K1 < 0)
                    {
                        dist = new double[2 + Room.ObjectMeshEdges[Sequence[order]][pathid].Polys.Count];
                        a = Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].a;
                        b = Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].b;
                        dist[0] += (PathVertices[pathid][0] - a).Length();
                        dist[1] += (PathVertices[pathid][0] - b).Length();
                        dist[2] += (PathVertices[pathid][0] - Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].Polys[0].Centroid).Length();
                        dist[3] += (PathVertices[pathid][0] - Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].Polys[1].Centroid).Length();
                    }
                    else
                    {
                        dist = new double[Room.ObjectMeshEdges[Sequence[order]][pathid].Polys.Count];
                        for (int i = 0; i < 2; i++)
                        {
                            int pid = pathid % Room.ObjectMeshEdges[Sequence[order]].Length;
                            Hare.Geometry.Edge e = Room.ObjectMeshEdges[Sequence[order]][pid];
                            Hare.Geometry.Vector p = Room.ObjectMeshEdges[Sequence[order]][pid].TributaryLength[i] * e.Tangents[i];
                            dist[i] += (PathVertices[pid][0] - (PathVertices[pid][1] + p)).Length();
                        }
                    }
                }

                Hare.Geometry.Vector dir1 = -1 * (PathVertices[pathid][order + 2] - PathVertices[pathid][order + 1]);
                double l = dir1.Length();
                dir1 /= l;
                a = Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].a;
                b = Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].b;

                //TODO: Adjust this and other statements to determine if the relationship to the incident sounds is indeed concave...
                if (K0 < 0 && K1 < 0)
                {
                    if (dist == null) dist = new double[2 + Room.ObjectMeshEdges[Sequence[order]][pathid].Polys.Count];
                    //doubly concave
                    dist[0] += (PathVertices[pathid][order + 2] - a).Length();
                    dist[1] += (PathVertices[pathid][order + 2] - b).Length();
                    dist[2] += (PathVertices[pathid][order + 2] - Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].Polys[0].Centroid).Length();
                    dist[3] += (PathVertices[pathid][order + 2] - Room.ObjectMeshEdges[Sequence[order]][pathid % Room.ObjectMeshEdges[Sequence[order]].Length].Polys[1].Centroid).Length();
                }
                else if (K0 < 0 || K1 < 0)
                {
                    if (dist == null) dist = new double[Room.ObjectMeshEdges[Sequence[order]][pathid].Polys.Count];

                    for (int i = 0; i < Room.ObjectMeshEdges[Sequence[order]][pathid].Polys.Count; i++)
                    {
                        int pid = pathid % Room.ObjectMeshEdges[Sequence[order]].Length;
                        Hare.Geometry.Edge e = Room.ObjectMeshEdges[Sequence[order]][pid];
                        Hare.Geometry.Vector p = Room.ObjectMeshEdges[Sequence[order]][pid].TributaryLength[i] * e.Tangents[i];
                        dist[i] += (PathVertices[pathid][order + 2] - (PathVertices[pid][order + 1] + p)).Length();
                    }
                }
                t_limits[0] += dist.Min() / c_sound;
                t_limits[1] += dist.Max() / c_sound;

                //Rotate surface frame to coincide with principle radii of ray.
                int element = pathid % Room.ObjectMeshEdges[Sequence[order]].Length;
                Hare.Geometry.Vector MinFrame = Room.Edge_Frames[Sequence[order]][element][0];
                Hare.Geometry.Vector MaxFrame = Room.Edge_Frames[Sequence[order]][element][1];
                MathNet.Numerics.LinearAlgebra.Double.DenseMatrix K = new MathNet.Numerics.LinearAlgebra.Double.DenseMatrix(2);
                K[0, 0] = K0;
                K[1, 1] = K1;
                Hare.Geometry.Vector Ti = Hare.Geometry.Hare_math.Cross(dir1, Normal);
                double cosphi = Hare_math.Dot(MinFrame, Ti); //TODO: figure out how to make sure that non-generic wavefront frames also conform to the surface frame.
                double sinphi = Math.Sqrt(1 - cosphi * cosphi); //TODO: is this the fastest way to convert to sin(phi)?

                if (W_Kurvature[0, 1] != 0)
                {
                    //Rotate Wavefront
                    double cosphi2 = Hare_math.Dot(Ti, W_Frame[0]);
                    double sinphi2 = Math.Sqrt(1 - cosphi2 * cosphi2); //TODO: is this the fastest way to convert to sin(phi)?
                    //rotate wavefont by phi.
                    DenseMatrix wmod1 = new DenseMatrix(2);
                    DenseMatrix wmod2 = new DenseMatrix(2);
                    wmod1[0, 0] = wmod1[1, 1] = wmod2[0, 0] = wmod2[1, 1] = cosphi2;
                    wmod1[0, 1] = wmod2[1, 0] = -sinphi2;
                    wmod2[0, 1] = wmod1[1, 0] = sinphi2;
                    W_Kurvature = wmod1 * W_Kurvature * wmod2;
                }

                W_Frame[0] = Ti;
                W_Frame[1] = Hare_math.Cross(Ti, dir1);

                //Perform operation to modify wavefront curvature per surface shape.
                DenseMatrix mod1 = new DenseMatrix(2);
                DenseMatrix mod2 = new DenseMatrix(2);
                mod1[0, 0] = (double)cosphi;
                mod1[1, 1] = (double)cosphi;
                mod2[0, 0] = (double)cosphi;
                mod2[1, 1] = (double)cosphi;
                mod1[0, 1] = -sinphi;
                mod2[1, 0] = -sinphi;
                mod2[0, 1] = (double)sinphi;
                mod1[1, 0] = (double)sinphi;
                K = mod1 * K * mod2;
                //TODO: Simplify this to a faster method not involving the matrix classes...
                double costheta = Hare.Geometry.Hare_math.Dot(dir1, Normal);

                K[0, 0] *= Math.Abs(1.0 / costheta);
                K[1, 1] *= Math.Abs(costheta);
                //if (costheta < 0) K *= -1;

                W_Kurvature[0, 1] *= -1;
                W_Kurvature[1, 0] *= -1;
                W_Kurvature += 2 * K;
                W_Kurvature = project_wavefront(l, W_Kurvature);

                if (K[0, 0] < 0 && K[1, 1] < 0)
                {
                    if (Room.Edge_Kurvatures[Sequence[order]][element][0] < 0 && Room.Edge_Kurvatures[Sequence[order]][element][1] < 0)
                    {
                        double A = 0;
                        for (int i = 0; i < Room.ObjectMeshEdges[Sequence[order]][element].TributaryArea.Count; i++)
                        {
                            A += Room.ObjectMeshEdges[Sequence[order]][element].TributaryArea[i];
                        }
                        h *= A;
                    }
                    else
                    {
                        double L = 0;
                        for (int i = 0; i < Room.ObjectMeshEdges[Sequence[order]][element].TributaryLength.Count; i++)
                        {
                            L += Room.ObjectMeshEdges[Sequence[order]][element].TributaryLength[i];
                        }
                        h *= L;
                    }
                }

                Hare.Geometry.Vector dir2 = (PathVertices[pathid][order + 1] - PathVertices[pathid][order + 2]);
                dir2.Normalize();
                Hare.Geometry.Vector local_N = dir1 + dir2;
                local_N.Normalize();
                W_Frame[0] -= local_N * Hare_math.Dot(W_Frame[0], local_N) * 2;
                W_Frame[1] -= local_N * Hare_math.Dot(W_Frame[1], local_N) * 2;
            }

            if (!hasedge) h *= Math.Sqrt(Math.Abs((1 / (4 * Math.PI)) * W_Kurvature.Determinant()));
            else h *= Math.Sqrt(Math.Abs(W_Kurvature.Determinant()));
            return h;
            //return h * Math.Sqrt(Math.Abs(W_Kurvature.Determinant()));
        }

        private MathNet.Numerics.LinearAlgebra.Double.DenseMatrix project_wavefront(double dist, MathNet.Numerics.LinearAlgebra.Double.DenseMatrix M)
        {
            //Do we need to rotate the frame so that it is oriented to principle axes before projecting?
            //Yes, we should rotate it. Not rotating can introduce several meters' worth in error...

            //bool rotate = false;
            //DenseMatrix wmod1 = new DenseMatrix(2);
            //DenseMatrix wmod2 = new DenseMatrix(2);
            //DenseMatrix M_p = new DenseMatrix(2);
            ////Rotate vectors (will be necessary to do more than first order).
            //if (M[0, 1] != 0)
            //{
            //    double phi = Math.Atan2(M[0, 0], M[1, 0]);
            //    double cosphi = Math.Cos(phi);
            //    double sinphi = Math.Sqrt(1 - cosphi * cosphi);// * Math.Sign(cosphi);
            //    rotate = true;
            //    //double phi = Math.Atan2(M[0, 0], M[1, 0]);
            //    //double cosphi = Math.Cos(phi);
            //    //double sinphi = Math.Sqrt(1 - cosphi * cosphi);// * Math.Sign(cosphi);
            //    //rotate wavefont by phi.
            //    wmod1[0, 0] = wmod1[1, 1] = wmod2[0, 0] = wmod2[1, 1] = cosphi;
            //    wmod1[0, 1] = wmod2[1, 0] = -sinphi;
            //    wmod2[0, 1] = wmod1[1, 0] = sinphi;
            //    M_p = wmod1 * M * wmod2;
            //}

            M[0, 0] = M[0, 0] == double.PositiveInfinity ? 0 : 1 / M[0, 0];
            M[1, 1] = M[1, 1] == double.PositiveInfinity ? 0 : 1 / M[1, 1];
            M[0, 1] = M[0, 1] == double.PositiveInfinity ? 0 : 1 / M[0, 1];
            M[1, 0] = M[1, 0] == double.PositiveInfinity ? 0 : 1 / M[1, 0];
            M[0, 0] += dist;
            M[1, 1] += dist;
            //M[0, 1] += dist;
            //M[1, 0] += dist;
            M[0, 0] = M[0, 0] == double.PositiveInfinity ? 0 : 1 / M[0, 0];
            M[1, 1] = M[1, 1] == double.PositiveInfinity ? 0 : 1 / M[1, 1];
            M[0, 1] = M[0, 1] == double.PositiveInfinity ? 0 : 1 / M[0, 1];
            M[1, 0] = M[1, 0] == double.PositiveInfinity ? 0 : 1 / M[1, 0];

            //M[0, 0] = 1 + dist * M[0, 0];
            //M[1, 1] = 1 + dist * M[1, 1];
            //M[0, 0] *= dist;
            //M[1, 1] *= dist;

            //M[0, 0] = dist + dist * M[0, 0];
            //M[1, 1] = dist + dist * M[1, 1];
            //M[0, 0] *= dist;
            //M[1, 1] *= dist;

            //if (rotate) M = wmod2 * M * wmod1;

            return M;
        }

        public double[][] B_Function_NoInterp(List<Point[]> PathVertices, int[] Sequence, double[][] intensity_mod, double c_sound, int Threadid, out double[][] Times, out int[] ElementIDs)
        {
            List<double[]> Bs = new List<double[]>();
            List<int> E_ids = new List<int>();

            //Calculate each leg of the compound path (if a composite of multiple compound paths...
            List<double> Tmin = new List<double>(), Tmax = new List<double>();
            for (int i = 0; i < PathVertices.Count; i++)
            {
                if (PathVertices[i] == null) continue;
                E_ids.Add(i);
                double[] t_limits = new double[2];
                DenseMatrix W_orig = new DenseMatrix(2); W_orig[0, 0] = double.Epsilon; W_orig[1, 1] = double.Epsilon;

                Hare.Geometry.Vector[] W_Frame = new Hare.Geometry.Vector[2];
                double p = Math.Abs(Power_Recursion(i, Sequence.Length - 1, PathVertices, Sequence, c_sound, t_limits, Threadid, ref W_orig, ref W_Frame));
                double[] B_oct = new double[BandCount];
                for (int oct = 0; oct < BandCount; oct++)
                {
                    double mod = Math.Sqrt(intensity_mod[i][oct] * Math.Pow(10, -.1 * Room.Attenuation(0, ThirdOctave)[oct] * Room.Sound_speed(0) * 0.5 * (t_limits[0] + t_limits[1])));
                    B_oct[oct] = p * mod;
                }

                Bs.Add(B_oct);
                Tmin.Add(Math.Min(t_limits[0], t_limits[1]));
                Tmax.Add(Math.Max(t_limits[0], t_limits[1]));
            }

            ElementIDs = E_ids.ToArray();
            Times = new double[2][] { Tmin.ToArray(), Tmax.ToArray() };
            return Bs.ToArray();
        }

        public double[][][] H_Function(double[][] Bs, List<Point[]> PathVertices, int[] Sequence, int[] elementIds, double[][] Times, out double mintime, int threadid)
        {
            mintime = Times[0].Min();
            double maxtime = Times[1].Max();
            double[][][] H_directional = new double[7][][];

            PathVertices.RemoveAll(item => item == null);
            int order = PathVertices[0].Length - 3;

            if (Sequence[0] > Room.ObjectCount - 1)
            {
                //Split energy over straddled bins.
                for (int d = 0; d < 7; d++)
                {
                    H_directional[d] = new double[BandCount][];
                    for (int oct = 0; oct < BandCount; oct++) H_directional[d][oct] = new double[(int)Math.Ceiling((maxtime - mintime) * 44100)];
                }

                for (int i = 0; i < Times[0].Length; i++)
                {
                    Hare.Geometry.Vector DIR = PathVertices[i][Sequence.Length] - PathVertices[i][Sequence.Length + 1];
                    DIR.Normalize();
                    Hare.Geometry.Vector emission = PathVertices[i][1] - PathVertices[i][0];
                    emission.Normalize();
                    double[] dmod = (ThreadSources == null ? Src : ThreadSources[threadid]).DirPower(threadid, Rnd[threadid].Next(), emission).Select(power => AcousticalMath.Pressure_Intensity(power, Room.Rho_C(0))).ToArray();

                    double delta = Times[1][i] - Times[0][i]; //duration of the power receipt.
                    double t0 = Times[0][i] - mintime; //start time of the power receipt.
                    double t1 = t0 + delta; //end time of the power receipt.
                    int s0 = (int)Math.Floor(t0 * 44100.00);//start index
                    int s1 = (int)Math.Floor(t1 * 44100.00);//start index
                    double frac0, frac1;
                    if (s0 == s1) { frac0 = 0; frac1 = 1; }
                    else { frac1 = (t1 - (double)s1 / 44100.00) / delta; frac0 = 1 - frac1; }

                    for (int oct = 0; oct < BandCount; oct++)
                    {
                        double p1 = frac1 * Bs[i][oct] * dmod[oct], p0 = frac0 * Bs[i][oct] * dmod[oct];
                        if (p1 > 0)
                        {
                            H_directional[0][oct][s1] += p1;
                            H_directional[DIR.dx > 0 ? 1 : 2][oct][s1] += Math.Sqrt(Math.Abs(DIR.dx)) * p1;
                            H_directional[DIR.dy > 0 ? 3 : 4][oct][s1] += Math.Sqrt(Math.Abs(DIR.dy)) * p1;
                            H_directional[DIR.dz > 0 ? 5 : 6][oct][s1] += Math.Sqrt(Math.Abs(DIR.dz)) * p1;
                        }
                        if (p0 > 0)
                        {
                            H_directional[0][oct][s0] += p0;
                            H_directional[DIR.dx > 0 ? 1 : 2][oct][s0] += Math.Sqrt(Math.Abs(DIR.dx)) * p0;
                            H_directional[DIR.dy > 0 ? 3 : 4][oct][s0] += Math.Sqrt(Math.Abs(DIR.dy)) * p0;
                            H_directional[DIR.dz > 0 ? 5 : 6][oct][s0] += Math.Sqrt(Math.Abs(DIR.dz)) * p0;
                        }
                    }
                }
                double[] H = H_directional[0][4];
            }
            else
            {
                //for curved and planar surfaces
                int l = Math.Max(1, (int)(Math.Ceiling((Times[1].Max() - Times[0].Min()) * 44100)));
                for (int j = 0; j < 7; j++)
                {
                    H_directional[j] = new double[BandCount][];
                    for(int oct = 0; oct < BandCount; oct++) H_directional[j][oct] = new double[l];
                }
                for (int i = 0; i < PathVertices.Count; i++)
                {
                    if (PathVertices[i] == null) continue;
                    bool k1, k2;
                    Hare.Geometry.Vector DIR = PathVertices[i][Sequence.Length] - PathVertices[i][Sequence.Length + 1];
                    DIR.Normalize();
                    Hare.Geometry.Vector emission = PathVertices[i][1] - PathVertices[i][0];
                    emission.Normalize();
                    double[] dmod = (ThreadSources == null ? Src : ThreadSources[threadid]).DirPower(threadid, Rnd[threadid].Next(), emission).Select(power => AcousticalMath.Pressure_Intensity(power, Room.Rho_C(0))).ToArray();

                    ////
                    if (Room.Edge_Normals[Sequence[PathVertices[i].Length - 3]].Length > 0)
                    {
                        Hare.Geometry.Vector d0 = PathVertices[i][order + 1] - PathVertices[i][order];
                        d0.Normalize();
                        double dir = Hare.Geometry.Hare_math.Dot(d0, Room.Edge_Normals[Sequence[order]][elementIds[i]]);
                        //
                        if (dir > 0) { k1 = Room.Edge_Kurvatures[Sequence[order]][elementIds[i]][0] * -1 < 0; k2 = Room.Edge_Kurvatures[Sequence[order]][elementIds[i]][1] * -1 < 0; }
                        else { k1 = Room.Edge_Kurvatures[Sequence[order]][elementIds[i]][0] * -1 < 0; k2 = Room.Edge_Kurvatures[Sequence[order]][elementIds[i]][1] * -1 < 0; }
                    }
                    else
                    {
                        k1 = false; k2 = false;
                    }

                    if (k1 || k2)
                    {
                        //Curved surface
                        int samplenbr = (int)Math.Ceiling(Math.Abs(Times[1][i] - Times[0][i]) * 44100);

                        double dt = 1.0 / 44100;
                        double[][] BFct = (k1 && k2) ? EdgeDistribution_2d(samplenbr, Room.ObjectMeshEdges[Sequence[order]][elementIds[i]], Bs[i], PathVertices[i][0], PathVertices[i][1], PathVertices[i][2]) : EdgeDistribution_1d(samplenbr, Room.ObjectMeshEdges[Sequence[order]][elementIds[i]], Bs[i], Times[0][i], Times[1][i], PathVertices[i][0], PathVertices[i][1], PathVertices[i][2]);

                        double x_2 = Math.Sqrt(Math.Abs(DIR.dx)); double y_2 = Math.Sqrt(Math.Abs(DIR.dy)); double z_2 = Math.Sqrt(Math.Abs(DIR.dz));
                        double[] x_out = new double[BFct.Length]; double[] y_out = new double[BFct.Length]; double[] z_out = new double[BFct.Length];
                        double[] T = new double[BFct.Length];

                        for (int j = 0; j < BFct.Length; j++)
                        {
                            for (int oct = 0; oct < BandCount; oct++)
                            {
                                double t = Times[0][i] + dt * j - mintime;
                                double pressure = BFct[j][oct] * dmod[oct];
                                H_directional[0][oct][(int)(t / dt)] += pressure;
                                if (DIR.dx > 0) H_directional[1][oct][(int)(t / dt)] += pressure * x_2;
                                else H_directional[2][oct][(int)(t / dt)] += pressure * x_2;
                                if (DIR.dy > 0) H_directional[3][oct][(int)(t / dt)] += pressure * y_2;
                                else H_directional[4][oct][(int)(t / dt)] += pressure * y_2;
                                if (DIR.dz > 0) H_directional[5][oct][(int)(t / dt)] += pressure * z_2;
                                else H_directional[6][oct][(int)(t / dt)] += pressure * z_2;
                            }
                        }
                    }
                    else
                    {
                        //Planar surface
                        double t = (Times[0][i] + Times[1][i]) * 0.5;

                        int s = (int)Math.Floor(44100 * (t - mintime));
                        for (int oct = 0; oct < BandCount; oct++)
                        {
                            double pressure = Bs[i][oct] * dmod[oct];
                            H_directional[0][oct][s] += pressure; //Omni Channel...
                            double x_2 = Math.Sqrt(Math.Abs(DIR.dx)); double y_2 = Math.Sqrt(Math.Abs(DIR.dy)); double z_2 = Math.Sqrt(Math.Abs(DIR.dz));
                            if (DIR.dx > 0) H_directional[1][oct][s] += pressure * x_2;
                            else H_directional[2][oct][s] += pressure * x_2;
                            if (DIR.dy > 0) H_directional[3][oct][s] += pressure * y_2;
                            else H_directional[4][oct][s] += pressure * y_2;
                            if (DIR.dz > 0) H_directional[5][oct][s] += pressure * z_2;
                            else H_directional[6][oct][s] += pressure * z_2;
                        }
                    }
                }
            }

            return H_directional;
        }

        public double[][] EdgeDistribution_1d(int samplect, Hare.Geometry.Edge e, double[] power, double Tmin, double Tmax, Hare.Geometry.Point p1, Hare.Geometry.Point p1_5, Hare.Geometry.Point p2)
        {
            int samples = Math.Max((int)Math.Floor((Tmax - Tmin) * 44100), 1);
            double[][] result = new double[samples][];

            for (int j = 0; j < samples; j++)
            {
                result[j] = new double[BandCount];
                for(int oct = 0; oct < BandCount; oct++) result[j][oct] += power[oct] / samples;
            }

            return result;
        }

        public double[][] EdgeDistribution_2d(int samplect, Hare.Geometry.Edge e, double[] power, Hare.Geometry.Point p1, Hare.Geometry.Point p1_5, Hare.Geometry.Point p2)
        {
            samplect = Math.Max(samplect, 1);
            double start = ((p1 - p1_5).Length() + (p1_5 - p2).Length()) / Speed_of_Sound;
            double e1 = ((p1 - e.a).Length() + (p2 - e.a).Length()) / Speed_of_Sound;
            double e2 = ((p1 - e.b).Length() + (p2 - e.b).Length()) / Speed_of_Sound;
            double minT = Math.Min(e1, e2);
            double maxT = Math.Max(e1, e2);
            double[] intermediates = new double[e.Polys.Count];
            int[] sample = new int[e.Polys.Count];
            int[] upperct = new int[e.Polys.Count];
            double[] result = new double[samplect];
            double max = 0;
            bool invert = false;
            double denominator = Math.Abs(e1 - e2);

            for (int i = 0; i < e.Polys.Count; i++)
            {
                intermediates[i] = ((p1 - e.Polys[i].Centroid).Length() + (p2 - e.Polys[i].Centroid).Length()) / Speed_of_Sound;
                double normalized = denominator > 1E-12 ? (intermediates[i] - minT) / denominator : 0.5;
                normalized = Math.Max(0, Math.Min(1, normalized));
                sample[i] = samplect == 1 ? 0 : (int)Math.Round(normalized * (samplect - 1));
                sample[i] = Math.Max(0, Math.Min(samplect - 1, sample[i]));
                upperct[i] = Math.Max(1, samplect - sample[i]);
                invert |= sample[i] < minT || sample[i] > maxT;
            }

            for (int i = 0; i < e.Polys.Count; i++)
            {
                for (int j = 0; j < samplect; j++)
                {
                    double s;
                    if (j < sample[i])
                    {
                        int lowerct = Math.Max(1, sample[i]);
                        s = e.TributaryLength[i] * j / lowerct;
                    }
                    else
                    {
                        s = e.TributaryLength[i] * (samplect - j) / upperct[i];
                    }
                    result[j] += s;
                    max += s;
                }
            }

            double[][] result2 = new double[result.Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result2[i] = new double[BandCount];
                if (max <= 0) continue;
                double shape = Math.Sqrt(Math.Max(0, result[i]) / max);
                for (int oct = 0; oct < BandCount; oct++) result2[i][oct] = Math.Abs(power[oct]) * shape;
            }

            return result2;
        }



        /// <summary>
        /// Processes image source paths with input of only a sequence of indices.
        /// </summary>
        /// <param name="Sequence">The input sequence of surface indices.</param>
        /// <param name="Threadid">The id of the calling thread.</param>
        /// <param name="rec_id">The id of the receiver.</param>
        /// <returns>True if a valid path, false if not.</returns>
        private bool ProcessPath(int[] Sequence, int Threadid, int rec_id)
        {
            Hare.Geometry.Point RefPoint;
            Hare.Geometry.Point NextPoint = new Hare.Geometry.Point();
            Hare.Geometry.Point[] Images = new Hare.Geometry.Point[Sequence.Length];
            RefPoint = ThreadSources[Threadid].Origin;
            int[] Seq_Polys = new int[Sequence.Length];
            double[] Trans_Mod = new double[BandCount];
            for (int t_oct = 0; t_oct < BandCount; t_oct++) Trans_Mod[t_oct] = 1;

            //Find all Source Images
            for (int q = 0; q < Sequence.Length; q++)
            {
                Images[q] = Room.Image(RefPoint, 0, Room.ObjectMembers[Sequence[q]][0]);
                RefPoint = new Point( Images[q].x, Images[q].y, Images[q].z);
            }

            int[] s = new int[Sequence.Length];
            List<int[]> sub = new List<int[]>();
            for (int i = 0; i < Sequence.Length; i++) s[i] = -1;
            //im[Order] = Room.Image(Order > 0 ? Images0[r][Images0[r].Length - 1] : this.Src.Origin, 0, Room.ObjectMembers[Sequence[Order]][0]);
            //s[Order] = -1;//no sub needed...
            //Images.Add(im);
            sub.Add(s);

            ProcessImages(new Point[1][] { Images }, Sequence, Threadid, sub, rec_id);
            //ProcessImages(Images, Sequence, rec_id, Threadid);

            return true;
        }

        /// <summary>
        /// Checks a path for occlusions in the model.
        /// </summary>
        /// <param name="Src"></param>
        /// <param name="EndPt"></param>
        /// <param name="Poly_X"></param>
        /// <param name="X_Point"></param>
        /// <param name="Thread_Id">The id of the calling thread.</param>
        /// <returns></returns>
        private bool FinalOcclusion(Hare.Geometry.Point Src, Hare.Geometry.Point EndPt, int Poly_X, ref double[] Trans_Mod, int Thread_Id)
        {
            Hare.Geometry.Vector D = (Src - EndPt);
            double L = D.Length();
            double L2 = 0;
            D.Normalize();
            Ray R = new Ray(EndPt, D, Thread_Id, Rnd[Thread_Id].Next());

            Hare.Geometry.X_Event X;
            do
            {
                if (Room.shoot(R, 0, out X) && Room.IsTransmissive[X.Poly_id] && Poly_X != Room.ObjectID(X.Poly_id) && X.t < L)
                {
                    ///The ray hit something transparent. Account for this...
                    double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad(ThirdOctave);
                    for (int oct = 0; oct < BandCount; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * Room.TransmissionSpectrum(X.Poly_id, ThirdOctave)[oct];// * (1-Scattering[oct]);
                    R.x = X.X_Point.x;
                    R.y = X.X_Point.y;
                    R.z = X.X_Point.z;
                    R.Ray_ID = Rnd[Thread_Id].Next();
                    L2 += X.t;
                    continue;
                }
                break;
            } while (true);

            ///If it hits nothing, then there is nothing occluding... IsOccluded = false
            if (!X.Hit) return false;

            ///If the thing it hit is closer than the source, then it is occluded... IsOccluded = true
            L2 += X.t;
            if (L2 < L) return true;

            ///If we got this far, then there is nothing occluding... 
            return false;
        }

        //Edge Diffraction Version of the Final Occlusion Check.
        private bool FinalOcclusion(Hare.Geometry.Point Src, Hare.Geometry.Point EndPt, List<int> polyids, ref double[] Trans_Mod, int Thread_Id)
        {
            Hare.Geometry.Vector D = (EndPt - Src);
            double L = D.Length();
            L *= .99; //For tolerance...
            double L2 = 0;
            D /= L;
            Ray R = new Ray(Src, D, Thread_Id, Rnd[Thread_Id].Next());

            Hare.Geometry.X_Event X;
            do
            {
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Change this so that it checks for the edge the BREPS the edge is on when x.t is less than tolerance.
                ////////////////////////////////////////////////////////////////////////////////////////////////////////

                if (Room.shoot(R, 0, out X, polyids[0], polyids.Count > 1 ? polyids[1] : -1))
                {
                    if (Room.IsTransmissive[X.Poly_id])//&& (Room.BrepID(X.Poly_id) == Poly_1 || Room.BrepID(X.Poly_id) == Poly_2)
                    {
                        ///The ray hit something transparent. Account for this...
                        double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad(ThirdOctave);
                        for (int oct = 0; oct < BandCount; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * Room.TransmissionSpectrum(X.Poly_id, ThirdOctave)[oct];
                        R.x = X.X_Point.x;
                        R.y = X.X_Point.y;
                        R.z = X.X_Point.z;
                        R.Ray_ID = Rnd[Thread_Id].Next();
                        L2 += X.t;
                        continue;
                    }
                    //else if (X.t < tol)
                    else if (polyids.Contains(Room.PlaneID(X.Poly_id)))
                    {
                        R.x = X.X_Point.x;
                        R.y = X.X_Point.y;
                        R.z = X.X_Point.z;
                        R.Ray_ID = Rnd[Thread_Id].Next();
                        L2 += X.t;
                        continue;
                    }
                }

                break;
            } while (true);

            ///If it hits nothing, then there is nothing occluding... IsOccluded = false
            if (!X.Hit) 
                return false;

            ///If the thing it hit is closer than the source, then it is occluded... IsOccluded = true
            L2 += X.t;
            if (L2 < L) 
                return true;

            ///If we got this far, then there is nothing occluding... 
            return false;
        }

        /// <summary>
        /// Checks a path for occlusions in the model.
        /// </summary>
        /// <param name="Origin"></param>
        /// <param name="EndPt"></param>
        /// <param name="Poly_X"></param>
        /// <param name="X_Point"></param>
        /// <param name="Thread_Id">The id of the calling thread.</param>
        /// <returns></returns>
        private bool OcclusionIntersect(Hare.Geometry.Point Origin, Hare.Geometry.Point EndPt, int Poly_X, ref double[] Trans_Mod, ref Hare.Geometry.Point X_Point, ref int Poly_Seq, int Thread_Id)
        {
            Hare.Geometry.Vector D = (EndPt - Origin);
            double L = D.Length();
            //L *= .99; //For tolerance...
            D /= L;
            Ray R = new Ray(Origin, D, Thread_Id, Rnd[Thread_Id].Next());

            Hare.Geometry.X_Event X = new Hare.Geometry.X_Event();
            do
            {
                if (Room.shoot(R, 0, out X) && Room.IsTransmissive[X.Poly_id] && Poly_X != Room.ObjectID(X.Poly_id))
                {
                    double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad(ThirdOctave);
                    for (int oct = 0; oct < BandCount; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * (Room.TransmissionSpectrum(X.Poly_id, ThirdOctave)[oct]);
                    R.x = X.X_Point.x;
                    R.y = X.X_Point.y;
                    R.z = X.X_Point.z;
                    R.Ray_ID = Rnd[Thread_Id].Next();
                    continue;
                }
                break;
            } while (true);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (X.t < 0.000001)
                return false;
            if (Poly_X != Room.ObjectID(X.Poly_id))
            {
                return false;
            }
            Poly_Seq = X.Poly_id;
            X_Point = X.X_Point;
            return true;
        }

        /// <summary>
        /// Checks a path for occlusions in the model.
        /// </summary>
        /// <param name="Origin"></param>
        /// <param name="EndPt"></param>
        /// <param name="Poly_X"></param>
        /// <param name="X_Point"></param>
        /// <param name="Thread_Id">The id of the calling thread.</param>
        /// <returns></returns>
        private bool OcclusionIntersectED(Hare.Geometry.Point Origin, Hare.Geometry.Point EndPt, int Poly_X, int sub, ref double[] Trans_Mod, int Thread_Id)
        {
            Hare.Geometry.Vector D = Origin - EndPt;
            double L = D.Length();
            //L *= .99; //For tolerance...
            D /= L;
            Ray R = new Ray(EndPt, D, Thread_Id, Rnd[Thread_Id].Next());
            Hare.Geometry.X_Event X = new Hare.Geometry.X_Event();

            List<int> planeids = new List<int>();
            if (Poly_X > Room.ObjectCount - 1)
            {
                //Edge Diffraction Case
                //for (int p = 0; p < Room.Raw_Edge[Poly_X - Room.ObjectCount].Polys.Count; p++) planeids.Add(Room.Edge_Nodes[Poly_X - Room.ObjectCount].Polys[p].Plane_ID);
                planeids = new List<int>() { Room.Edge_Nodes[Poly_X - Room.ObjectCount].EdgeSources[sub].Poly_1, Room.Edge_Nodes[Poly_X - Room.ObjectCount].EdgeSources[sub].Poly_2 };
            }
            else
            {
                //Curved Surface Source Case
                //throw new Exception("Used an edge diffraction occlusion check for something other than edge diffraction...");
                planeids = new List<int> { Room.ObjectMeshEdges[Poly_X][sub].Polys[0].Poly_ID, Room.ObjectMeshEdges[Poly_X][sub].Polys[1].Poly_ID};
            }

            do
            {
                if (Room.shoot(R, 0, out X, planeids[0], planeids.Count > 1 ? planeids[1] : -1) && Room.IsTransmissive[X.Poly_id])// && Poly_X != Room.PlaneID(X.Poly_id))
                {
                    double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad(ThirdOctave);
                    for (int oct = 0; oct < BandCount; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * Room.TransmissionSpectrum(X.Poly_id, ThirdOctave)[oct];
                    R.x = X.X_Point.x;
                    R.y = X.X_Point.y;
                    R.z = X.X_Point.z;
                    R.Ray_ID = Rnd[Thread_Id].Next();
                    L -= X.t;
                    continue;
                }
                break;
            } while (true);

            if (X.t < L - 0.0001 && X.t != 0) 
                return false;

            return true;
        }

        /// <summary>
        /// Write the image source paths to a file for storage.
        /// </summary>
        /// <param name="BW">The binary writer to which data will be written.</param>
        public void Write_Data(ref System.IO.BinaryWriter BW)
        {
            //1. Write an indicator to signify that there is image source data
            bool coherent = ValidPaths.Any(paths => paths.Any(path => path.CoherentSpeaker));
            BW.Write(coherent ? "Image-Source_Data coherent" : "Image-Source_Data");

            for (int q = 0; q < ValidPaths.Length; q++)
            {
                //2. Write the receiver number:int
                BW.Write(q);
                //3. Write number of paths:int
                BW.Write(ValidPaths[q].Count);
                //Write for all specular paths...
                //TODO: Sort for Specular Paths, filter out Compound Paths.
                for (int i = 0; i < ValidPaths[q].Count; i++)
                {

                    if (coherent)
                    {
                        BW.Write(ValidPaths[q][i].CoherentSpeaker);
                        double[] delays = ValidPaths[q][i].SteeringDelaySeconds ?? new double[] { 0 };
                        BW.Write(delays.Length);
                        foreach (double delay in delays) BW.Write(delay);
                    }
                    //V2 procedure: Introducing pressure based simulation...
                    //3a.1 Write an enumeration for the kind of reflection it is(0 for specular, 1 for compound). (int)
                    if (ValidPaths[q][i] is Specular_Path)
                    {
                        ///Specular Path:
                        BW.Write((short)0);

                        //4. Write the number of reflection path points
                        BW.Write(ValidPaths[q][i].Path[0].Length);
                        //5. Write the reflection path:double
                        for (int r = 0; r < ValidPaths[q][i].Path[0].Length; r++)
                        {
                            BW.Write(ValidPaths[q][i].Path[0][r].x);
                            BW.Write(ValidPaths[q][i].Path[0][r].y);
                            BW.Write(ValidPaths[q][i].Path[0][r].z);
                        }

                        //6a.1 Write the energy levels
                        for (int band = 0; band < BandCount; band++) BW.Write(ValidPaths[q][i].Energy(band, 44100)[0]);

                        //6a.2 Write a bool for whether it has a special materials filter...
                        if ((ValidPaths[q][i] as Specular_Path).Special_Filter != null)
                        {
                            BW.Write(true);
                            BW.Write((ValidPaths[q][i] as Specular_Path).Special_Filter.Length);
                            foreach (System.Numerics.Complex val in (ValidPaths[q][i] as Specular_Path).Special_Filter)
                            {
                                BW.Write(val.Real);
                                BW.Write(val.Imaginary);
                            }
                        }
                        else
                        {
                            BW.Write(false);
                        }
                        double[] prms = (ValidPaths[q][i] as Specular_Path).prms;
                        for (int j = 0; j < prms.Length; j++) BW.Write(prms[j]);
                        //7. Write the arrival time:double
                        BW.Write(ValidPaths[q][i].TravelTime);

                        //8. Write the Reflection Sequence:int
                        for (int r = 0; r < ValidPaths[q][i].Reflection_Sequence.Length; r++)
                        {
                            BW.Write(ValidPaths[q][i].Reflection_Sequence[r]);
                        }
                    }
                    else if (ValidPaths[q][i] is Compound_Path)
                    {
                        ///Specular Path:
                        BW.Write((short)1);
                        ///Compound Path:
                        //4. Write the number of reflection path points
                        //a. number of sub-paths
                        BW.Write(ValidPaths[q][i].Path.Length);
                        //Hare.Geometry.Point[][] PTS = new Hare.Geometry.Point[BR.ReadInt32()][];
                        //b. number of points in each sub-path
                        BW.Write(ValidPaths[q][i].Path[0].Length);
                        //int ptct = BR.ReadInt32();

                        //5. Write the reflection path:double
                        for (int r = 0; r < ValidPaths[q][i].Path.Length; r++)
                        {
                            for (int s = 0; s < ValidPaths[q][i].Path[r].Length; s++)
                            {
                                BW.Write(ValidPaths[q][i].Path[r][s].x);
                                BW.Write(ValidPaths[q][i].Path[r][s].y);
                                BW.Write(ValidPaths[q][i].Path[r][s].z);
                            }
                        }

                        //6a. Write the energy values
                        //for (int r = 0; r < ValidPaths[q][i].Path.Length; r++)
                        //{
                        //    for (int s = 0; s < ValidPaths[q][i].Path[r].Length; s++)
                        //    {
                        //        for (int oct = 0; oct < BandCount; oct++) BW.Write((ValidPaths[q][i] as Compound_Path).PathEnergy[s][oct]);
                        //    }
                        //}

                        //6a.2 Write a bool for whether it has a special material. s filter...
                        BW.Write(false);
                        //bool Special_Filter = BR.ReadBoolean();//Defunct for now - should be false.

                        //6a.3 Write a bool for whether it has a special materials filter...
                        //bool OctavePower = BR.ReadBoolean();
                        BW.Write(false);
                        //6a.3.1 Write Lx8 octave power mod factors...
                        //foreach (double[] powermod in (ValidPaths[q][i] as Compound_Path).Octave_Power)
                        //{
                        //    foreach (double w in powermod) BW.Write(w);
                        //}

                        //7. Write the arrival time:double
                        BW.Write(ValidPaths[q][i].TravelTime);

                        //8. Write the Reflection Sequence:int
                        for (int r = 0; r < ValidPaths[q][i].Reflection_Sequence.Length; r++)
                        {
                            BW.Write(ValidPaths[q][i].Reflection_Sequence[r]);
                        }

                        //9.1 Write the H function (omni first)
                        //9.1a - Write the number of points
                        Compound_Path compound = ValidPaths[q][i] as Compound_Path;
                        int filterlength = compound.H[0].Length;
                        BW.Write(filterlength);
                        for (int time = 0; time < filterlength; time++)
                            for (int band = 0; band < BandCount; band++) BW.Write(compound.H[band][time]);
                        for (int d = 0; d < 6; d++)
                            for (int time = 0; time < filterlength; time++)
                                for (int band = 0; band < BandCount; band++) BW.Write(compound.Hdir[d][band][time]);

                    }
                }
            }
        }

        /// <summary>
        /// Constructor which takes a Binary Reader at the appropriate point, from which calculated data will be extracted.
        /// </summary>
        /// <param name="BR"></param>
        /// <param name="Rec_CT"></param>
        /// <param name="Direct"></param>
        /// <returns></returns>
        public static ImageSourceData Read_Data(ref System.IO.BinaryReader BR, int Rec_CT, Direct_Sound Direct, bool Edges, double[] RhoC, int Src_ID, string version, IProgressFeedback VB)
        {
            return Read_Data(ref BR, Rec_CT, Direct, Edges, RhoC, Src_ID, version, VB, false);
        }

        public static ImageSourceData Read_Data(ref System.IO.BinaryReader BR, int Rec_CT, Direct_Sound Direct, bool Edges, double[] RhoC, int Src_ID, string version, IProgressFeedback VB, bool coherentSpeakers)
        {
            ImageSourceData IS = new ImageSourceData();
            int BandCount = Direct.SWL.Length;
            IS.BandCount = BandCount;
            IS.ValidPaths = new List<Deterministic_Reflection>[Rec_CT];
            IS.Direct_Time = new double[Rec_CT];
            double v = double.Parse(version.Substring(0, 3), System.Globalization.CultureInfo.InvariantCulture);

            for (int q = 0; q < Rec_CT; q++)
            {
                IS.ValidPaths[q] = new List<Deterministic_Reflection>();
                //2. Write the receiver number:int
                BR.ReadInt32();
                //3. Write number of paths:int
                int PathCt = BR.ReadInt32();
                for (int i = 0; i < PathCt; i++)
                {
                    bool coherent = coherentSpeakers && BR.ReadBoolean();
                    double[] delays = null;
                    if (coherentSpeakers)
                    {
                        int count = BR.ReadInt32();
                        if (count != 1 && count != 8 && count != 24) throw new System.IO.InvalidDataException("Invalid speaker-delay band count.");
                        delays = new double[count];
                        for (int band = 0; band < count; band++) delays[band] = BR.ReadDouble();
                    }
                    //3.1 Write the kind of reflection (0 = simple specular reflection, 1 = compound path (made from many paths)
                    int ReflectionType = BR.ReadInt16();
                    if (ReflectionType == 0)
                    {
                        //Specular Reflection
                        //4. Write the number of reflection path points
                        Hare.Geometry.Point[] PTS = new Hare.Geometry.Point[BR.ReadInt32()];

                        //5. Write the reflection path:double
                        for (int r = 0; r < PTS.Length; r++)
                        {
                            PTS[r] = new Hare.Geometry.Point(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble());
                        }

                        //Previously, Pachyderm performed the deterministic part in intensity only...                    
                        //6a. Write the energy values
                        double[] Energy = new double[BandCount];
                        for (int b = 0; b < BandCount; b++) Energy[b] = BR.ReadDouble();

                        bool Special_Filter = BR.ReadBoolean();
                        System.Numerics.Complex[] Filter = null;
                   
                        if (Special_Filter)
                        {
                            //6aa1. Write length of filter...
                            int Filter_Length = BR.ReadInt32();
                            Filter = new System.Numerics.Complex[Filter_Length];
                            //6aa2. Write filter...
                            for(int j = 0; j < Filter.Length; j++)
                            {
                                Filter[j] = new System.Numerics.Complex(BR.ReadDouble(), BR.ReadDouble());
                            }
                            //6aa3. Write octave band root mean square pressure...
                        }

                        double[] prms = new double[BandCount];
                        for (int j = 0; j < prms.Length; j++) prms[j] = BR.ReadDouble();

                        //7. Write the arrival time:double
                        double Time = BR.ReadDouble();

                        //8. Write the Reflection Sequence:int
                        int[] Sequence = new int[PTS.Length - 2];
                        for (int r = 0; r < Sequence.Length; r++)
                        {
                            Sequence[r] = BR.ReadInt32();
                        }

                        IS.ValidPaths[q].Add(new Specular_Path(PTS, Energy, prms, Filter, Time, Sequence, Direct.Min_Time(q), Src_ID)
                        { CoherentSpeaker = coherent, SteeringDelaySeconds = delays });
                    }
                    else if (ReflectionType == 1)
                    {
                        ///Compound Path:
                        //4. Write the number of reflection path points
                        //a. number of sub-paths
                        Hare.Geometry.Point[][] PTS = new Hare.Geometry.Point[BR.ReadInt32()][];
                        //b. number of points in each sub-path
                        int ptct = BR.ReadInt32();

                        //5. Write the reflection path:double
                        for (int r = 0; r < PTS.Length; r++)
                        {
                            PTS[r] = new Point[ptct];
                            for(int s = 0; s < ptct; s++) PTS[r][s] = new Hare.Geometry.Point(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble());
                        }

                        //6a. Write the energy values
                        //double[][] Energy = new double[PTS.Length][];
                        //for (int r = 0; r < PTS.Length; r++)
                        //{
                        //    Energy[r] = new double[BandCount];
                        //    Energy[r][0] = BR.ReadDouble();
                        //    Energy[r][1] = BR.ReadDouble();
                        //    Energy[r][2] = BR.ReadDouble();
                        //    Energy[r][3] = BR.ReadDouble();
                        //    Energy[r][4] = BR.ReadDouble();
                        //    Energy[r][5] = BR.ReadDouble();
                        //    Energy[r][6] = BR.ReadDouble();
                        //    Energy[r][7] = BR.ReadDouble();
                        //}

                        //6a.2 Write a bool for whether it has a special materials filter...
                        bool Special_Filter = BR.ReadBoolean();//Defunct for now - should be false.
                        if (Special_Filter) throw new Exception("Unsupported filter function in comound path... do you need to update your version of Pachyderm?");
                        //if (Special_Filter)
                        //{
                        //    //6aa1. Write length of filter...
                        //    int Filter_Length = BR.ReadInt32();
                        //    Filter = new System.Numerics.Complex[Filter_Length];
                        //    //6aa2. Write filter...
                        //    for (int j = 0; j < Filter.Length; j++)
                        //    {
                        //        Filter[j] = new System.Numerics.Complex(BR.ReadDouble(), BR.ReadDouble());
                        //    }
                        //    //6aa3. Write octave band root mean square pressure...
                        //}

                        //6a.3 Write a bool for whether it has a special materials filter...
                        bool OctavePower = BR.ReadBoolean();
                        double[][] PathPowerMod = new double[PTS.Length][];
                        if (OctavePower)
                        {
                            //6a.3.1 Write Lx8 octave power mod factors...
                            for (int j = 0; j < PTS.Length; j++)
                            {
                                PathPowerMod[j] = new double[BandCount];
                                for (int band = 0; band < BandCount; band++) PathPowerMod[j][band] = BR.ReadDouble();
                            }
                        }

                        //7. Write the arrival time:double
                        double T = BR.ReadDouble();

                        //8. Write the Reflection Sequence:int
                        int[] sequence = new int[ptct - 2];
                        for (int r = 0; r < sequence.Length; r++)
                        {
                            sequence[r] = BR.ReadInt32();
                        }

                        //9.1 Write the H function (omni first)
                        //9.1a - Write the number of points
                        int filterlength = BR.ReadInt32();
                        double[][] H_Omni = new double[filterlength][];
                        //9.1b - Write the H Function (omno, time, octave)
                        for (int k = 0; k < filterlength; k++)
                        {
                            H_Omni[k] = new double[BandCount];
                            for (int oct = 0; oct < BandCount; oct++) H_Omni[k][oct] = BR.ReadDouble();
                        }
                        //9.2 Write the H function (directional 6 channels, time, octave)
                        double[][][] H_Dir = new double[6][][];
                        for (int d = 0; d < H_Dir.Length; d++)
                        {
                            H_Dir[d] = new double[filterlength][];
                            for (int k = 0; k < filterlength; k++)
                            {
                                H_Dir[d][k] = new double[BandCount];
                                for (int oct = 0; oct < BandCount; oct++) H_Dir[d][k][oct] = BR.ReadDouble();
                            }
                        }
                        if (v >= 2.7)
                        {
                            double[][] omni = new double[BandCount][];
                            double[][][] directional = new double[6][][];
                            for (int band = 0; band < BandCount; band++)
                            {
                                omni[band] = new double[filterlength];
                                for (int time = 0; time < filterlength; time++) omni[band][time] = H_Omni[time][band];
                            }
                            for (int d = 0; d < 6; d++)
                            {
                                directional[d] = new double[BandCount][];
                                for (int band = 0; band < BandCount; band++)
                                {
                                    directional[d][band] = new double[filterlength];
                                    for (int time = 0; time < filterlength; time++) directional[d][band][time] = H_Dir[d][time][band];
                                }
                            }
                            H_Omni = omni;
                            H_Dir = directional;
                        }
                        double[][][] H = new double[7][][];
                        H[0] = H_Omni;
                        for(int j =  0; j < H_Dir.GetLength(0); j++)
                        {
                            H[j+1] = H_Dir[j];
                        }

                        IS.ValidPaths[q].Add(new Compound_Path(PTS, Src_ID, sequence, H, T, RhoC[q], ref Direct.Time_Pt[q], 0)
                        { CoherentSpeaker = coherent, SteeringDelaySeconds = delays });
                    }
                }
            }
            IS.Create_Filter(Direct.SWL, 16384, VB);
            return IS;
        }

        public void Set_Power(double[] Factor)
        {
            for(int i = 0; i < Paths.Length; i++) for(int j = 0; j < Paths[i].Count; j++)
            {
                Paths[i][j].Set_Power(Factor);
            }
        }

        public void Create_Filter(double[] SWL, int length, IProgressFeedback VB = null)
        {
            int reflections = 0;
            for (int i = 0; i < ValidPaths.Length; i++) reflections += ValidPaths[i].Count;
            if (reflections == 0) return;

            int completed = 0;
            // Observe worker exceptions on the calling thread instead of terminating the host.
            System.Threading.Tasks.Task work = System.Threading.Tasks.Task.Run(() =>
            {
                for (int i = 0; i < Paths.Length; i++)
                    for (int index = 0; index < Paths[i].Count; index++)
                    {
                        Deterministic_Reflection path = Paths[i][index];
                        try
                        {
                            path.Create_Filter(length, 0);
                        }
                        catch (Exception error)
                        {
                            // Keep the original exception and stack, with the path needed to reproduce it.
                            error.Data["Receiver"] = i;
                            error.Data["Image path"] = index;
                            error.Data["Path type"] = path.GetType().Name;
                            error.Data["Band mode"] = ThirdOctave ? "Third octave" : "Octave";
                            throw;
                        }
                        System.Threading.Interlocked.Increment(ref completed);
                    }
            });

            while (!work.IsCompleted)
            {
                if (VB != null) VB.Report((int)(100.0 * System.Threading.Volatile.Read(ref completed) / reflections));
                System.Threading.Thread.Sleep(50);
            }
            work.GetAwaiter().GetResult();
            if (VB != null) VB.Report(100);
        }

        public Hare.Geometry.Vector[] Dir_Energy(int rec_id, int index, int Octave, double alt, double azi, bool degrees, bool Figure8 = false)
        {
            return ValidPaths[rec_id][index].Dir_Energy(Octave, alt, azi, degrees);
        }

        public Hare.Geometry.Vector[] Dir_Energy(int rec_id, int index, int Octave, Hare.Geometry.Vector V)
        {
            return ValidPaths[rec_id][index].Dir_Energy(Octave, V);
        }

        /// <summary>
        /// returns the list of calculated paths.
        /// </summary>
        public List<Deterministic_Reflection>[] Paths
        {
            get { return ValidPaths; }
        }

        /// <summary>
        /// The number of image source paths for a given receiver.
        /// </summary>
        /// <param name="rec_id">The index of the chosen receiver.</param>
        /// <returns></returns>
        public int Count(int rec_id)
        {
            return ValidPaths[rec_id].Count;
        }
    }

    public abstract class Deterministic_Reflection
    {
        protected string Identifier;
        public bool CoherentSpeaker = false;
        public double[] SteeringDelaySeconds = null;
        public double FilterTime { get { return TravelTime + (CoherentSpeaker && SteeringDelaySeconds != null ? SteeringDelaySeconds.Min() : 0); } }

        public int FilterSample(int sampleFrequency)
        {
            return CoherentSpeaker ? (int)Math.Floor(FilterTime * sampleFrequency) : (int)Math.Ceiling(TravelTime * sampleFrequency);
        }

        protected double[] ApplySourceDelay(double[] filter, int sampleFrequency, int threadid)
        {
            if (!CoherentSpeaker) return filter;
            double minimum = SteeringDelaySeconds == null ? 0 : SteeringDelaySeconds.Min();
            double[] relative = SteeringDelaySeconds == null ? new double[] { 0 } : SteeringDelaySeconds.Select(delay => delay - minimum).ToArray();
            double fraction = FilterTime - Math.Floor(FilterTime * sampleFrequency) / sampleFrequency;
            return Audio.Pach_SP.Fractional_Delay(filter, relative, Math.Max(0, fraction), sampleFrequency, threadid);
        }

        public double[] CoherentBandPressure(int band, int sampleFrequency, double rhoC)
        {
            return CoherentBandPressure(band, sampleFrequency, rhoC, 0);
        }

        // Directional filters use the compound path convention: +X, -X, +Y, -Y, +Z, -Z.
        internal double[] CoherentBandPressure(int band, int sampleFrequency, double rhoC, int dimension)
        {
            if (!CoherentSpeaker) throw new InvalidOperationException("The path is not part of a coherent speaker group.");
            bool third = this is Specular_Path specular ? specular.prms.Length == 24 : (this as Compound_Path).H.Length == 24;
            double[] spectrum = Enumerable.Repeat(120.0, third ? 24 : 8).ToArray();
            Compound_Path compound = this as Compound_Path;
            bool cached = compound != null && sampleFrequency == 44100;
            if (cached && (compound.F == null || (dimension > 0 && compound.Fdir == null))) compound.Create_Filter(16384, 0);
            double[] filter = cached ? dimension == 0 ? compound.F : compound.Fdir[dimension - 1] : Create_Filter(spectrum, sampleFrequency, 16384, dimension, 0);
            double[] reference = this is Specular_Path spec ? Audio.Pach_SP.Filter.Transfer_Function(spec.prms, sampleFrequency, 16384, 0) : cached ? compound.F : dimension == 0 ? filter : Create_Filter(spectrum, sampleFrequency, 16384, 0, 0);
            reference = third ? Audio.Pach_SP.FIR_Third_Bandpass(reference, band, sampleFrequency, 0) : Audio.Pach_SP.FIR_Bandpass(reference, band, sampleFrequency, 0);
            double[] pressure = third ? Audio.Pach_SP.FIR_Third_Bandpass(filter, band, sampleFrequency, 0) : Audio.Pach_SP.FIR_Bandpass(filter, band, sampleFrequency, 0);
            double sum = reference.Sum(value => value * value);
            double target = Energy(band, sampleFrequency).Sum() * rhoC;
            double scale = sum > 0 ? Math.Sqrt(Math.Max(0, target) / sum) : 0;
            if (dimension > 0 && this is Specular_Path)
            {
                Point[] vertices = Path[0];
                Hare.Geometry.Vector arrival = vertices[vertices.Length - 2] - vertices[vertices.Length - 1];
                arrival.Normalize();
                double component = arrival.byint((dimension - 1) / 2);
                scale *= Math.Sqrt(Math.Max(0, dimension % 2 == 1 ? component : -component));
            }
            for (int sample = 0; sample < pressure.Length; sample++) pressure[sample] *= scale;
            return pressure;
        }

        public override string ToString()
        {
            return Identifier;
        }

        public abstract Hare.Geometry.Point[][] Path{get;}

        public abstract double TravelTime { get;}

        public abstract int[] Reflection_Sequence { get; }

        public abstract double[] Energy(int Octave, int Samplefrequency);
        public abstract Hare.Geometry.Vector[] Dir_Energy(int Octave);
        public abstract Hare.Geometry.Vector[] Dir_Energy(int Octave, Hare.Geometry.Vector V);
        public abstract double[] Dir_Energy(int Octave, int dir);
        public abstract Hare.Geometry.Vector[] Dir_Energy(int Octave, double alt, double azi, bool degrees, bool Figure8 = false);

        public abstract Hare.Geometry.Vector[] Dir_EnergySum(Hare.Geometry.Vector V);
        public abstract Hare.Geometry.Vector[] Dir_EnergySum(double alt, double azi, bool degrees);
        public abstract void Set_Power(double[] factor);

        //public abstract double[] Pressure { get; }
        //public abstract double[] Dir_Pressure(int Rec_ID, double alt, double azi, bool degrees, bool Figure8, int sampleFreq);
        //public abstract double[][] Dir_Pressure(int Rec_ID, double alt, double azi, bool degrees, int sampleFreq);
        public abstract double[] Filter { get; }
        public abstract double[] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, bool Figure8, int sampleFreq, bool flat);
        public abstract double[][] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, int sampleFreq, bool flat);

        public abstract void Create_Filter(int length, int threadid);
        public abstract double[] Create_Filter(double[] SWL, int SampleFrequency, int length, int dim, int threadid);
        public abstract double[][] Create_Filter(double[] SWL, int SampleFrequency, int length, int threadid);
        //public abstract void Create_Pressure(double[] SWL, int SampleFrequecny, int length, int threadid);
    }

    /// <summary>
    /// Class written to store screen paths (a simplified occulsion/diffraction calculation).
    /// </summary>
    public class Screen_Path: Specular_Path
    {
        public Screen_Path(Hare.Geometry.Point[] Path, int[] Sequence_polys, int[] Sequence_planes, Scene Room, Source Src, double C_Sound, double[] H, double[] Trans_Mod, ref double Direct_Time, int threadid, int Rnd)
        :base(Path, Sequence_polys, Sequence_planes, Room, Src, C_Sound , H, Trans_Mod, ref Direct_Time, threadid, Rnd)
        {
            base.Identifier = "Screen- " + base.Identifier;
        }
    }

    /// <summary>
    /// Class written to store the specular image source path.
    /// </summary>
    public class Specular_Path: Deterministic_Reflection
    {
        private Hare.Geometry.Point[] ValidPath;
        private double[] PathEnergy;
        private int BandCount { get { return PathEnergy.Length; } }
        private double Time;
        private double Length;
        private int[] Sequence;
        public double[] P;
        public double[] F;
        public double[] prms;//Octave band rms pressure.
        public System.Numerics.Complex[] Special_Filter;//Special circumstances filter (usually for detailed materials...)

        public Specular_Path(Hare.Geometry.Point[] Path, double[]Energy, double[] p, System.Numerics.Complex[] Filter, double T, int[] Seq, double Direct_Time, int SrcID)
        {
            ValidPath = Path;
            PathEnergy = Energy;
            prms = p;
            Special_Filter = Filter;
            Time = T;
            Sequence = Seq;
            Identify(SrcID, Direct_Time);
        }

        public Specular_Path(Hare.Geometry.Point[] Path, int[] Seq_planes, int[] Seq_Polys, Scene Room, Source Src, double C_Sound, double[] H, double[] Trans_Mod, ref double Direct_Time, int thread, int Rnd)
            : this(Path, Seq_planes, Seq_Polys, Room, Src, C_Sound, H, Trans_Mod, ref Direct_Time, thread, Rnd, null)
        { }

        internal Specular_Path(Hare.Geometry.Point[] Path, int[] Seq_planes, int[] Seq_Polys, Scene Room, Source Src, double C_Sound, double[] H, double[] Trans_Mod, ref double Direct_Time, int thread, int Rnd, int[] bands)
        {
            PathEnergy = new double[Src.BandCount];
            ValidPath = Path;
            //Build an Identifier
            Sequence = Seq_planes;
            double Ptx, Pty, Ptz;

            for (int q = 1; q < ValidPath.Length; q++)
            {
                Ptx = ValidPath[q].x - ValidPath[q - 1].x;
                Pty = ValidPath[q].y - ValidPath[q - 1].y;
                Ptz = ValidPath[q].z - ValidPath[q - 1].z;
                Length += Math.Sqrt(Ptx * Ptx + Pty * Pty + Ptz * Ptz);
            }

            Time = Length / C_Sound;
            Hare.Geometry.Vector DIR = ValidPath[1] - ValidPath[0];
            DIR.Normalize();
            Random rnd = new Random(Rnd);
            
            ///Energy based formulation
            Identify(Src.Source_ID(), Direct_Time);

            for (int oct = 0; oct < BandCount; oct++)
            {
                double H_rc = AcousticalMath.Intensity_Pressure(H[oct], Room.Rho_C(0));
                // H already includes source power, transmission and atmospheric attenuation.
                PathEnergy[oct] = H_rc;
            }

            foreach (int q in Seq_Polys)
            {
                if (!(Room.AbsorptionValue[q] is Basic_Material)) continue;
                double[] AbsorptionData = Room.AbsorptionValue[q].Coefficient_A_Broad(Src.ThirdOctave);
                double[] ScatteringData = Room.ScatteringValue[q].Coefficient();
                double[] TransmissionData = Room.TransmissionSpectrum(q, Src.ThirdOctave);
                for (int t = 0; t < BandCount; t++)
                {
                    PathEnergy[t] *= (1 - AbsorptionData[t]) * (1 - ScatteringData[Src.ThirdOctave ? t / 3 : t]) * (1 - TransmissionData[t]);
                    if (PathEnergy[t] == 0) PathEnergy[t] = 1E-22;  
                }
            }

            // The historical zero-energy floor applies only to calculated bands.
            if (bands != null)
                for (int band = 0; band < BandCount; band++)
                    if (!bands.Contains(band)) PathEnergy[band] = 0;

            prms = new double[BandCount];

            for (int i = 0; i < BandCount; i++) prms[i] = AcousticalMath.Pressure_Intensity(PathEnergy[i], Room.Rho_C(Path[0]));

            Special_Filter = new System.Numerics.Complex[16384];
            for (int i = 0; i < Special_Filter.Length; i++) Special_Filter[i] = 1;

            for (int i = 0; i < Seq_Polys.Length; i++)
            {
                if (Room.AbsorptionValue[Seq_Polys[i]] is Basic_Material) continue;
                Hare.Geometry.Vector d = Path[i + 1] - Path[i + 2];
                d.Normalize();
                System.Numerics.Complex[] Ref = Room.AbsorptionValue[Seq_Polys[i]].Reflection_Spectrum(44100, 16384, Room.Normal(Seq_Polys[i]), d, thread);
                int count = Math.Min(Special_Filter.Length, Ref.Length);
                for (int j = 0; j < count; j++) Special_Filter[j] *= Ref[j];
            }
        }

        public override void Set_Power(double[] factor)
        {
            //TODO: May not be setting power mod correctly here... Do we need this?
            //for (int i = 0; i < BandCount; i++) tf_spec[i] = prms[i] * Math.Pow(10, (120 - SWL[i]) / 20);
            for (int i = 0; i < PathEnergy.Length; i++) { PathEnergy[i] *= factor[i]; prms[i] *= Math.Sqrt(factor[i]); }
        }

        public override void Create_Filter(int length, int threadid)
        {
            double[] baseFilter = Audio.Pach_SP.Filter.Transfer_Function(prms, 44100, length, threadid);
            F = ApplySourceDelay(Apply_Special_Filter(baseFilter, 44100, threadid), 44100, threadid);
        }

        public override double[] Create_Filter(double[] SWL, int SampleFrequency, int length, int dim, int threadid)
        {
            double[] tf_spec = new double[BandCount];
            for (int i = 0; i < BandCount; i++) tf_spec[i] = prms[i] * Math.Pow(10, (120 - SWL[i]) / 20);

            double[] baseFilter = Audio.Pach_SP.Filter.Transfer_Function(tf_spec, SampleFrequency, length, threadid);
            return ApplySourceDelay(Apply_Special_Filter(baseFilter, SampleFrequency, threadid), SampleFrequency, threadid);
        }

        public override double[][] Create_Filter(double[] SWL, int SampleFrequency, int length, int threadid)
        {
            double[] tf_spec = new double[BandCount];
            for (int i = 0; i < BandCount; i++) tf_spec[i] = prms[i] * Math.Pow(10, (120 - SWL[i]) / 20);

            double[] baseFilter = Audio.Pach_SP.Filter.Transfer_Function(tf_spec, SampleFrequency, length, threadid);
            return new double[1][] { ApplySourceDelay(Apply_Special_Filter(baseFilter, SampleFrequency, threadid), SampleFrequency, threadid) };
        }

        private double[] Apply_Special_Filter(double[] baseFilter, int sampleFrequency, int threadid)
        {
            if (baseFilter == null || baseFilter.Length == 0 || Special_Filter == null || Special_Filter.Length == 0) return baseFilter;

            System.Numerics.Complex[] spectrum = Audio.Pach_SP.FFT_General(baseFilter, threadid);

            if (sampleFrequency == 44100 && spectrum.Length == Special_Filter.Length)
            {
                for (int i = 0; i < spectrum.Length; i++) spectrum[i] *= Special_Filter[i];
            }
            else
            {
                int N = spectrum.Length;
                int half = N / 2;
                int sourceHalf = Special_Filter.Length / 2;
                System.Numerics.Complex[] adapted = new System.Numerics.Complex[N];

                for (int k = 0; k <= half; k++)
                {
                    double frequency = k * (double)sampleFrequency / N;
                    double sourcePosition = frequency * Special_Filter.Length / 44100.0;
                    sourcePosition = Math.Max(0, Math.Min(sourceHalf, sourcePosition));
                    int i0 = (int)Math.Floor(sourcePosition);
                    int i1 = Math.Min(sourceHalf, i0 + 1);
                    double fraction = sourcePosition - i0;
                    System.Numerics.Complex value = Special_Filter[i0] * (1.0 - fraction) + Special_Filter[i1] * fraction;
                    adapted[k] = value;
                    if (k > 0 && k < N - k) adapted[N - k] = System.Numerics.Complex.Conjugate(value);
                }

                for (int i = 0; i < N; i++) spectrum[i] *= adapted[i];
            }

            double[] output = Audio.Pach_SP.IFFT_Real_General(spectrum, threadid);
            Audio.Pach_SP.Scale(ref output);
            return output;
        }

        private void Identify(int SrcID, double Direct_Time)
        {
            if (SrcID < 10)
            {
                Identifier = string.Format("S00{0}-", SrcID);
            }
            else if (SrcID < 100)
            {
                Identifier = string.Format("S0{0}-", SrcID);
            }
            else
            {
                Identifier = string.Format("S{0}-", SrcID);
            }

            if (Sequence.Length < 10)
            {
                Identifier = string.Concat(Identifier, string.Format("Order 00{0}: ", Sequence.Length));
            }
            else if (Sequence.Length < 100)
            {
                Identifier = string.Concat(Identifier, string.Format("Order 0{0}: ", Sequence.Length));
            }
            else
            {
                Identifier = string.Concat(Identifier, string.Format("Order {0}:", Sequence.Length));
            }

            Identifier = string.Concat(Identifier, string.Format("{0} ms. ", Math.Round((Time - Direct_Time) * 1000)));

            foreach (int Digit in Sequence)
            {
                Identifier = string.Concat(Identifier, Digit.ToString(), " ");
            }
        }

        public override Hare.Geometry.Point[][] Path
        {
            get { return new Hare.Geometry.Point[][] { ValidPath }; }
        }

        public override double TravelTime
        {
            get { return Time; }
        }

        public override int[] Reflection_Sequence
        {
            get { return Sequence; }
        }

        public override double[] Energy(int Octave, int SampleFrequency)
        {
            if (Octave < BandCount) return new double[] { PathEnergy[Octave] };
            else return new double[] { PathEnergy.Sum() };
        }

        public override double[] Dir_Energy(int Octave, int dir)
        {
            Hare.Geometry.Vector Dir = Path[0][Path.Length - 1] - Path[0][Path.Length - 2];
            Dir.Normalize();
            switch (dir)
            {
                case 0:
                    return new double[] { Dir.dx * PathEnergy[Octave] };
                case 1:
                    return new double[] { Dir.dy * PathEnergy[Octave] };
                case 2:
                    return new double[] { Dir.dz * PathEnergy[Octave] };
                default:
                    throw new Exception("indexed directions must conform to 0 = x, 1 = y and 2 = z") ;
            }
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave)
        {
            Hare.Geometry.Vector Dir = Path[0][Path[0].Length - 1] - Path[0][Path[0].Length - 2];
            Dir.Normalize();
            double I = 0;
            if (Octave == BandCount) for (int oct = 0; oct < BandCount; oct++) I += PathEnergy[oct];
            else I = PathEnergy[Octave];
            return new Hare.Geometry.Vector[] { Dir * I };
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, double alt, double azi, bool degrees, bool Figure8 = false)
        {
            Hare.Geometry.Vector[] V = Dir_Energy(Octave);
            return new Hare.Geometry.Vector[] { Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V[0], azi, 0, degrees), 0, alt, degrees)};
        }

        public override double[][] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, int SampleFreq, bool flat)
        {
            Hare.Geometry.Vector V = Path[0][Path[0].Length - 1] - Path[0][Path[0].Length - 2];
            V.Normalize();
            Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V, azi, 0, degrees), 0, alt, degrees);

            double[] F_Chosen = (SampleFreq == 44100 && flat) ? F : this.Create_Filter(SWL, SampleFreq, 16384, 0)[0];
            double[][] Fn = new double[F_Chosen.Length][];

            for (int i = 0; i < F_Chosen.Length; i++)
            {
                Fn[i] = new double[3] { Vn.dx * F_Chosen[i], Vn.dy * F_Chosen[i], Vn.dz * F_Chosen[i] };
            }
            return Fn;
        }

        public override double[] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, bool Figure8, int SampleFreq, bool flat)
        {
            Hare.Geometry.Vector V = -1 * (Path[0][Path[0].Length - 1] - Path[0][Path[0].Length - 2]);
            V.Normalize();
            Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V, azi, 0, degrees), 0, alt, degrees);

            double[] F_Chosen = (SampleFreq == 44100 && flat) ? F : this.Create_Filter(SWL, SampleFreq, 16384, 0)[0];
            double[] Fn = new double[F_Chosen.Length];
            if (Figure8)
            {
                for (int i = 0; i < F_Chosen.Length; i++)
                {
                    Fn[i] = Vn.dx * F_Chosen[i];
                }
            }
            else if (Vn.dx > 0)
            {
                for (int i = 0; i < F_Chosen.Length; i++)
                {
                    Fn[i] = Vn.dx * F_Chosen[i];
                }
            }
            return Fn;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, Hare.Geometry.Vector V)
        {
            double l = Math.Sqrt(V.dz * V.dz + V.dx * V.dx);
            double azi = Math.Asin(V.dy / l);
            double alt = Math.Atan2(V.dx, V.dz);
            return Dir_Energy(Octave, alt, azi, false);
        }

        public int ReflectionOrder
        {
            get { return Sequence.Length; }
        }

        public int SurfaceIndex(int x)
        {
            return Sequence[x];
        }

        public double EnergySum()
        {
            double sum = 0;
            foreach(double e in PathEnergy) sum += e;
            return sum;
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(Hare.Geometry.Vector V)
        {
            double l = Math.Sqrt(V.dz * V.dz + V.dx * V.dx);
            double azi = Math.Asin(V.dy / l);
            double alt = Math.Atan2(V.dx, V.dz);

            Hare.Geometry.Vector E = new Hare.Geometry.Vector();
            for (int oct = 0; oct < BandCount; oct++) E += Dir_Energy(oct, alt, azi, false)[0];
            return new Hare.Geometry.Vector[] { E };
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(double alt, double azi, bool degrees)
        {
            Hare.Geometry.Vector E = new Hare.Geometry.Vector();
            for (int oct = 0; oct < BandCount; oct++) E += Dir_Energy(oct, alt, azi, degrees)[0];
            return new Hare.Geometry.Vector[] { E };
        }

        public override double[] Filter
        {
            get
            {
                return F;
            }
        }

        public override string ToString()
        {
            return Identifier;
        }
    }

    public class Compound_Path : Deterministic_Reflection
    {
        private Hare.Geometry.Point[][] ValidPath; //Geometry of the reflection.
        public double[][] PathEnergy;
        private double Time; //Earliest time of the reflection.
        private int[] Sequence; //Unique identifying indices for each reflecting element.
        private double RhoC;
        //public double[][] Octave_Power; //Contains power modifying information for each sample in H (such as absorption or transmission coefficients, air attenuation, etc.)
        private int BandCount { get { return H.Length; } }
        public double[][] H; //Initial H-function of the reflection (contains info about diffraction, compression, etc. that is a consequence of the resulting wave-form.) [oct][t]
        public double[][][] Hdir; //Directional form of H.[d][oct][t]
        public double[] F; //Filter form of the reflection time signature.
        public double[][] Fdir; //Directional form of F.[d][t]


        public Compound_Path(Hare.Geometry.Point[][] PathVertices, int Src_id, int[] Seq_Planes, double[][][] _H, double mintime, double Rho_C, ref double Direct_Time, int Threadid)
        {
            //Reconcile Source Power as SWL with power with directivity. source directivity needs to be added to OctavePower.
            RhoC = Rho_C;
            Time = mintime;
            ///Here, a compound reflection collector, and any interpolation that must be done.
            ValidPath = PathVertices;
            List<Hare.Geometry.Point[]> Paths = new List<Hare.Geometry.Point[]>();
            Sequence = Seq_Planes;
            H = _H[0];
            Hdir = new double[6][][] { _H[1], _H[2], _H[3], _H[4], _H[5], _H[6] };

            double[] timeaxis = new double[H.Length];
            for (int t = 0; t < timeaxis.Length; t++) { timeaxis[t] = (double)t / 44100f;}

            //Build an Identifier
            Identify(Src_id, Direct_Time);
        }
        
        public Compound_Path(Hare.Geometry.Point[][] PathVertices, int Src_id, int[] Seq_Planes, double[][][] _H, double mintime, double Rho_C, ref double Direct_Time, int Threadid, double[][] Absorption)
            : this(PathVertices, Src_id, Seq_Planes, _H, mintime, Rho_C, ref Direct_Time, Threadid)
        {
            if (Absorption == null || Absorption.Length != Sequence.Length || Absorption.Any(a => a == null || a.Length != BandCount))
                throw new ArgumentException("Provide one absorption spectrum per path interaction, matching the response band count.", nameof(Absorption));

            H = (double[][])H.Clone();
            for (int d = 0; d < 6; d++) Hdir[d] = (double[][])Hdir[d].Clone();
            for (int oct = 0; oct < BandCount; oct++)
            {
                double Reflection = 1;
                for (int q = 0; q < Absorption.Length; q++) Reflection *= Math.Sqrt(1 - Absorption[q][oct]);
                H[oct] = (double[])H[oct].Clone();
                for (int t = 0; t < H[oct].Length; t++) H[oct][t] *= Reflection;
                for (int d = 0; d < 6; d++)
                {
                    Hdir[d][oct] = (double[])Hdir[d][oct].Clone();
                    for (int t = 0; t < Hdir[d][oct].Length; t++) Hdir[d][oct][t] *= Reflection;
                }
            }
        }

        public override double[][] Create_Filter(double[] SWL, int SampleFrequency, int LengthofPulse, int Threadid)
        {
            double[][] Fdir_out = new double[6][];
            double[] F_out = ApplySourceDelay(Audio.Pach_SP.ETCToFilter(H, SWL, 44100, SampleFrequency), SampleFrequency, Threadid);
            for (int i = 0; i < 6; i++)
            {
                Fdir_out[i] = ApplySourceDelay(Audio.Pach_SP.ETCToFilter(Hdir[i], SWL, 44100, SampleFrequency), SampleFrequency, Threadid);
            }
            return new double[7][] { F_out, Fdir_out[0], Fdir_out[1], Fdir_out[2], Fdir_out[3], Fdir_out[4], Fdir_out[5] };
        }

        public override double[] Create_Filter(double[] SWL, int SampleFrequency, int LengthofPulse, int dim, int Threadid)
        {
            if (dim < 1)
            {
                return ApplySourceDelay(Audio.Pach_SP.ETCToFilter(H, SWL, 44100, SampleFrequency), SampleFrequency, Threadid);
            }
            else
            {
                return ApplySourceDelay(Audio.Pach_SP.ETCToFilter(Hdir[dim-1], SWL, 44100, SampleFrequency), SampleFrequency, Threadid);
            }
        }

        public override void Create_Filter(int length, int threadid)
        {
            Fdir = new double[6][];

            F = ApplySourceDelay(Audio.Pach_SP.ETCToFilter(H, Enumerable.Repeat(120.0, BandCount).ToArray()), 44100, threadid);
            for (int i = 0; i < 6; i++)
            {
                Fdir[i] = ApplySourceDelay(Audio.Pach_SP.ETCToFilter(Hdir[i], Enumerable.Repeat(120.0, BandCount).ToArray()), 44100, threadid);
            }
        }

        private void Identify(int SrcID, double Direct_Time)
        {
            if (SrcID < 10)
            {
                Identifier = string.Format("S00{0}-", SrcID);
            }
            else if (SrcID < 100)
            {
                Identifier = string.Format("S0{0}-", SrcID);
            }
            else
            {
                Identifier = string.Format("S{0}-", SrcID);
            }

            if (Sequence.Length < 10)
            {
                Identifier = string.Concat(Identifier, string.Format("Order 00{0}: ", Sequence.Length));
            }
            else if (Sequence.Length < 100)
            {
                Identifier = string.Concat(Identifier, string.Format("Order 0{0}: ", Sequence.Length));
            }
            else
            {
                Identifier = string.Concat(Identifier, string.Format("Order {0}:", Sequence.Length));
            }

            Identifier = string.Concat(Identifier, string.Format("{0} ms. ", Math.Round((Time - Direct_Time) * 1000)));

            foreach (int Digit in Sequence)
            {
                Identifier = string.Concat(Identifier, Digit.ToString(), " ");
            }
        }

        public override double[] Energy(int Octave, int Sample_Frequency)
        {
            return Energy(Octave, Sample_Frequency, 0);
        }

        internal double[] Energy(int Octave, int Sample_Frequency, int dimension)
        {
            if (Octave >= BandCount)
            {
                double[] total = null;
                for (int band = 0; band < BandCount; band++)
                {
                    double[] component = Energy(band, Sample_Frequency, dimension);
                    if (total == null) total = new double[component.Length];
                    for (int sample = 0; sample < component.Length; sample++) total[sample] += component[sample];
                }
                return total;
            }
            if (F == null || (dimension > 0 && Fdir == null)) Create_Filter(16384, 0);
            double[] filter = dimension == 0 ? F : Fdir[dimension - 1];
            double[] H_FS = BandCount == 24 ? Audio.Pach_SP.FIR_Third_Bandpass(filter, Octave, 44100, 0) : Audio.Pach_SP.FIR_Bandpass(filter, Octave, 44100, 0);

            double[] energy = new double[(int)Math.Ceiling(H_FS.Length * (double)Sample_Frequency / 44100)];

            for (int i = 0; i < H_FS.Length; i++)
            {
                int sample = (int)(i * (double)Sample_Frequency / 44100);
                energy[sample] += AcousticalMath.Intensity_Pressure(H_FS[i], RhoC);
            }
            return energy;
        }

        public override double[] Dir_Energy(int Octave, int dir)
        {
            if (dir < 0 || dir > 2) throw new ArgumentOutOfRangeException("dir");
            double[] De = new double[H[0].Length];
            int first = Octave < BandCount ? Octave : 0, last = Octave < BandCount ? Octave + 1 : BandCount;
            for (int oct = first; oct < last; oct++) for (int t = 0; t < De.Length; t++)
                De[t] += AcousticalMath.Intensity_Pressure(Hdir[2 * dir][oct][t], RhoC) - AcousticalMath.Intensity_Pressure(Hdir[2 * dir + 1][oct][t], RhoC);
            return De;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave)
        {
            Hare.Geometry.Vector[] De = new Hare.Geometry.Vector[H[0].Length];
            for (int t = 0; t < De.Length; t++) De[t] = new Hare.Geometry.Vector();
            int first = Octave < BandCount ? Octave : 0, last = Octave < BandCount ? Octave + 1 : BandCount;
            for (int oct = first; oct < last; oct++) for (int t = 0; t < De.Length; t++)
            {
                double x = AcousticalMath.Intensity_Pressure(Hdir[0][oct][t], RhoC) - AcousticalMath.Intensity_Pressure(Hdir[1][oct][t], RhoC);
                double y = AcousticalMath.Intensity_Pressure(Hdir[2][oct][t], RhoC) - AcousticalMath.Intensity_Pressure(Hdir[3][oct][t], RhoC);
                double z = AcousticalMath.Intensity_Pressure(Hdir[4][oct][t], RhoC) - AcousticalMath.Intensity_Pressure(Hdir[5][oct][t], RhoC);
                De[t] += new Hare.Geometry.Vector(x, y, z);
            }
            return De;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, double alt, double azi, bool degrees, bool Figure8 = false)
        {
            Hare.Geometry.Vector[] Fn = new Hare.Geometry.Vector[H[0].Length];

            double[][] Hdir_out = new double[6][];
            for(int d = 0; d < 6; d++)
            {
                Hdir_out[d] = new double[H[0].Length];
                if (Octave < BandCount)
                {
                    for (int t = 0; t < Hdir[d][0].Length; t++)
                    {
                        Hdir_out[d][t] = AcousticalMath.Intensity_Pressure(Hdir[d][Octave][t], RhoC);
                    } 
                }
                else
                {
                    for (int o = 0; o < BandCount; o++)
                    {
                        for (int t = 0; t < Hdir[d][0].Length; t++)
                        {
                            Hdir_out[d][t] += AcousticalMath.Intensity_Pressure(Hdir[d][o][t], RhoC);
                        }
                    }
                }
            }

            if (Figure8)
            {
                for (int i = 0; i < Hdir_out[0].Length; i++)
                {
                    Fn[i] = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[0][i] - Hdir_out[1][i], Hdir_out[2][i] - Hdir_out[3][i], Hdir_out[4][i] - Hdir_out[5][i]), azi, 0, degrees), 0, alt, degrees);
                }
            }
            else
            {
                int[] ids = new int[3];
                ids[0] = (azi > 90 && azi < 270) ? 1 : 0;
                ids[1] = (azi <= 180) ? 2 : 3;
                ids[2] = (alt > 0) ? 4 : 5;
                for (int i = 0; i < Hdir_out[0].Length; i++)
                {
                    Fn[i] = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[ids[0]][i], Hdir_out[ids[1]][i], Hdir_out[ids[2]][i]), azi, 0, true), 0, alt, true);
                }
            }
            return Fn;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, Hare.Geometry.Vector V)
        {            double l = Math.Sqrt(V.dz * V.dz + V.dx * V.dx);
            double azi = Math.Asin(V.dy / l);
            double alt = Math.Atan2(V.dx, V.dz);
            return Dir_Energy(Octave, alt, azi, false);
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(Hare.Geometry.Vector V)
        {
            double l = Math.Sqrt(V.dz * V.dz + V.dx * V.dx);
            double azi = Math.Asin(V.dy / l);
            double alt = Math.Atan2(V.dx, V.dz);

            return Dir_EnergySum(alt, azi, false);
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(double alt, double azi, bool degrees)
        {
            Hare.Geometry.Vector[] E = new Hare.Geometry.Vector[H[0].Length];
            for(int i = 0; i < E.Length; i++) E[i] = new Hare.Geometry.Vector();
            
            for (int oct = 0; oct < BandCount; oct++)
            {
                Hare.Geometry.Vector[] d = Dir_Energy(oct, alt, azi, degrees);
                for (int i = 0; i < d.Length; i++) E[i] += d[i];
            }

            return E;
        }

        public override Hare.Geometry.Point[][] Path
        {
            get { return ValidPath; }
        }

        public override int[] Reflection_Sequence
        {
            get{ return Sequence;}
        }

        public override double TravelTime
        {
            get { return Time; }
        }

        public override void Set_Power(double[] factor)
        {
            if (factor.Length != BandCount) throw new ArgumentException("Power factors must match the simulation band count.");
            for (int band = 0; band < BandCount; band++)
            {
                double scale = Math.Sqrt(factor[band]);
                for (int t = 0; t < H[band].Length; t++) H[band][t] *= scale;
                for (int d = 0; d < 6; d++) for (int t = 0; t < Hdir[d][band].Length; t++) Hdir[d][band][t] *= scale;
            }
            if (F != null) Create_Filter(16384, 0);
            }

        public override double[] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, bool Figure8, int sampleFreq, bool flat)
        {
            double[][] Hdir_out = (sampleFreq == 44100 && flat) ? Fdir : Create_Filter(SWL, sampleFreq, 16384, 0);
            double[] Fn = new double[Hdir_out[0].Length];
            if (Figure8)
            {
                for (int i = 0; i < Hdir_out[0].Length; i++)
                {
                    Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[0][i] - Hdir_out[1][i], Hdir_out[2][i] - Hdir_out[3][i], Hdir_out[4][i] - Hdir_out[5][i]), azi, 0, degrees), 0, alt, degrees);
                    Fn[i] = Vn.dx;
                }
            }
            else
            {
                int[] ids = new int[3];
                ids[0] = (azi > 90 && azi < 270) ? 1 : 0;
                ids[1] = (azi <= 180) ? 2 : 3;
                ids[2] = (alt > 0) ? 4 : 5;
                for (int i = 0; i < Hdir_out[0].Length; i++)
                {
                    Hare.Geometry.Vector V = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[ids[0]][i], Hdir_out[ids[1]][i], Hdir_out[ids[2]][i]), azi, 0, true), 0, alt, true);
                    Fn[i] = V.dx;
                }
            }
            return Fn;
        }

        public override double[][] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, int sampleFreq, bool flat)
        {
            double[][] Hdir_out = (sampleFreq == 44100 && flat) ? Fdir : Create_Filter(SWL, sampleFreq, 16384, 0);
            double[][] Fn = new double[Hdir_out[0].Length][];

            for (int i = 0; i < Hdir_out[0].Length; i++)
            {
                Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[0][i] - Hdir_out[1][i], Hdir_out[2][i] - Hdir_out[3][i], Hdir_out[4][i] - Hdir_out[5][i]), azi, 0, degrees), 0, alt, degrees);
                Fn[i] = new double[3] { Vn.dx, Vn.dy, Vn.dz };
            }
            return Fn;
        }

        public override double[] Filter
        {
            get
            {
                return F;
            }
        }
    }
}