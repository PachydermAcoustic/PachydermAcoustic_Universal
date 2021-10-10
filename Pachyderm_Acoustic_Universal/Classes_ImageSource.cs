//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2020, Arthur van der Harten 
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
        private Source Src;
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
        : this(Source, Receiver, Direct, Rm, new int[2] { 0, 7 }, MaxOrder_in, false, SourceID_in)
        { }

        public ImageSourceData(Source Source, Receiver_Bank Receiver, Direct_Sound Direct, Polygon_Scene Rm, int MaxOrder_in, bool ED, int SourceID_in)
            : this(Source, Receiver, Direct, Rm, new int[2] { 0, 7 }, MaxOrder_in, ED, SourceID_in)
        { }

        /// <summary>
        /// Constructor prepares image source calculation to run.
        /// </summary>
        /// <param name="Source"></param>
        /// <param name="Receiver"></param>
        /// <param name="Direct"></param>
        /// <param name="Rm"></param>
        /// <param name="MaxOrder_in">The maximum order to be calculated for.</param>
        public ImageSourceData(Source Source, Receiver_Bank Receiver, Direct_Sound Direct, Polygon_Scene Rm, int[] Octaves, int MaxOrder_in, bool Edge_Diffraction, int SourceID_in)
        {
            IncludeEdges = Edge_Diffraction;
            Diffraction = Edge_Diffraction;
            Oct_choice = new int[Octaves[1] - Octaves[0] + 1];
            for (int i = 0; i < Octaves.Length; i++) Oct_choice[i] = i + Octaves[0];
            SrcNo = SourceID_in;
            ValidPaths = new List<Deterministic_Reflection>[Receiver.Count];
            Speed_of_Sound = Rm.Sound_speed(0);
            MaxOrder = MaxOrder_in;
            Src = Source;
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
                if (T.ThreadState == System.Threading.ThreadState.Running) return System.Threading.ThreadState.Running;
            }
            return System.Threading.ThreadState.Stopped;
        }

        /// <summary>
        /// Inherited member. Divides the simulation into threads, and begins.
        /// </summary>
        public override void Begin()
        {
            processorCT = Pach_Properties.Instance.ProcessorCount();
            ThreadPaths = new List<Deterministic_Reflection>[Rec.Length, processorCT];
            int SrfCT = Room.ObjectCount;
            //
            if (Diffraction) SrfCT += Room.EdgeCount;
            //
            CurrentSrf = new int[processorCT];
            CurrentEdge = new int[processorCT];
            Rnd = new Random[processorCT];
            T_List = new System.Threading.Thread[processorCT];
            elementCt = Room.ObjectCount + ((this.Diffraction) ? Room.EdgeCount : 0);

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
                Sequence[0] = CurrentSrf[Params.ThreadID] + Params.StartIndex;
                TraverseOrders(0, Sequence, Params.ThreadID, Images);
            }
        }

        /// <summary>
        /// Used by specular raytracer.
        /// </summary>
        /// <param name="Sequences">List of surface index sequences to try.</param>
        public void Lookup_Sequences(List<int[]>[] Sequences)
        {
            processorCT = Pach_Properties.Instance.ProcessorCount();
            ThreadPaths = new List<Deterministic_Reflection>[Rec.Length, processorCT];
            int SrfCT = Room.ObjectCount;
            CurrentSrf = new int[processorCT];
            Random R = new Random();
            T_List = new System.Threading.Thread[processorCT];
            for (int p = 0; p < Rec.Length; p++)
            {
                for (int P_I = 0; P_I < processorCT; P_I++)
                {
                    ThreadPaths[p, P_I] = new List<Deterministic_Reflection>();
                    Calc_Params T = new Calc_Params(P_I * Sequences[p].Count / processorCT, (P_I + 1) * Sequences[p].Count / processorCT, P_I, R.Next());
                    System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate
                    {
                        for (int i = T.StartIndex; i < T.EndIndex; i++)
                        {
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

            for (int i = 0; i < processorCT; i++) Srf_CT += CurrentSrf[i];

            if (Diffraction)
            {
                int Edge_CT = 0;
                for (int i = 0; i < processorCT; i++) Edge_CT += CurrentEdge[i];

                TS = new TimeSpan((long)(TS.Ticks * (((double)(Room.ObjectCount - Srf_CT) + (Room.EdgeCount - Edge_CT)) / (Srf_CT + Edge_CT + 2))));
                return string.Format("Calculating Image/Edge Set {0} of {1}. ({2} hours,{3} minutes,{4} seconds Remaining.) Press 'Esc' to Cancel...", Srf_CT + Edge_CT, Room.ObjectCount + Room.EdgeCount, TS.Hours, TS.Minutes, TS.Seconds);
            }
            else
            {
                TS = new TimeSpan((long)(TS.Ticks * (((double)(Room.ObjectCount - Srf_CT)) / (Srf_CT + 1))));
                return string.Format("Calculating Image Set {0} of {1}. ({2} hours,{3} minutes,{4} seconds Remaining.) Press 'Esc' to Cancel...", Srf_CT, Room.ObjectCount, TS.Hours, TS.Minutes, TS.Seconds);
            }
        }

        /// <summary>
        /// Aborts all threads, effectively ending the simulation.
        /// </summary>
        public override void Abort_Calculation()
        {
            foreach (System.Threading.Thread T in T_List) T.Abort();
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
        private void TraverseOrders(int Order, int[] Sequence, int ThreadId, List<Hare.Geometry.Point[]> Images0)
        {
            Array.Resize(ref Sequence, Order + 1);
            int ct = Order == 0 ? 1 : elementCt;

            for (int q = 0; q < ct; q++)
            {
                List<Hare.Geometry.Point[]> Images = new List<Hare.Geometry.Point[]>();
                if (Order < MaxOrder)
                {
                    if (Order > 0) Sequence[Order] = q;

                    if (Images0.Count == 0) Images0.Add(new Point[0]);
                    for (int r = 0; r < Images0.Count; r++)
                    {
                        Hare.Geometry.Point[] imnext = new Point[Images0[r].Length + 1];
                        for (int i = 0; i < Images0[r].Length; i++) imnext[i] = Images0[r][i];

                        if (Sequence[Order] > Room.ObjectCount - 1)
                        {
                            //Its an edge...
                            if (Order > 0 && Sequence[Order] != Sequence[Sequence.Length - 1]) continue;
                            int edge_id = Sequence[Order] - Room.ObjectCount;
                            if (Room.Edge_Nodes[edge_id].EdgeSources.Count < 5) continue;
                            //if (Sequence[Order] != Sequence[Sequence.Length - 2]) continue;
                            for (int i = 0; i < Room.Edge_Nodes[edge_id].EdgeSources.Count; i++)
                            {
                                Hare.Geometry.Point[] im = imnext.Clone() as Hare.Geometry.Point[];
                                im[Order] = Room.Edge_Nodes[edge_id].EdgeSources[i].Z_mid;
                                Images.Add(im.Clone() as Hare.Geometry.Point[]);
                            }
                        }
                        else if (Room.IsPlanar(Sequence[Order]))
                        {
                            //It's Planar...
                            if (Order > 0 && Sequence[Order] != Sequence[Sequence.Length - 1]) continue;
                            Hare.Geometry.Point[] im = imnext.Clone() as Hare.Geometry.Point[];
                            im[Order] = Room.Image(Order > 0 ? Images0[r][Images0[r].Length - 1] : this.Src.H_Origin(), 0, Room.ObjectMembers[Sequence[Order]][0]);
                            Images.Add(im);
                        }
                        else
                        {
                            //It's curved
                            //Repeats allowed for this case...
                            for (int i = 0; i < Room.ObjectMeshEdges[Sequence[Order]].Length; i++)
                            {
                                Hare.Geometry.Point[] im = imnext.Clone() as Hare.Geometry.Point[];
                                im[Order] = Room.ObjectMeshEdges[Sequence[Order]][i].mid;
                                Images.Add(im);
                            }
                        }
                        ProcessImages(Images.ToArray(), Sequence, ThreadId);
                        if (Order + 1 < MaxOrder)
                        {
                            TraverseOrders(Order + 1, Sequence, ThreadId, Images);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// This function calculates the actual path of the specular reflection.
        /// </summary>
        /// <param name="Images">The list of images.</param>
        /// <param name="Sequence">The list of surface indices for reflection.</param>
        /// <param name="Threadid">The id of the calling thread.</param>
        private void ProcessImages(Hare.Geometry.Point[][] Images, int[] Sequence, int Threadid)
        {
            for (int rec_id = 0; rec_id < Rec.Length; rec_id++)
            {
                double c_sound = Room.Sound_speed(Rec[rec_id]);

                double[][] Trans_Mod = new double[Images.Length][];
                int[] Seq_Polys = new int[Sequence.Length];
                List<Hare.Geometry.Point[]> PathVertices = new List<Hare.Geometry.Point[]>();
                Hare.Geometry.Point S = Src.H_Origin();
                Hare.Geometry.Point E = Rec[rec_id];
                double df = SampleRate * .5 / 4096;

                //Find all Path Legs from Receiver to Source
                for (int r = 0; r < Images.Length; r++)
                {
                    Trans_Mod[r] = new double[8];
                    for (int t_oct = 0; t_oct < 8; t_oct++) Trans_Mod[r][t_oct] = 1;

                    Hare.Geometry.Point[] path = new Hare.Geometry.Point[Sequence.Length + 2];
                    path[0] = S;
                    path[path.Length - 1] = E;

                    for (int q = Sequence.Length - 1; q >= 0; q--)
                    {
                        if (Sequence[q] > Room.ObjectCount - 1)
                        {
                            //It's an edge
                            if (!OcclusionIntersectED(path[q + 2], Images[r][q], Sequence[q] - Room.ObjectCount, ref Trans_Mod[r], ref path[q + 1], Threadid))
                            { path = null; break; }
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
                                if (!this.Check_Edge_For_Specular(path[q + 2], path[q] == null ? Images[r][q] : path[q], new Hare.Geometry.Vector(Room.ObjectMeshEdges[Sequence[q]][r].Polys[0].Normal * -1), new Hare.Geometry.Vector(Room.ObjectMeshEdges[Sequence[q]][r].Polys[1].Normal * -1), Room.ObjectMeshEdges[Sequence[q]][r].a, Room.ObjectMeshEdges[Sequence[q]][r].b, ref path[q + 1])) { path = null; break; }
                                if (!OcclusionIntersect(path[q + 2], path[q + 1], Sequence[q], ref Trans_Mod[r], ref path[q + 1], ref Seq_Polys[q], Threadid)) { path = null; break; }
                            }
                        }
                    }
                    PathVertices.Add(path);
                }

                //Check that any path was unoccluded... if so, then record this entry. If not, move on...
                if (PathVertices.Count(item => item != null) == 0) continue; //goto Next;

                //Final Occlusion Check:
                for (int r = 0; r < PathVertices.Count; r++)
                {
                    if (PathVertices[r] == null) continue;
                    if (Sequence[0] < Room.ObjectCount)
                    {
                        if (FinalOcclusion(PathVertices[r][0], PathVertices[r][1], Sequence[0], ref Trans_Mod[r], Threadid))
                            PathVertices[r] = null;
                    }
                    else
                    {
                        int edge_id = Sequence[0] - Room.ObjectCount;
                        if (FinalOcclusion(PathVertices[r][0], PathVertices[r][1], 0.00001, ref Trans_Mod[r], Threadid))
                            PathVertices[r] = null;
                    }
                }

                //Check again for null(occluded) paths...
                if (PathVertices.Count(item => item != null) == 0) continue; //goto Next;

                ///Process all paths for pulse entry...
                if (PathVertices.Count == 0) continue;

                for (int i = 0; i < PathVertices.Count; i++)
                {
                    if (PathVertices[i] == null) continue;
                    Hare.Geometry.Vector d = PathVertices[i][1] - PathVertices[i][0];
                    d.Normalize();
                    double[] En = Src.DirPower(Threadid, 0, d);
                    for (int oct = 0; oct < 8; oct++) Trans_Mod[i][oct] *= En[oct];
                }

                //////////////////////////////
                //Process Compound Path before storing it.
                double mintime;
                double[][] H_all = H_Function_NoInterp(PathVertices, Sequence, c_sound, Threadid, out mintime);

                if (H_all == null) return;

                //TODO: For portions of the reflection that have differential power characteristics (such as transmission through a material) collect bundles of paths and find the H_Function for each collection of path conditons individually. Archive as separate reflections.
                List<double[]> tm = Trans_Mod.ToList();
                tm.RemoveAll(item => item == null);
                Trans_Mod = tm.ToArray() as double[][];

                ///Enter the reflection
                PathVertices.RemoveAll(item => item == null);

                if (PathVertices.Count > 1) ThreadPaths[rec_id, Threadid].Add(new Compound_Path(Room, PathVertices.ToArray(), Sequence, Src, Trans_Mod, H_all[0], new double[][] { H_all[1], H_all[2], H_all[3], H_all[4], H_all[5], H_all[6] }, mintime, ref Direct_Time[rec_id], Threadid));
                else
                {
                    //double[] En = new double[8];
                    //for (int oct = 0; oct < 8; oct++) En[oct] *= H_all[0][0] * H_all[0][0];
                    ThreadPaths[rec_id, Threadid].Add(new Specular_Path(PathVertices[0], Sequence, Seq_Polys, Room, Src, c_sound, H_all[0], Trans_Mod[0], ref Direct_Time[rec_id], Threadid, Rnd[Threadid].Next()));
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
            //Determine closest ppint to a midpoint between reflection points before and after...
            Hare.Geometry.Point Cpt = (S + R) / 2;
            Hare.Geometry.Vector ZDir = Z2 - Z1;
            Hare.Geometry.Vector tDir = Cpt - Z1;
            double t = Hare_math.Dot(tDir, ZDir);
            double denom;
            if (t < 0)
                return false;
            else
            {
                denom = Hare_math.Dot(ZDir, ZDir);
                if (t >= denom)
                    return false;
            }

            Hare.Geometry.Point Z_pt = Z1 + ZDir * t / denom;
            //Hare.Geometry.Vector TheoreticalNormal = Cpt - Z_pt;
            Hare.Geometry.Vector v1 = Z_pt - S, v2 = Z_pt - R;
            v1.Normalize(); v2.Normalize();
            Hare.Geometry.Vector TheoreticalNormal = v1 + v2;
            TheoreticalNormal.Normalize();

            Hare.Geometry.Vector bsect = T1_Norm + T2_Norm;
            bsect.Normalize();

            double dot_Norm1 = Math.Abs(Hare_math.Dot(T1_Norm, bsect));
            double dot_Norm2 = Math.Abs(Hare_math.Dot(T2_Norm, bsect));
            double testdot = Math.Abs(Hare_math.Dot(bsect, TheoreticalNormal));
            //If it is between the two normals adjacent to these edges, then reflection is valid. Return true.
            if (testdot >= dot_Norm1 && testdot >= dot_Norm2)
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

                if (order == Sequence.Length-1)
                {
                    Hare.Geometry.Vector d1 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0] - Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[0]; // Direction Hare.Geometry.Vector of segment S1
                    Hare.Geometry.Vector d2 = Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1] - Room.Edge_Nodes[Sequence[order] - Room.ObjectCount].EdgeSources[pathid].Z_limits[1]; // Direction Hare.Geometry.Vector of segment S2
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
                W_Kurvature = project_wavefront(l, W_Kurvature);

                Hare.Geometry.Vector dir2 = (PathVertices[pathid][order + 1] - PathVertices[pathid][order ]);
                Hare.Geometry.Vector local_N = dir1 + dir2;
                local_N.Normalize();
                W_Frame[0] -= local_N * Hare_math.Dot(W_Frame[0], local_N) * 2;
                W_Frame[1] -= local_N * Hare_math.Dot(W_Frame[1], local_N) * 2;

                if (order < Sequence.Length -1) return h;
            }
            else
            {
                hasedge = false;
                Hare.Geometry.Point a;
                Hare.Geometry.Point b;
                double K0 = 0 + Room.Edge_Kurvatures[Sequence[order]][pathid % Room.ObjectMeshEdges.Length][0], K1 = 0 + Room.Edge_Kurvatures[Sequence[order]][pathid % Room.ObjectMeshEdges.Length][1];
                Hare.Geometry.Vector Normal = new Hare.Geometry.Vector(Room.Edge_Normals[Sequence[order]][pathid].x, Room.Edge_Normals[Sequence[order]][pathid].y, Room.Edge_Normals[Sequence[order]][pathid].z);
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
                            Hare.Geometry.Point p = Room.ObjectMeshEdges[Sequence[order]][pid].TributaryLength[i] * e.Tangents[i];
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
                        Hare.Geometry.Point p = Room.ObjectMeshEdges[Sequence[order]][pid].TributaryLength[i] * e.Tangents[i];
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

            //TODO: is there a W_Kurvature at the end of an edge leg? If so, get rid of the if statement, and always do this.
            if (!hasedge) h *= Math.Sqrt((1 / (4 * Math.PI)) * Math.Abs(W_Kurvature.Determinant()));

            return h;
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

        public double[][] H_Function_NoInterp(List<Point[]> PathVertices, int[] Sequence, double c_sound, int Threadid, out double mintime)
        {
            ///////////////////////////////////
            List<double> Time = new List<double>();
            List<double> Bs = new List<double>();
            List<DenseMatrix> Wavefront_Kurvatures = new List<DenseMatrix>();
            List<double> xdir = new List<double>();
            List<double> ydir = new List<double>();
            List<double> zdir = new List<double>();
            List<int> E_ids = new List<int>();
            //Q: Was there a way to tell when it makes large jumps?
            //A: Detect Inflections and split...

            //Calculate each leg of the compound path (if a composite of multiple compound paths...
            List<double> Tmin = new List<double>(), Tmax = new List<double>();
            for (int i = 0; i < PathVertices.Count; i++)
            {
                if (PathVertices[i] == null) continue;
                E_ids.Add(i);
                double[] t_limits = new double[2];
                DenseMatrix W_orig = new DenseMatrix(2); W_orig[0, 0] = double.Epsilon; W_orig[1, 1] = double.Epsilon;

                Hare.Geometry.Vector[] W_Frame = new Hare.Geometry.Vector[2];
                double p = Power_Recursion(i, Sequence.Length - 1, PathVertices, Sequence, c_sound, t_limits, Threadid, ref W_orig, ref W_Frame);

                //double duration_s = Math.Abs(t_limits[0] - t_limits[1]) * 44100.0 / 88200;
                //p *= Math.Max(duration_s , 1.0/88200.0);

                Hare.Geometry.Vector DIR;
                DIR = PathVertices[i][Sequence.Length - 1] - PathVertices[i][Sequence.Length];
                DIR.Normalize();
                Bs.Add(p);
                Wavefront_Kurvatures.Add(W_orig);
                Time.Add(0.5 * (t_limits[0] + t_limits[1]));
                Tmin.Add(Math.Min(t_limits[0], t_limits[1]));
                Tmax.Add(Math.Max(t_limits[0], t_limits[1]));
                xdir.Add(DIR.x);
                ydir.Add(DIR.y);
                zdir.Add(DIR.z);
            }

            mintime = double.PositiveInfinity;
            double maxtime = double.NegativeInfinity;

            double[][] H_directional = new double[7][];

            if (Bs.Count > 1)
            {
                mintime = Tmin.Min();
                maxtime = Tmax.Max();
                int samplect = (int)Math.Ceiling((maxtime - mintime) * 44100.0);
                for (int i = 0; i < 7; i++) H_directional[i] = new double[samplect];

                if (Sequence[0] > Room.ObjectCount - 1)
                {
                    for (int i = 0; i < Time.Count; i++)
                    {
                        //Split energy over straddled bins.
                        double delta = Tmax[i] - Tmin[i];
                        double t0 = Tmin[i] - mintime;
                        double t1 = t0 + delta;
                        int s1 = (int)Math.Ceiling(44100 * t0), s0 = s1-1;
                        double frac0 = ((double)s1/44100.0 - t0)/delta, frac1 = (t1 - (double)s1/44100.0)/delta;

                        double p1 = frac1 * Bs[s1], p0 = frac0 * Bs[s0];
                        H_directional[0][s1] += p1;
                        H_directional[xdir[s1]>0? 1 : 2][s1] += xdir[s1] * p1;
                        H_directional[ydir[s1] > 0 ? 3 : 4][s1] += ydir[s1] * p1;
                        H_directional[zdir[s1] > 0 ? 5 : 6][s1] += zdir[s1] * p1;
                        H_directional[0][s0] += p0;
                        H_directional[xdir[s0] > 0 ? 1 : 2][s0] += xdir[s0] * p0;
                        H_directional[ydir[s0] > 0 ? 3 : 4][s0] += ydir[s0] * p0;
                        H_directional[ydir[s0] > 0 ? 5 : 6][s0] += zdir[s0] * p0;
                    }
                }
                else
                {
                    //for curved surfaces
                    int totalsamplenbr = (int)(Math.Ceiling((Tmax.Max() - Tmin.Min()) * 88200));

                    for (int j = 0; j < 7; j++) H_directional[j] = new double[totalsamplenbr];

                    for (int i = 0; i < Time.Count; i++)
                    {
                        int samplenbr = (int)Math.Ceiling(Math.Abs(Tmax[i] - Tmin[i]) * 88200);
                        //if (Wavefront_Kurvatures[i][0, 0] < 0 || Wavefront_Kurvatures[i][1, 1] < 0)
                        int order = PathVertices[E_ids[i]].Length - 3;
                        Hare.Geometry.Vector d0 = PathVertices[E_ids[i]][order + 1] - PathVertices[E_ids[i]][order];
                        d0.Normalize();
                        double dir = Hare.Geometry.Hare_math.Dot(d0, Room.Edge_Normals[Sequence[order]][E_ids[i]]);
                        bool k1;
                        bool k2;

                        if (dir > 0) { k1 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][0] * -1 < 0; k2 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][1] * -1 < 0; }
                        else { k1 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][0] * -1 < 0; k2 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][1] * -1 < 0; }

                        //bool k1 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][0] * dir < 0;
                        //bool k2 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][1] * dir < 0;
                        if (k1 || k2)
                        {
                            double dt = 1.0 / 88200;// (Tmax[i] - Tmin[i]) / samplenbr;
                            //double[] BFct = Utilities.PachTools.NormalDistribution(samplenbr, Bs[i]);
                            double[] BFct = (k1 && k2) ? EdgeDistribution_2d(samplenbr, Room.ObjectMeshEdges[Sequence[0]][E_ids[i]], Bs[i], PathVertices[E_ids[i]][0], PathVertices[E_ids[i]][1], PathVertices[E_ids[i]][2]) : EdgeDistribution_1d(samplenbr, Room.ObjectMeshEdges[Sequence[0]][E_ids[i]], Bs[i], Tmin[i], Tmax[i], PathVertices[E_ids[i]][0], PathVertices[E_ids[i]][1], PathVertices[E_ids[i]][2]);

                            double x_2 = Math.Sqrt(Math.Abs(xdir[i])); double y_2 = Math.Sqrt(Math.Abs(ydir[i])); double z_2 = Math.Sqrt(Math.Abs(zdir[i]));
                            double[] x_out = new double[BFct.Length]; double[] y_out = new double[BFct.Length]; double[] z_out = new double[BFct.Length];
                            double[] Times = new double[BFct.Length];

                            for (int j = 0; j < BFct.Length; j++)
                            {
                                double t = Tmin[i] + dt * j - mintime;
                                double pressure = BFct[j];
                                H_directional[0][(int)(t / dt)] += pressure;
                                if (xdir[i] > 0) H_directional[1][(int)(t / dt)] += pressure * x_2;
                                else H_directional[2][(int)(t / dt)] += pressure * x_2;
                                if (ydir[i] > 0) H_directional[3][(int)(t / dt)] += pressure * y_2;
                                else H_directional[4][(int)(t / dt)] += pressure * y_2;
                                if (zdir[i] > 0) H_directional[5][(int)(t / dt)] += pressure * z_2;
                                else H_directional[6][(int)(t / dt)] += pressure * z_2;
                            }
                        }
                        else
                        {
                            double t = Time[i];
                            double pressure = Bs[i];
                            int s = (int)Math.Floor(88200 * (t - mintime));
                            H_directional[0][s] += pressure; //Omni Channel...
                            double x_2 = Math.Sqrt(Math.Abs(xdir[i])); double y_2 = Math.Sqrt(Math.Abs(ydir[0])); double z_2 = Math.Sqrt(Math.Abs(zdir[0]));
                            if (xdir[i] > 0) H_directional[1][s] += pressure * x_2;
                            else H_directional[2][s] += pressure * x_2;
                            if (ydir[i] > 0) H_directional[3][s] += y_2;
                            else H_directional[4][s] += y_2;
                            if (zdir[i] > 0) H_directional[5][s] += z_2;
                            else H_directional[6][s] += z_2;
                        }
                    }
                }
            }
            else
            {
                //Assumed either a convex, or planar reflection
                for (int i = 0; i < 7; i++) H_directional[i] = new double[1];
                double t = Time[0];
                double pressure = Bs[0];
                H_directional[0][0] += pressure; //Omni Channel...
                double x_2 = Math.Sqrt(Math.Abs(xdir[0])); double y_2 = Math.Sqrt(Math.Abs(ydir[0])); double z_2 = Math.Sqrt(Math.Abs(zdir[0]));
                if (xdir[0] > 0) H_directional[1][0] = pressure * x_2;
                else H_directional[2][0] = pressure * x_2;
                if (ydir[0] > 0) H_directional[3][0] = pressure * y_2;
                else H_directional[4][0] = pressure * y_2;
                if (zdir[0] > 0) H_directional[5][0] = pressure * z_2;
                else H_directional[6][0] = pressure * z_2;
            }

            return H_directional;
        }

        public double[][] H_Function(List<Point[]> PathVertices, int[] Sequence, double c_sound, int Threadid, out double mintime)
        {
            double p = 0;
            MathNet.Numerics.LinearAlgebra.Double.DenseMatrix Km = new MathNet.Numerics.LinearAlgebra.Double.DenseMatrix(2);

            ///////////////////////////////////
            List<double> Time = new List<double>();
            List<double> Bs = new List<double>();
            List<DenseMatrix> Wavefront_Kurvatures = new List<DenseMatrix>();
            List<double> xdir = new List<double>();
            List<double> ydir = new List<double>();
            List<double> zdir = new List<double>();
            List<int> E_ids = new List<int>();
            //Q: Was there a way to tell when it makes large jumps?
            //A: Detect Inflections and split...

            //Calculate each leg of the compound path (if a composite of multiple compound paths...
            List<double> Tmin = new List<double>(), Tmax = new List<double>();
            for (int i = 0; i < PathVertices.Count; i++)
            {
                if (PathVertices[i] == null) continue;
                E_ids.Add(i);
                double[] t_limits = new double[2];
                DenseMatrix W_orig = new DenseMatrix(2); W_orig[0, 0] = double.Epsilon; W_orig[1, 1] = double.Epsilon;

                Hare.Geometry.Vector[] W_Frame = new Hare.Geometry.Vector[2];
                p = Power_Recursion(i, Sequence.Length - 1, PathVertices, Sequence, c_sound, t_limits, Threadid, ref W_orig, ref W_Frame);

                //double duration_s = Math.Abs(t_limits[0] - t_limits[1]) * 44100.0 / 88200;
                //p *= Math.Max(duration_s , 1.0/88200.0);

                Hare.Geometry.Vector DIR;
                DIR = PathVertices[i][Sequence.Length - 1] - PathVertices[i][Sequence.Length];
                DIR.Normalize();
                Bs.Add(p);
                Wavefront_Kurvatures.Add(W_orig);
                Time.Add(0.5 * (t_limits[0] + t_limits[1]));
                Tmin.Add(Math.Min(t_limits[0], t_limits[1]));
                Tmax.Add(Math.Max(t_limits[0], t_limits[1]));
                xdir.Add(DIR.x);
                ydir.Add(DIR.y);
                zdir.Add(DIR.z);
            }

            mintime = double.PositiveInfinity;
            double maxtime = double.NegativeInfinity;

            double[][] H_directional = new double[7][];

            if (Bs.Count > 1)
            {
                mintime = Tmin.Min();
                maxtime = Tmax.Max();
                int samplect = (int)Math.Ceiling((maxtime - mintime) * 88200.0);
                for (int i = 0; i < 7; i++) H_directional[i] = new double[samplect];

                if (Sequence[0] > Room.ObjectCount - 1)
                {
                    List<MathNet.Numerics.Interpolation.CubicSpline> Pr_Spline = new List<MathNet.Numerics.Interpolation.CubicSpline>();
                    List<MathNet.Numerics.Interpolation.CubicSpline> x_Spline = new List<MathNet.Numerics.Interpolation.CubicSpline>();
                    List<MathNet.Numerics.Interpolation.CubicSpline> y_Spline = new List<MathNet.Numerics.Interpolation.CubicSpline>();
                    List<MathNet.Numerics.Interpolation.CubicSpline> z_Spline = new List<MathNet.Numerics.Interpolation.CubicSpline>();
                    List<List<double>> AllTimes = new List<List<double>>();

                    List<List<double>> Pr = new List<List<double>>();
                    List<List<double>> xout = new List<List<double>>();
                    List<List<double>> yout = new List<List<double>>();
                    List<List<double>> zout = new List<List<double>>();
                    List<List<double>> Times = new List<List<double>>();
                    if (Bs.Count > 2) Split_at_Inflection(Bs, Time, xdir, ydir, zdir, out Pr, out Times, out xout, out yout, out zout);
                    double[][] Spline_Times = new double[Pr.Count][];

                    for (int i = 0; i < Times.Count; i++)
                    {
                        for (int j = 0; j < Pr[i].Count; j++)
                        {
                            Pr[i][j] = Math.Sign(Pr[i][j]) * Math.Log10(Pr[i][j] * Pr[i][j]);
                            xout[i][j] = Math.Sign(Pr[i][j]) * Math.Log10(xout[i][j] * xout[i][j]);
                            yout[i][j] = Math.Sign(Pr[i][j]) * Math.Log10(yout[i][j] * yout[i][j]);
                            zout[i][j] = Math.Sign(Pr[i][j]) * Math.Log10(zout[i][j] * zout[i][j]);
                        }
                        Spline_Times[i] = new double[] { Times[i].Min(), Times[i].Max() };
                        if (Pr[i].Count < 5) continue;
                        Pr_Spline.Add(MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(Times[i], Pr[i]));
                        x_Spline.Add(MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(Times[i], xout[i]));
                        y_Spline.Add(MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(Times[i], yout[i]));
                        z_Spline.Add(MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(Times[i], zout[i]));
                    }
                    AllTimes.AddRange(Times);

                    if (Pr_Spline.Count == 0) return null;

                    for (int i = 0; i < Pr_Spline.Count; i++)
                    {
                        int start = (int)Math.Floor(88200 * (Spline_Times[i][0] - mintime));
                        int end = (int)Math.Floor(88200 * (Spline_Times[i][1] - mintime));
                        for (int j = start; j < end; j++)
                        {
                            double t = j / 88200.0 + mintime;
                            double pit = Pr_Spline[i].Interpolate(t);
                            double pressure = Math.Sign(pit) * Math.Sqrt(Math.Pow(10, Math.Abs(pit)));
                            H_directional[0][j] += pressure; //Omni Channel...
                            double xit = x_Spline[i].Interpolate(t);
                            double x = Math.Sign(xit) * Math.Sqrt(Math.Pow(10, Math.Abs(xit))); //TODO: Coordinate approach to directionality with other aspects of software.
                            if (x > 0) H_directional[1][j] = x;
                            else H_directional[2][j] = x;
                            double yit = y_Spline[i].Interpolate(t);
                            double y = Math.Sign(yit) * Math.Sqrt(Math.Pow(10, Math.Abs(yit)));
                            if (y > 0) H_directional[3][j] = y;
                            else H_directional[4][j] = y;
                            double zit = z_Spline[i].Interpolate(t);
                            double z = Math.Sign(zit) * Math.Sqrt(Math.Pow(10, Math.Abs(zit)));
                            if (z > 0) H_directional[5][j] = z;
                            else H_directional[6][j] = z;
                        }
                    }
                }
                else
                {
                    //for curved surfaces
                    int totalsamplenbr = (int)(Math.Ceiling((Tmax.Max() - Tmin.Min()) * 88200));

                    for (int j = 0; j < 7; j++) H_directional[j] = new double[totalsamplenbr];

                    for (int i = 0; i < Time.Count; i++)
                    {
                        int samplenbr = (int)Math.Ceiling(Math.Abs(Tmax[i] - Tmin[i]) * 88200);
                        //if (Wavefront_Kurvatures[i][0, 0] < 0 || Wavefront_Kurvatures[i][1, 1] < 0)
                        int order = PathVertices[E_ids[i]].Length - 3;
                        Hare.Geometry.Vector d0 = PathVertices[E_ids[i]][order + 1] - PathVertices[E_ids[i]][order];
                        d0.Normalize();
                        double dir = Hare.Geometry.Hare_math.Dot(d0, Room.Edge_Normals[Sequence[order]][E_ids[i]]);
                        bool k1;
                        bool k2;

                        if (dir > 0) { k1 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][0] * -1 < 0; k2 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][1] * -1 < 0; }
                        else { k1 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][0] * -1 < 0; k2 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][1] * -1 < 0; }

                        //bool k1 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][0] * dir < 0;
                        //bool k2 = Room.Edge_Kurvatures[Sequence[0]][E_ids[i]][1] * dir < 0;
                        if (k1 || k2)
                        {
                            double dt = 1.0 / 88200;// (Tmax[i] - Tmin[i]) / samplenbr;
                            //double[] BFct = Utilities.PachTools.NormalDistribution(samplenbr, Bs[i]);
                            double[] BFct = (k1 && k2) ? EdgeDistribution_2d(samplenbr, Room.ObjectMeshEdges[Sequence[0]][E_ids[i]], Bs[i], PathVertices[E_ids[i]][0], PathVertices[E_ids[i]][1], PathVertices[E_ids[i]][2]) : EdgeDistribution_1d(samplenbr, Room.ObjectMeshEdges[Sequence[0]][E_ids[i]], Bs[i], Tmin[i], Tmax[i], PathVertices[E_ids[i]][0], PathVertices[E_ids[i]][1], PathVertices[E_ids[i]][2]);

                            double x_2 = Math.Sqrt(Math.Abs(xdir[i])); double y_2 = Math.Sqrt(Math.Abs(ydir[i])); double z_2 = Math.Sqrt(Math.Abs(zdir[i]));
                            double[] x_out = new double[BFct.Length]; double[] y_out = new double[BFct.Length]; double[] z_out = new double[BFct.Length];
                            double[] Times = new double[BFct.Length];

                            for (int j = 0; j < BFct.Length; j++)
                            {
                                double t = Tmin[i] + dt * j - mintime;
                                double pressure = BFct[j];
                                H_directional[0][(int)(t / dt)] += pressure;
                                if (xdir[i] > 0) H_directional[1][(int)(t / dt)] += pressure * x_2;
                                else H_directional[2][(int)(t / dt)] += pressure * x_2;
                                if (ydir[i] > 0) H_directional[3][(int)(t / dt)] += pressure * y_2;
                                else H_directional[4][(int)(t / dt)] += pressure * y_2;
                                if (zdir[i] > 0) H_directional[5][(int)(t / dt)] += pressure * z_2;
                                else H_directional[6][(int)(t / dt)] += pressure * z_2;
                            }
                        }
                        else
                        {
                            double t = Time[i];
                            double pressure = Bs[i];
                            int s = (int)Math.Floor(88200 * (t - mintime));
                            H_directional[0][s] += pressure; //Omni Channel...
                            double x_2 = Math.Sqrt(Math.Abs(xdir[i])); double y_2 = Math.Sqrt(Math.Abs(ydir[0])); double z_2 = Math.Sqrt(Math.Abs(zdir[0]));
                            if (xdir[i] > 0) H_directional[1][s] += pressure * x_2;
                            else H_directional[2][s] += pressure * x_2;
                            if (ydir[i] > 0) H_directional[3][s] += y_2;
                            else H_directional[4][s] += y_2;
                            if (zdir[i] > 0) H_directional[5][s] += z_2;
                            else H_directional[6][s] += z_2;
                        }
                    }
                }
            }
            else
            {
                //Assumed either a convex, or planar reflection
                for (int i = 0; i < 7; i++) H_directional[i] = new double[1];
                double t = Time[0];
                double pressure = Bs[0];
                H_directional[0][0] += pressure; //Omni Channel...
                double x_2 = Math.Sqrt(Math.Abs(xdir[0])); double y_2 = Math.Sqrt(Math.Abs(ydir[0])); double z_2 = Math.Sqrt(Math.Abs(zdir[0]));
                if (xdir[0] > 0) H_directional[1][0] = pressure * x_2;
                else H_directional[2][0] = pressure * x_2;
                if (ydir[0] > 0) H_directional[3][0] = pressure * y_2;
                else H_directional[4][0] = pressure * y_2;
                if (zdir[0] > 0) H_directional[5][0] = pressure * z_2;
                else H_directional[6][0] = pressure * z_2;
            }

            return H_directional;
        }

        public double[] EdgeDistribution_1d(int samplect, Hare.Geometry.Edge e, double power, double Tmin, double Tmax, Hare.Geometry.Point p1, Hare.Geometry.Point p1_5, Hare.Geometry.Point p2)
        {
            int samples = Math.Max((int)Math.Floor((Tmax - Tmin) * 88200), 1);
            double[] result = new double[samples];

            for (int j = 0; j < samples; j++)
            {
                result[j] += power / samples;
            }

            return result;
        }

        public double[] EdgeDistribution_2d(int samplect, Hare.Geometry.Edge e, double power, Hare.Geometry.Point p1, Hare.Geometry.Point p1_5, Hare.Geometry.Point p2)
        {
            double start = ((p1 - p1_5).Length() + (p1_5 - p2).Length()) / Speed_of_Sound;
            double e1 = ((p1 - e.a).Length() + (p2 - e.a).Length()) / Speed_of_Sound;
            double e2 = ((p1 - e.b).Length() + (p2 - e.b).Length()) / Speed_of_Sound;
            double min = Math.Min(e1, e2);
            double maxt = Math.Max(e1, e2);
            double[] intermediates = new double[e.Polys.Count];
            int[] sample = new int[e.Polys.Count];
            int[] upperct = new int[e.Polys.Count];
            double[] result = new double[samplect];
            double max = 0;
            bool invert = false;
            for (int i = 0; i < e.Polys.Count; i++)
            {
                intermediates[i] = ((p1 - e.Polys[i].Centroid).Length() + (p2 - e.Polys[i].Centroid).Length()) / Speed_of_Sound;
                sample[i] = (int)((intermediates[i] - min) / Math.Abs(e1 - e2)) * samplect;
                upperct[i] = (samplect - sample[i]);
                invert |= sample[i] < min || sample[i] > maxt;
            }

            //if (invert)
            //{
            //    //The adjacent centroids are located outside of the time window. Rotate the traversal in order to ensure that energy is recorded at the correct times.
            //    for (int i = 0; i < e.Polys.Count; i++)
            //    {
            //        Hare.Geometry.Point c = e.Polys[i].Centroid;
            //        double triblength = Hare.Geometry.Hare_math.Cross((c - e.a), (c - e.b)).Length() / (e.a - e.b).Length(); ;
            //        for (int j = 0; j < samplect; j++)
            //        {


            //        }
            //    }
            //}
            //else
            //{
            for (int i = 0; i < e.Polys.Count; i++)
            {
                for (int j = 0; j < samplect; j++)
                {
                    if (j < sample[i])
                    {
                        double s = e.TributaryLength[i] * j / sample[i];
                        result[j] += s;
                        max += s;
                    }
                    else
                    {
                        double s = e.TributaryLength[i] * (samplect - j) / upperct[i];
                        result[j] += s;
                        max += s;
                    }
                }
            }
            //}
            for (int i = 0; i < result.Length; i++) result[i] *= power * power / max;
            for (int i = 0; i < result.Length; i++) result[i] = Math.Sqrt(result[i]);
            //for (int i = 0; i < result.Length; i++) result[i] /= max;
            return result;
        }

        private void Split_at_Inflection(List<double> BS, List<double> Time, List<double> xdir, List<double> ydir, List<double> zdir, out List<List<double>> BS_out, out List<List<double>> T_out, out List<List<double>> xout, out List<List<double>> yout, out List<List<double>> zout)
        {
            BS_out = new List<List<double>>();
            T_out = new List<List<double>>();
            xout = new List<List<double>>();
            yout = new List<List<double>>();
            zout = new List<List<double>>();

            double t_last = Time[0];
            int start = 0;
            List<double> B = new List<double>();
            for (int i = 1; i < Time.Count - 1; i++)
            {
                double dt1 = (Time[i] - t_last);
                double dt2 = (Time[i + 1] - Time[i]);
                B.Add(BS[i]);
                if (Math.Sign(dt1) != Math.Sign(dt2) || i == Time.Count-2)  //Double check that this is the right way to find extreme inflections.
                {
                    BS_out.Add(B);
                    T_out.Add(Time.GetRange(start, i - start));
                    List<double> x = xdir.GetRange(start, i - start);
                    List<double> y = ydir.GetRange(start, i - start);
                    List<double> z = zdir.GetRange(start, i - start);
                    for (int j = 0; j < B.Count; j++)
                    {
                        x[j] = Math.Sqrt(x[j]) * B[j];
                        y[j] = Math.Sqrt(y[j]) * B[j];
                        z[j] = Math.Sqrt(z[j]) * B[j];
                    }
                    xout.Add(x);
                    yout.Add(y);
                    zout.Add(z);
                    start = i;
                    B = new List<double>();
                }
                t_last = Time[i];
            }
            //BS_out.Add(BS.GetRange(start, BS.Count - 1));
            //T_out.Add(Time.GetRange(start, Time.Count - 1));
            //xout.Add(xdir.GetRange(start, xdir.Count - 1));
            //yout.Add(ydir.GetRange(start, ydir.Count - 1));
            //zout.Add(zdir.GetRange(start, zdir.Count - 1));
        }

        /// <summary>
        /// This function calculates the actual path of the specular reflection.
        /// </summary>
        /// <param name="Images">The list of images.</param>
        /// <param name="Sequence">The list of surface indices for reflection.</param>
        /// <param name="Threadid">The id of the calling thread.</param>
        private void ProcessImages(Hare.Geometry.Point[] Images, int[] Sequence, int rec_id, int Threadid)
        {
            double[] Trans_Mod = new double[8];
            int[] Seq_Polys = new int[Sequence.Length];
            for (int t_oct = 0; t_oct < 8; t_oct++) Trans_Mod[t_oct] = 1;
            Hare.Geometry.Point[] PathVertices = new Hare.Geometry.Point[Sequence.Length + 2];
            PathVertices[0] = Src.H_Origin();
            PathVertices[PathVertices.Length - 1] = Rec[rec_id];

            //Find all Path Legs from Receiver to Source
            for (int q = Sequence.Length; q > 0; q--) if (!OcclusionIntersect(PathVertices[q + 1], Images[q - 1], Sequence[q - 1], ref Trans_Mod, ref PathVertices[q], ref Seq_Polys[q - 1], Threadid)) return;

            //Final Occlusion Check:
            if (FinalOcclusion(PathVertices[0], PathVertices[1], Sequence[0], ref Trans_Mod, Threadid)) return;
            //Specular_Path SP = new Specular_Path(PathVertices, Sequence, Seq_Polys, Room, Src, Speed_of_Sound, Trans_Mod, ref Direct_Time[rec_id], Threadid, Rnd[Threadid].Next());
            //ThreadPaths[rec_id, Threadid].Add(SP);
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
            RefPoint = Src.H_Origin();
            int[] Seq_Polys = new int[Sequence.Length];
            double[] Trans_Mod = new double[8];
            for (int t_oct = 0; t_oct < 8; t_oct++) Trans_Mod[t_oct] = 1;

            //Find all Source Images
            for (int q = 0; q < Sequence.Length; q++)
            {
                //RefPM = Room.PlaneMembers[Sequence[q]];
                Images[q] = Room.Image(RefPoint, 0, Room.ObjectMembers[Sequence[q]][0]);
                RefPoint = Images[q];
            }

            ProcessImages(Images, Sequence, rec_id, Threadid);

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
                    double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad();
                    for (int oct = 0; oct < 8; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * Room.TransmissionValue[X.Poly_id][oct];
                    R.origin = X.X_Point;
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
        private bool FinalOcclusion(Hare.Geometry.Point Src, Hare.Geometry.Point EndPt, double tol, ref double[] Trans_Mod, int Thread_Id) //int Poly_1, int Poly_2,
        {
            Hare.Geometry.Vector D = (Src - EndPt);
            double L = D.Length();
            double L2 = 0;
            D.Normalize();
            Ray R = new Ray(EndPt, D, Thread_Id, Rnd[Thread_Id].Next());

            Hare.Geometry.X_Event X;
            do
            {
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Change this so that it checks for the edge the BREPS the edge is on when x.t is less than tolerance.
                ////////////////////////////////////////////////////////////////////////////////////////////////////////

                if (Room.shoot(R, 0, out X))
                {
                    if (Room.IsTransmissive[X.Poly_id])//&& (Room.BrepID(X.Poly_id) == Poly_1 || Room.BrepID(X.Poly_id) == Poly_2)
                    {
                        ///The ray hit something transparent. Account for this...
                        double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad();
                        for (int oct = 0; oct < 8; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * Room.TransmissionValue[X.Poly_id][oct];
                        R.origin = X.X_Point;
                        R.Ray_ID = Rnd[Thread_Id].Next();
                        L2 += X.t;
                        continue;
                    }
                    else if (X.t < tol)
                    {
                        R.origin = X.X_Point;
                        R.Ray_ID = Rnd[Thread_Id].Next();
                        L2 += X.t;
                        continue;
                    }
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
            D /= L;
            L *= .99; //For tolerance...
            Ray R = new Ray(Origin, D, Thread_Id, Rnd[Thread_Id].Next());

            Hare.Geometry.X_Event X = new Hare.Geometry.X_Event();
            do
            {
                if (Room.shoot(R, 0, out X) && Room.IsTransmissive[X.Poly_id] && Poly_X != Room.ObjectID(X.Poly_id))
                {
                    double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad();
                    for (int oct = 0; oct < 8; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * (Room.TransmissionValue[X.Poly_id][oct]);
                    R.origin = X.X_Point;
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
        private bool OcclusionIntersectED(Hare.Geometry.Point Origin, Hare.Geometry.Point EndPt, int Poly_X, ref double[] Trans_Mod, ref Hare.Geometry.Point X_Point, int Thread_Id)//ref double[] Trans_Mod
        {
            Hare.Geometry.Vector D = (EndPt - Origin);
            double L = D.Length();
            D /= L;
            Ray R = new Ray(Origin, D, Thread_Id, Rnd[Thread_Id].Next());
            Hare.Geometry.X_Event X = new Hare.Geometry.X_Event();

            do
            {
                if (Room.shoot(R, 0, out X) && Room.IsTransmissive[X.Poly_id])// && Poly_X != Room.PlaneID(X.Poly_id))
                {
                    double[] Absorption = Room.AbsorptionValue[X.Poly_id].Coefficient_A_Broad();
                    for (int oct = 0; oct < 8; oct++) Trans_Mod[oct] *= (1 - Absorption[oct]) * Room.TransmissionValue[X.Poly_id][oct];
                    R.origin = X.X_Point;
                    R.Ray_ID = Rnd[Thread_Id].Next();
                    L -= X.t;
                    continue;
                }
                break;
            } while (true);

            if (X.t < L - 0.0001 && X.t != 0) return false;
            //if (Poly_X != Room.PlaneID(X.Poly_id))
            //{
            //    //Guid G1 = Rhino.RhinoDoc.ActiveDoc.Objects.Add(new Rhino.Geometry.LineCurve(Utilities.PachTools.HPttoRPt(X.X_Point), Utilities.PachTools.HPttoRPt(EndPt)));
            //    //Guid G2 = Rhino.RhinoDoc.ActiveDoc.Objects.Add(new TextDot(Room.PlaneID(X.Poly_id).ToString(), Utilities.PachTools.HPttoRPt(EndPt)));
            //    //Rhino.RhinoDoc.ActiveDoc.Groups.Add(new Guid[2] { G1, G2 });
            //    return false;
            //}
            //Poly_Seq = X.Poly_id;
            //X_Point = X.X_Point;
            X_Point = EndPt;
            return true;
        }

        /// <summary>
        /// Write the image source paths to a file for storage.
        /// </summary>
        /// <param name="BW">The binary writer to which data will be written.</param>
        public void Write_Data(ref System.IO.BinaryWriter BW)
        {
            //1. Write an indicator to signify that there is image source data
            BW.Write("Image-Source_Data");

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
                        BW.Write(ValidPaths[q][i].Energy(0, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(1, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(2, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(3, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(4, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(5, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(6, 44100)[0]);
                        BW.Write(ValidPaths[q][i].Energy(7, 44100)[0]);

                        //6a.2 Write a bool for whether it has a special materials filter...
                        if ((ValidPaths[q][i] as Specular_Path).Special_Filter != null)
                        {
                            BW.Write(true);
                            BW.Write((ValidPaths[q][i] as Specular_Path).Special_Filter.Length);
                            foreach(System.Numerics.Complex val in (ValidPaths[q][i] as Specular_Path).Special_Filter)
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
                    }
                    else if (ValidPaths[q][i] is Compound_Path)
                    {
                        ////TODO: Find a robust format for compound reflection paths...
                        /////Compound Path:
                        //BW.Write((short)1);
                        ////Write the number of samples and the pressure signal down:
                        //BW.Write((ValidPaths[q][i] as Compound_Path).P.Length);
                        //foreach (double val in (ValidPaths[q][i] as Compound_Path).P)
                        //{
                        //    BW.Write(val);
                        //}
                        ////Write the number of samples and the pressure signal down:
                        //foreach(Hare.Geometry.Vector vct in (ValidPaths[q][i] as Compound_Path).Directions)
                        //{
                        //    BW.Write(vct.x);
                        //    BW.Write(vct.y);
                        //    BW.Write(vct.z);
                        //}
                    }

                    //7. Write the arrival time:double
                    BW.Write(ValidPaths[q][i].TravelTime);

                    //8. Write the Reflection Sequence:int
                    for (int r = 0; r < ValidPaths[q][i].Reflection_Sequence.Length; r++)
                    {
                        BW.Write(ValidPaths[q][i].Reflection_Sequence[r]);
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
        public static ImageSourceData Read_Data(ref System.IO.BinaryReader BR, int Rec_CT, Direct_Sound Direct, bool Edges, int Src_ID, string version)
        {
            ImageSourceData IS = new ImageSourceData();
            IS.ValidPaths = new List<Deterministic_Reflection>[Rec_CT];
            IS.Direct_Time = new double[Rec_CT];
            double v = double.Parse(version.Substring(0, 3));

            for (int q = 0; q < Rec_CT; q++)
            {
                IS.ValidPaths[q] = new List<Deterministic_Reflection>();
                //2. Write the receiver number:int
                BR.ReadInt32();
                //3. Write number of paths:int
                int PathCt = BR.ReadInt32();
                for (int i = 0; i < PathCt; i++)
                {
                    int ReflectionType = BR.ReadInt16();
                    if (ReflectionType == 0)
                    {
                        //Speculare Reflection
                        //4. Write the number of reflection path points
                        Hare.Geometry.Point[] PTS = new Hare.Geometry.Point[BR.ReadInt32()];

                        //5. Write the reflection path:double
                        for (int r = 0; r < PTS.Length; r++)
                        {
                            PTS[r] = new Hare.Geometry.Point(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble());
                        }

                        //Previously, Pachyderm performed the deterministic part in intensity only...                    
                        //6a. Write the energy values
                        double[] Energy = new double[8];
                        Energy[0] = BR.ReadDouble();
                        Energy[1] = BR.ReadDouble();
                        Energy[2] = BR.ReadDouble();
                        Energy[3] = BR.ReadDouble();
                        Energy[4] = BR.ReadDouble();
                        Energy[5] = BR.ReadDouble();
                        Energy[6] = BR.ReadDouble();
                        Energy[7] = BR.ReadDouble();

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

                        double[] prms = new double[8];
                        for (int j = 0; j < prms.Length; j++) prms[j] = BR.ReadDouble();

                        //7. Write the arrival time:double
                        double Time = BR.ReadDouble();

                        //8. Write the Reflection Sequence:int
                        int[] Sequence = new int[PTS.Length - 2];
                        for (int r = 0; r < Sequence.Length; r++)
                        {
                            Sequence[r] = BR.ReadInt32();
                        }

                        IS.ValidPaths[q].Add(new Specular_Path(PTS, Energy, prms, Filter, Time, Sequence, Direct.Min_Time(q), Src_ID));
                    }
                    else if (ReflectionType == 1)
                    {
                        //TODO: Find a robust format for compound reflection paths...
                        ///Specular Path:
                        //BW.Write((short)1);
                        //Write the number of samples and the pressure signal down:
                        //Write the number of samples and the pressure signal down:

                        //6a.2. Write the number of samples in the pressure signal.(int)
                        //6b. Write the pressure values
                    }
                }
            }
            IS.Create_Filter(Direct.SWL, 4096);
            return IS;
        }

        public void Set_Power(double[] Factor)
        {
            for(int i = 0; i < Paths.Length; i++) for(int j = 0; j < Paths[i].Count; j++)
            {
                Paths[i][j].Set_Power(Factor);
            }
        }

        public void Create_Filter(double[] SWL, int length)
        {
            for (int i = 0; i < Paths.Length; i++)
            {
                foreach (Deterministic_Reflection P in Paths[i])
                {
                    P.Create_Filter(length, 0);
                }
            }
        }

        public Hare.Geometry.Vector[] Dir_Energy(int rec_id, int index, int Octave, double alt, double azi, bool degrees)
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
        public abstract Hare.Geometry.Vector[] Dir_Energy(int Octave, double alt, double azi, bool degrees);

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
        {
            PathEnergy = new double[8];
            ValidPath = Path;
            //Build an Identifier
            Sequence = Seq_planes;
            Hare.Geometry.Point Pt;

            for (int q = 1; q < ValidPath.Length; q++)
            {
                Pt = ValidPath[q] - ValidPath[q - 1];
                Length += Math.Sqrt(Pt.x * Pt.x + Pt.y * Pt.y + Pt.z * Pt.z);
            }

            Time = Length / C_Sound;
            Hare.Geometry.Vector DIR = ValidPath[1] - ValidPath[0];
            DIR.Normalize();
            Random rnd = new Random(Rnd);
            
            ///Energy based formulation
            Identify(Src.Source_ID(), Direct_Time);

            for (int oct = 0; oct < 8; oct++)
            {
                PathEnergy[oct] = Math.Pow(10, -.1 * Room.Attenuation(0)[oct] * Length) * H[0] * H[0];
                PathEnergy[oct] *= Trans_Mod[oct];
            }

            foreach (int q in Seq_Polys)
            {
                if (!(Room.AbsorptionValue[q] is Basic_Material)) continue;
                double[] AbsorptionData = Room.AbsorptionValue[q].Coefficient_A_Broad();
                double[] ScatteringData = Room.ScatteringValue[q].Coefficient();
                double[] TransmissionData = Room.TransmissionValue[q];
                for (int t = 0; t <= 7; t++)
                {
                    PathEnergy[t] *= (1 - AbsorptionData[t]) * (1 - ScatteringData[t]) * (1 - TransmissionData[t]);
                }
            }

            prms = new double[8];

            for (int i = 0; i < 8; i++) prms[i] = Math.Sqrt(PathEnergy[i] * Room.Rho_C(Path[0]));

            Special_Filter = new System.Numerics.Complex[4096];
            for (int i = 0; i < Special_Filter.Length; i++) Special_Filter[i] = 1;

            foreach (int q in Seq_Polys)
            {
                if (Room.AbsorptionValue[q] is Basic_Material) continue;
                //Pressure based formulation of materials
                //TODO: Coordinate with intensity based absorption above.
                for (int i = 0; i < Seq_Polys.Length; i++)
                {
                    Hare.Geometry.Vector d = Path[i + 1] - Path[i + 2]; d.Normalize();
                    if (!(Room.AbsorptionValue[Seq_Polys[i]] is Basic_Material))
                    {
                        System.Numerics.Complex[] Ref = Room.AbsorptionValue[Seq_Polys[i]].Reflection_Spectrum(44100, 4096, Room.Normal(Seq_Polys[i]), d, thread);
                        for (int j = 0; j < Special_Filter.Length; j++) Special_Filter[j] *= Ref[j];
                    }
                }
            }
        }

        public override void Set_Power(double[] factor)
        {
            for (int i = 0; i < PathEnergy.Length; i++) PathEnergy[i] *= factor[i];
        }

        public override void Create_Filter(int length, int threadid)
        {
            F = Audio.Pach_SP.Filter.Transfer_Function(prms, 44100, length, threadid);
        }

        public override double[] Create_Filter(double[] SWL, int SampleFrequency, int length, int dim, int threadid)
        {
            double[] tf_spec = new double[8];
            for (int i = 0; i < 8; i++) tf_spec[i] = prms[i] * Math.Pow(10, (120 - SWL[i]) / 20);

            return Audio.Pach_SP.Filter.Transfer_Function(tf_spec, length, SampleFrequency, threadid);
        }

        public override double[][] Create_Filter(double[] SWL, int length, int Sample_Frequnency, int threadid)
        {
            double[] tf_spec = new double[8];
            for (int i = 0; i < 8; i++) tf_spec[i] = prms[i] * Math.Pow(10, (120 - SWL[i]) / 20);

            return new double[1][] { Audio.Pach_SP.Filter.Transfer_Function(tf_spec, Sample_Frequnency, length, threadid) };
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
            if (Octave < 8) return new double[] { PathEnergy[Octave] };
            else return new double[] { PathEnergy.Sum() };
        }

        public override double[] Dir_Energy(int Octave, int dir)
        {
            Hare.Geometry.Vector Dir = Path[0][Path.Length - 1] - Path[0][Path.Length - 2];
            Dir.Normalize();
            switch (dir)
            {
                case 0:
                    return new double[] { Dir.x * PathEnergy[Octave] };
                case 1:
                    return new double[] { Dir.y * PathEnergy[Octave] };
                case 2:
                    return new double[] { Dir.z * PathEnergy[Octave] };
                default:
                    throw new Exception("indexed directions must conform to 0 = x, 1 = y and 2 = z") ;
            }
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave)
        {
            Hare.Geometry.Vector Dir = Path[0][Path[0].Length - 1] - Path[0][Path[0].Length - 2];
            Dir.Normalize();
            return new Hare.Geometry.Vector[] { Dir * PathEnergy[Octave] };
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, double alt, double azi, bool degrees)
        {
            Hare.Geometry.Vector[] V = Dir_Energy(Octave);
            return new Hare.Geometry.Vector[] { Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V[0], azi, 0, degrees), 0, alt, degrees)};
        }

        public override double[][] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, int SampleFreq, bool flat)
        {
            Hare.Geometry.Vector V = Path[0][Path[0].Length - 1] - Path[0][Path[0].Length - 2];
            V.Normalize();
            Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V, azi, 0, degrees), 0, alt, degrees);
            double[] F_Chosen = (SampleFreq == 44100 && flat) ? F : this.Create_Filter(SWL, 4096, SampleFreq, 0)[0];
            double[][] Fn = new double[F_Chosen.Length][];

            for (int i = 0; i < F_Chosen.Length; i++)
            {
                Fn[i] = new double[3] { Vn.x * F_Chosen[i], Vn.y * F_Chosen[i], Vn.z * F_Chosen[i] };
            }
            return Fn;
        }

        public override double[] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, bool Figure8, int SampleFreq, bool flat)
        {
            Hare.Geometry.Vector V = Path[0][Path[0].Length - 1] - Path[0][Path[0].Length - 2];
            V.Normalize();
            Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V, azi, 0, degrees), 0, alt, degrees);
            double[] F_Chosen = (SampleFreq == 44100 && flat) ? F : this.Create_Filter(SWL, 4096, SampleFreq, 0)[0];
            double[] Fn = new double[F_Chosen.Length];
            if (Figure8)
            {
                for (int i = 0; i < F_Chosen.Length; i++)
                {
                    Fn[i] = Vn.x * F_Chosen[i];
                }
            }
            if (Vn.x > 0)
            {
                for (int i = 0; i < F_Chosen.Length; i++)
                {
                    Fn[i] = Vn.x * F_Chosen[i];
                }
            }
            return Fn;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, Hare.Geometry.Vector V)
        {
            double l = Math.Sqrt(V.z * V.z + V.x * V.x);
            double azi = Math.Asin(V.y / l);
            double alt = Math.Atan2(V.x, V.z);
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
            double l = Math.Sqrt(V.z * V.z + V.x * V.x);
            double azi = Math.Asin(V.y / l);
            double alt = Math.Atan2(V.x, V.z);

            Hare.Geometry.Vector E = new Hare.Geometry.Vector();
            for (int oct = 0; oct < 8; oct++) E += Dir_Energy(oct, alt, azi, false)[0];
            return new Hare.Geometry.Vector[] { E };
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(double alt, double azi, bool degrees)
        {
            Hare.Geometry.Vector E = new Hare.Geometry.Vector();
            for (int oct = 0; oct < 8; oct++) E += Dir_Energy(oct, alt, azi, degrees)[0];
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
        private double[] SWL; //Source sound power.
        private Hare.Geometry.Point[][] ValidPath; //Geometry of the reflection.
        private double[][] PathEnergy;
        private double Time; //Earliest time of the reflection.
        private int[] Sequence; //Unique identifying indices for each reflecting element.
        public double[] F; //Filter form of the reflection time signature.
        public double[][] Fdir; //Directional form of F.
        public double[] H; //Initual H-function of the reflection (contains info about diffraction, compression, etc. that is a consequence of the resulting wave-form.)
        public double[][] Hdir; //Directional form of H.
        public double[][] Octave_Power; //Contains power modifying information for each sample in H (such as absorption or transmission coefficients, air attenuation, etc.)
        private double speed_of_sound;

        public Compound_Path(Scene Room, Hare.Geometry.Point[][] PathVertices, int[] Seq_Planes, Source Src, double[][] Power_mods, double[] _H, double[][] _Hdir, double T0, ref double Direct_Time, int Threadid)
        {
            //Reconcile Source Power as SWL with power with directivity. source directivity needs to be added to OctavePower.
            SWL = Src.SWL().Clone() as double[];
            Time = T0;
            speed_of_sound = Room.Sound_speed(0);
            ///Here, a compound reflection collector, and any interpolation that must be done.
            ValidPath = PathVertices;
            List<Hare.Geometry.Point[]> Paths = new List<Hare.Geometry.Point[]>();
            Sequence = Seq_Planes;
            H = _H;
            Hdir = _Hdir;

            double[] pd = new double[H.Length];
            double[] timeaxis = new double[H.Length];
            for (int t = 0; t < timeaxis.Length; t++) { timeaxis[t] = (double)t / 44100f; pd[t] = (double)H[t]; }
            Octave_Power = new double[H.Length + 2048][];

            List<double[]> pm = Power_mods.ToList();
            pm.RemoveAll(item => item == null);
            Power_mods = pm.ToArray();

            for (int i = 0; i < PathVertices.Length; i++)
            {
                Hare.Geometry.Vector dir = PathVertices[i][1] - PathVertices[i][0];
                dir.Normalize();
                double[] power = Src.DirPower(Threadid, 0, dir);

                Octave_Power[i] = new double[8];
                for (int oct = 0; oct < 8; oct++)
                {
                    Octave_Power[i][oct] = Math.Pow(10, -.1 * Room.Attenuation(0)[oct] * Room.Sound_speed(0) * ((float)i/44100f));
                    Octave_Power[i][oct] *= power[oct] * Power_mods[i][oct];
                }
            }

            //Build an Identifier
            Identify(Src.Source_ID(), Direct_Time);
        }
        
        public override double[][] Create_Filter(double[] SWL, int SampleFrequency, int LengthofPulse, int Threadid)
        {
            double[] tf_spec = new double[8];
            for (int i = 0; i < 8; i++) tf_spec[i] = Math.Pow(10, (120 - SWL[i]) / 20); //When this becomes multi-order, we will need - SWL * OctavePower[i]

            double[][] Fdir_out = new double[6][];
            double[] pulse = Audio.Pach_SP.Filter.Transfer_Function(tf_spec, SampleFrequency, LengthofPulse, Threadid);
            double[] H_FS = Pachyderm_Acoustic.Audio.Pach_SP.Resample(H, 88200, SampleFrequency, Threadid);
            double[] F_out = Audio.Pach_SP.FFT_Convolution_double(H_FS, pulse, Threadid);
            for (int i = 0; i < 6; i++)
            {
                double[] HDir_FS = Pachyderm_Acoustic.Audio.Pach_SP.Resample(Hdir[i], 88200, SampleFrequency, Threadid);
                Fdir_out[i] = Audio.Pach_SP.FFT_Convolution_double(HDir_FS, pulse, Threadid);
            }
            return new double[7][] { F_out, Fdir_out[0], Fdir_out[1], Fdir_out[2], Fdir_out[3], Fdir_out[4], Fdir_out[5] };
        }

        public override double[] Create_Filter(double[] SWL, int SampleFrequency, int LengthofPulse, int dim, int Threadid)
        {
            double[] tf_spec = new double[8];
            for (int i = 0; i < 8; i++) tf_spec[i] = Math.Pow(10, (120 - SWL[i]) / 20);
            double[] pulse = Audio.Pach_SP.Filter.Transfer_Function(tf_spec, SampleFrequency, LengthofPulse, Threadid);
            if (dim < 1)
            {
                double[] H_FS = Pachyderm_Acoustic.Audio.Pach_SP.Resample(H, 88200, SampleFrequency, Threadid);
                return Audio.Pach_SP.FFT_Convolution_double(H_FS, pulse, Threadid);
            }
            else
            {
                double[] HDir_FS = Pachyderm_Acoustic.Audio.Pach_SP.Resample(Hdir[dim - 1], 88200, SampleFrequency, Threadid);
                return Audio.Pach_SP.FFT_Convolution_double(HDir_FS, pulse, Threadid);
            }
        }

        public override void Create_Filter(int length, int threadid)
        {
            double[] tf_spec = new double[8];
            for (int i = 0; i < 8; i++) tf_spec[i] = Octave_Power[0][i] * 20; //When this becomes multi-order, we will need - SWL * OctavePower[i]

            Fdir = new double[6][];
            double[] pulse = Audio.Pach_SP.Filter.Transfer_Function(tf_spec, 88200, length, 0);
            double[] H_FS = (H.Length < 5) ? new double[1]{ H.Sum() } : Pachyderm_Acoustic.Audio.Pach_SP.Resample_Cubic(H, 88200, 44100, threadid);

            F = Audio.Pach_SP.FFT_Convolution_double(H_FS, pulse, 0);
            for (int i = 0; i < 6; i++)
            {
                double[] HDir_FS = (Hdir[i].Length < 5) ? new double[1] { Hdir[i].Sum() } : Pachyderm_Acoustic.Audio.Pach_SP.Resample(Hdir[i], 88200, 44100, threadid);
                Fdir[i] = Audio.Pach_SP.FFT_Convolution_double(HDir_FS, pulse, threadid);
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
            double[] H_FS;
            double[] I = new double[Octave_Power.Length];
            if (Octave < 8)
            {
                //TODO: convolve each sample individually (may want to precompute...
                for(int t = 0; t < Octave_Power.Length; t++) if (Octave_Power[t] != null) I[t] = Octave_Power[t][Octave];
                H_FS = Pachyderm_Acoustic.Audio.Pach_SP.FIR_Bandpass(H, Octave, 88200, 0);
            }
            else
            {
                for (int t = 0; t < Octave_Power.Length; t++) if (Octave_Power[t] != null) I[t] = Octave_Power[t].Sum();
                H_FS = H;
            }
            //TODO - won't interpolate over dirac pulses... Fix it...
            H_FS = (H_FS.Length < 5) ? new double[1] { H_FS.Sum() } : Pachyderm_Acoustic.Audio.Pach_SP.Resample_Cubic(H_FS, 88200, Sample_Frequency, 0);
            double[] energy = new double[H_FS.Length];

            for (int i = 0; i < H_FS.Length; i++)
            {
                energy[i] +=  I[i] * (H_FS[i] * H_FS[i]);//I[i > I.Length-1? I.Length - 1: i] * 
            }
            return energy;
        }

        public override double[] Dir_Energy(int Octave, int dir)
        {
            double[] De = new double[PathEnergy.Length];
            for (int i = 0; i < PathEnergy.Length; i++)
            {
                Hare.Geometry.Vector Dir = Path[i][Path.Length - 1] - Path[i][Path.Length - 2];
                Dir.Normalize();
                switch (dir)
                {
                    case 0:
                        De[i] = Dir.x * PathEnergy[i][Octave]; break;
                    case 1:
                        De[i] = Dir.y * PathEnergy[i][Octave]; break;
                    case 2:
                        De[i] = Dir.z * PathEnergy[i][Octave]; break;
                    default:
                        throw new Exception("indexed directions must conform to 0 = x, 1 = y and 2 = z");
                }
            }
            return De;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave)
        {
            Hare.Geometry.Vector[] De = new Hare.Geometry.Vector[PathEnergy.Length];
            for (int i = 0; i < PathEnergy.Length; i++)
            {
                Hare.Geometry.Vector Dir = Path[i][Path[i].Length - 1] - Path[i][Path[i].Length - 2];
                Dir.Normalize();
                De[i] = Dir * PathEnergy[i][Octave];
            }
            return De;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, double alt, double azi, bool degrees)
        {
            Hare.Geometry.Vector[] V = Dir_Energy(Octave);
            for (int i = 0; i < V.Length; i++)
            {
                 Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V[i], azi, 0, degrees), 0, alt, degrees);
            }
            return V;
        }

        public override Hare.Geometry.Vector[] Dir_Energy(int Octave, Hare.Geometry.Vector V)
        {
            double l = Math.Sqrt(V.z * V.z + V.x * V.x);
            double azi = Math.Asin(V.y / l);
            double alt = Math.Atan2(V.x, V.z);
            return Dir_Energy(Octave, alt, azi, false);
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(Hare.Geometry.Vector V)
        {
            double l = Math.Sqrt(V.z * V.z + V.x * V.x);
            double azi = Math.Asin(V.y / l);
            double alt = Math.Atan2(V.x, V.z);

            return Dir_EnergySum(alt, azi, false);
        }

        public override Hare.Geometry.Vector[] Dir_EnergySum(double alt, double azi, bool degrees)
        {
            Hare.Geometry.Vector[] E = new Hare.Geometry.Vector[PathEnergy.Length];
            for(int i = 0; i < E.Length; i++) E[i] = new Hare.Geometry.Vector();
            
            for (int oct = 0; oct < 8; oct++)
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
            for (int i = 0; i < PathEnergy.Length; i++) for (int j = 0; j < PathEnergy[i].Length; j++) PathEnergy[i][j] *= factor[j];
        }

        public override double[] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, bool Figure8, int sampleFreq, bool flat)
        {
            double[] Fn = new double[H.Length];

            double[][] Hdir_out = (sampleFreq == 44100 && flat) ? Hdir : Create_Filter(SWL, sampleFreq, 4096, 0);

            if (Figure8)
            {
                for (int i = 0; i < Hdir_out.Length; i++)
                {
                    Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[i][0] - Hdir_out[i][1], Hdir_out[i][2] - Hdir_out[i][3], Hdir_out[i][4] - Hdir_out[i][5]), azi, 0, degrees), 0, alt, degrees);
                    Fn[i] = Vn.x;
                }
            }
            else
            {
                int[] ids = new int[3];
                ids[0] = (azi > 90 && azi < 270) ? 1 : 0;
                ids[0] = (azi <= 180) ? 3 : 4;
                ids[0] = (alt > 0) ? 4 : 5;
                for (int i = 0; i < Hdir_out[0].Length; i++)
                {
                    Hare.Geometry.Vector V = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[ids[0]][i], Hdir_out[ids[1]][i], Hdir_out[ids[2]][i]), azi, 0, true), 0, alt, true);
                    Fn[i] = V.x;
                }
            }
            return Fn;
        }

        public override double[][] Dir_Filter(double[] SWL, double alt, double azi, bool degrees, int sampleFreq, bool flat)
        {
            double[][] Hdir_out = (sampleFreq == 44100 && flat) ? Hdir : Create_Filter(SWL, sampleFreq, 4096, 0);
            double[][] Fn = new double[H.Length][];

            for (int i = 0; i < Hdir.Length; i++)
            {
                Hare.Geometry.Vector Vn = Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(new Hare.Geometry.Vector(Hdir_out[i][0] - Hdir_out[i][1], Hdir_out[i][2] - Hdir_out[i][3], Hdir_out[i][4] - Hdir_out[i][5]), azi, 0, degrees), 0, alt, degrees);
                Fn[i] = new double[3] { Vn.x, Vn.y, Vn.z };
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