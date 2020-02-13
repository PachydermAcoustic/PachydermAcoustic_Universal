//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2019, Arthur van der Harten 
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

using System;
using System.Collections.Generic;
using Hare.Geometry;
using Pachyderm_Acoustic.Environment;
using Pachyderm_Acoustic.Utilities;
using System.Linq;

namespace Pachyderm_Acoustic
{
    /// <summary>
    /// A receiver model which optimizes for very large numbers of receivers, as exist in mapping simulations.
    /// </summary>
    [Serializable]
    public class PachMapReceiver : Receiver_Bank
    {
        public string SrcType = "";
        public double[] SWL;
        public Hare.Geometry.Topology Map_Mesh;
        public AABB OBox;
        public int VoxelCtX, VoxelCtY, VoxelCtZ;
        public AABB[,,] Voxels;
        public List<int>[,,] Voxel_Inv;
        public Hare.Geometry.Point BoxDims;
        public Hare.Geometry.Point BoxDims_Inv;
        public Hare.Geometry.Point VoxelDims;
        public Hare.Geometry.Point VoxelDims_Inv;
        public double increment;
        public bool Z_Displacement;
        public bool Directional;
        public bool Rec_Vertex;
        public bool Mesh_Offset;
        public bool Time1Pt;
        public Hare.Geometry.Point Src;
        public double delay_ms;
        Scene _Sc;

        /// <summary>
        /// Use this constructor when reading in a file.
        /// </summary>
        public PachMapReceiver()
        {

        }

        /// <summary>
        /// General use constructor. Use this prior to starting a mapping calculation. It works seamlessly with ray-tracing operations.
        /// </summary>
        public PachMapReceiver(Hare.Geometry.Topology Map_MeshIn, Source Src_Pt, int SampleRate_in, double Increment_in, Scene Sc, int RCT, double Cutoff_time, bool Time1Pt, bool Z_displacementIn, bool DirectionalIn, bool RecOnVertex, bool Offset_Mesh)
        {
            Map_Mesh = Map_MeshIn;
            Fill_in(Src_Pt, SampleRate_in, Increment_in, Sc, RCT, Cutoff_time, Time1Pt, Z_displacementIn, DirectionalIn, RecOnVertex, Offset_Mesh);
            Partition(Src_Pt);
        }

        private void Fill_in(Source Src_Pt, int SampleRate_in, double Increment_in, Scene Sc, int RCT, double Cutoff_time, bool Time_1Pt, bool Z_displacementIn, bool DirectionalIn, bool RecOnVertex, bool Offset_Mesh)
        {
            _Sc = Sc;
            SWL = new double[8] { Src_Pt.SWL(0), Src_Pt.SWL(1), Src_Pt.SWL(2), Src_Pt.SWL(3), Src_Pt.SWL(4), Src_Pt.SWL(5), Src_Pt.SWL(6), Src_Pt.SWL(7) };
            Src = Src_Pt.Origin();
            SrcType = Src_Pt.Type();
            CutOffTime = Cutoff_time;
            SampleRate = SampleRate_in;
            increment = Increment_in;
            Rec_Vertex = RecOnVertex;
            Mesh_Offset = Offset_Mesh;
            Z_Displacement = Z_displacementIn;
            Directional = DirectionalIn;
            Time1Pt = Time_1Pt;

            SampleCT = Time1Pt ? 1 : (int)Math.Floor(CutOffTime * SampleRate / 1000);
            Max = new Hare.Geometry.Point(Double.NegativeInfinity, Double.NegativeInfinity, Double.NegativeInfinity);
            Min = new Hare.Geometry.Point(Double.PositiveInfinity, Double.PositiveInfinity, Double.PositiveInfinity);
            Rec_List = Rec_Vertex ? new Map_Receiver[Map_Mesh.Vertex_Count] : new Map_Receiver[Map_Mesh.Polygon_Count];
        }

        private void Partition(Source Src_Pt)//, bool ProcessMesh)
        {
            int processorCt = Pach_Properties.Instance.ProcessorCount();

            //if (ProcessMesh)
            //{
                Point p;
            if (Z_Displacement)
            {
                Hare.Geometry.Point V = new Hare.Geometry.Point(0, 0, increment * .5);
                if (Rec_Vertex)
                {
                    for (int i = 0; i < Map_Mesh.Vertex_Count; i++)
                    {
                        p = new Hare.Geometry.Point(Map_Mesh[i].x, Map_Mesh[i].y, Map_Mesh[i].z) + V;
                        Rec_List[i] = new Map_Receiver(increment, i, p, Src_Pt, _Sc.Sound_speed(p), _Sc.Rho(p), _Sc.Attenuation(p), 1000, CutOffTime, processorCt, Time1Pt, Directional);
                    }
                }
                else
                {
                    for (int i = 0; i < Map_Mesh.Polygon_Count; i++)
                    {
                        Point center = new Point();
                        center = Map_Mesh.Polygon_Vertices(i)[0] + Map_Mesh.Polygon_Vertices(i)[1] + Map_Mesh.Polygon_Vertices(i)[2];

                        if (Map_Mesh.Polygon_Vertices(i).Count() == 4)
                        {
                            center += Map_Mesh.Polygon_Vertices(i)[3];
                            center /= 4;
                        }
                        else { center /= 3; }

                        p = center + V;
                        Rec_List[i] = new Map_Receiver(increment, i, p, Src_Pt, _Sc.Sound_speed(p), _Sc.Rho(p), _Sc.Attenuation(p), 1000, CutOffTime, processorCt, Time1Pt, Directional);
                    }
                }
            }
            else
            {
                if (Rec_Vertex)
                {
                    for (int i = 0; i < Map_Mesh.Vertex_Count; i++)
                    {
                        Vector V = Map_Mesh.Vertex_Normals[i] * increment;
                        p = Map_Mesh[i] + V;
                        Rec_List[i] = new Map_Receiver(increment, i, p, Src_Pt, _Sc.Sound_speed(p), _Sc.Rho(p), _Sc.Attenuation(p), 1000, CutOffTime, processorCt, Time1Pt, Directional);
                    }
                }
                else
                {
                    for (int i = 0; i < Map_Mesh.Polys.Count; i++)
                    {
                        Point center = Map_Mesh[i, 0] + Map_Mesh[i, 1] + Map_Mesh[i, 2];

                        if (Map_Mesh.Polys[i].VertextCT == 4)
                        {
                            center += Map_Mesh[i, 3];
                            center /= 4;
                        }
                        else { center /= 3; }

                        Vector V = Map_Mesh.Normal(i);

                        p = center + V;
                        Rec_List[i] = new Map_Receiver(increment, i, p, Src_Pt, _Sc.Sound_speed(p), _Sc.Rho(p), _Sc.Attenuation(p), 1000, CutOffTime, processorCt, Time1Pt, Directional);
                    }
                }

                if (Mesh_Offset)
                {

                    for (int i = 0; i < Map_Mesh.Vertex_Count; i++)
                    {
                        Vector V = Map_Mesh.Vertex_Normals[i] * increment;
                        Map_Mesh[i].x = Map_Mesh[i].x + V.x;
                        Map_Mesh[i].y = Map_Mesh[i].y + V.y;
                        Map_Mesh[i].z = Map_Mesh[i].z + V.z;
                    }
                }
            }

            foreach (Spherical_Receiver R in Rec_List)
            {
                //Find max and min bounds of all spheres...
                if ((R.Origin.x + (increment * .5)) > Max.x) Max.x = R.Origin.x + (increment * .5);
                if ((R.Origin.y + (increment * .5)) > Max.y) Max.y = R.Origin.y + (increment * .5);
                if ((R.Origin.z + (increment * .5)) > Max.z) Max.z = R.Origin.z + (increment * .5);
                if ((R.Origin.x - (increment * .5)) < Min.x) Min.x = R.Origin.x - (increment * .5);
                if ((R.Origin.y - (increment * .5)) < Min.y) Min.y = R.Origin.y - (increment * .5);
                if ((R.Origin.z - (increment * .5)) < Min.z) Min.z = R.Origin.z - (increment * .5);
            }
            OBox = new AABB(Min, Max);

            //Divide the Min and Max into separate sub-voxels...
            VoxelCtX = (int)Math.Ceiling((OBox.Max_PT.x - OBox.Min_PT.x) / (increment * 4));
            VoxelCtY = (int)Math.Ceiling((OBox.Max_PT.y - OBox.Min_PT.y) / (increment * 4));
            VoxelCtZ = (int)Math.Ceiling((OBox.Max_PT.z - OBox.Min_PT.z) / (increment * 2.5));

            Voxels = new AABB[VoxelCtX, VoxelCtY, VoxelCtZ];
            Voxel_Inv = new List<int>[VoxelCtX, VoxelCtY, VoxelCtZ];

            BoxDims = (Max - Min);
            VoxelDims = new Hare.Geometry.Point(BoxDims.x / VoxelCtX, BoxDims.y / VoxelCtY, BoxDims.z / VoxelCtZ);
            VoxelDims_Inv = new Hare.Geometry.Point(1 / VoxelDims.x, 1 / VoxelDims.y, 1 / VoxelDims.z);
            BoxDims_Inv = new Hare.Geometry.Point(1 / BoxDims.x, 1 / BoxDims.y, 1 / BoxDims.z);

            for (int x = 0; x < VoxelCtX; x++)
            {
                for (int y = 0; y < VoxelCtY; y++)
                {
                    for (int z = 0; z < VoxelCtZ; z++)
                    {
                        Voxel_Inv[x, y, z] = new List<int>();
                        Hare.Geometry.Point voxelmin = new Hare.Geometry.Point(x * VoxelDims.x, y * VoxelDims.y, z * VoxelDims.z);
                        Hare.Geometry.Point voxelmax = new Hare.Geometry.Point((x + 1) * VoxelDims.x, (y + 1) * VoxelDims.y, (z + 1) * VoxelDims.z);
                        AABB Box = new AABB(voxelmin + OBox.Min_PT, voxelmax + OBox.Min_PT);
                        Voxels[x, y, z] = Box;
                        for (int i = 0; i < Rec_List.Length; i++)
                        {
                            //Check for intersection between voxel x,y,z with Receiver i...
                            Hare.Geometry.Point PT = Box.ClosestPt(Rec_List[i].Origin);
                            PT -= Rec_List[i].Origin;
                            if ((PT.x * PT.x + PT.y * PT.y + PT.z * PT.z) < Rec_List[i].Radius2)
                            {
                                Voxel_Inv[x, y, z].Add(i);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Adds the direct sound to the collection of histograms.
        /// </summary>
        /// <param name="D">a completed direct sound simulation</param>
        /// <param name="S">the source object</param>
        public void AddDirect(Direct_Sound D, Source S)
        {
            for (int i = 0; i < Rec_List.Length; i++)
            {
                Vector dir = S.H_Origin() - Rec_List[i].Origin;
                dir.Normalize();
                for (int oct = 0; oct < 8; oct++)
                {
                    int tsample = 0;
                    int length = D.Io[i][oct].Length;
                    for (int s = 0; s < length; s++)
                    {
                        if (D.EnergyValue(oct, s, i) > 0 && tsample == 0) tsample = s;
                        if (s < Rec_List[i].Recs.SampleCT && D.Validity[i]) Rec_List[i].Combine_Sample(s - tsample, D.EnergyValue(oct, s, i), dir, dir, oct);
                    }
                }
            }
        }

        /// <summary>
        /// Checks receivers for intersections with a Broadband Ray (specular only). Records detections. Use after shooting a ray, but before adding the distance of ray travel to the total of distance traveled.
        /// </summary>
        /// <param name="R">a broadband ray.</param>
        /// <param name="Sumlength">the length of the ray before detection. Adds on the distance from the last reflection to the receiver.</param>
        /// <param name="EndPt">The point of the d</param>
        /// <param name="X">The current voxel index 'X'</param>
        /// <param name="Y">The current voxel index 'Y'</param>
        /// <param name="Z">The current voxel index 'Z'</param>
        private void CheckReceivers_Broadband(BroadRay R, Hare.Geometry.Point EndPt, int X, int Y, int Z)
        {
            //PachTools.AddBox(Voxels[X, Y, 0].Min(), Voxels[X, Y, 0].Max());//Debug tool...
            for (int i = 0; i < Voxel_Inv[X, Y, Z].Count; i++)
            {
                if (!(Rec_List[Voxel_Inv[X, Y, Z][i]].Ray_ID[R.ThreadID] == R.Ray_ID))
                {
                    Rec_List[Voxel_Inv[X, Y, Z][i]].CheckBroadbandRay(R, EndPt);
                    //    //Debug Code
                    //    ///////////////////////////////////////////////////////////////////
                    //    //OnSphere sphere = new OnSphere(Rec_List[Voxel_Inv[X, Y, Z][i]].Origin,Rec_List[Voxel_Inv[X, Y, Z][i]].Radius);
                    //    //OnRevSurface sphere_srf = sphere.RevSurfaceForm();
                    //    //MRhinoSurfaceObject sphere_obj = new MRhinoSurfaceObject();
                    //    //sphere_obj.SetSurface(sphere_srf);
                    //    //RhUtil.RhinoApp().ActiveDoc().AddObject(sphere_obj);
                    //    ///////////////////////////////////////////////////////////////////
                    //}
                    Rec_List[Voxel_Inv[X, Y, Z][i]].Ray_ID[R.ThreadID] = R.Ray_ID;
                }
            }
        }

        /// <summary>
        /// Checks receivers for intersections with a single octave ray (specular only). Records detections. Use after shooting a ray, but before adding the distance of ray travel to the total of distance traveled.
        /// </summary>
        /// <param name="R">a single octave ray.</param>
        /// <param name="Sumlength">the length of the ray before detection. Adds on the distance from the last reflection to the receiver.</param>
        /// <param name="EndPt">The point of the d</param>
        /// <param name="X">The current voxel index 'X'</param>
        /// <param name="Y">The current voxel index 'Y'</param>
        /// <param name="Z">The current voxel index 'Z'</param>
        private void CheckReceivers_Octave(OctaveRay R, Hare.Geometry.Point EndPt, int X, int Y, int Z)
        {
            //PachTools.AddBox(Voxels[X, Y, 0].Min(), Voxels[X, Y, 0].Max());//Debug tool...
            for (int i = 0; i < Voxel_Inv[X, Y, Z].Count; i++)
            {
                if (!(Rec_List[Voxel_Inv[X, Y, Z][i]].Ray_ID[R.ThreadID] == R.Ray_ID))
                {
                    Rec_List[Voxel_Inv[X, Y, Z][i]].CheckRay(R, EndPt);
                    //    //Debug Code
                    //    ///////////////////////////////////////////////////////////////////
                    //    //OnSphere sphere = new OnSphere(Rec_List[Voxel_Inv[X, Y, Z][i]].Origin,Rec_List[Voxel_Inv[X, Y, Z][i]].Radius);
                    //    //OnRevSurface sphere_srf = sphere.RevSurfaceForm();
                    //    //MRhinoSurfaceObject sphere_obj = new MRhinoSurfaceObject();
                    //    //sphere_obj.SetSurface(sphere_srf);
                    //    //RhUtil.RhinoApp().ActiveDoc().AddObject(sphere_obj);
                    //    ///////////////////////////////////////////////////////////////////
                    //}
                    Rec_List[Voxel_Inv[X, Y, Z][i]].Ray_ID[R.ThreadID] = R.Ray_ID;
                }
            }
        }

        /// <summary>
        /// Create an integrated sound pressure level map based on the receiver positions extracted from the audience map surfaces.
        /// </summary>
        /// <param name="SPL_Bounds">The boundaries of SPL domain to be used.</param>
        /// <param name="T_Bounds">The boundaries of time domain to be used.</param>
        /// <param name="c">The Color Scale.</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <param name="SrcID">The IDs of the sources.</param>
        /// <erturns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_SPL_Map(PachMapReceiver[] Rec_List, double[] T_Bounds, int Octave, List<int> SrcID, bool Coherent_Superposition, bool ZeroAtDirect)
        {
            //T in ms.
            //Calculate SPL values...

            double[] SPL_Values = new double[Rec_List[0].Rec_List.Length];
            Topology MM = Rec_List[0].MapMesh();

            for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
            {
                double[] hist;

                if (Coherent_Superposition)
                {
                    hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];
                    double E_Sum = 0;
                    int T_Max = (int)Math.Floor((T_Bounds[0]) * 44.1);
                    int T_Min = (int)Math.Floor((T_Bounds[1]) * 44.1);

                    int zero = 0;
                    if (ZeroAtDirect)
                    {
                        double z = double.MaxValue;
                        foreach (int S_ID in SrcID)
                        {
                            double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                            if (z > t) z = t;
                        }
                        zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);
                    }

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    if (Octave < 8)
                    {
                        hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);
                    }
                    else
                    {
                        double[] histtemp;
                        double[] histtotal = new double[hist.Length];
                        for (int oct = 0; oct < 8; oct++)
                        {
                            histtemp = Audio.Pach_SP.FIR_Bandpass(hist, oct, 44100, 0);
                            for (int j = 0; j < histtemp.Length; j++) histtotal[j] += histtemp[j];
                        }
                        hist = histtotal;
                    }

                    for (int t = T_Min; t < T_Max; t++)
                    {
                        if (t >= hist.Length) break;
                        E_Sum += hist[t] * hist[t];
                    }
                    E_Sum /= Rec_List[0].Rec_List[i].Rho_C;

                    SPL_Values[i] = 10 * Math.Log10(E_Sum / 1E-12);
                }
                else
                {
                    hist = new double[Rec_List[0].SampleCT];
                    double E_Sum = 0;
                    double[] temp;
                    int T_Max = (int)Math.Floor(T_Bounds[0]);
                    int T_Min = (int)Math.Floor(T_Bounds[1]);

                    int zero = 0;
                    if (ZeroAtDirect)
                    {
                        double z = double.MaxValue;
                        foreach (int S_ID in SrcID)
                        {
                            double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                            if (z > t) z = t;
                        }
                        zero = (int)(z * Rec_List[0].SampleRate);
                    }

                    foreach (int S_ID in SrcID)
                    {
                        temp = Rec_List[S_ID].GetEnergyHistogram(Octave, Rec_List[S_ID].delay_ms, i);
                        int arrival = (int)((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate) - zero;
                        for (int j = 0; j < temp.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;
                            hist[t] += temp[j];
                        }
                    }

                    for (int t = T_Min; t < T_Max; t++)
                    {
                        if (t >= hist.Length) break;
                        E_Sum += hist[t];
                    }

                    SPL_Values[i] = AcousticalMath.SPL_Intensity(E_Sum);
                }
            }
            return SPL_Values;
        }

        /// <summary>
        /// Create an integrated sound pressure level map based on the receiver positions extracted from the audience map surfaces.
        /// </summary>
        /// <param name="SPL_Bounds">The boundaries of SPL domain to be used.</param>
        /// <param name="T_Bounds">The boundaries of time domain to be used.</param>
        /// <param name="c">The Color Scale.</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <param name="SrcID">The IDs of the sources.</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_SPLA_Map(PachMapReceiver[] Rec_List, double[] T_Bounds, List<int> SrcID, bool Coherent_Superposition)
        {
            //T in ms.
            //Calculate SPL values...
            int T_Max = (int)Math.Min(Rec_List[SrcID[0]].SampleCT, Math.Floor((T_Bounds[0])));
            int T_Min = (int)Math.Floor((T_Bounds[1]));
            double[] SPL_Values = new double[Rec_List[0].Rec_List.Length];

            double[] AFactors = new double[8] { Math.Pow(10, (-26.2 / 10)), Math.Pow(10, (-16.1 / 10)), Math.Pow(10, (-8.6 / 10)), Math.Pow(10, (-3.2 / 10)), 1, Math.Pow(10, (1.2 / 10)), Math.Pow(10, (1 / 10)), Math.Pow(10, (-1.1 / 10)) };
            for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
            {
                if (Coherent_Superposition)
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1);

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    double[] histtemp;
                    double histtotal = 0;
                    for (int oct = 0; oct < 8; oct++)
                    {
                        histtemp = Audio.Pach_SP.FIR_Bandpass(hist, oct, 44100, 0);
                        for (int j = 0; j < histtemp.Length; j++) histtotal += histtemp[j] * AFactors[oct];
                    }

                    SPL_Values[i] = 10 * Math.Log10(histtotal / 1E-12);
                }
                else
                {
                    double E_Sum = 0;
                    foreach (int S_ID in SrcID)
                    {
                        for (int oct = 0; oct < 8; oct++)
                        {
                            double E_oct_Sum = 0;
                            double[] hist = Rec_List[S_ID].GetEnergyHistogram(oct, Rec_List[S_ID].delay_ms, i);
                            for (int t = T_Min; t < T_Max; t++)
                            {
                                E_oct_Sum += hist[t];
                            }
                            E_Sum += E_oct_Sum * AFactors[oct];
                        }
                    }

                    SPL_Values[i] = AcousticalMath.SPL_Intensity(Math.Abs(E_Sum));
                }
            }

            return SPL_Values;
        }

        /// <summary>
        /// Plot energy directions based on the time range specified.
        /// </summary>
        /// <param name="SPL_Bounds">The boundaries of SPL domain to be used.</param>
        /// <param name="T_Bounds">The boundaries of time domain to be used.</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static Vector[] Get_Directional_Map(PachMapReceiver[] Rec_List, double[] T_Bounds, int Octave, List<int> SrcID)
        {
            //T in ms.
            //Calculate SPL values...
            int T_Max = (int)Math.Floor((T_Bounds[0]));// / (double)ms_per_bin);
            int T_Min = (int)Math.Floor((T_Bounds[1]));// / (double)ms_per_bin);
            Vector[] dir = new Vector[Rec_List[0].Rec_List.Length];

            for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
            {
                Vector V = new Vector();
                foreach (int S_ID in SrcID)
                {
                    double[] Hist = Rec_List[S_ID].GetEnergyHistogram(Octave, Rec_List[S_ID].delay_ms, i);
                    for (int t = T_Min; t < T_Max; t++)
                    {
                        V += Rec_List[S_ID].Directions_Pos(Octave, t, i);
                        V += Rec_List[S_ID].Directions_Neg(Octave, t, i);
                    }
                }
                V.Normalize();
                dir[i] = V;
            }

            return dir;
        }

        /// <summary>
        /// Assigns colors to parameter values on a Rhinoceros mesh.
        /// </summary>
        /// <param name="Values">The values of the data points</param>
        /// <param name="Bounds">Parameter bounds</param>
        /// <param name="H_OFFSET">H value offset (HSV Color system.)</param>
        /// <param name="H_BREADTH">H value breadth (HSV Color system.)</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <returns>the array of color assignments which fit neatly into the mesh </returns>
        public static System.Drawing.Color[] SetColors(double[] Values, double[] Bounds, Pach_Graphics.colorscale c)
        {
            System.Drawing.Color[] Mesh_Colors = new System.Drawing.Color[Values.Length];
            double Scale_Breadth = Bounds[1] - Bounds[0];

            for (int i = 0; i < Values.Length; i++)
            {
                System.Drawing.Color color = c.GetValue(Values[i], Bounds[0], Bounds[1]);
                Mesh_Colors[i] = color;
            }
            return Mesh_Colors;
        }

        public static double[] Get_EchoCritPercent_Map(PachMapReceiver[] Rec_List, int Octave, List<int> SrcID)
        {
            //T in ms.
            //Calculate C values...
            double[] E_Values = new double[Rec_List[0].Rec_List.Length];

            bool Echo10, Echo50;
            double[] EKG, PercEcho;

            for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
            {
                double[] hist = null;
                //double E_Sum = 0;

                int zero = 0;

                double z = double.MaxValue;
                foreach (int S_ID in SrcID)
                {
                    double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                    if (z > t) z = t;
                }
                zero = (int)(z * Rec_List[0].SampleRate); //*44.1

                hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1) + 4096];

                foreach (int S_ID in SrcID)
                {
                    double[] P;
                    Rec_List[S_ID].GetFilter(i, out P);
                    P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                    int arrival = (int)((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate) - zero;

                    for (int j = 0; j < P.Length; j++)
                    {
                        int t = j + arrival;
                        if (t >= hist.Length) break;

                        hist[t] += P[j];
                    }
                }

                if (Octave < 8) hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);

                AcousticalMath.EchoCriterion(hist, Rec_List[0].SampleRate, 0, true, out EKG, out PercEcho, out Echo10, out Echo50);
                E_Values[i] = PercEcho.Max();
            }

            return E_Values;
        }

        /// <summary>
        /// Creates a map of the C-80 values of the room
        /// </summary>
        /// <param name="C_Bounds"></param>
        /// <param name="C_Cutoff"></param>
        /// <param name="H_OFFSET">H value offset (HSV Color system.)</param>
        /// <param name="H_BREADTH">H value breadth (HSV Color system.)</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_Clarity_Map(PachMapReceiver[] Rec_List, int C_Cutoff, int Octave, List<int> SrcID, bool Coherent_Superposition)
        {
            //T in ms.
            //Calculate C values...
            double[] C_Values = new double[Rec_List[0].Rec_List.Length];

            if (Coherent_Superposition)
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    int zero = 0;
                    double z = double.MaxValue;
                    foreach (int S_ID in SrcID)
                    {
                        double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                        if (z > t) z = t;
                    }
                    zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    if (Octave < 8) hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);

                    for (int t = 0; t < hist.Length; t++) hist[t] *= hist[t];

                    C_Values[i] = AcousticalMath.Clarity(hist, Rec_List[0].SampleRate, (double)C_Cutoff/1000d, 0, false);
                }
            }
            else
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    C_Values[i] = AcousticalMath.Clarity(IR_Construction.ETCurve(null, null, Rec_List, Rec_List[0].CutOffTime, Rec_List[0].SampleRate, Octave, i, SrcID, true), Rec_List[0].SampleRate, (double)C_Cutoff / 1000d, 0, false);
                }
            }

            return C_Values;
        }

        /// <summary>
        /// Creates a map of the Definition value of the room.
        /// </summary>
        /// <param name="D_Bounds"></param>
        /// <param name="H_OFFSET">H value offset (HSV Color system.)</param>
        /// <param name="H_BREADTH">H value breadth (HSV Color system.)</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_Definition_Map(PachMapReceiver[] Rec_List, int Octave, List<int> SrcID, bool Coherent_Superposition)
        {
            //T in ms.
            //Calculate C values...
            double[] D_Values = new double[Rec_List[0].Rec_List.Length];

            if (Coherent_Superposition)
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    int zero = 0;
                    double z = double.MaxValue;
                    foreach (int S_ID in SrcID)
                    {
                        double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                        if (z > t) z = t;
                    }
                    zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    if (Octave < 8) hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);

                    for (int t = 0; t < hist.Length; t++) hist[t] *= hist[t];

                    D_Values[i] = AcousticalMath.Definition(hist, Rec_List[0].SampleRate, 0.05, 0, false);
                }
            }
            else
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    D_Values[i] = AcousticalMath.Definition(IR_Construction.ETCurve(null, null, Rec_List, Rec_List[0].CutOffTime, Rec_List[0].SampleRate, Octave, i, SrcID, true), Rec_List[0].SampleRate, 0.05, 0, false);
                }
            }

            return D_Values;
        }

        /// <summary>
        /// Creates a map of the Reverberation Time in the room.
        /// </summary>
        /// <param name="STI_Bounds"></param>
        /// <param name="H_OFFSET">H value offset (HSV Color system.)</param>
        /// <param name="H_BREADTH">H value breadth (HSV Color system.)</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_RT_Map(PachMapReceiver[] Rec_List, int Decay_Depth, int Octave, List<int> SrcID, bool Coherent_Superposition)
        {
            //T in ms.
            //Calculate T-X values...

            double[] RT_Values = new double[Rec_List[0].Rec_List.Length];

            if (Coherent_Superposition)
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    int zero = 0;
                    double z = double.MaxValue;
                    foreach (int S_ID in SrcID)
                    {
                        double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                        if (z > t) z = t;
                    }
                    zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    if (Octave < 8) hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);

                    for (int t = 0; t < hist.Length; t++) hist[t] *= hist[t];

                    double[] SI = AcousticalMath.Schroeder_Integral(hist);
                    RT_Values[i] = AcousticalMath.T_X(SI, Decay_Depth, Rec_List[0].SampleRate);
                }
            }
            else
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] SI = AcousticalMath.Schroeder_Integral(IR_Construction.ETCurve(null, null, Rec_List, Rec_List[0].CutOffTime, Rec_List[0].SampleRate, Octave, i, SrcID, true));
                    RT_Values[i] = AcousticalMath.T_X(SI, Decay_Depth, Rec_List[0].SampleRate);
                }
            }

            return RT_Values;
        }

        /// <summary>
        /// Creates a map of the Speech Transmission Index in the room.
        /// </summary>
        /// <param name="Rec_List"></param>
        /// <param name="STI_Bounds"></param>
        /// <param name="c"></param>
        /// <param name="NoiseSPL"></param>
        /// <param name="SrcID"></param>
        /// <param name="type"> 0 for 2003 general, 1 for Male, 2 for Female</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_STI_Map(PachMapReceiver[] Rec_List, double[] NoiseSPL, List<int> SrcID, int type, bool Coherent_Superposition)
        {
            double[] STI_Values = new double[Rec_List[0].Rec_List.Length];
            if (type < 0 || type > 2) return STI_Values;
            if (Coherent_Superposition)
            {
                System.Threading.Tasks.Parallel.For(0, Rec_List[0].Rec_List.Length, i =>
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    int zero = 0;
                    double z = double.MaxValue;
                    foreach (int S_ID in SrcID)
                    {
                        double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                        if (z > t) z = t;
                    }
                    zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    double[][] ETC = new double[8][];
                    for (int oct = 0; oct < 8; oct++)
                    {
                        double[] h = Audio.Pach_SP.FIR_Bandpass(hist, oct, 44100, 0);
                        for (int t = 0; t < hist.Length; t++) ETC[oct][t] = h[t] * h[t] / Rec_List[0].Rec_List[i].Rho_C;
                    }
                    STI_Values[i] = AcousticalMath.Speech_Transmission_Index(ETC, 1.2 * 343, NoiseSPL, Rec_List[0].SampleRate)[type];
                });
            }
            else
            {
                System.Threading.Tasks.Parallel.For(0, Rec_List[0].Rec_List.Length, i =>
                {
                    double[][] ETC = new double[8][];
                    for (int oct = 0; oct < 8; oct++) ETC[oct] = IR_Construction.ETCurve(null, null, Rec_List, Rec_List[0].CutOffTime, Rec_List[0].SampleRate, oct, i, SrcID, true);
                    STI_Values[i] = AcousticalMath.Speech_Transmission_Index(ETC, 1.2 * 343, NoiseSPL, Rec_List[0].SampleRate)[type];
                });
            }

            return STI_Values;
        }

        /// <summary>
        /// Creates an Early Decay Time Map
        /// </summary>
        /// <param name="EDT_Bounds"></param>
        /// <param name="H_OFFSET">H value offset (HSV Color system.)</param>
        /// <param name="H_BREADTH">H value breadth (HSV Color system.)</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_EDT_Map(PachMapReceiver[] Rec_List, int Octave, List<int> SrcID, bool Coherent_Superposition)
        {
            //T in ms.
            //Calculate C values...
            double[] EDT_Values = new double[Rec_List[0].Rec_List.Length];

            if (Coherent_Superposition)
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    int zero = 0;
                    double z = double.MaxValue;
                    foreach (int S_ID in SrcID)
                    {
                        double t = (Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time;
                        if (z > t) z = t;
                    }
                    zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);

                    foreach (int S_ID in SrcID)
                    {
                        double[] P;
                        Rec_List[S_ID].GetFilter(i, out P);
                        P = Audio.Pach_SP.Filter2Signal(P, Rec_List[S_ID].SWL, Rec_List[S_ID].SampleRate, 0);

                        int arrival = (int)Math.Floor((Rec_List[S_ID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                        for (int j = 0; j < P.Length; j++)
                        {
                            int t = j + arrival;
                            if (t >= hist.Length) break;

                            hist[t] += P[j];
                        }
                    }

                    if (Octave < 8) hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);

                    for (int t = 0; t < hist.Length; t++) hist[t] *= hist[t];

                    double[] SI = AcousticalMath.Schroeder_Integral(hist);
                    EDT_Values[i] = AcousticalMath.EarlyDecayTime(SI, Rec_List[0].SampleRate);
                }
            }
            else
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] SI = AcousticalMath.Schroeder_Integral(IR_Construction.ETCurve(null, null, Rec_List, Rec_List[0].CutOffTime, Rec_List[0].SampleRate, Octave, i, SrcID, true));
                    EDT_Values[i] = AcousticalMath.EarlyDecayTime(SI, Rec_List[0].SampleRate);
                }
            }

            return EDT_Values;
        }

        /// <summary>
        /// Creates a strength map
        /// </summary>
        /// <param name="G_Bounds"></param>
        /// <param name="Ref_Hist"></param>
        /// <param name="H_OFFSET">H value offset (HSV Color system.)</param>
        /// <param name="H_BREADTH">H value breadth (HSV Color system.)</param>
        /// <param name="S_OFFSET">S value offset (HSV Color system.)</param>
        /// <param name="S_BREADTH">S value breadth (HSV Color system.)</param>
        /// <param name="V_OFFSET">V value offset (HSV Color system.)</param>
        /// <param name="V_BREADTH">V value breadth (HSV Color system.)</param>
        /// <param name="Octave">The octave band to plot...</param>
        /// <returns>On Mesh with color assignments matching the input parameters and output variables.</returns>
        public static double[] Get_G_Map(PachMapReceiver[] Rec_List, int Octave, double SWL, int SrcID, bool Coherent_Superposition)
        {
            //T in ms.
            //Calculate C values...
            double[] G_Values = new double[Rec_List[0].Rec_List.Length];

            if (Coherent_Superposition)
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    double[] hist = new double[(int)Math.Ceiling(Rec_List[0].SampleCT * 44.1)];

                    int zero = 0;
                    double z = double.MaxValue;
                    //foreach (int S_ID in SrcID)
                    //{
                    double t = (Rec_List[SrcID].Rec_List[i] as Map_Receiver).Direct_Time;
                    if (z > t) z = t;
                    //}
                    zero = (int)Math.Floor(z * Rec_List[0].SampleRate * 44.1);
                    //}

                    //foreach (int S_ID in SrcID)
                    //{
                    double[] P;
                    Rec_List[SrcID].GetFilter(i, out P);
                    P = Audio.Pach_SP.Filter2Signal(P, Rec_List[SrcID].SWL, Rec_List[SrcID].SampleRate, 0);

                    int arrival = (int)Math.Floor((Rec_List[SrcID].Rec_List[i] as Map_Receiver).Direct_Time * Rec_List[0].SampleRate * 44.1) - zero;

                    for (int j = 0; j < P.Length; j++)
                    {
                        int ti = j + arrival;
                        if (ti >= hist.Length) break;

                        hist[ti] += P[j];
                    }
                    //}

                    if (Octave < 8) hist = Audio.Pach_SP.FIR_Bandpass(hist, Octave, 44100, 0);

                    for (int j = 0; j < hist.Length; j++) hist[j] *= hist[j];

                    G_Values[i] = AcousticalMath.Strength(hist, SWL, Coherent_Superposition);
                }
            }
            else
            {
                for (int i = 0; i < Rec_List[0].Rec_List.Length; i++)
                {
                    G_Values[i] = AcousticalMath.Strength(Rec_List[SrcID].Rec_List[i].GetEnergyHistogram(Octave), SWL, Coherent_Superposition);
                }
            }

            return G_Values;
        }

        /// <summary>
        /// The number of receivers in the map receiver.
        /// </summary>
        public override int Count
        {
            get
            {
                return Rec_List.Length;
            }
        }

        /// <summary>
        /// Checks receivers for intersections with a broadband ray (specular only). Records detections. Use after shooting a ray, but before adding the distance of ray travel to the total of distance traveled.
        /// </summary>
        /// <param name="R">a broadband ray.</param>
        /// <param name="Sumlength">the length of the ray before detection. Adds on the distance from the last reflection to the receiver.</param>
        /// <param name="EndPt">The point of the d</param>
        public override void CheckBroadbandRay(BroadRay R, Hare.Geometry.Point EndPt)
        {
            // R.ThreadID = this.Get_shotID();
            int X, Y, Z, Xend, Yend, Zend;
            //Identify whether the Origin is inside the voxel grid...
            //double Sumlength = R.t_sum;

            if (!OBox.IsPointInBox(R.origin))
            {
                double t0 = 0;
                if (!OBox.Intersect(R, ref t0, ref R.origin)) return;
                //Sumlength += t0;
            }

            //Identify where the ray enters the voxel grid...
            X = (int)Math.Floor((R.origin.x - OBox.Min_PT.x) / VoxelDims.x);
            Y = (int)Math.Floor((R.origin.y - OBox.Min_PT.y) / VoxelDims.y);
            Z = (int)Math.Floor((R.origin.z - OBox.Min_PT.z) / VoxelDims.z);

            Xend = (int)Math.Floor((EndPt.x - OBox.Min_PT.x) / VoxelDims.x);
            Yend = (int)Math.Floor((EndPt.y - OBox.Min_PT.y) / VoxelDims.y);
            Zend = (int)Math.Floor((EndPt.z - OBox.Min_PT.z) / VoxelDims.z);


            double tDeltaX, tDeltaY, tDeltaZ;
            double tMaxX = 0, tMaxY = 0, tMaxZ = 0;

            int stepX, stepY, stepZ, OutX, OutY, OutZ;

            if (X < 0) X = 0;
            if (X >= VoxelCtX) X = VoxelCtX - 1;
            if (Y < 0) Y = 0;
            if (Y >= VoxelCtY) Y = VoxelCtY - 1;
            if (Z < 0) Z = 0;
            if (Z >= VoxelCtZ) Z = VoxelCtZ - 1;

            ///
            if ((R.direction.x < 0 && X < Xend) || (R.direction.x >= 0 && X > Xend)) return;
            if ((R.direction.y < 0 && Y < Yend) || (R.direction.y >= 0 && Y > Yend)) return;
            if ((R.direction.z < 0 && Z < Zend) || (R.direction.z >= 0 && Z > Zend)) return;
            ///

            if (R.direction.x < 0)
            {
                OutX = -1;
                stepX = -1;
                tMaxX = (Voxels[X, Y, Z].Min_PT.x - R.origin.x) / R.direction.x;
                tDeltaX = VoxelDims.x / R.direction.x * stepX;
            }
            else
            {
                OutX = VoxelCtX;
                stepX = 1;
                tMaxX = (Voxels[X, Y, Z].Max_PT.x - R.origin.x) / R.direction.x;
                tDeltaX = VoxelDims.x / R.direction.x * stepX;
            }

            if (R.direction.y < 0)
            {
                OutY = -1;
                stepY = -1;
                tMaxY = (Voxels[X, Y, Z].Min_PT.y - R.origin.y) / R.direction.y;
                tDeltaY = VoxelDims.y / R.direction.y * stepY;
            }
            else
            {
                OutY = VoxelCtY;
                stepY = 1;
                tMaxY = (Voxels[X, Y, Z].Max_PT.y - R.origin.y) / R.direction.y;
                tDeltaY = VoxelDims.y / R.direction.y * stepY;
            }

            if (R.direction.z < 0)
            {
                OutZ = -1;
                stepZ = -1;
                tMaxZ = (Voxels[X, Y, Z].Min_PT.z - R.origin.z) / R.direction.z;
                tDeltaZ = VoxelDims.z / R.direction.z * stepZ;
            }
            else
            {
                OutZ = VoxelCtZ;
                stepZ = 1;
                tMaxZ = (Voxels[X, Y, Z].Max_PT.z - R.origin.z) / R.direction.z;
                tDeltaZ = VoxelDims.z / R.direction.z * stepZ;
            }

            do
            {
                CheckReceivers_Broadband(R, EndPt, X, Y, Z);

                if (tMaxX < tMaxY)
                {
                    if (tMaxX < tMaxZ)
                    {
                        X += stepX;
                        if (X < 0 || X >= VoxelCtX)
                            return; /* outside grid */
                        tMaxX = tMaxX + tDeltaX;
                    }
                    else
                    {
                        Z += stepZ;
                        if (Z < 0 || Z >= VoxelCtZ)
                            return; /* outside grid */
                        tMaxZ = tMaxZ + tDeltaZ;
                    }
                }

                else
                {
                    if (tMaxY < tMaxZ)
                    {
                        Y += stepY;
                        if (Y < 0 || Y >= VoxelCtY)
                            return; /* outside grid */
                        tMaxY = tMaxY + tDeltaY;
                    }
                    else
                    {
                        Z += stepZ;
                        if (Z < 0 || Z >= VoxelCtZ)
                            return; /* outside grid */
                        tMaxZ = tMaxZ + tDeltaZ;
                    }
                }
            } while (X != Xend && Y != Yend && Z != Zend);
        }

        /// <summary>
        /// Checks receivers for intersections with a single octave ray (specular only). Records detections. Use after shooting a ray, but before adding the distance of ray travel to the total of distance traveled.
        /// </summary>
        /// <param name="R">a broadband ray.</param>
        /// <param name="Sumlength">the length of the ray before detection. Adds on the distance from the last reflection to the receiver.</param>
        /// <param name="EndPt">The point of the d</param>
        public override void CheckRay(OctaveRay R, Hare.Geometry.Point EndPt)
        {
            //R.ThreadID = this.Get_shotID();
            int X, Y, Z, Xend, Yend, Zend;
            //Identify whether the Origin is inside the voxel grid...
            //double Sumlength = R.t_sum;

            if (!OBox.IsPointInBox(R.origin))
            {
                double t0 = 0;
                if (!OBox.Intersect(R, ref t0, ref R.origin)) return;
                //Sumlength += t0;
            }

            //Identify which voxel the Origin point is located in...

            //Identify where the ray enters the voxel grid...
            X = (int)Math.Floor((R.origin.x - OBox.Min_PT.x) / VoxelDims.x);
            Y = (int)Math.Floor((R.origin.y - OBox.Min_PT.y) / VoxelDims.y);
            Z = (int)Math.Floor((R.origin.z - OBox.Min_PT.z) / VoxelDims.z);

            Xend = (int)Math.Floor((EndPt.x - OBox.Min_PT.x) / VoxelDims.x);
            Yend = (int)Math.Floor((EndPt.y - OBox.Min_PT.y) / VoxelDims.y);
            Zend = (int)Math.Floor((EndPt.z - OBox.Min_PT.z) / VoxelDims.z);

            double tDeltaX, tDeltaY, tDeltaZ;
            double tMaxX = 0, tMaxY = 0, tMaxZ = 0;

            int stepX, stepY, stepZ, OutX, OutY, OutZ;

            if (X < 0) X = 0;
            if (X >= VoxelCtX) X = VoxelCtX - 1;
            if (Y < 0) Y = 0;
            if (Y >= VoxelCtY) Y = VoxelCtY - 1;
            if (Z < 0) Z = 0;
            if (Z >= VoxelCtZ) Z = VoxelCtZ - 1;

            ///
            if ((R.direction.x < 0 && X < Xend) || (R.direction.x >= 0 && X > Xend)) return;
            if ((R.direction.y < 0 && Y < Yend) || (R.direction.y >= 0 && Y > Yend)) return;
            if ((R.direction.z < 0 && Z < Zend) || (R.direction.z >= 0 && Z > Zend)) return;
            ///

            if (R.direction.x < 0)
            {
                OutX = -1;
                stepX = -1;
                tMaxX = (Voxels[X, Y, Z].Min_PT.x - R.origin.x) / R.direction.x;
                tDeltaX = VoxelDims.x / R.direction.x * stepX;
            }
            else
            {
                OutX = VoxelCtX;
                stepX = 1;
                tMaxX = (Voxels[X, Y, Z].Max_PT.x - R.origin.x) / R.direction.x;
                tDeltaX = VoxelDims.x / R.direction.x * stepX;
            }

            if (R.direction.y < 0)
            {
                OutY = -1;
                stepY = -1;
                tMaxY = (Voxels[X, Y, Z].Min_PT.y - R.origin.y) / R.direction.y;
                tDeltaY = VoxelDims.y / R.direction.y * stepY;
            }
            else
            {
                OutY = VoxelCtY;
                stepY = 1;
                tMaxY = (Voxels[X, Y, Z].Max_PT.y - R.origin.y) / R.direction.y;
                tDeltaY = VoxelDims.y / R.direction.y * stepY;
            }

            if (R.direction.z < 0)
            {
                OutZ = -1;
                stepZ = -1;
                tMaxZ = (Voxels[X, Y, Z].Min_PT.z - R.origin.z) / R.direction.z;
                tDeltaZ = VoxelDims.z / R.direction.z * stepZ;
            }
            else
            {
                OutZ = VoxelCtZ;
                stepZ = 1;
                tMaxZ = (Voxels[X, Y, Z].Max_PT.z - R.origin.z) / R.direction.z;
                tDeltaZ = VoxelDims.z / R.direction.z * stepZ;
            }

            do
            {
                //CheckReceivers_Octave(R, Sumlength, EndPt, X, Y, Z);
                CheckReceivers_Octave(R, EndPt, X, Y, Z);

                if (tMaxX < tMaxY)
                {
                    if (tMaxX < tMaxZ)
                    {
                        X += stepX;
                        if (X < 0 || X >= VoxelCtX)
                            return; /* outside grid */
                        tMaxX = tMaxX + tDeltaX;
                    }
                    else
                    {
                        Z += stepZ;
                        if (Z < 0 || Z >= VoxelCtZ)
                            return; /* outside grid */
                        tMaxZ = tMaxZ + tDeltaZ;
                    }
                }

                else
                {
                    if (tMaxY < tMaxZ)
                    {
                        Y += stepY;
                        if (Y < 0 || Y >= VoxelCtY)
                            return; /* outside grid */
                        tMaxY = tMaxY + tDeltaY;
                    }
                    else
                    {
                        Z += stepZ;
                        if (Z < 0 || Z >= VoxelCtZ)
                            return; /* outside grid */
                        tMaxZ = tMaxZ + tDeltaZ;
                    }
                }
            //} while (X != Xend && Y != Yend && Z != Zend);
            } while (X != Xend || Y != Yend || Z != Zend);
        }

        /// <summary>
        /// Directionality reporter.
        /// </summary>
        /// <param name="Octave">the octave band</param>
        /// <param name="q">the </param>
        /// <param name="Rec_Index"></param>
        /// <returns></returns>
        public override Vector Directions_Pos(int Octave, int q, int Rec_Index)
        {
            return Rec_List[Rec_Index].Recs.Directions_Pos(Octave, q);
        }

        /// <summary>
        /// Directionality reporter.
        /// </summary>
        /// <param name="Octave">the octave band</param>
        /// <param name="q">the </param>
        /// <param name="Rec_Index"></param>
        /// <returns></returns>
        public override Vector Directions_Neg(int Octave, int q, int Rec_Index)
        {
            return Rec_List[Rec_Index].Recs.Directions_Neg(Octave, q);
        }

        /// <summary>
        /// Get the mesh on which the map is based.
        /// </summary>
        /// <returns>The Mesh object.</returns>
        public Topology MapMesh()
        {
            return Map_Mesh;
        }

        /// <summary>
        /// Used after a multithreaded calculation to combine thread-local results.
        /// </summary>
        /// <param name="Rec_in">Thread-local results to be combined with this receiver.</param>
        public override void Combine_Clones(Receiver_Bank Rec_in)
        {
            PachMapReceiver Rec = (PachMapReceiver)Rec_in;
            //if (Rec.Rec_List[0].Recs.Directions(8,0).z = -1)
            for (int R = 0; R < Rec.Count; R++)
            {
                for (int oct = 0; oct < 8; oct++)
                {
                    for (int T = 0; T < Rec.Duration(); T++)
                    {
                        Rec_List[R].Combine_Sample(T, Rec.Rec_List[R].Energy(T, oct), Rec.Rec_List[R].Directions_Pos(oct, T), Rec.Rec_List[R].Directions_Neg(oct, T), oct);
                    }
                }
            }
        }

        /// <summary>
        /// Returns the kind of data stored.
        /// </summary>
        /// <returns>string indicating what kind of data</returns>
        public string Data_Type()
        {
            return Rec_List[0].Recs.Type();
        }

        /// <summary>
        /// The direct sound outputs a power mod factor. This takes that function's place, since there is no direct sound for this simulation type.
        /// </summary>
        /// <param name="new_SWL"></param>
        /// <returns></returns>
        public double[] PowerModFactor(double[] new_SWL)
        {
            double[] factor = new double[8];
            for (int i = 0; i < 8; i++) factor[i] = Utilities.AcousticalMath.Intensity_SPL(new_SWL[i]) / Utilities.AcousticalMath.Intensity_SPL(this.SWL[i]);
            return factor;
        }

        /// <summary>
        /// Specialized spherical receiver.
        /// </summary>
        [Serializable]
        public class Map_Receiver : Spherical_Receiver
        {
            public double Direct_Time;
            /// <summary>
            /// Use this constructor only for reading in data.
            /// </summary>
            /// <param name="Rec_Origin"></param>
            /// <param name="SourceCT"></param>
            /// <param name="Index"></param>
            /// <param name="SampleRate"></param>
            /// <param name="SampleCT"></param>
            /// <param name="Directional"></param>
            public Map_Receiver(Point Rec_Origin, Point Src_Origin, double Direct_Arrival_Time, double rho_c, int Index, int SampleRate, int SampleCT, bool Directional)
            {
                Rho_C = rho_c;
                Direct_Time = Direct_Arrival_Time;
                base.Origin = new Point(Rec_Origin.x, Rec_Origin.y, Rec_Origin.z);
                this.CO_Time = (double)SampleCT / (double)SampleRate;
                this.SampleRate = SampleRate;
                if (Directional)
                {
                    base.Recs = new Directional_Histogram(SampleRate, SampleCT);
                }
                else
                {
                    base.Recs = new Histogram(SampleRate, SampleCT);
                }
            }

            public Map_Receiver(double Diameter, int i, Hare.Geometry.Point Point, Source Src, double SoundSpeed, double rho, double[] Attenuation, int SampleRate_in, double COTime_in, int ProcessorSpec, bool Time1Pt, bool Directional)
            {
                Ray_ID = new int[ProcessorSpec];
                Radius = 1;
                Radius2 = Radius * Radius;
                CO_Time = COTime_in;
                SampleRate = SampleRate_in;
                Origin = Point;
                C_Sound = SoundSpeed;
                Atten = Attenuation;
                SizeMod = 1 / Math.PI;
                Origin = Point;
                Radius = Diameter * 0.5;
                Radius2 = Radius * Radius;
                this.Rho_C = rho * Sound_Speed;
                this.SizeMod = 1 / (Math.PI * Radius2);
                Point L1 = Src.Origin() - Origin;
                Direct_Time = Math.Sqrt(L1.x * L1.x + L1.y * L1.y + L1.z * L1.z) / C_Sound;

                if (!Time1Pt)
                {
                    if (!Directional)
                    {
                        base.Recs = new Histogram(SampleRate_in, COTime_in);
                    }
                    else
                    {
                        base.Recs = new Directional_Histogram(SampleRate_in, COTime_in);
                    }
                }
                else
                {
                    base.Recs = new Histogram_1Pt(SampleRate_in, COTime_in);
                }
            }

            public override double Energy(int Time, int oct)
            {
                return this.Recs.Energy[oct][Time];
            }

            /// <summary>
            /// This method checks a receiver for a ray using a Broadband Ray (these typically occur before rays are split in the Raytracing simulation.
            /// </summary>
            /// <param name="r">the ray.</param>
            /// <param name="endPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            public override void CheckBroadbandRay(BroadRay r, Hare.Geometry.Point endPt)
            {
                Vector m = r.origin - Origin;
                double b = Hare_math.Dot(m, r.direction);
                double c = Hare_math.Dot(m, m) - Radius2;
                if (c > 0 && b > 0) return;
                double discr = b * b - c;
                if (discr < 0) return;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);
                double t = (t1 + t2) / 2;
                if (t > 0 && t * t < SqDistance(endPt, r.origin))
                {
                    double raydist = t / C_Sound + r.t_sum - Direct_Time;
                    Vector dir = r.direction * -1;
                    Recs.Add(raydist, r.Energy[0] * Math.Pow(10, -.1 * Atten[0] * raydist) * SizeMod, dir, Rho_C, 0);
                    Recs.Add(raydist, r.Energy[1] * Math.Pow(10, -.1 * Atten[1] * raydist) * SizeMod, dir, Rho_C, 1);
                    Recs.Add(raydist, r.Energy[2] * Math.Pow(10, -.1 * Atten[2] * raydist) * SizeMod, dir, Rho_C, 2);
                    Recs.Add(raydist, r.Energy[3] * Math.Pow(10, -.1 * Atten[3] * raydist) * SizeMod, dir, Rho_C, 3);
                    Recs.Add(raydist, r.Energy[4] * Math.Pow(10, -.1 * Atten[4] * raydist) * SizeMod, dir, Rho_C, 4);
                    Recs.Add(raydist, r.Energy[5] * Math.Pow(10, -.1 * Atten[5] * raydist) * SizeMod, dir, Rho_C, 5);
                    Recs.Add(raydist, r.Energy[6] * Math.Pow(10, -.1 * Atten[6] * raydist) * SizeMod, dir, Rho_C, 6);
                    Recs.Add(raydist, r.Energy[7] * Math.Pow(10, -.1 * Atten[7] * raydist) * SizeMod, dir, Rho_C, 7);
                }
            }

            /// <summary>
            /// Checks receiver for an intersection with a ray 
            /// </summary>
            /// <param name="R">the ray.</param>
            /// <param name="endPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            public override void CheckRay(OctaveRay R, Hare.Geometry.Point endPt)
            {
                Vector m = R.origin - Origin;
                double b = Hare_math.Dot(m, R.direction);
                double c = Hare_math.Dot(m, m) - Radius2;
                if (c > 0 && b > 0) return;
                double discr = b * b - c;
                if (discr < 0) return;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);
                double t = (t1 + t2) / 2;
                if (t > 0 && t * t < SqDistance(endPt, R.origin))
                {
                    double raydist = t / C_Sound + R.t_sum;
                    Recs.Add(raydist - Direct_Time, R.Intensity * Math.Pow(10, -.1 * Atten[R.Octave] * raydist) * SizeMod, R.direction * -1, Rho_C, R.Octave);// / (1.33333333333333 * Math.PI * Min_Radius2 * Min_Radius), EndPt - R.origin, R.Octave);//R.Scat_Mod * 
                }
            }
        
            //public void Create_Pressure(double[] SWL, int ThreadID)
            //{
            //    Recs.P = Audio.Pach_SP.Filter_Interpolation(SWL, Recs.Energy, 1000, 44100, this.Rho_C);
            //    Recs.P = Audio.Pach_SP.Filter2Signal(Recs.P, SWL, SampleRate, ThreadID);
            //}
        }
    }
}