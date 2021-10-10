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

using System;
using System.Collections.Generic;
using Hare.Geometry;
using System.Linq;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        /// <summary>
        /// Scene base class. Do not use.
        /// </summary>C:\Users\User\Desktop\Pachyderm_Acoustic_2.0\Pachyderm_Acoustic\Classes_Scenes.cs
        [Serializable]
        public abstract class Scene
        {
            protected List<Material> AbsorptionData = new List<Material>();
            protected List<Scattering> ScatteringData = new List<Scattering>();
            protected List<double[]> TransmissionData = new List<double[]>();
            protected List<bool> Transmissive = new List<bool>();
            protected double TempC_S;
            protected double Pa_S;
            protected double hr_S;
            protected bool EdgeFC;
            protected int AC_S;
            protected double[] Area;
            protected Medium_Properties Env_Prop;
            protected double SplitRatio = 0.25;
            public System.Windows.Forms.DialogResult Status = System.Windows.Forms.DialogResult.OK;
            public Random R_Seed;
            public bool Valid;
            public bool Custom_Method;
            public bool Complete = true;
            public bool Partitioned = false;
            public List<Edge> Edge_Nodes = new List<Edge>();
            public bool IsHomogeneous = true;
            public bool hasnulllayers = false;
            //Raw edges
            protected List<MathNet.Numerics.Interpolation.CubicSpline[]> Edges;
            protected List<MathNet.Numerics.Interpolation.CubicSpline[]> Edge_Tangents;
            protected List<bool> Edge_isSoft;
            protected List<double> EdgeLength;
            #region Inheritables

            /// <summary>
            /// 
            /// </summary>
            /// <param name="Temp">temperature in Celsius</param>
            /// <param name="hr">relative humidity in percent</param>
            /// <param name="Pa">Pressure in Pascals</param>
            /// <param name="Air_Choice"></param>
            /// <param name="EdgeCorrection"></param>
            /// <param name="IsAcoustic"></param>
            public Scene(double Temp, double hr, double Pa, int Air_Choice, bool EdgeCorrection, bool IsAcoustic)
            {
                Custom_Method = IsAcoustic;
                Valid = false;
                TempC_S = Temp;
                Pa_S = Pa;          //TODO: Check that the pressure is in the correct units... hpa, kpa, pa.
                hr_S = hr;
                AC_S = Air_Choice;
                EdgeFC = EdgeCorrection;
                R_Seed = new Random();

                double TK = Temp + 273.15;
                //convert to Kelvins
                //Psat = Pr * 10 ^ (-6.8346 * (To1 / T) ^ 1.261 + 4.6151)
                double Psat = 101.325 * Math.Pow(10, (-6.8346 * Math.Pow((273.16 / TK), 1.261) + 4.6151));
                //6.1121 * Math.Exp((18.678 - T / 234.5) * TC / (257.14 + T))
                double h = hr * (Psat / Pa);

                Env_Prop = new Uniform_Medium(Air_Choice, Pa, TK, hr, EdgeCorrection);

                // Saturation water vapor pressure:
                // psv = exp(aa(1,1)*T^2 + aa(2,1)*T + aa(3,1) + aa(4,1)/T); % Formula Giacomo
                double psv = 101325 * Math.Pow(10, (4.6151 - 6.8346 * Math.Pow(((273.15 + 0.01) / (Temp + 273.15)), 1.261)));     // Formula ISO 9613-1:1993

                // Enhancement factor:
                double fpst = 1.00062 + 3.14e-8 * Atmospheric_Pressure + 5.6e-7 * Temp * Temp;

                // Mole fraction of water vapor in air:
                double xw = Relative_Humidity * psv * fpst / (100 * Atmospheric_Pressure);

                // Compressibility factor:
                double Z = 1 - Atmospheric_Pressure / (Temp + 273.15) * (1.58123e-6 + -2.9331e-8 * +1.1043e-10 * Temp * Temp + (5.707e-6 + -2.051e-8) * xw + (1.9898e-4 + -2.376e-6) * xw * xw) + Math.Pow((Atmospheric_Pressure / (Temp + 273.15)), 2) * (1.83e-11 + -0.765e-8 * xw * xw);

                //// Density of air:
                //rho = 3.48349 * 1e-3 * Atmospheric_Pressure / (Z * (Temp + 273.15)) * (1 - 0.3780 * xw);
            }

            public int EdgeCount
            {
                get
                {
                    return this.Edge_Nodes.Count;
                }
            }

            /// <summary>
            /// cast a ray within the model.
            /// </summary>
            /// <param name="R">A ray, complete with origin point and direction...</param>
            /// <param name="u">optional surface coordinate</param>
            /// <param name="v">optional surface coordinate</param>
            /// <param name="Poly_ID">the polygon of intersection</param>
            /// <param name="X_PT">the point of intersection</param>
            /// <param name="t">the distance traveled by the ray</param>
            /// <returns>true if successful, false if no hit</returns>
            public abstract bool shoot(Hare.Geometry.Ray R, out double u, out double v, out int Poly_ID, out Hare.Geometry.Point X_PT, out double t);
            /// <summary>
            /// cast a ray within the model.
            /// </summary>
            /// <param name="R">A ray, complete with origin point and direction...</param>
            /// <param name="u">optional surface coordinate</param>
            /// <param name="v">optional surface coordinate</param>
            /// <param name="Poly_ID">the polygon of intersection</param>
            /// <param name="X_PT">the point of intersection</param>
            /// <param name="t">the distance traveled by the ray</param>
            /// <param name="code">returns code indicating where the ray has been (for polyrefractive)</param>
            /// <returns>true if successful, false if no hit</returns>
            public abstract bool shoot(Ray R, out double u, out double v, out int Poly_ID, out List<Hare.Geometry.Point> X_PT, out List<double> t, out List<int> code);
            /// <summary>
            /// cast a ray within the model.
            /// </summary>
            /// <param name="R">A ray, complete with origin point and direction...</param>
            /// <param name="u">optional surface coordinate</param>
            /// <param name="v">optional surface coordinate</param>
            /// <param name="Poly_ID">the polygon of intersection</param>
            /// <param name="X_PT">the point of intersection</param>
            /// <param name="t">the distance traveled by the ray</param>
            /// <param name="code">returns code indicating where the ray has been (for polyrefractive)</param>
            /// <returns>true if successful, false if no hit</returns>
            public abstract bool shoot(Ray R, int topo, out X_Event Xpt);
            /// <summary>
            /// The local normal of a surface.
            /// </summary>
            /// <param name="i">surface index</param>
            /// <param name="u">u coordinate for NURBS tracing</param>
            /// <param name="v">v coordinate for NURBS tracing</param>
            /// <returns></returns>
            public abstract Hare.Geometry.Vector Normal(int i, double u, double v);
            /// <summary>
            /// The local normal of a surface.
            /// </summary>
            /// <param name="i">surface index</param>
            /// <returns></returns>
            public abstract Hare.Geometry.Vector Normal(int i);
            ///// <summary>
            ///// The local normal of a surface.
            ///// </summary>
            ///// <param name="u">u coordinate for NURBS tracing</param>
            ///// <param name="v">v coordinate for NURBS tracing</param>
            ///// <param name="i">surface index</param>
            ///// <returns></returns>
            //public abstract Vector3d Normal(double u, double v, int i);
            /// <summary>
            /// Returns the surface area of a surface in the model.
            /// </summary>
            /// <param name="x">the index of the surface.</param>
            /// <returns></returns>
            public abstract double SurfaceArea(int x);
            /// <summary>
            /// Returns the number of surface objects in the model.
            /// </summary>
            /// <returns></returns>
            public abstract int Count();
            /// <summary>
            /// Optimizes the model using the chosen spatial partition system.
            /// </summary>
            public abstract void partition();
            /// <summary>
            /// Optimizes the model using the chosen spatial partition system.
            /// </summary>
            /// <param name="SP_Param">a parameter describing some aspect of the partition</param>
            public abstract void partition(int SP_Param);
            ///// <summary>
            ///// Optimizes the model using the chosen spatial partition system.
            ///// </summary>
            ///// <param name="P">points to add to the spatial partition system.</param>
            ///// <param name="SP_PARAM">a parameter describing some aspect of the partition</param>
            //public abstract void partition(List<Point3d> P, int SP_PARAM);
            /// <summary>
            ///a parameter describing some aspect of the partition
            /// </summary>
            /// <param name="P">points to add to the spatial partition system.</param>
            /// <param name="SP_PARAM">a parameter describing some aspect of the partition</param>
            public abstract void partition(List<Hare.Geometry.Point> P, int SP_PARAM);
            ///// <summary>
            /////a parameter describing some aspect of the partition
            ///// </summary>
            ///// <param name="P">points to add to the spatial partition system.</param>
            ///// <param name="SP_PARAM">a parameter describing some aspect of the partition</param>
            //public abstract void partition(Point3d[] P, int SP_PARAM);
            /// <summary>
            ///a parameter describing some aspect of the partition
            /// </summary>
            /// <param name="P">points to add to the spatial partition system.</param>
            /// <param name="SP_PARAM">a parameter describing some aspect of the partition</param>
            public abstract void partition(Hare.Geometry.Point[] P, int SP_PARAM);
            /// <summary>
            /// Checks whether the specified surface is planar.
            /// </summary>
            /// <param name="i">the index of the surface.</param>
            /// <returns>true if planar, false if not.</returns>
            public abstract bool IsPlanar(int i);
            /// <summary>
            /// 
            /// </summary>
            /// <param name="P"></param>
            /// <param name="D"></param>
            /// <returns></returns>
            public abstract Hare.Geometry.Point ClosestPt(Hare.Geometry.Point P, ref double D);
            /// <summary>
            /// returns all points in the model.
            /// </summary>
            /// <param name="PTS"></param>
            /// <returns></returns>
            public abstract bool PointsInScene(List<Hare.Geometry.Point> PTS);
            public abstract string Scene_Type();
            public abstract void EdgeFrame_Tangents(Hare.Geometry.Point Origin, Hare.Geometry.Vector Normal, int[] SrfIDs, ref List<double> dist2, List<Vector> Dir, List<int> IDs);
            public abstract void Register_Edges(IEnumerable<Hare.Geometry.Point> S, IEnumerable<Hare.Geometry.Point> R);
            
            public void Register_Edges(IEnumerable<Source> S, Receiver_Bank R)
            {
                List<Hare.Geometry.Point> HS = new List<Hare.Geometry.Point>();
                List<Hare.Geometry.Point> HR = new List<Hare.Geometry.Point>();

                foreach (Source SPT in S) HS.Add(SPT.Origin());
                foreach (Point RPT in R.Origins()) HR.Add(RPT);

                Register_Edges(HS, HR);
            }
            #endregion

            public abstract double Sound_speed(int arg);
            public abstract double Sound_speed(Hare.Geometry.Point pt);
            
            public bool Edge_Frequency
            {
                get
                {
                    return EdgeFC;
                }
            }

            public double Temperature
            {
                get
                {
                    return TempC_S;
                }
            }

            public double Relative_Humidity
            {
                get
                {
                    return hr_S;
                }
            }

            public double Atmospheric_Pressure
            {
                get
                {
                    return Pa_S;
                }
            }

            public int Attenuation_Method
            {
                get
                {
                    return AC_S;
                }
            }

            public void AttenuationFilter(int no_of_elements, int sampleFrequency, double dist, ref double[] Freq, ref double[] Atten, Hare.Geometry.Point pt)
            {
                Env_Prop.AttenuationFilter(no_of_elements, sampleFrequency, dist, ref Freq, ref Atten, pt);
            }

            public double AttenuationPureTone(Hare.Geometry.Point pt, double frequency)
            {
                return this.Env_Prop.AttenuationPureTone(pt, frequency);
            }

            public double[] Attenuation(int arg)
            {
                return this.Env_Prop.Attenuation_Coef(arg);
            }

            public double[] Attenuation(Hare.Geometry.Point pt)
            {
                return this.Env_Prop.Attenuation_Coef(pt);
            }

            public double Rho(int arg)
            {
                return this.Env_Prop.Rho(arg);
            }

            public double Rho(Hare.Geometry.Point pt)
            {
                return this.Env_Prop.Rho(pt);
            }

            public double Rho_C(int arg)
            {
                return this.Env_Prop.Rho_C(arg);
            }

            public double Rho_C(Hare.Geometry.Point pt)
            {
                return this.Env_Prop.Rho_C(pt);
            }

            public Material Surface_Material(int id)
            {
                return AbsorptionData[id];
            }

            public abstract void Absorb(ref BroadRay Ray, out double cos_theta, double u, double v);
            public abstract void Absorb(ref OctaveRay Ray, out double cos_theta, double u, double v);
            public abstract void Scatter_Early(ref BroadRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, double cos_theta, double u, double v);
            public abstract void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, double cos_theta, double u, double v);
            public abstract void Scatter_Simple(ref OctaveRay Ray, ref Random rand, double cos_theta, double u, double v);
            
            public List<Material> AbsorptionValue
            {
                get
                {
                    return AbsorptionData;
                }
            }

            public List<Scattering> ScatteringValue
            {
                get
                {
                    return ScatteringData;
                }
            }

            public List<double[]> TransmissionValue
            {
                get
                {
                    return TransmissionData;
                }
            }

            public List<bool> IsTransmissive
            {
                get
                {
                    return Transmissive;
                }
            }

            /// <summary>
            /// Calculates the direction of a specular reflection.
            /// </summary>
            /// <param name="RayDirect"></param>
            /// <param name="u"></param>
            /// <param name="v"></param>
            /// <param name="x"></param>
            public virtual void ReflectRay(ref OctaveRay RayDirect, ref double u, ref double v, ref int x)
            {
                Vector local_N = Normal(x, u, v);
                RayDirect.direction -= local_N * Hare_math.Dot(RayDirect.direction, local_N) * 2;
            }

            public virtual void ReflectRay(ref BroadRay RayDirect, ref double u, ref double v, ref int x)
            {
                Vector local_N = Normal(x, u, v);
                RayDirect.direction -= local_N * Hare_math.Dot(RayDirect.direction, local_N) * 2;
            }

            public abstract Hare.Geometry.Point Min();
            public abstract Hare.Geometry.Point Max();
        }

        [Serializable]
        public class Polygon_Scene : Scene
        {
            protected Hare.Geometry.Topology[] Topo;
            private Hare.Geometry.Spatial_Partition SP;
            //Data kept to identify polygons as being part of larger planar entities...//
            //private List<Vector[]> PDir = new List<Vector[]>();
            //private List<double[]> Kurvatures = new List<double[]>();
            /// <summary>
            /// List of polygons with plane ids they are attributed.
            /// </summary>
            //private List<int> Plane_ID = new List<int>();
            /// <summary>
            /// List of planes with list of polygons attributed to each plane.
            /// </summary>
            //private List<List<int>> Plane_Members = new List<List<int>>();
            public double[] Plane_Area;
            private double[][] PolyPlaneFract;
            private int[][] Object_Members;
            private List<int> Object_ID;
            private double[][][] kurvature;
            private Vector[][][] frame_axes;
            public Vector[][][] Edge_Frames;
            public double[][][] Edge_Kurvatures;
            public Vector[][] Edge_Normals;
            private Hare.Geometry.Edge[][] Object_Edges;
            private bool[] iscurved;
            ////////////////////////////////////////////////////////////////////

            public Polygon_Scene(double Temp, double hr, double Pa, int Air_Choice, bool EdgeCorrection, bool IsAcoustic)
                : base(Temp, hr, Pa, Air_Choice, EdgeCorrection, IsAcoustic)
            {
            }

            public void Construct(Point[][][] Model, List<Material> Mat, List<Scattering> Scat, List<double[]> Trans, bool[] IsCurved = null, double[][][] Kurvatures = null, Vector[][][] Frame_Axes = null)
            {
                if (Mat.Count != Model.Length && Scat.Count != Model.Length) throw new Exception("The number of material codes must match the numer of objects in the model...");
                if (IsCurved == null || Kurvatures == null || Frame_Axes == null)
                {
                    iscurved = new bool[Model.Length];
                    Kurvatures = new double[Model.Length][][];
                    Frame_Axes = new Vector[Model.Length][][];
                    for (int i = 0; i < Model.Length; i++)
                    {
                        Kurvatures[i] = new double[Model[i].Length][];
                        Frame_Axes[i] = new Vector[Model[i].Length][];
                        for (int j = 0; j < Model[i].Length; j++)
                        {
                            Kurvatures[i][j] = new double[2];
                            Frame_Axes[i][j] = new Vector[2] { new Vector(), new Vector() };
                        }
                    }
                }

                kurvature = Kurvatures;
                frame_axes = Frame_Axes;
                iscurved = IsCurved;

                List<Point[]> PTS = new List<Point[]>();

                AbsorptionData = new List<Material>();
                ScatteringData = new List<Scattering>();
                TransmissionData = new List<double[]>();

                ///Get Bounds of Model
                double MinX = double.PositiveInfinity, MinY = double.PositiveInfinity, MinZ = double.PositiveInfinity;
                double MaxX = double.NegativeInfinity, MaxY = double.NegativeInfinity, MaxZ = double.NegativeInfinity;
                for (int i = 0; i < Model.Length; i++)
                {
                    for (int j = 0; j < Model[i].Length; j++)
                    {
                        PTS.Add(Model[i][j]);
                        AbsorptionData.Add(Mat[i]);
                        ScatteringData.Add(Scat[i]);
                        TransmissionData.Add(Trans[i]);
                        for (int k = 0; k < Model[i][j].Length; k++)
                        {
                            if (Model[i][j][k].x < MinX) MinX = Model[i][j][k].x;
                            if (Model[i][j][k].y < MinX) MinY = Model[i][j][k].y;
                            if (Model[i][j][k].z < MinX) MinZ = Model[i][j][k].z;
                            if (Model[i][j][k].x > MaxX) MaxX = Model[i][j][k].x;
                            if (Model[i][j][k].y > MaxY) MaxY = Model[i][j][k].y;
                            if (Model[i][j][k].z > MaxZ) MaxZ = Model[i][j][k].z;
                        }
                    }
                }
                ////////////////////////////////////////
                Topo = new Hare.Geometry.Topology[1];
                Topo[0] = new Topology(PTS.ToArray(), 2);// Utilities.Pach_Tools.RPttoHPt(Box.Min), Utilities.PachTools.RPttoHPt(Box.Max));
                ////////////////////////////////////////
                
                int done = 0;
                Object_ID = new List<int>();
                Object_Members = new int[Model.Length][];
                Object_Edges = new Hare.Geometry.Edge[Model.Length][];
                Edge_Frames = new Vector[Model.Length][][];
                Edge_Kurvatures = new double[Model.Length][][];
                Edge_Normals = new Vector[Model.Length][];
                for(int i = 0; i < Model.Length; i++)
                {
                    Object_Members[i] = new int[Model[i].Length];
                    Dictionary<int, Hare.Geometry.Edge> O_Edges = new Dictionary<int, Hare.Geometry.Edge>();
                    Dictionary<int, Vector[]> E_Frames = new Dictionary<int, Vector[]>();
                    Dictionary<int, List<Vector>> E_Norm = new Dictionary<int, List<Vector>>();
                    Dictionary<int, double[]> E_Kurvatures = new Dictionary<int, double[]>();
                    Dictionary<int, int> e_ct = new Dictionary<int, int>();
                    for (int j = 0; j < Model[i].Length; j++)
                    {
                        Object_Members[i][j] = done + j;
                        Object_ID.Add(i);

                        foreach (Hare.Geometry.Edge e in Topo[0].Polys[done+j].Edges)
                        {
                            int hash = e.GetHashCode();
                            if (!O_Edges.ContainsKey(e.GetHashCode()))
                            {
                                O_Edges.Add(hash, e);
                                E_Frames.Add(hash, Frame_Axes[i][j]);
                                E_Norm.Add(hash, new List<Vector> { Topo[0].Polys[i].Normal });
                                E_Kurvatures.Add(hash, Kurvatures[i][j].Clone() as double[]);
                                e_ct.Add(hash, 1);
                            }
                            else
                            {
                                E_Frames[hash][0] += Frame_Axes[i][j][0];
                                E_Frames[hash][1] += Frame_Axes[i][j][1];
                                E_Kurvatures[hash][0] += Kurvatures[i][j][0];
                                E_Kurvatures[hash][1] += Kurvatures[i][j][1];
                                E_Norm[hash].Add(Topo[0].Polys[i].Normal);
                                e_ct[hash]++;
                            }
                        }
                    }

                    Dictionary<int, Vector> e_norm = new Dictionary<int, Vector>();
                    int dir = 1;
                    foreach (int h in E_Frames.Keys)
                    {
                        //Hare.Geometry.Vector d = E_Norm[h][0] + E_Norm[h][1];
                        //d.Normalize();
                        if (E_Norm[h].Count > 1)
                        {
                            Hare.Geometry.Vector n = new Vector();
                            for (int j = 0; j < O_Edges[h].Polys.Count; j++) n += O_Edges[h].Polys[j].Normal;//E_Norm[h][j];
                            n.Normalize();

                            dir = (Hare.Geometry.Hare_math.Dot(E_Frames[h][0], n) < 0) ? 1 : -1;
                            break;
                        }
                    }

                    List<int> omissions = new List<int>();

                    foreach (int h in E_Frames.Keys)
                    {
                        Hare.Geometry.Vector n = new Vector();
                        Hare.Geometry.Vector abis = new Vector();

                        if (O_Edges[h].Polys.Count < 2 || Math.Abs(Hare.Geometry.Hare_math.Dot(O_Edges[h].Tangents[0], O_Edges[h].Tangents[1])) == 1) { omissions.Add(h); continue; } //IF not curved across edfe, omit this one...

                        for (int j = 0; j < O_Edges[h].Polys.Count; j++) { n += O_Edges[h].Polys[j].Normal; abis += O_Edges[h].Tangents[j]; }//E_Norm[h][j];
                        n.Normalize();
                        abis.Normalize();
                        
                        E_Frames[h][0] = E_Frames[h][0] / E_Frames[h][0].Length();
                        E_Frames[h][1] = E_Frames[h][1] / E_Frames[h][1].Length();

                        double dot0 = Math.Abs(Hare_math.Dot(O_Edges[h].Tangents[0], E_Frames[h][0]));
                        double dot1 = Math.Abs(Hare_math.Dot(O_Edges[h].Tangents[0], E_Frames[h][1]));

                        int EFC = (dot0 < dot1) ? 1 : 0;

                        dir = (Hare.Geometry.Hare_math.Dot(abis, n) < 0) ? 1 : -1;
                        e_norm.Add(h, n * dir);

                        if (E_Kurvatures[h][EFC] < 0)
                        {
                            E_Kurvatures[h][0] *= -1;
                            E_Kurvatures[h][1] *= -1;
                        }

                        E_Kurvatures[h][0] /= e_ct[h];
                        E_Kurvatures[h][1] /= e_ct[h];
                    }

                    foreach(int h in omissions) { O_Edges.Remove(h); E_Frames.Remove(h); E_Kurvatures.Remove(h); }

                    Object_Edges[i] = O_Edges.Values.ToArray();
                    Edge_Frames[i] = E_Frames.Values.ToArray();
                    Edge_Kurvatures[i] = E_Kurvatures.Values.ToArray();
                    Edge_Normals[i] = e_norm.Values.ToArray();
                    
                    done += Model[i].Length;
                }

                //Set up a system to find random points on planes.//
                Plane_Area = new double[Topo[0].Plane_Members.Count];
                PolyPlaneFract = new double[Topo[0].Plane_Members.Count][];

                for (int q = 0; q < Topo[0].Plane_Members.Count; q++)
                {
                    foreach (int t in Topo[0].Plane_Members[q])
                    {
                        Plane_Area[q] += Topo[0].Polygon_Area(t);
                    }
                }

                for (int q = 0; q < Topo[0].Plane_Members.Count; q++)
                {
                    PolyPlaneFract[q] = new double[Topo[0].Plane_Members[q].Count];
                    PolyPlaneFract[q][0] = Topo[0].Polygon_Area(Topo[0].Plane_Members[q][0]) / Plane_Area[q];
                    for (int t = 1; t < Topo[0].Plane_Members[q].Count; t++)
                    {
                        PolyPlaneFract[q][t] += PolyPlaneFract[q][t - 1] + Topo[0].Polygon_Area(Topo[0].Plane_Members[q][t]) / Plane_Area[q];
                    }
                }
                Valid = true;
            }

            public Hare.Geometry.Edge[] Raw_Edge
            {
                get
                {
                    return Topo[0].planeEdges;
                }
            }

            public Hare.Geometry.Topology Hare_Data
            {
                get
                {
                    return Topo[0];
                }
            }

            public override void Absorb(ref OctaveRay Ray, out double cos_theta, double u, double v)
            {
                AbsorptionData[Ray.Surf_ID].Absorb(ref Ray, out cos_theta, Normal(Ray.Surf_ID, u, v));
            }

            public override void Absorb(ref BroadRay Ray, out double cos_theta, double u, double v)
            {
                AbsorptionData[Ray.Surf_ID].Absorb(ref Ray, out cos_theta, Normal(Ray.Surf_ID, u, v));                
            }

            public override void Scatter_Early(ref BroadRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, double cos_theta, double u, double v)
            {
                ScatteringData[Ray.Surf_ID].Scatter_Early(ref Ray, ref Rays, ref rand, Normal(Ray.Surf_ID, u, v), cos_theta);
            }

            public override void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, double cos_theta, double u, double v)
            {
                ScatteringData[Ray.Surf_ID].Scatter_Late(ref Ray, ref Rays, ref rand, Normal(Ray.Surf_ID, u, v), cos_theta);
            }

            public override void Scatter_Simple(ref OctaveRay Ray, ref Random rand, double cos_theta, double u, double v)
            {
                ScatteringData[Ray.Surf_ID].Scatter_VeryLate(ref Ray, ref rand, Normal(Ray.Surf_ID, u, v), cos_theta);
            }

            public Point[] polygon(int id)
            {
                if (id >= this.Count()) return new Point[0];
                return Topo[0].Polygon_Vertices(id);
            }

            public Point polygon_centroid(int id)
            {
                if (id >= this.Count()) return null;
                return Topo[0].Polygon_Centroid(id);
            }

            //private bool coplanarmesh(Point[][] M)
            //{
            //    Vector meanN = new Vector();
            //    if (M.FaceNormals.Count == 0) M.FaceNormals.ComputeFaceNormals();
            //    meanN = M.FaceNormals[0];
            //    meanN.Normalize();
            //    double total = 0;
            //    foreach (Vector n in M.FaceNormals) total += Math.Sqrt(meanN.X * n.X + meanN.Y * n.Y + meanN.Z * n.Z);
            //    return total == 0;
            //}

            public override void Register_Edges(IEnumerable<Hare.Geometry.Point> S, IEnumerable<Hare.Geometry.Point> R)
            {
                //For Curved Edge Analysis using the universal Mathnet NURBS structures:
                for (int i = 0; i < Edges.Count; i++)
                {
                    Edge_Nodes.Add(new Edge_Curved(ref S, ref R, Edge_isSoft[i], EdgeLength[i], Env_Prop, Edges[i], Edge_Tangents[0], Edge_Tangents[1]));
                }
            }
 
            /// <summary>
            /// used by the image source method to mirror sources over a face.
            /// </summary>
            /// <param name="Point">the point to mirror</param>
            /// <param name="Top_Id">the id of the topology to reference (for multi-resolution)</param>
            /// <param name="q">the index of the surface to use</param>
            /// <returns>the mirrored point</returns>
            public Hare.Geometry.Point Image(Hare.Geometry.Point Point, int Top_Id, int q)
            {
                //Mirror the point along the plane of the polygon...
                double Dist = Topo[0].DistToPlane(Point, q);
                Hare.Geometry.Point PX = Point - 2 * Topo[Top_Id].Normal(q) * Dist;
                return PX;
            }

            /// <summary>
            /// Generates a random point on a polygon in the model
            /// </summary>
            /// <param name="Plane_ID">the plane index</param>
            /// <param name="Polygon">a number from 0 to 1, which aids selection of a polygon on the plane.</param>
            /// <param name="Rnd1">random number 1</param>
            /// <param name="Rnd2">random number 2</param>
            /// <param name="Rnd3">random number 3</param>
            /// <returns>the random point</returns>
            //public Hare.Geometry.Point RandomPoint(int Plane_ID, double Polygon, double Rnd1, double Rnd2, double Rnd3)
            //{
            //    int i = 0;
            //    for (i = 0; i < PolyPlaneFract[Plane_ID].Length; i++)
            //    {
            //        if (PolyPlaneFract[Plane_ID][i] > Polygon) break;
            //    }
            //    return Topo[0].Polys[Topo[0].Plane_Members[Plane_ID][i]].GetRandomPoint(Rnd1, Rnd2, Rnd3);
            //}

            public override void partition()
            {
                Partitioned = true;
                partition(new List<Point>());
            }

            public override void partition(int SP_Param)
            {
                Partitioned = true;
                partition(new List<Point>(), SP_Param);
            }

            public override void partition(Point[] P, int SP_PARAM)
            {
                Partitioned = true;
                partition(new List<Point>(P), SP_PARAM);
            }

            public override void partition(List<Hare.Geometry.Point> P, int SP_PARAM)
            {
                Partitioned = true;
                for (int i = 0; i < Topo.Length; i++)
                {
                    Topo[i].Finish_Topology(P);
                }

                int spspec = Pach_Properties.Instance.SP_Spec();
                if (spspec == 0)
                {
                    SP = new Hare.Geometry.Voxel_Grid(Topo, SP_PARAM, 3);
                }
                else if (spspec == 1)
                {
                    //TODO: implement an Octree...
                    throw new NotImplementedException();
                }
            }

            public void partition(List<Hare.Geometry.Point> P)
            {
                Partitioned = true;
                partition(P, Pach_Properties.Instance.VG_Domain());
            }
            
            public override bool shoot(Hare.Geometry.Ray R, out double u, out double v, out int Poly_ID, out Hare.Geometry.Point X_PT, out double t)
            {
                Hare.Geometry.X_Event X = new Hare.Geometry.X_Event();
                if (SP.Shoot(R, 0, out X))
                {
                    Poly_ID = X.Poly_id;
                    X_PT = X.X_Point;
                    t = X.t;
                    u = 0;
                    v = 0;
                    return true;
                }
                Poly_ID = 0;
                X_PT = new Hare.Geometry.Point();
                //Rhino.RhinoDoc.ActiveDoc.Objects.Add(new Rhino.Geometry.LineCurve(Utilities.PachTools.HPttoRPt(R.origin), Utilities.PachTools.HPttoRPt(R.origin + R.direction)));
                t = 0;
                u = 0;
                v = 0;
                return false;
            }

            public override bool shoot(Hare.Geometry.Ray R, out double u, out double v, out int Poly_ID, out List<Hare.Geometry.Point> X_PT, out List<double> t, out List<int> code)
            {
                ///////////////////////
                //double L2 = (R.direction.x * R.direction.x + R.direction.y * R.direction.y + R.direction.z *R.direction.z);
                //if (L2 > 1.05 || L2 < 0.95)
                //{
                //    Rhino.RhinoApp.Write("Vectors have lost normalization...");
                //}
                ///////////////////////
                //R.direction.Normalize();
                ///////////////////////
                Hare.Geometry.X_Event X = new Hare.Geometry.X_Event();
                if (SP.Shoot(R, 0, out X))
                {
                    Poly_ID = X.Poly_id;
                    X_Event_NH XNH = X as X_Event_NH;
                    if (XNH != null)
                    {
                        X_PT = XNH.P_Points;
                        t = XNH.t_trav;
                        code = XNH.SPCode;
                    }
                    else
                    {
                        X_PT = new List<Hare.Geometry.Point>{X.X_Point};
                        t = new List<double>{X.t};
                        code = new List<int>{0};
                    }

                    u = 0;
                    v = 0;
                    return true;
                }
                Poly_ID = 0;
                X_PT = new List<Hare.Geometry.Point>();
                //Rhino.RhinoDoc.ActiveDoc.Objects.Add(new Rhino.Geometry.LineCurve(Utilities.PachTools.HPttoRPt(R.origin), Utilities.PachTools.HPttoRPt(R.origin + R.direction)));
                t = new List<double>();
                u = 0;
                v = 0;
                code = new List<int>() { 0 };
                return false;
            }

            public override bool shoot(Hare.Geometry.Ray R, int top_id, out Hare.Geometry.X_Event X)
            {
                if (SP.Shoot(R, top_id, out X)) return true;
                return false;
            }

            public override Hare.Geometry.Point ClosestPt(Hare.Geometry.Point P, ref double Dist)
            {
                double Max = double.MaxValue;
                Hare.Geometry.Point RP = new Hare.Geometry.Point();
                Hare.Geometry.Point TP;
                for (int i = 0; i < Topo[0].Polygon_Count; i++)
                {
                    TP = Topo[0].Closest_Point(P, i);
                    Dist = TP.x * TP.x + TP.y * TP.y + TP.z * TP.z;
                    if (Dist < Max)
                    {
                        RP = TP;
                        Max = Dist;
                    }
                }
                Dist = Math.Sqrt(Max);
                return RP;
            }

            /// <summary>
            /// Identifies the scene type.
            /// </summary>
            /// <returns></returns>
            public override string Scene_Type()
            {
                return "Polygon_Scene";
            }

            public override bool PointsInScene(List<Hare.Geometry.Point> PTS)
            {
                throw new NotImplementedException();
            }

            public override double Sound_speed(Hare.Geometry.Point pt)
            {
                return Env_Prop.Sound_Speed(pt);
            }

            public override double Sound_speed(int arg)
            {
                return Env_Prop.Sound_Speed(arg);
            }

            public override Hare.Geometry.Vector Normal(int i, double u, double v)
            {
                return Topo[0].Normal(i);
            }

            public override Hare.Geometry.Vector Normal(int i)
            {
                return Topo[0].Normal(i);
            }

            public override int Count()
            {
                return Topo[0].Polygon_Count;
            }

            public override bool IsPlanar(int object_id)
            {
                return !iscurved[object_id];
            }

            public double[] Kurvature(int object_id, int polyid)
            {
                return kurvature[object_id][polyid];
            }

            public Vector[] Frame_Axes(int object_id, int polyid)
            {
                return frame_axes[object_id][polyid];
            }

            public override double SurfaceArea(int x)
            {
                return Topo[0].Polygon_Area(x);
            }

            public override void EdgeFrame_Tangents(Hare.Geometry.Point Origin, Vector Normal, int[] PlaneIDs, ref List<double> dist2, List<Vector> Dir, List<int> EdgeIDs)
            {
                double d = Hare_math.Dot(Normal, Origin);
                //for (int i = 0; i < PlaneCount; i++)
                foreach(int i in EdgeIDs)
                {
                    Hare.Geometry.Point[] Pts = new Hare.Geometry.Point[2];
                    double d_min = double.MaxValue;
                    double d_max = double.MinValue;
                    foreach (int j in Topo[0].Plane_Members[i])
                    {
                        //Do the polygon/plane intersection for each member 'j' of i.
                        Hare.Geometry.Point[] vtx = Topo[0].Polygon_Vertices(j);
                        Hare.Geometry.Point[] temp = new Hare.Geometry.Point[1];
                        uint tmpcount = 0;

                        for (int k = 0, h = 1; k < vtx.Length; k++, h++)
                        {
                            Vector ab = vtx[h % vtx.Length] - vtx[k];
                            double t = (d - Hare_math.Dot(Normal, vtx[k])) / Hare_math.Dot(Normal, ab);

                            // If t in [0..1] compute and return intersection point
                            if (t >= 0.0f && t <= 1.0f)
                            {
                                temp[tmpcount] = vtx[k] + t * ab;
                                tmpcount++;
                            }
                            if (h == 0 && tmpcount == 0) break;
                            if (tmpcount > 1) break;
                        }
                        foreach (Hare.Geometry.Point p in temp)
                        {
                            Hare.Geometry.Point v = Origin - p;
                            double tmp = v.x * v.x + v.y * v.y + v.z * v.z;
                            if (tmp > d_max)
                            {
                                Pts[1] = p;
                                d_max = tmp;
                            }
                            if (tmp < d_min)
                            {
                                Pts[0] = p;
                                d_min = tmp;
                            }
                        }
                    }
                    dist2.Add(d_min);
                    EdgeIDs.Add(i);
                    Vector direction = (Pts[1] - Pts[0]);
                    direction.Normalize();
                    Dir.Add(direction);
                }
            }

            /// <summary>
            /// returns the number of planes. (not polygons)
            /// </summary>
            public int PlaneCount
            {
                get
                {
                    return Topo[0].Plane_Members.Count;
                }
            }

            /// <summary>
            /// returns the number of objects. (not polygons)
            /// </summary>
            public int ObjectCount
            {
                get
                {
                    return this.Object_Members.Length;
                }
            }

            /// <summary>
            /// gets the list of plane ids by polygon index.
            /// </summary>
            public int PlaneID(int i)
            {
                   return Topo[0].Polys[i].Plane_ID;
            }

            /// <summary>
            /// gets the list of plane ids by polygon index.
            /// </summary>
            public int ObjectID(int i)
            {
                return Object_ID[i];
            }

            public bool Box_Intersect(AABB box, out double abs, out Vector V)//, out int[] PolyIds, out double[] Abs, out double[] Trans, out double[] Scat)
            {
                List<int> isect;
                abs = 0;
                V = new Vector();
                SP.Box_Intersect(box, out isect);
                if (isect.Count == 0) return false;

                //if (isect.Count > 1)
                //{
                //    Rhino.RhinoApp.WriteLine(isect.Count.ToString());
                //}

                foreach(int i in isect)
                {
                    V += this.Normal(i);
                    abs += (1 - this.AbsorptionData[i].Reflection_Narrow(0)).Magnitude;
                    ///TODO:Find some intelligent way of applying absorption.
                }

                //if (double.IsNaN(V.y))
                //{
                //    Rhino.RhinoApp.WriteLine("Doh!");
                //}

                V.Normalize();
                abs /= isect.Count;

                return true;
            }

            /// <summary>
            /// gets the list of polygons on a each plane.
            /// </summary>
            public List<List<int>> PlaneMembers
            {
                get
                {
                    return Topo[0].Plane_Members;
                }
            }

            /// <summary>
            /// gets the list of polygons on each object.
            /// </summary>
            public int[][] ObjectMembers
            {
                get
                {
                    return this.Object_Members;
                }
            }

            public Hare.Geometry.Edge[][] ObjectMeshEdges
            {
                get
                {
                    return this.Object_Edges;
                }
            }

            public override Hare.Geometry.Point Max()
            {
                return this.Topo[0].Max;
            }

            public override Hare.Geometry.Point Min()
            {
                return this.Topo[0].Min;
            }
        }

        public class Empty_Scene : Scene
        {
            Hare.Geometry.Point minpt;
            Hare.Geometry.Point maxpt;

            public Empty_Scene(double Temp, double hr, double Pa, int Air_Choice, bool EdgeCorrection, bool IsAcoustic, Hare.Geometry.Point min, Hare.Geometry.Point max)
                : base(Temp, hr, Pa, Air_Choice, EdgeCorrection, IsAcoustic)
            {
                minpt = min;
                maxpt = max;
            }

            public override void Absorb(ref BroadRay Ray, out double cos_theta, double u, double v)
            {
                cos_theta = 0;
            }

            public override void Absorb(ref OctaveRay Ray, out double cos_theta, double u, double v)
            {
                cos_theta = 0;
            }

            public override Hare.Geometry.Point ClosestPt(Hare.Geometry.Point P, ref double D)
            {
                return new Hare.Geometry.Point(0, 0, 0);
            }

            public override int Count()
            {
                return 0;
            }

            public override void EdgeFrame_Tangents(Hare.Geometry.Point Origin, Vector Normal, int[] SrfIDs, ref List<double> dist2, List<Vector> Dir, List<int> IDs)
            {
                return;
            }

            public override bool IsPlanar(int i)
            {
                return false;
            }

            public override Hare.Geometry.Point Max()
            {
                return maxpt;
            }

            public override Hare.Geometry.Point Min()
            {
                return minpt;
            }

            public override Vector Normal(int i)
            {
                return new Vector(0, 0, 0);
            }

            public override void partition()
            {
                return;
            }

            public override void partition(Hare.Geometry.Point[] P, int SP_PARAM)
            {
                return;
            }

            public override void partition(int SP_Param)
            {
                return;
            }

            public override void partition(List<Hare.Geometry.Point> P, int SP_PARAM)
            {
                return;
            }

            public override Vector Normal(int i, double u, double v)
            {
                return new Vector(0, 0, 0);
            }

            public override bool PointsInScene(List<Hare.Geometry.Point> PTS)
            {
                return false;
            }

            public override void Register_Edges(IEnumerable<Hare.Geometry.Point> S, IEnumerable<Hare.Geometry.Point> R)
            {
                return;
            }

            public override void Scatter_Early(ref BroadRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, double cos_theta, double u, double v)
            {
                return;
            }

            public override void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, double cos_theta, double u, double v)
            {
                return;
            }

            public override void Scatter_Simple(ref OctaveRay Ray, ref Random rand, double cos_theta, double u, double v)
            {
                return;
            }

            public override string Scene_Type()
            {
                return "Empty Scene";
            }

            public override bool shoot(Ray R, out double u, out double v, out int Poly_ID, out Hare.Geometry.Point X_PT, out double t)
            {
                u = -1;
                v = -1;
                Poly_ID = -1;
                X_PT = new Hare.Geometry.Point();
                t = -1;
                return false;
            }

            public override bool shoot(Ray R, int topo, out X_Event Xpt)
            {
                Xpt = new X_Event(new Hare.Geometry.Point(), -1, -1, -1, -1);
                return false;
            }

            public override bool shoot(Ray R, out double u, out double v, out int Poly_ID, out List<Hare.Geometry.Point> X_PT, out List<double> t, out List<int> code)
            {
                u = -1;
                v = -1;
                Poly_ID = -1;
                X_PT = new List<Hare.Geometry.Point>();
                t = new List<double>();
                code = new List<int>();
                return false;
            }

            public override double Sound_speed(Hare.Geometry.Point pt)
            {
                return Env_Prop.Sound_Speed(pt);
            }

            public override double Sound_speed(int arg)
            {
                return Env_Prop.Sound_Speed(arg);
            }

            public override double SurfaceArea(int x)
            {
                return 0;
            }
        }
    }
}