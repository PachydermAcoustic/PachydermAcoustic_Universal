//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2023, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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
using System.ComponentModel.Design.Serialization;
using System.Linq;
using Hare.Geometry;
using Pachyderm_Acoustic.Utilities;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        /// <summary>
        /// Source object base class. Do not use this object directly.
        /// </summary>
        [Serializable]
        public abstract class Source
        {
            protected Point Center;
            protected double[] SPL = new double[8];
            protected double[] SourcePower = new double[8];
            protected string type = "";
            protected int S_ID = 0;
            public double Rho_C = 411.6;
            /// <summary>
            /// Most explicit constructor.
            /// </summary>
            /// <param name="power_In_db">Sound Power of the source</param>
            /// <param name="Source">Origin of the source</param>
            /// <param name="TotalRays_in">Number of rays </param>
            /// <param name="Broadband_Time"></param>
            /// <param name="ID">The identifier of the source</param>
            public Source(double[] power_In_db, Point Source, int ID, bool Third_Octave)
            {
                S_ID = ID;
                Center = Source;
                if (Third_Octave)
                {
                    if (power_In_db.Length == 24) 
                    {
                        SourcePower = new double[24];
                        SPL = new double[24];
                        for (int i = 0; i < 24; i++)
                        {
                            SPL[i] = power_In_db[i];
                            SourcePower[i] = 1E-12 * Math.Pow(10, .1 * SPL[i]);
                        }
                    }
                    else
                    {
                        //...an interpolation algorithm for octave to third octave 
                        double[] freq = new double[9] {16, 63, 125, 250, 500, 1000, 2000, 4000, 8000 };
                        power_In_db = power_In_db.Reverse().ToArray();
                        Array.Resize(ref power_In_db, 9);
                        power_In_db[8] = double.Epsilon;
                        power_In_db = power_In_db.Reverse().ToArray();
                        MathNet.Numerics.Interpolation.CubicSpline powerspline = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkimaInplace(freq, power_In_db);
                        SourcePower = new double[24];
                        for (int oct = 1; oct < 9; oct++)
                        {
                            double total = 0;
                            for (int third = 1; third < 4; third++)
                            {
                                double f = 20 * Math.Pow(2, oct + (double)third / 3.0);
                                SourcePower[3 * oct + third - 4] = 1E-12 * Math.Pow(10, .1 * powerspline.Interpolate(f)) / 3.0;
                                total += SourcePower[3 * oct + third - 4] = 1E-12 * Math.Pow(10, .1 * powerspline.Interpolate(f)) / 3.0;
                            }
                            double ideal = 1E-12 * Math.Pow(10, .1 * power_In_db[1]);
                            for (int third = 1; third < 4; third++)
                            {
                                SourcePower[3 * oct + third - 4] *= ideal / total;
                                SPL[3 * oct + third - 4] = 10 * Math.Log10(SourcePower[3 * oct + third - 4] / 1E-12);
                            }
                        }
                    }
                }
                else
                {
                    SourcePower = new double[8];

                    for (int i = 0; i < 8; i++)
                    {
                        SPL[i] = power_In_db[i];
                        SourcePower[i] = 1E-12 * Math.Pow(10, .1 * SPL[i]);
                    }
                }
            }

            /// <summary>
            /// Constructor assigns a default sound power level of 120 dB per octave.
            /// </summary>
            /// <param name="Source">Origin of the source</param>
            /// <param name="ID">The identifier of the source</param>
            /// <param name="Third_Octave">Is the simulation to be performed in third octaves?</param>
            public Source(Point Source, int ID, bool Third_Octave)
            {
                S_ID = ID;
                Center = Source;

                int no_of_bands = Third_Octave ? 24 : 8;
                double oct_mod = Third_Octave ? Math.Log(3) : 0;

                for(int band = 0; band < no_of_bands; band++)
                {
                    SPL[band] = 120 - oct_mod;
                    SourcePower[band] = 1E-12 * Math.Pow(10, .1 * SPL[band]);
                }
            }

            public virtual void AppendPts(ref List<Point> SPT)
            {
                SPT.Add(Origin);
            }

            /// <summary>
            /// get the sound power of the source by octave band.
            /// </summary>
            public virtual double[] SoundPower
            {
                get
                {
                    return SourcePower;
                }
                set
                {
                    SourcePower = value;
                }
            }

            public virtual double[] Dir_Filter(int threadid, int random, Vector Direction, int sample_frequency, int length_starttofinish)
            {
                double[] power = DirPower(threadid, random, Direction);
                double[] OUT = new double[8];
                for (int oct = 0; oct < 8; oct++) OUT[oct] = Math.Sqrt(Math.Abs(power[oct]));
                return Audio.Pach_SP.Magnitude_Filter(OUT, sample_frequency, length_starttofinish, threadid);
            }

            /// <summary>
            /// Provides source power. (Deterministic Only)
            /// </summary>
            /// <param name="threadid"></param>
            /// <param name="random"></param>
            /// <param name="Direction"></param>
            /// <returns></returns>
            public virtual double[] DirPower(int threadid, int random, Vector Direction)
            {
                return SourcePower.Clone() as double[];
            }

            /// <summary>
            /// Provides source power. (Deterministic Only)
            /// </summary>
            /// <param name="threadid"></param>
            /// <param name="random"></param>
            /// <param name="Direction"></param>
            /// <returns></returns>
            public virtual double[] DirPressure(int threadid, int random, Vector Direction)
            {
                double[] H = SourcePower.Clone() as double[];
                for (int i = 0; i < H.Length; i++) H[i] = AcousticalMath.Pressure_Intensity(H[i], this.Rho_C);
                return H;
            }

            /// <summary>
            /// Sound Power Level of the Source. Octave 8 yields the logarithmic sum of all sound power levels.
            /// </summary>
            /// <param name="octave"></param>
            /// <returns></returns>
            public virtual double SWL(int octave)
            {
                if (octave < 8) return SPL[octave];
                return 10 * Math.Log10(Math.Pow(10, .1 * SPL[0]) + Math.Pow(10, .1 * SPL[1]) + Math.Pow(10, .1 * SPL[2]) + Math.Pow(10, .1 * SPL[3]) + Math.Pow(10, .1 * SPL[4]) + Math.Pow(10, .1 * SPL[5]) + Math.Pow(10, .1 * SPL[6]) + Math.Pow(10, .1 * SPL[7]));
            }

            public virtual double[] SWL()
            {
                return SPL.Clone() as double[];
            }

            public virtual int Revolution_Period()
            {
                return 80; //Default value for the revolution period - rays per revolution.
            }

            /// <summary>
            /// The list of directions. For stochastic calculations only.
            /// </summary>
            public abstract BroadRay Directions(int thread, ref Random random);
            public abstract BroadRay Directions(int thread, ref Random random, int[] Octaves);
            
            /// <summary>
            /// The Origin of the source.
            /// </summary>
            /// <returns></returns>
            public Point Origin
            {
                get
                {
                    return Center;
                }
            }

            public virtual string Type()
            {
                return type;
            }   
        
            public int Source_ID()
            {
                return S_ID;
            }

            public virtual void Lighten()
            { }
        }

        /// <summary>
        /// source object based on a random number generator.
        /// </summary>
        [Serializable]
        public class RandomSource: Source
        {
            public RandomSource(double[] power_in_db, Point Source, int ID, bool third_octave)
                :base(power_in_db, Source, ID, third_octave)
            {
                type = "PseudoRandom";
            }

            public override BroadRay Directions(int thread, ref Random random)
            {
                double Theta = random.NextDouble() * 2 * System.Math.PI;
                double Phi = random.NextDouble() * 2 * System.Math.PI;
                Vector Direction = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));

                BroadRay B = BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, Direction.dx, Direction.dy, Direction.dz, random.Next(), thread, DirPower(thread, random.Next(), Direction), 0, Source_ID());
                return B;
            }

            public override BroadRay Directions(int thread, ref Random random, int[] Octaves)
            {
                double Theta = random.NextDouble() * 2 * System.Math.PI;
                double Phi = random.NextDouble() * 2 * System.Math.PI;
                Vector Direction = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));

                //BroadRay B = new BroadRay(Origin.x, Origin.y, Origin.z, Direction.dx, Direction.dy, Direction.dz, random.Next(), thread, DirPower(thread, random.Next(), Direction), 0, Source_ID(), Octaves);
                BroadRay B = BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, Direction.dx, Direction.dy, Direction.dz, random.Next(), thread, DirPower(thread, random.Next(), Direction), 0, Source_ID(), Octaves);
                return B;
            }
        }

        /// <summary>
        /// source object based on the equal area faces of a geodesic sphere.
        /// </summary>
        [Serializable]
        public class GeodesicSource : Source
        {
            protected Topology T;
            int Fnum = 0;
            protected int rayct = -1;
            
            public GeodesicSource(double[] power_in_db, Point Source, int ID, bool Third_Octave)
                :base(power_in_db, Source, ID, Third_Octave)
            {
                Random RAND = new Random();
                type = "Geodesic";
                GeoSphere(1);
            }

            /// <summary>
            /// creates a geodesic sphere.
            /// </summary>
            /// <param name="order"></param>
            public void GeoSphere(int order)
            {
                double sqr5 = System.Math.Sqrt(5.0);
                double phi = (1.0 + sqr5) * 0.5; // golden ratio
                double ratio = System.Math.Sqrt(10.0 + (2.0 * sqr5)) / (4.0 * phi);
                double a = (.25 / ratio) * 0.5;
                double b = (.25 / ratio) / (2.0 * phi);

                // Define the icosahedron's 12 vertices
                Vector P0 = new Vector(0, b, -a);
                Vector P1 = new Vector(b, a, 0);
                Vector P2 = new Vector(-b, a, 0);
                Vector P3 = new Vector(0, b, a);
                Vector P4 = new Vector(0, -b, a);
                Vector P5 = new Vector(-a, 0, b);
                Vector P6 = new Vector(0, -b, -a);
                Vector P7 = new Vector(a, 0, -b);
                Vector P8 = new Vector(a, 0, b);
                Vector P9 = new Vector(-a, 0, -b);
                Vector P10 = new Vector(b, -a, 0);
                Vector P11 = new Vector(-b, -a, 0);

                P0.Normalize();
                P1.Normalize();
                P2.Normalize();
                P3.Normalize();
                P4.Normalize();
                P5.Normalize();
                P6.Normalize();
                P7.Normalize();
                P8.Normalize();
                P9.Normalize();
                P10.Normalize();
                P11.Normalize();

                Hare.Geometry.Point[][] P = new Hare.Geometry.Point[20 * (int)Math.Pow(4, order)][];

                //Create the icosahedron's 20 triangular faces
                triangle(P0, P1, P2, 0, order, ref P);
                triangle(P3, P2, P1, 0, order, ref P);
                triangle(P3, P4, P5, 0, order, ref P);
                triangle(P3, P8, P4, 0, order, ref P);
                triangle(P0, P6, P7, 0, order, ref P);
                triangle(P0, P9, P6, 0, order, ref P);
                triangle(P4, P10, P11, 0, order, ref P);
                triangle(P6, P11, P10, 0, order, ref P);
                triangle(P2, P5, P9, 0, order, ref P);
                triangle(P11, P9, P5, 0, order, ref P);
                triangle(P1, P7, P8, 0, order, ref P);
                triangle(P10, P8, P7, 0, order, ref P);
                triangle(P3, P5, P2, 0, order, ref P);
                triangle(P3, P1, P8, 0, order, ref P);
                triangle(P0, P2, P9, 0, order, ref P);
                triangle(P0, P7, P1, 0, order, ref P);
                triangle(P6, P9, P11, 0, order, ref P);
                triangle(P6, P10, P7, 0, order, ref P);
                triangle(P4, P11, P5, 0, order, ref P);
                triangle(P4, P8, P10, 0, order, ref P);

                T = new Topology(P);
            }
            
            private void triangle(Vector P0, Vector P1, Vector P2, int Ord, int max, ref Hare.Geometry.Point[][] P)
            {
                if (Ord < max)
                {
                    Vector P3 = (P0 + P1) / 2;
                    P3.Normalize();
                    
                    Vector P4 = (P1 + P2) / 2;
                    P4.Normalize();
                    
                    Vector P5 = (P2 + P0) / 2;
                    P5.Normalize();
                    
                    this.triangle(P0, P3, P5, Ord + 1, max, ref P);
                    this.triangle(P1, P4, P3, Ord + 1, max, ref P);
                    this.triangle(P2, P5, P4, Ord + 1, max, ref P);
                    this.triangle(P3, P4, P5, Ord + 1, max, ref P);
                }
                else
                {
                    P[Fnum] = new Hare.Geometry.Point[3]{new Point(P0.dx, P0.dy, P0.dz), new Point(P1.dx, P1.dy, P1.dz), new Point(P2.dx, P2.dy, P2.dz)};
                    Fnum++;
                }
            }

            public override BroadRay Directions(int thread, ref Random random)
            {
                lock (ctlock)
                {
                    rayct++;
                    Hare.Geometry.Point Pt = T.Polys[rayct % T.Polygon_Count].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
                    Hare.Geometry.Vector P = new Vector(Pt.x, Pt.y, Pt.z);
                    P.Normalize();
                    return BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, P.dx, P.dy, P.dz, random.Next(), thread, DirPower(thread, random.Next(), P), 0, Source_ID()); //Provides divided Power[stochastic]
                }
            }

            object ctlock = new object();
            public override BroadRay Directions(int thread, ref Random random, int[] Octaves)
            {
                Hare.Geometry.Point Pt;
                lock (ctlock)
                {
                    rayct++;
                    Pt = T.Polys[rayct % T.Polygon_Count].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
                }
                Hare.Geometry.Vector P = new Vector(Pt.x, Pt.y, Pt.z);
                P.Normalize();

                //return new BroadRay(Origin.x, Origin.y, Origin.z, P.dx, P.dy, P.dz, random.Next(), thread, DirPower(thread, random.Next(), P), 0, Source_ID(), Octaves); //Provides divided Power[stochastic]
                return BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, P.dx, P.dy, P.dz, random.Next(), thread, DirPower(thread, random.Next(), P), 0, Source_ID(), Octaves); //Provides divided Power[stochastic]
            }
        }

        [Serializable]
        public class DirectionalSource : GeodesicSource
        {
            Balloon _S;
            Hare.Geometry.Voxel_Grid Balloon;

            public DirectionalSource(Balloon S, double[] power_in_db, Point Source, int[] Bands, int ID, bool Third_Oct)
                :base(power_in_db, Source, ID, Third_Oct)
            {
                _S = S;
                for (int oct = 0; oct < 8; oct++)
                {
                    if (oct < Bands[0] || oct > Bands[1])
                    {
                        base.SourcePower[oct] = 0;
                        base.SPL[oct] = 0;
                    }
                }

                type = "Directional";
                //Balloon = new Voxel_Grid(S.Balloons(new double[8] {120,120,120,120,120,120,120,120}), 1);
                ///Testing///
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[0]);
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[1]);
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[2]);
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[3]); 
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[4]);
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[5]);
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[6]);
                //Utilities.PachTools.Plot_Hare_Topology(Balloon.Model[7]);
                /////////////
            }

            /// <summary>
            /// Provides a ray from the source. (Stochastic only)
            /// </summary>
            /// <param name="index"></param>
            /// <param name="thread"></param>
            /// <param name="random"></param>
            /// <returns></returns>
            public override BroadRay Directions(int thread, ref Random random)
            {
                if (Balloon == null) Balloon = new Voxel_Grid(_S.Balloons(new double[8] { 120, 120, 120, 120, 120, 120, 120, 120 }), 1);

                X_Event X = new X_Event();

                double[] RayPower = new double[8];

                Hare.Geometry.Point Pt = T.Polys[base.rayct % T.Polygon_Count].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
                Hare.Geometry.Vector P = new Vector(Pt.x, Pt.y, Pt.z);
                P.Normalize();

                for (int oct = 0; oct < 8; oct++)
                {
                    if (base.SPL[oct] == 0)
                    {
                        RayPower[oct] = 0;
                    }
                    else
                    {
                        Balloon.Shoot(new Ray(new Hare.Geometry.Point(0, 0, 0), P, thread, random.Next()), oct, out X);
                        RayPower[oct] = 1E-12 * Math.Pow(10, .1 * X.t);
                    }
                }

                base.rayct++;

                //return new BroadRay(Origin.x, Origin.y, Origin.z, P.dx, P.dy, P.dz, random.Next(), thread, RayPower, 0, Source_ID());
                return BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, P.dx, P.dy, P.dz, random.Next(), thread, RayPower, 0, Source_ID());
            }

            /// <summary>
            /// Deterministic Source Power
            /// </summary>
            /// <param name="thread"></param>
            /// <param name="random"></param>
            /// <param name="DIR"></param>
            /// <returns></returns>
            public override double[] DirPower(int thread, int random, Vector DIR)
            {
                if (Balloon == null) Balloon = new Voxel_Grid(_S.Balloons(new double[8] { 120, 120, 120, 120, 120, 120, 120, 120 }), 1);

                X_Event X = new X_Event();
                double[] RayPower = new double[8];
                for (int oct = 0; oct < 8; oct++)
                {
                    if (base.SPL[oct] == 0)
                    {
                        RayPower[oct] = 0;
                    }
                    else
                    {
                        Balloon.Shoot(new Ray(new Hare.Geometry.Point(0, 0, 0), DIR, thread, random), oct, out X);
                        RayPower[oct] = 1E-12 * Math.Pow(10, .1 * X.t);
                    }
                }
                return RayPower;
            }

            public override double[] DirPressure(int threadid, int random, Vector Direction)
            {
                double[] H = DirPower (threadid, random, Direction);
                for (int i = 0; i < H.Length; i++) H[i] = AcousticalMath.Pressure_Intensity(H[i], this.Rho_C);
                return H;
            }

            public override void Lighten()
            {
                base.Lighten();
                _S = null;
                GC.Collect();
            }
        }

        /// <summary>
        /// source object based on the vertices of a geodesic sphere.
        /// </summary>
        [Serializable]
        public class GeodesicMeshSource : Source
        {
            public Topology T;
            Hare.Geometry.Point[][] P;
            int Fnum = 0;
            int ct = 0;

            public GeodesicMeshSource(double[] power_in_db, Point Source, int MinRays, int ID, bool Third_Octave)
                : base(power_in_db, Source, ID, Third_Octave)
            {
                Random RAND = new Random();
                type = "Geodesic";
                int j = 2;
                while(true)
                {
                    Fnum = 0;
                    j++;
                    GeoSphere(j);
                    if (T.Vertex_Count > MinRays) break;
                }
            }

            /// <summary>
            /// creates a geodesic sphere.
            /// </summary>
            /// <param name="order"></param>
            public void GeoSphere(int order)
            {
                double sqr5 = System.Math.Sqrt(5.0);
                double phi = (1.0 + sqr5) * 0.5; // golden ratio
                double ratio = System.Math.Sqrt(10.0 + (2.0 * sqr5)) / (4.0 * phi);
                double a = (.25 / ratio) * 0.5;
                double b = (.25 / ratio) / (2.0 * phi);

                // Define the icosahedron's 12 vertices
                Vector P0 = new Vector(0, b, -a);
                Vector P1 = new Vector(b, a, 0);
                Vector P2 = new Vector(-b, a, 0);
                Vector P3 = new Vector(0, b, a);
                Vector P4 = new Vector(0, -b, a);
                Vector P5 = new Vector(-a, 0, b);
                Vector P6 = new Vector(0, -b, -a);
                Vector P7 = new Vector(a, 0, -b);
                Vector P8 = new Vector(a, 0, b);
                Vector P9 = new Vector(-a, 0, -b);
                Vector P10 = new Vector(b, -a, 0);
                Vector P11 = new Vector(-b, -a, 0);

                P0.Normalize();
                P1.Normalize();
                P2.Normalize();
                P3.Normalize();
                P4.Normalize();
                P5.Normalize();
                P6.Normalize();
                P7.Normalize();
                P8.Normalize();
                P9.Normalize();
                P10.Normalize();
                P11.Normalize();

                P = new Hare.Geometry.Point[20 * (int)Math.Pow(4, order)][];

                //Create the icosahedron's 20 triangular faces
                triangle(P0, P1, P2, 0, order);
                triangle(P3, P2, P1, 0, order);
                triangle(P3, P4, P5, 0, order);
                triangle(P3, P8, P4, 0, order);
                triangle(P0, P6, P7, 0, order);
                triangle(P0, P9, P6, 0, order);
                triangle(P4, P10, P11, 0, order);
                triangle(P6, P11, P10, 0, order);
                triangle(P2, P5, P9, 0, order);
                triangle(P11, P9, P5, 0, order);
                triangle(P1, P7, P8, 0, order);
                triangle(P10, P8, P7, 0, order);
                triangle(P3, P5, P2, 0, order);
                triangle(P3, P1, P8, 0, order);
                triangle(P0, P2, P9, 0, order);
                triangle(P0, P7, P1, 0, order);
                triangle(P6, P9, P11, 0, order);
                triangle(P6, P10, P7, 0, order);
                triangle(P4, P11, P5, 0, order);
                triangle(P4, P8, P10, 0, order);

                T = new Topology(P);
            }

            private void triangle(Vector P0, Vector P1, Vector P2, int Ord, int max)
            {
                if (Ord < max)
                {
                    //Vector P3 = (P0 + P1) / 2;
                    Vector P3 = new Vector(P0.dx + P1.dx, P0.dy + P1.dy, P0.dz + P1.dz);
                    P3.Normalize();

                    //Vector P4 = (P1 + P2) / 2;
                    Vector P4 = new Vector(P1.dx + P2.dx, P1.dy + P2.dy, P1.dz + P2.dz);
                    P4.Normalize();

                    //Vector P5 = (P2 + P0) / 2;
                    Vector P5 = new Vector(P2.dx + P0.dx, P2.dy + P0.dy, P2.dz + P0.dz);
                    P5.Normalize();

                    this.triangle(P0, P3, P5, Ord + 1, max);
                    this.triangle(P1, P4, P3, Ord + 1, max);
                    this.triangle(P2, P5, P4, Ord + 1, max);
                    this.triangle(P3, P4, P5, Ord + 1, max);
                }
                else
                {
                    P[Fnum] = new Hare.Geometry.Point[3] { new Point(P0.dx, P0.dy, P0.dz), new Point( P1.dx, P1.dy, P1.dz), new Point(P2.dx, P2.dy,P2.dz) };
                    Fnum++;
                }
            }

            private Vector Dir_Random(ref Random rnd)
            {
                Vector d = T.Polys[ct%T.Polygon_Count].GetRandomPoint(rnd.Next(), rnd.Next(), rnd.Next()) - Center;
                ct++;
                d.Normalize();
                return d;
            }

            public override BroadRay Directions(int thread, ref Random random)
            {
                Vector d = Dir_Random(ref random);
                return BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, d.dx, d.dy, d.dz, random.Next(), thread, new double[8] { 1, 1, 1, 1, 1, 1, 1, 1 }, 0, Source_ID());
            }

            public override BroadRay Directions(int thread, ref Random random, int[] Octaves)
            {
                Vector d = Dir_Random(ref random);
                return BroadRayPool.Instance.new_BroadRay(Origin.x, Origin.y, Origin.z, d.dx, d.dy, d.dy, random.Next(), thread, new double[8] { 1, 1, 1, 1, 1, 1, 1, 1 }, 0, Source_ID(), Octaves);
            }
        }

        public class SourceCluster: Source
        {
            public List<Source> Sources = new List<Source>();
            public bool Third_Octave = false;
            int revolution;
            protected int rayct;

            public SourceCluster(List<Source> sources, int id, bool third_octave = false)
                :base(new double[8] { 0, 0, 0, 0, 0, 0, 0, 0 }, new Point(0,0,0), id, third_octave)
            {
                Third_Octave = third_octave;
                Sources = sources;
                for(int i = 0; i < Sources.Count; i++)
                {
                    for (int o = 0; o < 8; o++) this.SPL[o] = Math.Max(Sources[i].SWL()[o], SPL[o]);
                    for (int o = 0; o < 8; o++) this.SourcePower[o] = Math.Max(Sources[i].SoundPower[o], SourcePower[o]);
                    revolution += Sources[i].Revolution_Period();
                }
                this.type = "Cluster";
            }

            public override double[] DirPower(int threadid, int random, Vector Direction)
            {
                throw new InvalidOperationException();
            }  

            public override double[] DirPressure(int threadid, int random, Vector Direction)
            {
                throw new InvalidOperationException();
            }

            public override void Lighten()
            {
                foreach (Source S in Sources)
                {
                    S.Lighten();
                }
            }

            public override BroadRay Directions(int thread, ref Random random)
            {
                rayct++;
                return Sources[rayct%Sources.Count].Directions(thread, ref random);
            }

            public override BroadRay Directions(int thread, ref Random random, int[] Octaves)
            {
                rayct++;
                return Sources[rayct % Sources.Count].Directions(thread, ref random, Octaves);
            }

            public override void AppendPts(ref List<Point> SPT)
            {
                foreach (Source S in Sources)
                {
                    S.AppendPts(ref SPT);
                }   
            }

            public override double[] Dir_Filter(int threadid, int random, Vector Direction, int sample_frequency, int length_starttofinish)
            {
                return base.Dir_Filter(threadid, random, Direction, sample_frequency, length_starttofinish);
            }
        }
    }
}