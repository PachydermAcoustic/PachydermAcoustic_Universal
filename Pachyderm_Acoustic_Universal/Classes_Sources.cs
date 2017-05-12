//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2015, Arthur van der Harten 
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
            protected Hare.Geometry.Point H_Center;
            protected string type = "";
            protected int S_ID = 0;
            protected double delay = 0;
            protected Phase_Regime ph;

            /// <summary>
            /// Most explicit constructor.
            /// </summary>
            /// <param name="power_In_db">Sound Power of the source</param>
            /// <param name="Source">Origin of the source</param>
            /// <param name="TotalRays_in">Number of rays </param>
            /// <param name="Broadband_Time"></param>
            /// <param name="ID">The identifier of the source</param>
            public Source(double[] power_In_db, Point Source, Phase_Regime p, int ID)
            {
                ph = p;
                S_ID = ID;
                H_Center = new Hare.Geometry.Point(Source.x, Source.y, Source.z);
                Center = Source;                

                SPL[0] = power_In_db[0];
                SPL[1] = power_In_db[1];
                SPL[2] = power_In_db[2];
                SPL[3] = power_In_db[3];
                SPL[4] = power_In_db[4];
                SPL[5] = power_In_db[5];
                SPL[6] = power_In_db[6];
                SPL[7] = power_In_db[7];

                SourcePower[0] = 1E-12 * Math.Pow(10, .1 * SPL[0]);
                SourcePower[1] = 1E-12 * Math.Pow(10, .1 * SPL[1]);
                SourcePower[2] = 1E-12 * Math.Pow(10, .1 * SPL[2]);
                SourcePower[3] = 1E-12 * Math.Pow(10, .1 * SPL[3]);
                SourcePower[4] = 1E-12 * Math.Pow(10, .1 * SPL[4]);
                SourcePower[5] = 1E-12 * Math.Pow(10, .1 * SPL[5]);
                SourcePower[6] = 1E-12 * Math.Pow(10, .1 * SPL[6]);
                SourcePower[7] = 1E-12 * Math.Pow(10, .1 * SPL[7]);
            }

            /// <summary>
            /// Constructor assigns a default sound power level of 120 dB per octave.
            /// </summary>
            /// <param name="Source">Origin of the source</param>
            /// <param name="TotalRays_in">Number of rays </param>
            /// <param name="Broadband_Time"></param>
            /// <param name="ID">The identifier of the source</param>
            public Source(Point Source, int ID)
            {
                S_ID = ID;
                H_Center = new Hare.Geometry.Point(Source.x, Source.y, Source.z);
                Center = Source;

                SPL[0] = 120;
                SPL[1] = 120;
                SPL[2] = 120;
                SPL[3] = 120;
                SPL[4] = 120;
                SPL[5] = 120;
                SPL[6] = 120;
                SPL[7] = 120;

                SourcePower[0] = 1E-12 * Math.Pow(10, .1 * SPL[0]);
                SourcePower[1] = 1E-12 * Math.Pow(10, .1 * SPL[1]);
                SourcePower[2] = 1E-12 * Math.Pow(10, .1 * SPL[2]);
                SourcePower[3] = 1E-12 * Math.Pow(10, .1 * SPL[3]);
                SourcePower[4] = 1E-12 * Math.Pow(10, .1 * SPL[4]);
                SourcePower[5] = 1E-12 * Math.Pow(10, .1 * SPL[5]);
                SourcePower[6] = 1E-12 * Math.Pow(10, .1 * SPL[6]);
                SourcePower[7] = 1E-12 * Math.Pow(10, .1 * SPL[7]);
            }

            public virtual void AppendPts(ref List<Point> SPT)
            {
                SPT.Add(Origin());
            }

            /// <summary>
            /// get the sound power of the source by octave band.
            /// </summary>
            public virtual double[] SoundPower
            {
                get
                {
                    return SourcePower;// new double[] { SourcePower[0], SourcePower[1], SourcePower[2], SourcePower[3], SourcePower[4], SourcePower[5], SourcePower[6], SourcePower[7] };
                }
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
                return SoundPower;
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

            /// <summary>
            /// The list of directions. For stochastic calculations only.
            /// </summary>
            public abstract BroadRay Directions(int index, int thread, ref Random random);
            public abstract BroadRay Directions(int index, int thread, ref Random random, int[] Octaves);

            /// <summary>
            /// return the sound power of the source.
            /// </summary>
            /// <param name="Octave"></param>
            /// <returns></returns>
            public double TotalPower(int Octave)
            {
                if (Octave < 8) return SourcePower[Octave];
                return SourcePowerSum();
            }

            public double SourcePowerSum()
            {
                return SourcePower[0] + SourcePower[1] + SourcePower[2] + SourcePower[3] + SourcePower[4] + SourcePower[5] + SourcePower[6] + SourcePower[7];
            }
            
            /// <summary>
            /// The origin of the source.
            /// </summary>
            /// <returns></returns>
            public Point Origin()
            {
                return Center;
            }

            /// <summary>
            /// The origin of the source.
            /// </summary>
            /// <returns></returns>
            public Hare.Geometry.Point H_Origin()
            {
                return H_Center;
            }

            public virtual string Type()
            {
                return type;
            }   
        
            public int Source_ID()
            {
                return S_ID;
            }

            public double Delay
            {
                get
                {
                    return this.delay;
                }
                set
                {
                    delay = value;
                }
            }

            public enum Phase_Regime
            {
                PhaseMatched,
                PhaseRandom
            }
        }

        /// <summary>
        /// source object based on a random number generator.
        /// </summary>
        [Serializable]
        public class RandomSource: Source
        {
            public RandomSource(double[] power_in_db, Point Source, double delay_in, int ID)
            :this(power_in_db, Source, ID)
            {
                delay = delay_in;
            }

            public RandomSource(double[] power_in_db, Point Source, int ID)
                :base(power_in_db, Source, Phase_Regime.PhaseMatched, ID)
            {
                type = "PseudoRandom";
            }

            public RandomSource(double[] power_in_db, Point Source, Phase_Regime ph, int ID)
                : base(power_in_db, Source, ph, ID)
            {
                type = "PseudoRandom";
            }

            public override BroadRay Directions(int index, int thread, ref Random random)
            {
                double Theta = random.NextDouble() * 2 * System.Math.PI;
                double Phi = random.NextDouble() * 2 * System.Math.PI;
                Vector Direction = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));
                
                BroadRay B = new BroadRay(H_Center, Direction, random.Next(), thread, SourcePower, delay, Source_ID()); //Provides divided Power[stochastic]
                return B;
            }

            public override BroadRay Directions(int index, int thread, ref Random random, int[] Octaves)
            {
                double Theta = random.NextDouble() * 2 * System.Math.PI;
                double Phi = random.NextDouble() * 2 * System.Math.PI;
                Vector Direction = new Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));
                
                BroadRay B = new BroadRay(H_Center, Direction, random.Next(), thread, SourcePower, delay, Source_ID(), Octaves); //Provides divided Power[stochastic]
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
            protected int rayct = 0;
            
            public GeodesicSource(double[] power_in_db, Point Source, double delay_in, int ID)
            :this(power_in_db, Source, ID)
            {
                delay = delay_in;
            }

            public GeodesicSource(double[] power_in_db, Point Source, int ID)
                :base(power_in_db, Source, Phase_Regime.PhaseMatched, ID)
            {
                Random RAND = new Random();
                type = "Geodesic";
                GeoSphere(3);
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
                    P[Fnum] = new Hare.Geometry.Point[3]{P0, P1, P2};
                    Fnum++;
                }
            }

            public override BroadRay Directions(int index, int thread, ref Random random)
            {
                Hare.Geometry.Point Pt = T.Polys[rayct%T.Polygon_Count].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
                Hare.Geometry.Vector P = new Vector(Pt.x, Pt.y, Pt.z);
                P.Normalize();
                rayct++;

                return new BroadRay(H_Center, P, random.Next(), thread, SourcePower, delay, Source_ID()); //Provides divided Power[stochastic]
            }

            public override BroadRay Directions(int index, int thread, ref Random random, int[] Octaves)
            {
                Hare.Geometry.Point Pt = T.Polys[rayct % T.Polygon_Count].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
                Hare.Geometry.Vector P = new Vector(Pt.x, Pt.y, Pt.z);
                P.Normalize();
                rayct++;

                return new BroadRay(H_Center, P, random.Next(), thread, SourcePower, delay, Source_ID(), Octaves); //Provides divided Power[stochastic]
            }
        }

        [Serializable]
        public class SpeakerSource : GeodesicSource
        {
            Hare.Geometry.Voxel_Grid Balloon;

            public SpeakerSource(Speaker_Balloon S, double[] power_in_db, Point Source, int[] Bands, double delay_in, int ID)
            :this(S, power_in_db, Source, Bands, ID)
            {
                delay = delay_in;
            }

            public SpeakerSource(Speaker_Balloon S, double[] power_in_db, Point Source, int[] Bands, int ID)
                :base(power_in_db, Source, 0, ID)
            {
                for (int oct = 0; oct < 8; oct++)
                {
                    if (oct < Bands[0] || oct > Bands[1])
                    {
                        base.SourcePower[oct] = 0;
                        base.SPL[oct] = 0;
                    }
                }

                type = "Loudspeaker";
                Balloon = new Voxel_Grid(S.Balloons(power_in_db), 1);
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
            public override BroadRay Directions(int index, int thread, ref Random random)
            {
                X_Event X = new X_Event();
                
                double[] RayPower = new double[8];

                Hare.Geometry.Point Pt = T.Polys[base.rayct%T.Polygon_Count].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
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
                        Balloon.Shoot(new Ray(new Hare.Geometry.Point(0,0,0), P, thread, random.Next()), oct, out X);
                        RayPower[oct] = 1E-12 * Math.Pow(10, .1 * X.t);
                    }
                }

                double[] phtemp = new double[8];
                if (ph == Phase_Regime.PhaseRandom) for (int o = 0; o < 8; o++) phtemp[o] = random.Next() * 2 * Math.PI;
                else for (int o = 0; o < 8; o++) phtemp[o] = 0 - Delay * Utilities.Numerics.angularFrequency[o];
                base.rayct++;

                return new BroadRay(H_Center, P, random.Next(), thread, RayPower, delay, Source_ID());
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
                        Balloon.Shoot(new Ray(new Hare.Geometry.Point(0,0,0), DIR, thread, random), oct, out X);
                        RayPower[oct] = 1E-12 * Math.Pow(10, .1 * X.t);
                    }
                }
                return RayPower;
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

            public GeodesicMeshSource(double[] power_in_db, Point Source, int MinRays, int ID)
                : base(power_in_db, Source, Phase_Regime.PhaseMatched, ID)
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
                    Vector P3 = (P0 + P1) / 2;
                    P3.Normalize();

                    Vector P4 = (P1 + P2) / 2;
                    P4.Normalize();

                    Vector P5 = (P2 + P0) / 2;
                    P5.Normalize();

                    this.triangle(P0, P3, P5, Ord + 1, max);
                    this.triangle(P1, P4, P3, Ord + 1, max);
                    this.triangle(P2, P5, P4, Ord + 1, max);
                    this.triangle(P3, P4, P5, Ord + 1, max);
                }
                else
                {
                    P[Fnum] = new Hare.Geometry.Point[3] { P0, P1, P2 };
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

            public override BroadRay Directions(int index, int thread, ref Random random)
            {

                return new BroadRay(H_Center, Dir_Random(ref random), random.Next(), thread, SourcePower, delay, Source_ID());
            }

            public override BroadRay Directions(int index, int thread, ref Random random, int[] Octaves)
            {
                return new BroadRay(H_Center, Dir_Random(ref random), random.Next(), thread, SourcePower, delay, Source_ID(), Octaves);
            }
        }
    }
}