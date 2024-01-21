//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2023, Arthur van der Harten 
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
using System.Linq;
using Hare.Geometry;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        public class LineSource : Source
        {
            /// <summary>
            /// Sample points on curves.
            /// </summary>
            public Point[] Samples;
            /// <summary>
            /// Sound Power of each sample on the curve. (10^Lp/10 * L / no_of_samples)
            /// </summary>m
            public double[] Power;
            /// <summary>
            /// User designated sound power of each line source.
            /// </summary>
            double[] Level;
            /// <summary>
            /// number of samples per meter.
            /// </summary>
            double samplespermeter = 16;
            /// <summary>
            /// Private member controlling the directional characteristics of the line source.
            /// </summary>
            public Directionality D;
            /// <summary>
            /// The bounding box of line/curve source.
            /// </summary>
            public AABB bounds;
            /// <summary>
            /// Describes construction for a simple (generic) line source, with no unusual directionality characteristics.
            /// </summary>
            /// <param name="samples"></param>
            /// <param name="length"></param>
            /// <param name="Code"></param>
            /// <param name="el_m"></param>
            /// <param name="SrcID"></param>
            /// <param name="ph"></param>
            public LineSource(Hare.Geometry.Point[] samples, double length, string Code, double el_m, int SrcID, bool Third_Octave)
                : base(new double[8] { 60, 49, 41, 35, 31, 28, 26, 24 }, new Point(0, 0, 0), SrcID, Third_Octave)
            {
                //TODO: Accommodate third octave
                samplespermeter = el_m;

                //Divide curve up in ~equal length segments.
                Samples = samples;
                D = new Simple();

                Level = Utilities.PachTools.DecodeSourcePower(Code);
                Power = new double[8];

                double minx = double.PositiveInfinity, maxx = double.NegativeInfinity;
                for (int x = 0; x < Samples.Length; x++)
                {
                    minx = Math.Min(minx, Samples[x].x);
                    maxx = Math.Max(maxx, Samples[x].x);
                }
                double miny = double.PositiveInfinity, maxy = double.NegativeInfinity;
                for (int x = 0; x < Samples.Length; x++)
                {
                    miny = Math.Min(miny, Samples[x].y);
                    maxy = Math.Max(maxy, Samples[x].y);
                }
                double minz = double.PositiveInfinity, maxz = double.NegativeInfinity;
                for (int x = 0; x < Samples.Length; x++)
                {
                    minz = Math.Min(minz, Samples[x].z);
                    maxz = Math.Max(maxz, Samples[x].z);
                }

                bounds = new AABB(new Point(minx, miny, minz), new Point(maxx, maxy, maxz));

                double PowerMod = length / (double)Samples.Length;
                for (int oct = 0; oct < 8; oct++) Power[oct] = 1E-12 * Math.Pow(10, .1 * Level[oct]) * PowerMod;
            }

            public Hare.Geometry.Point ClosestPoint(Hare.Geometry.Point p, out int sel)
            {
                double min = double.PositiveInfinity;
                sel = -1;
                for (int i = 0; i < Samples.Length; i++)
                {
                    Hare.Geometry.Vector v = p - Samples[i];
                    double lsq = Hare_math.Dot(v, v);
                    if (lsq < min)
                    {
                        min = lsq;
                        sel = i;
                    }
                }
                return Samples[sel];
            }

            /// <summary>
            /// Describes construction for an aircraft runway (ANCON method).
            /// </summary>
            /// <param name="samples"></param>
            /// <param name="length"></param>
            /// <param name="_velocity"></param>
            /// <param name="_delta"></param>
            /// <param name="Code"></param>
            /// <param name="el_m"></param>
            /// <param name="SrcID"></param>
            /// <param name="ph"></param>
            public LineSource(Hare.Geometry.Point[] samples, double length, double _velocity, double _delta, string Code, int el_m, int SrcID, bool Third_Octave)
            : base(new double[8] { 60, 49, 41, 35, 31, 28, 26, 24 }, new Point(0, 0, 0), SrcID, Third_Octave)
            {
                //TODO: accommodate third octave...
                samplespermeter = el_m;

                Samples = samples;
                double velocity = _velocity;
                double delta = _delta;
                D = new ANCON(delta, velocity);

                Level = Utilities.PachTools.DecodeSourcePower(Code);
                Power = new double[8];

                double PowerMod = length / (double)Samples.Length;
                for (int oct = 0; oct < 8; oct++) Power[oct] = 1E-12 * Math.Pow(10, .1 * Level[oct]) * PowerMod;
            }

            public override string Type()
            {
                return "Line Source";
            }

            public override void AppendPts(ref List<Hare.Geometry.Point> SPT)
            {
                for (int j = 0; j < Samples.Length; j++)
                {
                    SPT.Add(Samples[j]);
                }
            }

            public override BroadRay Directions(int thread, ref Random random, int[] Octaves)
            {
                BroadRay B = Directions(thread, ref random);
                //B.Octaves = Octaves;
                return B;
            }

            public override BroadRay Directions(int thread, ref Random random)
            {
                return D.Directions(thread, ref random, ref Samples, ref Power, ref S_ID);
            }

            public double[] DirPower(Vector Direction, int sample_id)
            {
                //int sel;
                //ClosestPoint(pt, out sel);
                return D.DirPower(Direction, ref Samples, sample_id, ref Power);
            }

            public abstract class Directionality
            {
                protected int ct = -1;
                public abstract BroadRay Directions(int thread, ref Random random, ref Hare.Geometry.Point[] samples, ref double[] DomainPower, ref int S_ID);
                public abstract double[] DirPower(Vector Direction, ref Hare.Geometry.Point[] samples, int i, ref double[] DomainPower);
            }

            public class Simple : Directionality
            {
                public override BroadRay Directions(int thread, ref Random random, ref Hare.Geometry.Point[] samples, ref double[] DomainPower, ref int S_ID)
                {
                    ct++;
                    double pos = random.NextDouble();
                    Point P = samples[(int)Math.Floor(pos * samples.Length)];//ct%samples.Length];

                    double Theta = random.NextDouble() * 2 * System.Math.PI;
                    double Phi = random.NextDouble() * 2 * System.Math.PI;
                    Vector Direction = new Hare.Geometry.Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));

                    return new BroadRay(P, Direction, random.Next(), thread, DomainPower, 0, S_ID);
                }

                public override double[] DirPower(Vector Direction, ref Hare.Geometry.Point[] samples, int i, ref double[] DomainPower)
                {
                    return DomainPower.Clone() as double[];
                }
            }

            public double ElementsPerMeter
            { 
                get 
                {
                    return samplespermeter; 
                } 
            }

            public class ANCON : Directionality
            {
                double delta; //slant angle in radians
                double reciprocal_velocity; //m/s
                double dLinf = 10 * Math.Pow(10, 0.8 / 10);

                public ANCON(double delta_in, double velocity)
                {
                    delta = delta_in * Math.PI / 180; //degrees to radians
                    reciprocal_velocity = 1 / velocity;
                }

                public override double[] DirPower( Vector direction, ref Hare.Geometry.Point[] samples, int i, ref double[] domainpower)
                {
                    double[] raypower = new double[8];

                    Hare.Geometry.Point P = samples[i];
                    Hare.Geometry.Vector fore;
                    if (i == 0)
                    {
                        fore = P - samples[1];
                    }
                    else if (i == samples.Length - 1)
                    {
                        fore = P - samples[i - 1];
                    }
                    else
                    {
                        fore = ((P - samples[i - 1]) + (samples[i + 1] - P)) / 2;
                    }

                    double cosphi = Hare.Geometry.Hare_math.Dot(direction, fore);
                    double sinphi = Math.Sqrt(1 - cosphi * cosphi);
                    double tanphi = sinphi / cosphi;
                    double f_r = sinphi * sinphi * Math.Pow((tanphi * tanphi + 1) / (tanphi * tanphi + (1 + tanphi * Math.Tan(delta))), 1.5);

                    for (int oct = 0; oct < 8; oct++)
                    {
                        if (domainpower[oct] == 0)
                        {
                            raypower[oct] = 0;
                        }
                        else
                        {
                            raypower[oct] = domainpower[oct] * reciprocal_velocity * dLinf * f_r;
                        }
                    }
                    return raypower;
                }

                public override BroadRay Directions(int thread, ref Random random, ref Hare.Geometry.Point[] samples, ref double[] DomainPower, ref int S_ID)
                {
                    ct++;
                    double pos = random.NextDouble();
                    int i = ct % samples.Length;
                    Point P = samples[i];
                    Hare.Geometry.Vector fore;
                    if (i == 0)
                    {
                        fore = P - samples[1];
                    }
                    else if (i == samples.Length-1)
                    {
                        fore = P - samples[i-1];
                    }
                    else
                    {
                        fore = ((P - samples[i - 1]) + (samples[i + 1] - P)) / 2;
                    }

                    double Theta = random.NextDouble() * 2 * System.Math.PI;
                    double Phi = random.NextDouble() * 2 * System.Math.PI;
                    Vector Direction = new Hare.Geometry.Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));

                    double cosphi = Hare.Geometry.Hare_math.Dot(Direction, fore);
                    double sinphi = Math.Sqrt(1 - cosphi * cosphi);
                    double tanphi = sinphi / cosphi;
                    double F_r = sinphi * sinphi * Math.Pow((tanphi * tanphi + 1) / (tanphi * tanphi + (1 + tanphi * Math.Tan(delta))), 1.5);

                    double[] power = new double[8];
                    for (int oct = 0; oct < 8; oct++) power[oct] = DomainPower[oct] * reciprocal_velocity * dLinf * F_r;

                    return new BroadRay(P, Direction, random.Next(), thread, DomainPower, 0, S_ID);
                }
            }
        }
    }
}