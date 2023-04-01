﻿//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
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
using Hare.Geometry;
using System.Linq;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        public class SurfaceSource: Source
        {
            /// <summary>
            /// Sample points on curves.
            /// </summary>
            public Hare.Geometry.Point[] Samples;
            /// <summary>
            /// Sound Power of each sample on the surface. (10^Lp/10 * L / no_of_samples)
            /// </summary>m
            double[] DomainPower;
            /// <summary>
            /// User designated sound power of surface source.
            /// </summary>
            public double[] DomainLevel;
            /// <summary>
            /// number of samples per meter.
            /// </summary>
            int samplespermeter = 16;
            /// <summary>
            /// Topology of the meshed surface.
            /// </summary>
            Topology T;

            public SurfaceSource(Hare.Geometry.Point[] samples, Topology T_in, String CodeList, double area, int el_m, int SrcID)
            : base(new double[8] { 0, 0, 0, 0, 0, 0, 0, 0 }, new Hare.Geometry.Point(0, 0, 0), SrcID)
            {
                T = T_in;

                samplespermeter = el_m;
                Samples = samples;

                for (int i = 0; i < Samples.Length; i++)
                {
                    DomainLevel = Utilities.PachTools.DecodeSourcePower(CodeList);
                    DomainPower = new double[8];
                    double PowerMod = area;
                    for (int oct = 0; oct < 8; oct++) DomainPower[oct] = 1E-12 * Math.Pow(10, .1 * DomainLevel[oct]) / PowerMod;
                }
            }

            public override string Type()
            {
                return "Special - Surface Source";
            }

            public override void AppendPts(ref List<Hare.Geometry.Point> SPT)
            {                
                    for (int j = 0; j < Samples.Length; j++)
                    {
                        SPT.Add(Samples[j]);
                    }
            }

            public override BroadRay Directions(int index, int thread, ref Random random, int[] Octaves)
            {
                BroadRay Ray = Directions(index, thread, ref random);
                Ray.Octaves = Octaves;
                return Ray;
            }

            public override BroadRay Directions(int index, int thread, ref Random random)
            {
                int i = (int)(random.Next() * (double)T.Polygon_Count);
                Point P = T.Polys[i].GetRandomPoint(random.NextDouble(), random.NextDouble(), 0);
                double Theta = random.NextDouble() * 2 * System.Math.PI;
                double Phi = random.NextDouble() * 2 * System.Math.PI;
                Hare.Geometry.Vector Direction = new Hare.Geometry.Vector(Math.Sin(Theta) * Math.Cos(Phi), Math.Sin(Theta) * Math.Sin(Phi), Math.Cos(Theta));
                
                return new BroadRay(P, Direction, random.Next(), thread, DomainPower, 0, S_ID);
            }
        }
    }
}