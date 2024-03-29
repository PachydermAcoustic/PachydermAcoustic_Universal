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

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        /// <summary>
        /// Ray for a single octave band.
        /// </summary>
        public class OctaveRay: Ray
        {
            public double Intensity;
            public int Octave;
            public double t_sum;
            public int Surf_ID;
            public int Source_ID;
            public double Decimation_threshold;
            public double Decim_Inv;

            public OctaveRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, int srf_id = -1)
            : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                Source_ID = SrcID;
                Surf_ID = srf_id;
            }

            public OctaveRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave,  int SrcID, double Decimation_in, int srf_id = -1)
                : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                Source_ID = SrcID;
                Surf_ID = srf_id;
                Decimation = Decimation_in;
            }

            public OctaveRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, double Decimation_Inv, int srf_id = -1)
            : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                Source_ID = SrcID;
                Surf_ID = srf_id;
                Decimation_threshold = Decimation_in;
                Decim_Inv = Decimation_Inv;
            }


            /// <summary>
            /// Appends the distance of the last ray traversal to the total distance the ray has traveled.
            /// </summary>
            /// <param name="Length"></param>
            public void AddLeg(double t)
            {
                t_sum += t;
            }

            /// <summary>
            /// Creates an octave ray with the percentage of the parent ray's energy specified by the scattering coeficient. Duely modifies the parent ray.
            /// </summary>
            /// <param name="E_Mod_Coef"></param>
            /// <returns></returns>
            public OctaveRay SplitRay( double E_Mod_Coef)
            {
                OctaveRay O = new OctaveRay(origin, direction, Ray_ID, ThreadID, Intensity * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation_threshold, Decim_Inv, Surf_ID);
                this.Intensity *= (1 - E_Mod_Coef);
                return O;
            }

            /// <summary>
            /// Creats an octave ray with the percentage of the parent ray's energy specified by the scattering coeficient. Duely modifies the parent ray.
            /// </summary>
            /// <param name="E_Mod_Coef"></param>
            /// <param name="phase_delay"></param>
            /// <returns></returns>
            public OctaveRay SplitRay(double E_Mod_Coef, double phase_delay)
            {
                OctaveRay O = new OctaveRay(origin, direction, Ray_ID, ThreadID, Intensity * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation_threshold, Decim_Inv, Surf_ID);
                this.Intensity *= (1 - E_Mod_Coef);
                return O;
            }

            public OctaveRay Clone()
            {
                OctaveRay O = new OctaveRay(this.origin, this.direction, this.Ray_ID, this.ThreadID, this.Intensity, this.t_sum, this.Octave, this.Source_ID, this.Decimation_threshold, this.Decim_Inv, this.Surf_ID);
                return O;
            }

            public double Decimation
            {
                set 
                {
                    Decimation_threshold = value;
                    Decim_Inv = 1.0 / value;
                }
            }
        }

        /// <summary>
        /// Ray for all octave bands.
        /// </summary>
        public class BroadRay : Ray
        {
            public double[] Energy;
            public double t_sum;
            public int Surf_ID;
            public int Source_ID;
            public int[] Octaves = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
            public double[] Decimation = new double[8];
            public int[] Freq_Bands = new int[] { };
            //public int[] Octaves = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
            //public int[] Third_Octaves = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };

            public BroadRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID)
                : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum += time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy,0);
                Source_ID = SrcID;
                Surf_ID = -1;
            }

            public BroadRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID, int[] _Octaves)
            : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum += time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Surf_ID = -1;
                Octaves = _Octaves;
            }

            public BroadRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, double[] Decim, int SrcID)
            : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum += time;
                Energy = new double[8];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Decimation = Decim;
                Surf_ID = -1;
            }

            public BroadRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, double[] Decim, int SrcID, int[] _Octaves)
            : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum += time;
                Energy = new double[8];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Decimation = Decim;
                Surf_ID = -1;
                Freq_Bands = _Octaves;
            }


            /// <summary>
            /// Appends the distnace of the last ray traversal to the total distance the ray has traveled.
            /// </summary>
            /// <param name="Length"></param>
            public void AddLeg(double time)
            {
                t_sum += time;
            }

            public BroadRay Clone()
            {
                return new BroadRay(origin, direction, Ray_ID, ThreadID, Energy, t_sum, Source_ID);
            }

            /// <summary>
            /// Creates an octave ray for the specified octave. After this operation, the Broadray will no longer trace this octave.
            /// </summary>
            /// <param name="Octave"></param>
            /// <returns></returns>
            public OctaveRay SplitRay(int Octave)
            {
                OctaveRay O = new OctaveRay(origin, direction, Ray_ID, ThreadID, Energy[Octave], t_sum, Octave, Source_ID, Decimation[Octave], Surf_ID);
                this.Energy[Octave] = 0;
                return O;
            }

            /// <summary>
            /// Creats an octave ray for the specified octave with the percentage of the Broadray's energy specified by the scattering coeficient.
            /// </summary>
            /// <param name="Octave"></param>
            /// <param name="Scattering_Coef"></param>
            /// <returns></returns>
            public OctaveRay SplitRay(int Octave, double E_Mod_Coef)
            {
                OctaveRay O = new OctaveRay(origin, direction, Ray_ID, ThreadID, Energy[Octave] * E_Mod_Coef,  t_sum, Octave, Source_ID, Decimation[Octave], Surf_ID);
                this.Energy[Octave] *= (1 - E_Mod_Coef);
                return O;
            }
        }
    }
}