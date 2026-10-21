//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2026, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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
using System.Collections.Concurrent;
using System.Threading;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        /// <summary>
        /// One geometric ray carrying the energies in a contiguous range of spectral bands.
        /// </summary>
        public class OctaveRay: Ray
        {
            public double Intensity;
            public int BandStart;
            public int BandCount = 1;
            public bool ThirdOctave { get { return BandCount == 3; } }
            public double[] Fractions;
            // Intensity is the total; Fractions holds the pooled spectral shares.
            public int Octave;
            public double t_sum;
            public int Surf_ID;
            public int Source_ID;
            public double Decimation_threshold;
            public double Decim_Inv;
            public short Scatter_Mask;
            public OctaveRay()
            : base(new Point(0,0,0), new Vector(0,0,0), 0, 0)
            { }

            public OctaveRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, int srf_id = -1)
            : base(x, y, z, dx, dy, dz, ThreadID_IN, ID)
            {
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                BandStart = octave; BandCount = 1;
                Scatter_Mask = 0;
                Source_ID = SrcID;
                Surf_ID = srf_id;
                Decimation = Decimation_in;
            }

            public OctaveRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, double Decimation_Inv, int srf_id = -1)
            : base(x, y, z, dx, dy, dz, ThreadID_IN, ID)
            {
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                BandStart = octave; BandCount = 1;
                Scatter_Mask = 0;
                Source_ID = SrcID;
                Surf_ID = srf_id;
                Decimation_threshold = Decimation_in;
                Decim_Inv = Decimation_Inv;
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, int srf_id = -1)
            {
                poly_origin1 = poly_origin2 = 0;
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                ThreadID = ThreadID_IN;
                Ray_ID = ID;
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                BandStart = octave; BandCount = 1;
                Scatter_Mask = 0;
                Source_ID = SrcID;
                Surf_ID = srf_id;
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, int srf_id = -1)
            {
                poly_origin1 = poly_origin2 = 0;
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                ThreadID = ThreadID_IN;
                Ray_ID = ID;
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                BandStart = octave; BandCount = 1;
                Scatter_Mask = 0;
                Source_ID = SrcID;
                Surf_ID = srf_id;
                Decimation = Decimation_in;
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, double Decimation_Inv, int srf_id = -1)
            {
                poly_origin1 = poly_origin2 = 0;
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                ThreadID = ThreadID_IN;
                Ray_ID = ID;
                t_sum = Time;
                Octave = octave;
                Intensity = Intensity_in;
                BandStart = octave; BandCount = 1;
                Scatter_Mask = 0;
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
                //OctaveRay O = new OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, Intensity * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation_threshold, Decim_Inv, Surf_ID);
                OctaveRay O = OctaveRayPool.Instance.new_OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, Intensity * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation_threshold, Decim_Inv, Scatter_Mask, Surf_ID);
                O.BandStart = BandStart; O.BandCount = BandCount;
                if (BandCount > 1)
                {
                    if (O.Fractions == null || O.Fractions.Length < BandCount) O.Fractions = new double[BandCount];
                    System.Array.Copy(Fractions, O.Fractions, BandCount);
                }
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
                //OctaveRay O = new OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, Intensity * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation_threshold, Decim_Inv, Surf_ID);
                OctaveRay O = OctaveRayPool.Instance.new_OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, Intensity * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation_threshold, Decim_Inv, Scatter_Mask, Surf_ID);
                O.BandStart = BandStart; O.BandCount = BandCount;
                if (BandCount > 1)
                {
                    if (O.Fractions == null || O.Fractions.Length < BandCount) O.Fractions = new double[BandCount];
                    System.Array.Copy(Fractions, O.Fractions, BandCount);
                }
                this.Intensity *= (1 - E_Mod_Coef);
                return O;
            }

            public OctaveRay SplitRay(double[] coefficients)
            {
                OctaveRay branch = Clone();
                if (branch.BandCount == 1) branch.Intensity *= coefficients[branch.BandStart];
                else
                {
                    double total = 0;
                    for (int i = 0; i < branch.BandCount; i++) { branch.Fractions[i] = branch.Intensity * branch.Fractions[i] * coefficients[branch.BandStart + i]; total += branch.Fractions[i]; }
                    branch.Intensity = total;
                    for (int i = 0; i < branch.BandCount; i++) branch.Fractions[i] = total > 0 ? branch.Fractions[i] / total : 1.0 / branch.BandCount;
                }
                if (BandCount == 1) Intensity *= (1 - coefficients[BandStart]);
                else
                {
                    double total = 0;
                    for (int i = 0; i < BandCount; i++) { Fractions[i] = Intensity * Fractions[i] * (1 - coefficients[BandStart + i]); total += Fractions[i]; }
                    Intensity = total;
                    for (int i = 0; i < BandCount; i++) Fractions[i] = total > 0 ? Fractions[i] / total : 1.0 / BandCount;
                }
                return branch;
            }
            public OctaveRay SplitRay(double fraction0, double fraction1, double fraction2)
            {
                if (BandCount == 1) return SplitRay(fraction0);
                OctaveRay branch = Clone();
                if (BandCount != 3) throw new System.InvalidOperationException("Use the indexed coefficient spectrum for this band range.");
                double branchTotal = 0, remainingTotal = 0;
                for (int i = 0; i < BandCount; i++)
                {
                    double fraction = i == 0 ? fraction0 : i == 1 ? fraction1 : fraction2;
                    branch.Fractions[i] = branch.Intensity * branch.Fractions[i] * fraction; branchTotal += branch.Fractions[i];
                    Fractions[i] = Intensity * Fractions[i] * (1 - fraction); remainingTotal += Fractions[i];
                }
                branch.Intensity = branchTotal; Intensity = remainingTotal;
                for (int i = 0; i < BandCount; i++)
                {
                    branch.Fractions[i] = branchTotal > 0 ? branch.Fractions[i] / branchTotal : 1.0 / BandCount;
                    Fractions[i] = remainingTotal > 0 ? Fractions[i] / remainingTotal : 1.0 / BandCount;
                }
                return branch;
            }

            public OctaveRay Clone()
            {
                //OctaveRay O = new OctaveRay(x, y, z, dx, dy, dz, this.Ray_ID, this.ThreadID, this.Intensity, this.t_sum, this.Octave, this.Source_ID, this.Decimation_threshold, this.Decim_Inv, this.Surf_ID);
                OctaveRay O = OctaveRayPool.Instance.new_OctaveRay(x, y, z, dx, dy, dz, this.Ray_ID, this.ThreadID, this.Intensity, this.t_sum, this.Octave, this.Source_ID, this.Decimation_threshold, this.Decim_Inv, this.Scatter_Mask, this.Surf_ID);
                O.BandStart = BandStart; O.BandCount = BandCount;
                if (BandCount > 1)
                {
                    if (O.Fractions == null || O.Fractions.Length < BandCount) O.Fractions = new double[BandCount];
                    System.Array.Copy(Fractions, O.Fractions, BandCount);
                }
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

            public void SetScattered()
            {
                Scatter_Mask |= 1 << 0;    // Set first bit to 1 to indicate current state is Scattered
                Scatter_Mask |= 1 << 1;    // Set second bit to indicate it was ever Scattered
            }

            public void SetSpecular()
            {
                Scatter_Mask &= ~(1 << 0); // Clear first bit to 0 to indicate current state is Specular
            }

            public bool IsSpecular
            {
                get { return (Scatter_Mask & 1) == 0; } // Check if first bit is 0
            }

            public bool WasScattered
            {
                get { return (Scatter_Mask & 2) == 2; }// Check if second bit is 1
            }
        }

        /// <summary>
        /// Broadband ray carrying indexed spectral energies.
        /// </summary>
        public class BroadRay : Ray
        {
            public double[] Energy = new double[8];
            public double t_sum;
            public int Surf_ID;
            public int Source_ID;
            public double[] Decimation = new double[8];
            // Geometric octave groups (0 through 7), independent of spectral resolution.
            public int[] Freq_Bands = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
            //public int[] Octaves = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
            //public int[] Third_Octaves = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };

            public BroadRay()
                :base(new Point(0,0,0),new Vector(0,0,0), 0, 0)
            {
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID)
            {
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                this.ThreadID = ThreadID_IN;
                this.Ray_ID = ID;
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Surf_ID = -1;
                Decimation = new double[8];
                Freq_Bands = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
                poly_origin1 = poly_origin2 = 0;
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID, int[] _Octaves)
            {
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                this.ThreadID = ThreadID_IN;
                this.Ray_ID = ID;
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Surf_ID = -1;
                Decimation = new double[8];
                poly_origin1 = poly_origin2 = 0;
                Freq_Bands = _Octaves;
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, double[] Decim, int SrcID)
            {
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                this.ThreadID = ThreadID_IN;
                this.Ray_ID = ID;
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Decimation = Decim;
                Surf_ID = -1;
                Freq_Bands = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
                poly_origin1 = poly_origin2 = 0;
            }

            public void reassign(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, double[] Decim, int SrcID, int[] _Octaves)
            {
                x = StartPt.x; y = StartPt.y; z = StartPt.z;
                dx = Direction.dx; dy = Direction.dy; dz = Direction.dz;
                this.ThreadID = ThreadID_IN;
                this.Ray_ID = ID;
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Decimation = Decim;
                Surf_ID = -1;
                poly_origin1 = poly_origin2 = 0;
                Freq_Bands = _Octaves;
            }

            public BroadRay(Point StartPt, Vector Direction, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID)
                : base(StartPt, Direction, ThreadID_IN, ID)
            {
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy,0);
                Source_ID = SrcID;
                Surf_ID = -1;
            }

            public BroadRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID)
            : base(x, y, z, dx, dy, dz, ThreadID_IN, ID)
            {
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Surf_ID = -1;
            }

            public BroadRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double[] energy_in, double time, int SrcID, int[] _Octaves)
            : base(x, y, z, dx, dy, dz, ThreadID_IN, ID)
            {
                t_sum = time;
                Energy = new double[energy_in.Length];
                energy_in.CopyTo(Energy, 0);
                Source_ID = SrcID;
                Surf_ID = -1;
                Freq_Bands = _Octaves;
            }

            //public BroadRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double[] energy_in, double time, double[] Decim, int SrcID)
            //: base(x, y, z, dx, dy, dz, ThreadID_IN, ID)
            //{
            //    t_sum = time;
            //    Energy = new double[8];
            //    energy_in.CopyTo(Energy, 0);
            //    Source_ID = SrcID;
            //    Decimation = Decim;
            //    Surf_ID = -1;
            //}

            //public BroadRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double[] energy_in, double time, double[] Decim, int SrcID, int[] _Octaves)
            //: base(x,y,z,dx,dy,dz, ThreadID_IN, ID)
            //{
            //    t_sum = time;
            //    Energy = new double[8];
            //    energy_in.CopyTo(Energy, 0);
            //    Source_ID = SrcID;
            //    Decimation = Decim;
            //    Surf_ID = -1;
            //    Freq_Bands = _Octaves;
            //}

            /// <summary>
            /// Appends the distnace of the last ray traversal to the total distance the ray has traveled.
            /// </summary>
            /// <param name="Length"></param>
            public void AddLeg(double time)
            {
                t_sum += time;
            }

            public int BandCount { get { return Energy.Length; } }
            public int BandsPerGroup { get { return BandCount / Decimation.Length; } }
            public int BandsPerOctave { get { return BandsPerGroup; } }
            public bool ThirdOctave { get { return BandsPerOctave == 3; } }
            public double OctaveEnergy(int octave)
            {
                int bandStart = octave * BandsPerGroup;
                double energy = Energy[bandStart];
                for (int component = 1; component < BandsPerGroup; component++) energy += Energy[bandStart + component];
                return energy;
            }

            public BroadRay Clone()
            {
                BroadRay ray = BroadRayPool.Instance.new_BroadRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, (double[])Energy.Clone(), t_sum, Source_ID, (int[])Freq_Bands.Clone(), Surf_ID);
                ray.Decimation = (double[])Decimation.Clone();
                return ray;
            }

            /// <summary>
            /// Creates an octave ray for the specified octave. After this operation, the Broadray will no longer trace this octave.
            /// </summary>
            /// <param name="Octave"></param>
            /// <returns></returns>
            public OctaveRay SplitRay(int Octave)
            {
                int bandStart = Octave * BandsPerGroup;
                OctaveRay O = OctaveRayPool.Instance.new_OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, OctaveEnergy(Octave), t_sum, Octave, Source_ID, Decimation[Octave], 0, Surf_ID);
                O.BandStart = bandStart; O.BandCount = BandsPerGroup;
                if (O.BandCount > 1)
                {
                    if (O.Fractions == null || O.Fractions.Length < O.BandCount) O.Fractions = new double[O.BandCount];
                    for (int i = 0; i < O.BandCount; i++) O.Fractions[i] = O.Intensity > 0 ? Energy[bandStart + i] / O.Intensity : 1.0 / O.BandCount;
                }
                for (int band = bandStart; band < bandStart + BandsPerGroup; band++) Energy[band] = 0;
                return O;
            }

            /// <summary>
            /// Creats an octave ray for the specified octave with the percentage of the Broadray's energy specified by the scattering coeficient.
            /// </summary>
            /// <param name="Octave"></param>
            /// <param name="Scattering_Coef"></param>
            /// <returns></returns>
            public OctaveRay SplitRay(int group, double[] coefficients)
            {
                int bandStart = group * BandsPerGroup;
                OctaveRay branch = OctaveRayPool.Instance.new_OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, 0, t_sum, group, Source_ID, Decimation[group], 0, Surf_ID);
                branch.BandStart = bandStart; branch.BandCount = BandsPerGroup;
                if (branch.BandCount == 1) branch.Intensity = Energy[bandStart] * coefficients[bandStart];
                else
                {
                    if (branch.Fractions == null || branch.Fractions.Length < branch.BandCount) branch.Fractions = new double[branch.BandCount];
                    branch.Intensity = 0;
                    for (int i = 0; i < branch.BandCount; i++) { branch.Fractions[i] = Energy[bandStart + i] * coefficients[bandStart + i]; branch.Intensity += branch.Fractions[i]; }
                    for (int i = 0; i < branch.BandCount; i++) branch.Fractions[i] = branch.Intensity > 0 ? branch.Fractions[i] / branch.Intensity : 1.0 / branch.BandCount;
                }
                for (int band = bandStart; band < bandStart + BandsPerGroup; band++) Energy[band] *= (1 - coefficients[band]);
                return branch;
            }

            public OctaveRay SplitRay(int Octave, double E_Mod_Coef)
            {
                int bandStart = Octave * BandsPerGroup;
                OctaveRay ray = OctaveRayPool.Instance.new_OctaveRay(x, y, z, dx, dy, dz, Ray_ID, ThreadID, OctaveEnergy(Octave) * E_Mod_Coef, t_sum, Octave, Source_ID, Decimation[Octave], 0, Surf_ID);
                ray.BandStart = bandStart; ray.BandCount = BandsPerGroup;
                if (ray.BandCount > 1)
                {
                    if (ray.Fractions == null || ray.Fractions.Length < ray.BandCount) ray.Fractions = new double[ray.BandCount];
                    ray.Intensity = 0;
                    for (int i = 0; i < ray.BandCount; i++) { ray.Fractions[i] = Energy[bandStart + i] * E_Mod_Coef; ray.Intensity += ray.Fractions[i]; }
                    for (int i = 0; i < ray.BandCount; i++) ray.Fractions[i] = ray.Intensity > 0 ? ray.Fractions[i] / ray.Intensity : 1.0 / ray.BandCount;
                }
                for (int band = bandStart; band < bandStart + BandsPerGroup; band++) Energy[band] *= (1 - E_Mod_Coef);
                return ray;
            }
        }

        public class BroadRayPool
        {
            protected int count;
            protected int position = 0;
            private SemaphoreSlim ctr;
            private static BroadRayPool instance = null;
            protected object insLock = new object();
            protected bool initialized;
            public static void Initialize()
            {
                instance = new BroadRayPool();
                instance.initialized = true;
            }
            public static BroadRayPool Instance
            {
                get
                {
                    if (instance == null)
                    {
                        instance = new BroadRayPool();
                    }
                    return instance;
                }
            }

            private new BroadRay[] Pool;
            object lockpos = new object();

            private BroadRayPool(int capacity = 200000) 
            {
                instance = this;
                count = capacity; 
                ctr = new SemaphoreSlim(0, count);
                Pool = new BroadRay[capacity];
                for(int i = 0; i < Pool.Length; i++) { Pool[i] = new BroadRay(); }
                ctr.Release(count);
            }

            public BroadRay new_BroadRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double[] Intensity_in, double Time, int SrcID, int srf_id = -1)
            {
                BroadRay BR = pull;
                BR.x = x; BR.y = y; BR.z = z;
                BR.dx = dx; BR.dy = dy; BR.dz = dz;
                BR.Ray_ID = ID;
                BR.ThreadID = ThreadID_IN;
                BR.Energy = Intensity_in;
                BR.t_sum = Time;
                BR.Freq_Bands = new int[8] {0,1,2,3,4,5,6,7};
                BR.Source_ID = SrcID;
                BR.Surf_ID = srf_id;
                return BR;
            }

            public BroadRay new_BroadRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double[] Intensity_in, double Time, int SrcID, int[] _Octaves, int srf_id = -1)
            {
                BroadRay BR = pull;
                BR.x = x; BR.y = y; BR.z = z;
                BR.dx = dx; BR.dy = dy; BR.dz = dz;
                BR.Ray_ID = ID;
                BR.ThreadID = ThreadID_IN;
                BR.Energy = Intensity_in;
                BR.t_sum = Time;
                BR.Freq_Bands = _Octaves;
                BR.Source_ID = SrcID;
                BR.Surf_ID = srf_id; 
                return BR;
            }

            private int increment
            {
                get 
                {
                    int p;
                    lock(lockpos)
                    {
                        position++;
                        if (position >= count) position = 0;
                        p = position;
                    }
                    return p;
                }
            }

            public BroadRay pull
            {
                get
                {
                    ctr.Wait();
                    BroadRay R;
                    int p = increment;
                    //R = Pool[p] == null ? new BroadRay() : Pool[p];
                    R = Pool[p];
                    R.dx = 0;
                    R.dy = 0;
                    R.dz = 0;
                    R.x = 0;
                    R.y = 0;
                    R.z = 0;
                    R.Energy = new double[8];
                    R.Decimation = new double[8]; R.Freq_Bands = new int[0];
                    R.t_sum = 0;
                    R.Ray_ID = 0;
                    R.Surf_ID = -1;
                    R.Source_ID = 0;
                    R.poly_origin1 = 0;
                    R.poly_origin2 = 0;
                    R.ThreadID = 0;
                    return R;
                }
            }

            public void release()
            {
                ctr.Release();
            }
        }

        public static class GeometryPool
        {
            public static BroadRayPool MainPool = BroadRayPool.Instance;
            public static OctaveRayPool OctPool = OctaveRayPool.Instance;
        }

        public class OctaveRayPool
        {
            protected int count;
            private SemaphoreSlim ctr;
            private SemaphoreSlim memmod;
            private static OctaveRayPool instance = null;
            protected bool initialized;
            private ConcurrentBag<OctaveRay> Pool;
            public static void Initialize()
            {
                instance = new OctaveRayPool();
                instance.initialized = true;
            }

            private OctaveRayPool(int capacity = 5000)
            {
                instance = this;
                count = capacity;
                ctr = new SemaphoreSlim(0,10000000);
                memmod = new SemaphoreSlim(1, 1);
                Pool = new ConcurrentBag<OctaveRay>();
                for (int i = 0; i < capacity; i++) { Pool.Add(new OctaveRay()); }
                ctr.Release(count);
            }
            public static OctaveRayPool Instance
            {
                get
                {
                    if (instance == null)
                    {
                        instance = new OctaveRayPool();
                    }
                    return instance;
                }
            }            

            public OctaveRay new_OctaveRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, double Decimation_Inv, short ScatterMask, int srf_id = -1)
            {
                OctaveRay OR = pull;
                OR.x = x; OR.y = y; OR.z = z;
                OR.dx = dx; OR.dy = dy; OR.dz = dz;
                OR.Ray_ID = ID;
                OR.ThreadID = ThreadID_IN;
                OR.Intensity = Intensity_in;
                OR.t_sum = Time;
                OR.Octave = octave;
                OR.BandStart = octave; OR.BandCount = 1;
                OR.Source_ID = SrcID;
                OR.Decimation = Decimation_in;
                OR.Decim_Inv = Decimation_Inv;
                OR.Scatter_Mask = ScatterMask;
                OR.Surf_ID = srf_id;
                return OR;
            }

            public OctaveRay new_OctaveRay(double x, double y, double z, double dx, double dy, double dz, int ID, int ThreadID_IN, double Intensity_in, double Time, int octave, int SrcID, double Decimation_in, short ScatterMask, int srf_id = -1)
            {
                OctaveRay OR = pull;
                OR.x = x; OR.y = y; OR.z = z;
                OR.dx = dx; OR.dy = dy; OR.dz = dz;
                OR.Ray_ID = ID;
                OR.ThreadID = ThreadID_IN;
                OR.Intensity = Intensity_in;
                OR.t_sum = Time;
                OR.Octave = octave;
                OR.BandStart = octave; OR.BandCount = 1;
                OR.Source_ID = SrcID;
                OR.Decimation = Decimation_in;
                OR.Decim_Inv = 1.0 / Decimation_in;
                OR.Scatter_Mask = ScatterMask;
                OR.Surf_ID = srf_id;
                return OR;
            }

            public OctaveRay pull
            {
                get
                {
                    if (!ctr.Wait(0))
                    {
                        // If the pool is empty, expand the pool...                        
                        if (memmod.Wait(0))
                        {
                            do
                            {
                                if (count + 100 < 10000000)
                                {
                                    for (int i = 0; i < 100; i++) Pool.Add(new OctaveRay());
                                    count += 100;
                                    ctr.Release(100);
                                }
                            } while (ctr.CurrentCount == 0);
                            memmod.Release();
                        }
                        ctr.Wait();
                    }

                    OctaveRay R;
                    Pool.TryTake(out R);
                    R.x = 0;
                    R.y = 0;
                    R.z = 0;
                    R.dx = 0;
                    R.dy = 0;
                    R.dz = 0;
                    R.Octave = 0;
                    R.Intensity = 0;
                    R.BandStart = 0; R.BandCount = 1;
                    R.Scatter_Mask = 0; R.Decimation_threshold = R.Decim_Inv = 0;
                    R.poly_origin1 = 0;
                    R.poly_origin2 = 0;
                    R.Ray_ID = 0;
                    R.Surf_ID = -1;
                    R.t_sum = 0;
                    R.Source_ID = 0;
                    R.ThreadID = 0;
                    return R;
                }
            }

            public void release(OctaveRay R)
            {
                Pool.Add(R);
                ctr.Release();
            }
        }
    }
}