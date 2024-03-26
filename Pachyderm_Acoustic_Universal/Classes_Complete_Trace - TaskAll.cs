//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2024, Arthur van der Harten 
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
using System.Linq;
using System.Threading;
using System.Collections.Concurrent;
using Pachyderm_Acoustic.Utilities;
using System.Threading.Tasks;
using Pachyderm_Acoustic.Pach_Graphics;
using System.Threading.Tasks.Sources;

namespace Pachyderm_Acoustic
{
    /// <summary>
    /// This simulation type calculates the energy time curve of a room by ray tracing.
    /// </summary>
    public class SplitRayTracer : Simulation_Type
    {
        protected Scene Room;
        protected Source Source;
        protected Receiver_Bank RecMain;
        protected int Raycount;
        protected double COTime;
        public int[] _currentRay = new int[1];
        private DateTime _st;
        private double _lost;
        private double _eInit;
        private double _lost_percent;
        private int IS_Order;
        private int _processorCt;
        public int[] _rayTotal = new int[1];
        public TimeSpan _ts;
        private int h_oct;
        private int[] _octaves;
        private Convergence_Check[] check;
        private bool conclude = false;
        I_Conv_Progress Vis_Feedback;

        /// <summary>
        /// Constructor for the general case ray tracer.
        /// </summary>
        /// <param name="sourceIn"></param>
        /// <param name="receiverIn"></param>
        /// <param name="roomIn"></param>
        /// <param name="cutOffLengthIn"></param>
        /// <param name="rayCountIn"></param>
        /// <param name="octaveRange">Two values - the lowest octave to be calculated and the highest octave to be calculated - 0 being 62.5, and 7 being 8k.</param>
        /// <param name="isOrderIn">The highest order for which image source was calcualted. If no Image Source, then enter 0.</param>
        /// <param name="partitionedReceiver">Is the receiver partitioned... i.e., did you use a mapping receiver bank?</param>
        public SplitRayTracer(Source sourceIn, Receiver_Bank receiverIn, Scene roomIn, double cutoffTime, int[] octaveRange, int isOrderIn, int rayCountIn, I_Conv_Progress Vis_Feedback)
        {
            BroadRayPool.Initialize();
            OctaveRayPool.Initialize();
            IS_Order = isOrderIn;
            Room = roomIn;
            RecMain = receiverIn;
            Raycount = rayCountIn;
            _octaves = new int[octaveRange[1] - octaveRange[0] + 1];
            for (int o = octaveRange[0]; o <= octaveRange[1]; o++) _octaves[o - octaveRange[0]] = o;
            COTime = cutoffTime;
            Source = sourceIn;
            double[] totalAbs = new double[8];
            this.Vis_Feedback = Vis_Feedback;

            for (int s = 0; s < Room.Count(); s++)
            {
                double area = Room.SurfaceArea(s);
                foreach (int o in _octaves)
                {
                    totalAbs[o] += area * Room.AbsorptionValue[s].Coefficient_A_Broad(o);
                }
            }

            double min = double.PositiveInfinity;
            double octpower = 0;
            foreach (int o in _octaves) if (Source.SoundPower[o] > octpower) { octpower = Source.SoundPower[o]; }

            foreach (int o in _octaves)
            {
                if (Source.SoundPower[o] < octpower * 1E-6) continue;
                if (min > totalAbs[o])
                {
                    min = totalAbs[o];
                    h_oct = o;
                }
            }

            if (Raycount < 1)
            {
                List<Point> spt = new List<Point>();
                if (!(Source is LineSource) && !(Source is SurfaceSource))
                {
                    spt.Add(Source.Origin);
                }
                else if (Source is LineSource)
                {
                    foreach (Point p in (Source as LineSource).Samples)
                    {
                        spt.Add(p);
                    }
                }
                else
                {
                    foreach (Point p in (Source as SurfaceSource).Samples)
                    {
                        spt.Add(p);
                    }
                }
                Random r = new Random();

                double maxT_T = 0;
                int T_id = -1;
                double maxT_F = 0;
                int F_id = -1;
                double t;
                for (int i = 0; i < RecMain.Rec_List.Length; i++)
                {
                    Hare.Geometry.Point s = new Point();
                    double d = double.MaxValue;
                    foreach (Point p in spt)
                    {
                        double dt = (p - RecMain.Rec_List[i].Origin).Length();
                        if (dt < d)
                        {
                            s = p;
                            d = dt;
                        }
                    }
                    if (Check_Validity(s, i, r.Next(), out t))
                    {
                        if (maxT_T < t)
                        {
                            T_id = i;
                            maxT_T = t;
                        }
                    }
                    else
                    {
                        if (maxT_F < t)
                        {
                            F_id = i;
                            maxT_F = t;
                        }
                    }
                }
                check = (rayCountIn < 0) ? new Convergence_Check[2] { T_id < 0 ? null : new Minimum_Convergence_Check(this.Source, this.Room, receiverIn, T_id, h_oct, 0, Vis_Feedback), F_id < 0 ? null : new Minimum_Convergence_Check(this.Source, this.Room, receiverIn, F_id, h_oct, 1, Vis_Feedback) }
                : new Convergence_Check[2] { T_id < 0 ? null : new Detailed_Convergence_Check(receiverIn, T_id, h_oct, 0, Vis_Feedback), F_id < 0 ? null : new Detailed_Convergence_Check(receiverIn, F_id, h_oct, 1, Vis_Feedback) };
            }
            else check = null;
        }

        private bool Check_Validity(Point p, int rec_id, int rnd, out double dist)
        {
            Hare.Geometry.Vector d = RecMain.Origin(rec_id) - p;
            dist = d.Length();
            d.Normalize();
            Ray R = new Ray(Source.Origin, d, 0, rnd);
            double x1 = 0, x2 = 0;
            double t;
            int x3 = 0;
            Point x4;

            while (true)
            {
                if (Room.shoot(R, out x1, out x2, out x3, out x4, out t))
                {
                    return (t >= dist);
                }
                else
                {
                    return true;
                }
            }
        }

        ValueTask Prime;
        List<ValueTask> MainRays = new List<ValueTask>();

        /// <summary>
        /// Inherited member. Divides the simulation into tasks, and begins.
        /// </summary>
        public override void Begin()
        {
            _st = DateTime.Now;
            Random Rnd = new Random((int)DateTime.Now.ToFileTimeUtc());
            _processorCt = Pach_Properties.Instance.ProcessorCount();

            if (Raycount > 0)
            {
                for (int P_I = 0; P_I < _processorCt; P_I++)
                {
                    int start = (int)Math.Floor((double)P_I * Raycount / _processorCt);
                    int end;
                    if (P_I == _processorCt - 1) end = Raycount;
                    else end = (int)Math.Floor((double)(P_I + 1) * Raycount / _processorCt);
                    Calc_Params T = new Calc_Params(Room, start, end, P_I, Rnd.Next());
                    ValueTask t = new ValueTask(Task.Run(() => Calculate(T)));
                    MainRays.Add(t);
                }

                Task.Run(async () => 
                {   
                    while (_currentRay[0] < Raycount) { await Task.Delay(1500); }
                    for (int i = 0; i < MainRays.Count; i++) if (MainRays[i].IsCompletedSuccessfully);
                    conclude = true;
                    _ts = DateTime.Now - _st;
                });

            }
            else
            {
                Prime = new ValueTask(Task.Run(() => Calculate_Conv(_processorCt, Rnd.Next())));
            }
        }

        /// <summary>
        /// Called by Pach_RunSim_Command. Indicates whether or not the simulation has completed.
        /// </summary>
        /// <returns>Returns running if any threads in this simulation are still running. Returns stopped if all have stopped.</returns>
        public override System.Threading.ThreadState ThreadState()
        {
            if (!conclude) return System.Threading.ThreadState.Running;
            return System.Threading.ThreadState.Stopped;
        }

        public void Conclude_Simulation()
        {
            conclude = true;
        }

        object ctlock = new object();

        /// <summary>
        /// Called by each thread from the begin method for calculations with a user-defined length.
        /// </summary>
        /// <param name="i">the object is type "Calc_Params" which holds all the necessary information to run a portion of the simulation.</param>
        public async void Calculate(object i)
        {
            ///Homogeneous media only...
            Calc_Params Params = (Calc_Params)i;
            Random Rnd = new Random(Params.RandomSeed);
            //LastRays = new List<Task>[_processorCt];

            for (int ray = 0; ray < Params.EndIndex - Params.StartIndex; ray++)
            {
                await cast_ray_distributed(Params.Room, Params.ThreadID, Rnd);
                lock (ctlock)
                {
                    _currentRay[0]++;
                }
            }
        }

        /// <summary>
        /// Called by each thread from the begin method.
        /// </summary>
        /// <param name="i">the object is type "Calc_Params" which holds all the necessary information to run a portion of the simulation.</param>
        public async void Calculate_Conv(object i, object seed)
        {
            ///Homogeneous media only...
            Random rnd = new Random((int)seed);
            _st = DateTime.Now;
            bool conv = true;

            int rev = 0;

            do
            {
                rev++;
                if (conclude)
                    break;

                for (int P_I = 0; P_I < 80; P_I++)
                {    
                    Calc_Params T = new Calc_Params(Room, P_I, P_I + 1, P_I, rnd.Next());
                    ValueTask t = new ValueTask(Task.Run(async () =>
                    {
                        Random Rnd = new Random(T.RandomSeed);
                        await cast_ray_distributed(T.Room, T.ThreadID, Rnd);
                        lock (ctlock)
                        {
                            _currentRay[0]++;
                        }                
                    }));
                    MainRays.Add(t);
                }

                for (var j = 0; j < MainRays.Count; j++) await MainRays[j];

                GC.Collect();
                conv = true;
                foreach (Convergence_Check c in check) if (c != null) conv &= c.Check();
            } while (!conv);

            Conclude_Simulation();

            _ts = DateTime.Now - _st;
        }

        double thr2_template = Math.Pow(10, -4.5);

        object lostlock = new object();
        object initlock = new object();
        private async Task cast_ray_distributed(Scene Room, int ThreadID, Random Rnd)
        {
            List<ValueTask> LastRays = new List<ValueTask>();
            List<Point> Start;
            List<double> leg;
            List<int> code;
            double Threshold_Power_1 = 0;
            double[] Threshold_Power_2 = new double[8];
            double[] Threshold_Power_3 = new double[8];

            Queue<OctaveRay> Rays = new Queue<OctaveRay>();
            BroadRay R;

            R = Source.Directions(ThreadID, ref Rnd, _octaves);
            //for (int oct = 0; oct < 8; oct++) R.Energy[oct] /= Raycount;
            lock (initlock) { for (int j = 0; j < 8; j++) _eInit += R.Energy[j]; }
            Threshold_Power_1 = R.Energy[h_oct] * 1E-2; //1E-3
            foreach (int o in _octaves)
            {
                Threshold_Power_2[o] = R.Energy[o] * thr2_template;//Cease splitting at 45 dB down.
                Threshold_Power_3[o] = R.Energy[o] * 1E-9;//-10;//Finish 3 digits (more than generous) under 60 dB of decay.
            }
            R.Decimation = Threshold_Power_3;
            int order = 0;
            double u = 0, v = 0;
            do
            {
                if (conclude) return;
                R.Ray_ID = Rnd.Next();

                if (!Room.shoot(R, out u, out v, out R.Surf_ID, out Start, out leg, out code))
                {
                    //Ray is lost... move on...
                    if (order > IS_Order)
                        RecMain.CheckBroadbandRay(R, new Point(R.x + R.dx * 1000000, R.y + R.dy * 1000000, R.z + R.dz * 1000000));
                    lock (lostlock)
                    {
                        for (int j = 0; j < 8; j++) _lost += R.Energy[j];
                    }
                    BroadRayPool.Instance.release();
                    goto Do_Scattered;
                }

                if (order > IS_Order)
                {
                    //Specular ray order exceeds order of image source calculation. Start logging specular part.
                    RecMain.CheckBroadbandRay(R, Start[0]);
                }
                R.Energy[0] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[0] * leg[0]);
                R.Energy[1] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[1] * leg[0]);
                R.Energy[2] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[2] * leg[0]);
                R.Energy[3] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[3] * leg[0]);
                R.Energy[4] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[4] * leg[0]);
                R.Energy[5] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[5] * leg[0]);
                R.Energy[6] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[6] * leg[0]);
                R.Energy[7] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[7] * leg[0]);
                R.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                R.x = Start[0].x; R.y = Start[0].y; R.z = Start[0].z;

                order++;

                ///////////////////////////////////////////////////////////////
                //The Split Case : Specular and Diffuse Explicit
                ///Standard Energy Conservation Routine-
                // 1. Multiply by Reflection energy. (1 - Absorption)
                double cos_theta;
                Room.Absorb(ref R, out cos_theta, u, v);
                // 2. Apply Transmission (if any). 
                //// a. Create new source for transmitted energy (E * Transmission).
                //// b. Modify E (E * 1 - Transmission).
                //foreach (int oct in _octaves)
                //{
                //    if (Room.TransmissionValue[R.Surf_ID][oct] > 0.0)
                //    {
                //        Rays.Enqueue(R.SplitRay(oct, Room.TransmissionValue[R.Surf_ID][oct]));
                //    }
                //}
                //3. Apply Scattering
                Room.Scatter_Early(ref R, ref Rays, ref Rnd, cos_theta, u, v, Room.TransmissionValue[R.Surf_ID]);
                ///////////////////////////////////////////////////////////////

                //Utilities.PachTools.Ray_Acoustics.SpecularReflection(ref R.direction, ref Room, ref _u[Params.ThreadID], ref _v[Params.ThreadID], ref R.Surf_ID);
            } while (R.Energy[h_oct] > Threshold_Power_1);

            foreach (int o in _octaves)
            {
                OctaveRay OR = R.SplitRay(o);
                //OR.Decimation = Threshold_Power_3[o];
                Rays.Enqueue(OR);
            };//Split all rays into individual octaves.

            BroadRayPool.Instance.release();
        Do_Scattered:

            if (Rays.Count > 0)
            {
                do
                {
                    OctaveRay OR = Rays.Dequeue();
                    do
                    {
                        OR.Ray_ID = Rnd.Next();
                        if (!Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                        {
                            RecMain.CheckRay(OR, new Point(OR.x + OR.dx * 1000000, OR.y + OR.dy * 1000000, OR.z + OR.dz * 1000000));
                            lock (lostlock) { _lost += OR.Intensity; }
                            OctaveRayPool.Instance.release();
                            goto Do_Scattered;
                        }

                        RecMain.CheckRay(OR, Start[0]);

                        OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                        OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                        OR.x = Start[0].x; OR.y = Start[0].y; OR.z = Start[0].z;

                        ///Standard Energy Conservation Routine-
                        // 1. Multiply by Reflection energy. (1 - Absorption)
                        double cos_theta;
                        Room.Absorb(ref OR, out cos_theta, u, v);

                        // 2. Apply Transmission (if any).
                        //if (Room.TransmissionValue[OR.Surf_ID][OR.Octave] > 0.0)
                        //{
                        //    //_rayTotal[Params.ThreadID]++;
                        //    //// a. Create new source for transmitted energy (E * Transmission).
                        //    //// b. Modify E (E * 1 - Transmission).
                        //    if (Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave])
                        //    {
                        //        OR.direction *= -1;
                        //        Room.ReflectRay(ref OR, ref _u[OR.ThreadID], ref _v[OR.ThreadID], ref OR.Surf_ID);
                        //    }
                        //}

                        // 2,3. Apply Scattering & Transmission.
                        Room.Scatter_Late(ref OR, ref Rays, ref Rnd, cos_theta, u, v, Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave]);
                    }
                    while (OR.t_sum < COTime && OR.Intensity > Threshold_Power_2[OR.Octave]);

                    int seed = Rnd.Next();
                    OR.Ray_ID = seed;

                    ValueTask t = new ValueTask(Task.Run(() => End_Trace(OR)));
                    LastRays.Add(t);
                }
                while (Rays.Count > 0);
            }

            for (int i = 0; i < LastRays.Count; i++)
                await LastRays[i];
        }

        object rt_Lock = new object();
        private void End_Trace(OctaveRay OR)
        {
            lock (rt_Lock) { _rayTotal[0]++; }
            double u, v;
            List<Point> Start;
            List<double> leg;
            List<int> code;
            Random Rnd = new Random(OR.Ray_ID);

            do
            {
                OR.Ray_ID = Rnd.Next();

                if (!Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                {
                    RecMain.CheckRay(OR, new Point(OR.x + OR.dx * 1000000, OR.y + OR.dy * 1000000, OR.z + OR.dz * 1000000));
                    lock (lostlock) { _lost += OR.Intensity; }
                    break;
                }

                RecMain.CheckRay(OR, Start[0]);
                OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                OR.x = Start[0].x; OR.y = Start[0].y; OR.z = Start[0].z;

                ///Standard Energy Conservation Routine-
                // 1. Multiply by Reflection energy. (1 - Absorption)
                double cos_theta;
                Room.Absorb(ref OR, out cos_theta, u, v);

                // 2. Apply Transmission (if any).
                //if (Room.TransmissionValue[OR.Surf_ID][OR.Octave] > 0.0)
                //{
                //    //// a. Create new source for transmitted energy (E * Transmission).
                //    //// b. Modify E (E * 1 - Transmission).
                //    if (Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave])
                //    {
                //        OR.direction *= -1;
                //        Room.ReflectRay(ref OR, ref _u[OR.ThreadID], ref _v[OR.ThreadID], ref OR.Surf_ID);
                //    }
                //}

                // 2,3. Apply Scattering & Transmission.
                Room.Scatter_Simple(ref OR, ref Rnd, cos_theta, u, v, Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave]);
            }
            while (OR.t_sum < COTime && OR.Intensity > OR.Decimation_threshold);
            OctaveRayPool.Instance.release();
        }

        private void cast_ray(Scene Room, int ThreadID, ref Random Rnd)
        {
            List<Point> Start;
            List<double> leg;
            List<int> code;
            double Threshold_Power_1 = 0;
            double[] Threshold_Power_2 = new double[8];
            double[] Threshold_Power_3 = new double[8];
            Queue<OctaveRay> Rays = new Queue<OctaveRay>();
            BroadRay R;

            R = Source.Directions(ThreadID, ref Rnd, _octaves);

            lock (initlock)
            {
                for (int j = 0; j < 8; j++) _eInit += R.Energy[j];
            }
            Threshold_Power_1 = R.Energy[h_oct] * 1E-2; //1E-3
            foreach (int o in _octaves)
            {
                Threshold_Power_2[o] = R.Energy[o] * Math.Pow(10, -4.5);//Cease splitting at 0.0001 sound intensity, or 40 dB down.
                Threshold_Power_3[o] = R.Energy[o] * 1E-9;//Finish 4 digits (more than generous) under 60 dB of decay. Full double Precision would be 15 digits.
            }
            int order = 0;
            double u = 0, v = 0;
            do
            {
                R.Ray_ID = Rnd.Next();

                if (!Room.shoot(R, out u, out v, out R.Surf_ID, out Start, out leg, out code))
                {
                    //Ray is lost... move on...
                    if (order > IS_Order)
                        RecMain.CheckBroadbandRay(R, new Point(R.x + R.dx * 1000000, R.y + R.dy * 1000000, R.z + R.dz * 1000000));
                    lock (lostlock) { for (int j = 0; j < 8; j++) _lost += R.Energy[j]; }
                    BroadRayPool.Instance.release();
                    goto Do_Scattered;
                }

                if (order > IS_Order)
                {
                    //Specular ray order exceeds order of image source calculation. Start logging specular part.
                    RecMain.CheckBroadbandRay(R, Start[0]);
                }
                R.Energy[0] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[0] * leg[0]);
                R.Energy[1] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[1] * leg[0]);
                R.Energy[2] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[2] * leg[0]);
                R.Energy[3] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[3] * leg[0]);
                R.Energy[4] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[4] * leg[0]);
                R.Energy[5] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[5] * leg[0]);
                R.Energy[6] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[6] * leg[0]);
                R.Energy[7] *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[7] * leg[0]);
                R.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                R.x = Start[0].x; R.y = Start[0].y; R.z = Start[0].z;

                order++;

                ///////////////////////////////////////////////////////////////
                //The Split Case : Specular and Diffuse Explicit
                ///Standard Energy Conservation Routine-
                // 1. Multiply by Reflection energy. (1 - Absorption)
                double cos_theta;
                Room.Absorb(ref R, out cos_theta, u, v);
                //// 2. Apply Transmission (if any). 
                ////// a. Create new source for transmitted energy (E * Transmission).
                ////// b. Modify E (E * 1 - Transmission).
                ////foreach (int oct in _octaves)
                //{
                //    if (Room.TransmissionValue[R.Surf_ID][oct] > 0.0)
                //    {
                //        Rays.Enqueue(R.SplitRay(oct, Room.TransmissionValue[R.Surf_ID][oct]));
                //    }
                //}
                //2, 3. Apply Scattering & Transmission
                Room.Scatter_Early(ref R, ref Rays, ref Rnd, cos_theta, u, v, Room.TransmissionValue[R.Surf_ID]);
                ///////////////////////////////////////////////////////////////

                //Utilities.PachTools.Ray_Acoustics.SpecularReflection(ref R.direction, ref Room, ref _u[Params.ThreadID], ref _v[Params.ThreadID], ref R.Surf_ID);
            } while (R.Energy[h_oct] > Threshold_Power_1);

            foreach (int o in _octaves) Rays.Enqueue(R.SplitRay(o));//Split all rays into individual octaves.
            BroadRayPool.Instance.release();
        Do_Scattered:

            if (Rays.Count > 0)
            {
                do
                {
                    OctaveRay OR = Rays.Dequeue();
                    lock (rt_Lock) { _rayTotal[0]++; }
                    do
                    {
                        OR.Ray_ID = Rnd.Next();
                        if (!Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                        {
                            RecMain.CheckRay(OR, new Point(OR.x + OR.dx * 1000000, OR.y + OR.dy * 1000000, OR.z + OR.dz * 1000000));
                            lock (lostlock) { _lost += OR.Intensity; }
                            OctaveRayPool.Instance.release();
                            goto Do_Scattered;
                        }

                        RecMain.CheckRay(OR, Start[0]);

                        OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                        OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                        OR.x = Start[0].x; OR.y = Start[0].y; OR.z = Start[0].z;

                        ///Standard Energy Conservation Routine-
                        // 1. Multiply by Reflection energy. (1 - Absorption)
                        double cos_theta;
                        Room.Absorb(ref OR, out cos_theta, u, v);

                        //// 2. Apply Transmission (if any).
                        //if (Room.TransmissionValue[OR.Surf_ID][OR.Octave] > 0.0)
                        //{
                        //    //_rayTotal[Params.ThreadID]++;
                        //    //// a. Create new source for transmitted energy (E * Transmission).
                        //    //// b. Modify E (E * 1 - Transmission).
                        //    if (Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave])
                        //    {
                        //        OR.direction *= -1;
                        //        Room.ReflectRay(ref OR, ref _u[OR.ThreadID], ref _v[OR.ThreadID], ref OR.Surf_ID);
                        //    }
                        //}

                        // 2,3. Apply Scattering and transmission
                        Room.Scatter_Late(ref OR, ref Rays, ref Rnd, cos_theta, u, v, Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave]);
                    }
                    while (OR.t_sum < COTime && OR.Intensity > Threshold_Power_2[OR.Octave]);

                    do
                    {
                        OR.Ray_ID = Rnd.Next();

                        if (!Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                        {
                            RecMain.CheckRay(OR, new Point(OR.x + OR.dx * 1000000, OR.y + OR.dy * 1000000, OR.z + OR.dz * 1000000));
                            lock (lostlock) { _lost += OR.Intensity; }
                            OctaveRayPool.Instance.release();
                            goto Do_Scattered;
                        }

                        RecMain.CheckRay(OR, Start[0]);
                        //OctChecks.Enqueue(new object[] { OR.Clone(), new Hare.Geometry.Point(Start[0].x, Start[0].y, Start[0].z) });
                        OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                        OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                        OR.x = Start[0].x; OR.y = Start[0].y; OR.z = Start[0].z;

                        ///Standard Energy Conservation Routine-
                        // 1. Multiply by Reflection energy. (1 - Absorption)
                        double cos_theta;
                        Room.Absorb(ref OR, out cos_theta, u, v);

                        //// 2. Apply Transmission (if any).
                        //if (Room.TransmissionValue[OR.Surf_ID][OR.Octave] > 0.0)
                        //{
                        //    //// a. Create new source for transmitted energy (E * Transmission).
                        //    //// b. Modify E (E * 1 - Transmission).
                        //    if (Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave])
                        //    {
                        //        OR.direction *= -1;
                        //        Room.ReflectRay(ref OR, ref _u[OR.ThreadID], ref _v[OR.ThreadID], ref OR.Surf_ID);
                        //    }
                        //}

                        // 2, 3. Apply Scattering & Transmission.
                        Room.Scatter_Simple(ref OR, ref Rnd, cos_theta, u, v, Rnd.NextDouble() < Room.TransmissionValue[OR.Surf_ID][OR.Octave]);
                    }
                    while (OR.t_sum < COTime && OR.Intensity > Threshold_Power_3[OR.Octave]);
                    OctaveRayPool.Instance.release();
                }
                while (Rays.Count > 0);
            }
        }

        /// <summary>
        /// Returns the completed receiver object. This is used to extract the result histogram.
        /// </summary>
        public Receiver_Bank GetReceiver
        {
            get
            {
                return RecMain;
            }
        }

        public Convergence_Check[] Convergence_Report
        {
            get { return check; }
        }

        /// <summary>
        /// A string to identify the type of simulation being run.
        /// </summary>
        /// <returns></returns>
        public override string Sim_Type()
        {
            return "Stochastic Ray Tracing";
        }

        /// <summary>
        /// Called by Pach_RunSim_Command. Get a string describing the status of the simulation for display.
        /// </summary>
        /// <returns></returns>
        public override string ProgressMsg()
        {
            if (check == null)
            {
                _ts = DateTime.Now - _st;
                //int Ray_CT = 0;
                //int subrays = 0;
                //for (int i = 0; i < _processorCt; i++)
                //{
                //Ray_CT += _currentRay;
                //subrays += _rayTotal;
                //}
                if (_currentRay[0] == 0) return "Preparing Threads";
                _ts = new TimeSpan((long)(_ts.Ticks * (((double)(Raycount - _currentRay[0])) / (_currentRay[0] + 1))));
                return string.Format("Calculating Ray {0} of {1} (sub-ray {2}). ({3} hours,{4} minutes,{5} seconds Remaining.) Press 'Esc' to Cancel...", _currentRay[0], Raycount, _rayTotal[0], _ts.Hours, _ts.Minutes, _ts.Seconds);
            }
            else
            {
                if (this.Vis_Feedback != null) Vis_Feedback.Populate();
                foreach (Convergence_Check c in check)
                {
                    if (c == null) break;
                    c.Update();
                    //if (!c.Check(this._currentRay.Sum())||conclude) continue;
                    if (!conclude) break;
                    //this.Abort_Calculation();
                    return "Concluding Simulation...";
                }

                //Convergence_Progress_WinForms.Instance.Refresh();

                //int Ray_CT = 0;
                //int subrays = 0;
                //for (int i = 0; i < _processorCt; i++)
                //{
                //    Ray_CT += _currentRay;
                //    subrays += _rayTotal[i];
                //}
                return string.Format("Revolution {0} - Ray {1} : sub {2} - Convergence not yet reached...", (int)Math.Ceiling((double)_currentRay[0] / 80), _currentRay[0], _rayTotal[0]);
            }
        }

        /// <summary>
        /// Consolidates output from all threads into a single set of output.
        /// </summary>
        public override void Combine_ThreadLocal_Results()
        {
            RecMain.Scale(_currentRay[0]);

            this._ts = DateTime.Now - _st;

            double LTotal = 0;
            double RTotal = 0;
            //foreach (double L in _lost)
            //{
            //    LTotal += L;
            //}
            //foreach (double T in _eInit)
            //{
            //    RTotal += T;
            //}
            _lost_percent = 100 * _lost / _eInit;
        }

        public double PercentLost
        {
            get
            {
                return _lost_percent;
            }
        }

        public abstract class Convergence_Check
        {
            protected double[] RunningSim;
            protected int oct;
            public int check_no;

            public delegate void PlotHandler(double[] Conv1, double Conv2, double ConvInf, int ID, int count, double corr);
            public event PlotHandler On_Convergence_Check;

            public Convergence_Check(Receiver_Bank R, int id, int _oct, int check_id, I_Conv_Progress Vis_Feedback)
            {
                oct = _oct;
                RunningSim = R.Rec_List[id].Recs.Energy[oct];
                check_no = check_id;
                On_Convergence_Check += Vis_Feedback.Fill;
            }

            protected async void Plot_Feedback(double[] diff, double conv, double ConvInf, int ID, int count, double corr)
            {
                if (On_Convergence_Check != null) On_Convergence_Check(diff, conv, ConvInf, ID, count, corr);
            }

            public abstract bool Check();
            public abstract bool Check(int No_of_Rays);
            public abstract void Update();
        }

        public class Detailed_Convergence_Check : Convergence_Check
        {
            double[] Snapshot;
            int step;
            int binct;
            int count = 0;
            int RayNo = 0;

            double[] diff;
            double conv;

            public Detailed_Convergence_Check(Receiver_Bank R, int id, int _oct, int check_id, I_Conv_Progress Vis_Feedback)
                : base(R, id, _oct, check_id, Vis_Feedback)
            {
                step = (R.SampleRate / 1000);
                binct = RunningSim.Length / step;
                Snapshot = new double[binct];
            }

            public override bool Check()
            {
                throw new NotImplementedException();
            }

            public override bool Check(int No_of_Rays)
            {
                double delta_rays = No_of_Rays - RayNo;
                if (delta_rays < 100) return false;
                RayNo = No_of_Rays;

                double[] sn = new double[Snapshot.Length];
                for (int i = 0; i < Snapshot.Length; i++)
                {
                    for (int j = 0; j < step; j++)
                    {
                        sn[i] += RunningSim[i * step + j];
                    }
                }

                double[] SI = Pachyderm_Acoustic.Utilities.AcousticalMath.Schroeder_Integral(sn);
                double inf60 = SI[0] * 1E-6;

                int stop = SI.Length - 1;
                while (sn[stop] < inf60) stop--;
                double conv = 0;

                double[] diff = new double[sn.Length];

                for (int i = 0; i < stop; i++)
                {
                    if (sn[i] == 0 || Snapshot[i] == 0) continue;
                    double d = Math.Abs(sn[i] - Snapshot[i]) / Snapshot[i];
                    diff[i] = d;
                    conv = Math.Max(conv, d);
                }
                Snapshot = sn;

                if (conv == 0) return false;

                if (conv > 0.1)
                {
                    count = 0;
                }
                count++;


                if (count > 10) return true;
                return false;
            }

            public async override void Update()
            {
                //await Task.Run(() => { Plot_Feedback(diff, conv, 1000, 0, count, check_no);});
                Plot_Feedback(diff, conv, 1000, 0, count, check_no);
            }
        }

        public class Minimum_Convergence_Check : Convergence_Check
        {
            double snapshot50, snapshot80, snapshotinf;
            int SampleStart, Sample50, Sample80, SampleInf;
            int count = 0;
            int RayNo = 0;
            double[] Schr_old;
            double conv1, conv2, convinf, r;

            public Minimum_Convergence_Check(Source S, Scene Sc, Receiver_Bank R, int id, int _oct, int check_id, I_Conv_Progress Vis_Feedback)
                : base(R, id, _oct, check_id, Vis_Feedback)
            {
                SampleStart = (int)Math.Floor((S.Origin - R.Origin(id)).Length() / Sc.Sound_speed(R.Origin(id)) * R.SampleRate);
                Sample50 = (int)Math.Floor(50.0 * R.SampleRate / 1000) + SampleStart;
                Sample80 = (int)Math.Floor(80.0 * R.SampleRate / 1000) + SampleStart;
                SampleInf = R.SampleCT;
            }

            public override bool Check()
            {
                double sn50 = 0;
                double sn80 = 0;
                double sninf = 0;
                for (int i = SampleStart; i < Sample50; i++)
                {
                    sn50 += RunningSim[i];
                }
                sn80 += sn50;
                for (int i = Sample50; i < Sample80; i++)
                {
                    if (i < RunningSim.Length) sn80 += RunningSim[i];
                }
                for (int i = Sample80; i < SampleInf; i++)
                {
                    sninf += RunningSim[i];
                }

                double[] Schr_new = AcousticalMath.Log10Data(AcousticalMath.Schroeder_Integral(RunningSim), -70);

                r = Schr_old == null ? 0 : MathNet.Numerics.Statistics.Correlation.Spearman(Schr_new, Schr_old);
                //MathNet.Numerics.Statistics.Correlation.Pearson(Schr_new, Schr_old);

                Schr_old = Schr_new;

                conv1 = Math.Abs(sn50 - snapshot50) / (snapshot50);
                conv2 = Math.Abs(sn80 - snapshot80) / (snapshot80);
                convinf = Math.Abs(sninf - snapshotinf) / (snapshotinf);

                if (sn50 == 0
                    || conv1 > 0.02
                    || sn80 == 0
                    || conv2 > 0.02
                    || sninf == 0
                    || convinf > 0.1)
                {
                    count = 0;
                }
                else
                {
                    count++;
                }

                snapshot50 = sn50;
                snapshot80 = sn80;
                snapshotinf = sninf;

                if (count > 5) return true;
                return false;
            }

            public override void Update()
            {
                Plot_Feedback(new double[1] { conv1 }, conv2, convinf, 0, check_no, r);
                //return true;
            }

            public override bool Check(int No_of_Rays)
            {
                double delta_rays = No_of_Rays - RayNo;
                if (delta_rays < 80) return false;
                RayNo = No_of_Rays;
                delta_rays /= 80;

                double sn50 = 0;
                double sn80 = 0;
                double sninf = 0;
                for (int i = SampleStart; i < Sample50; i++)
                {
                    sn50 += RunningSim[i];
                }
                sn80 += sn50;
                for (int i = Sample50; i < Sample80; i++)
                {
                    if (i < RunningSim.Length) sn80 += RunningSim[i];
                }
                for (int i = Sample80; i < SampleInf; i++)
                {
                    sninf += RunningSim[i];
                }

                double conv1 = Math.Abs(sn50 - snapshot50) / (snapshot50 * delta_rays), conv2 = Math.Abs(sn80 - snapshot80) / (snapshot80 * delta_rays), convinf = Math.Abs(sninf - snapshotinf) / (snapshotinf * delta_rays);

                if (sn50 == 0
                    || conv1 > 0.02
                    || sn80 == 0
                    || conv2 > 0.02
                    || sninf == 0
                    || convinf > 0.1)
                {
                    count = 0;
                }
                else
                {
                    count++;
                }

                snapshot50 = sn50;
                snapshot80 = sn80;
                snapshotinf = sninf;

                Plot_Feedback(new double[1] { conv1 }, conv2, convinf, check_no, count, 1);

                if (count > 10) return true;
                return false;
            }
        }
    }

    public interface I_Conv_Progress
    {
        void Populate();
        void Fill(double[] Conv1, double Conv2, double ConvInf, int ID, int count, double corr);
    }

    public class RayQueue<T> : ConcurrentQueue<T>
    {
        int ProcessorCt;
        public RayQueue(int no_of_processors)
        : base()
        {
            ProcessorCt = no_of_processors;
        }

        public AutoResetEvent Ray_Loaded = new AutoResetEvent(false);
        public AutoResetEvent Rays_dwindling = new AutoResetEvent(false);

        public new void Enqueue(T item)
        {
            base.Enqueue(item);
            Ray_Loaded.Set();
        }

        public new bool TryDequeue(out T item)
        {
            if (Count < ProcessorCt) Rays_dwindling.Set();
            Rays_dwindling.Set();
            return base.TryDequeue(out item);
        }
    }
}