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
using System.ComponentModel;
using System.Windows.Forms;

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
        public int[] _currentRay;
        private DateTime _st;
        private double[] _u, _v;
        private double[] _lost;
        private double[] _eInit;
        private double _lost_percent;
        private int IS_Order;
        private System.Threading.Thread[] _tlist;
        private int _processorCt;
        public int[] _rayTotal;
        public TimeSpan _ts;
        private int h_oct;
        private int[] _octaves;
        private Convergence_Check[] check;
        private bool conclude = false;

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
        public SplitRayTracer(Source sourceIn, Receiver_Bank receiverIn, Scene roomIn, double cutoffTime, int[] octaveRange, int isOrderIn, int rayCountIn)
        {
            IS_Order = isOrderIn;
            Room = roomIn;
            RecMain = receiverIn;
            Raycount = rayCountIn;
            _octaves = new int[octaveRange[1] - octaveRange[0] + 1];
            for (int o = octaveRange[0]; o <= octaveRange[1]; o++) _octaves[o - octaveRange[0]] = o;
            COTime = cutoffTime;
            Source = sourceIn;
            double[] totalAbs = new double[8];

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
                    spt.Add(Source.Origin());
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
                check = (rayCountIn < 0) ? new Convergence_Check[2] { T_id < 0 ? null : new Minimum_Convergence_Check(this.Source, this.Room, receiverIn, T_id, h_oct, 0), F_id < 0 ? null : new Minimum_Convergence_Check(this.Source, this.Room, receiverIn, F_id, h_oct, 1) }
                : new Convergence_Check[2] { T_id < 0 ? null : new Detailed_Convergence_Check(receiverIn, T_id, h_oct, 0), F_id < 0 ? null : new Detailed_Convergence_Check(receiverIn, F_id, h_oct, 1) };
                //Pachyderm_Acoustic.UI.Convergence_Progress.Instance.Show();
            }
            else check = null;
        }

        private bool Check_Validity(Point p, int rec_id, int rnd, out double dist)
        {
            Hare.Geometry.Vector d = RecMain.Origin(rec_id) - p;
            dist = d.Length();
            d.Normalize();
            Ray R = new Ray(Source.H_Origin(), d, 0, rnd);
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

        /// <summary>
        /// Inherited member. Divides the simulation into threads, and begins.
        /// </summary>
        public override void Begin()
        {
            Random Rnd = new Random((int)DateTime.Now.ToFileTimeUtc());
            _processorCt = Pach_Properties.Instance.ProcessorCount();
            _currentRay = new int[_processorCt];
            _rayTotal = new int[_processorCt];
            _lost = new double[_processorCt];
            _eInit = new double[_processorCt];
            _u = new double[_processorCt];
            _v = new double[_processorCt];

            //System.Runtime.GCSettings.LatencyMode = System.Runtime.GCLatencyMode.SustainedLowLatency;

            if (Raycount > 0)
            {
                _tlist = new System.Threading.Thread[_processorCt];
                for (int P_I = 0; P_I < _processorCt; P_I++)
                {
                    int start = (int)Math.Floor((double)P_I * Raycount / _processorCt);
                    int end;
                    if (P_I == _processorCt - 1) end = Raycount;
                    else end = (int)Math.Floor((double)(P_I + 1) * Raycount / _processorCt);
                    Calc_Params T = new Calc_Params(Room, start, end, P_I, Rnd.Next());
                    System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { Calculate(T); });
                    _tlist[P_I] = new System.Threading.Thread(TS);
                    _tlist[P_I].Start();
                }
            }
            else
            {
                _tlist = new System.Threading.Thread[_processorCt + 1];
                System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { this.Calculate_Conv(_processorCt, Rnd.Next()); });
                _tlist[0] = new System.Threading.Thread(TS);
                _tlist[0].Start();
            }
        }

        /// <summary>
        /// Called by Pach_RunSim_Command. Indicates whether or not the simulation has completed.
        /// </summary>
        /// <returns>Returns running if any threads in this simulation are still running. Returns stopped if all have stopped.</returns>
        public override System.Threading.ThreadState ThreadState()
        {
            foreach (System.Threading.Thread T in _tlist)
            {
                if (!(T == null) && (T.ThreadState == System.Threading.ThreadState.Running || T.ThreadState == System.Threading.ThreadState.WaitSleepJoin)) 
                    return System.Threading.ThreadState.Running;
            }
            return System.Threading.ThreadState.Stopped;
        }

        /// <summary>
        /// Aborts all threads, effectively ending the simulation.
        /// </summary>
        public override void Abort_Calculation()
        {
            foreach (System.Threading.Thread T in _tlist)
            {
                T.Abort();
            }
        }

        public void Conclude_Simulation()
        {
            conclude = true;
        }

        /// <summary>
        /// Called by each thread from the begin method for calculations with a user-defined length.
        /// </summary>
        /// <param name="i">the object is type "Calc_Params" which holds all the necessary information to run a portion of the simulation.</param>
        public void Calculate(object i)
        {
            ///Homogeneous media only...
            Calc_Params Params = (Calc_Params)i;
            Random Rnd = new Random(Params.RandomSeed);
            _st = DateTime.Now;

            LastRays = new BlockingCollection<OctaveRay>[_processorCt];
            _st = DateTime.Now;
            //Thread[] lastraytracers = new Thread[_processorCt];
            //Thread lastraytracers;

            for (int P_I = 0; P_I < _processorCt; P_I++)
            {
                LastRays[P_I] = new BlockingCollection<OctaveRay>();
            }

            //for (int P_I = 0; P_I < _processorCt; P_I++)
            //{
            //    System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { TraceWorker(Rnd.Next(), P_I); });
            //    lastraytracers[P_I] = new System.Threading.Thread(TS);
            //    lastraytracers[P_I].Name = "Last Trace " + P_I;
            //    lastraytracers[P_I].Start();
            //}

            //System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { TraceWorker(Rnd.Next(), Params.ThreadID); });
            //lastraytracers = new Thread(TS);
            //lastraytracers.Name = "Last Trace " + Params.ThreadID;
            //lastraytracers.Start();

            for (int ray = 0; ray < Params.EndIndex - Params.StartIndex; ray++)
            {
                cast_ray_distributed(Params.Room, Params.ThreadID, ref Rnd);
                _currentRay[Params.ThreadID]++;
            }

            //lastraytracers.Join();

            int ct;
            do
            {
                ct = 0;
                for (int j = 0; j < _processorCt; j++) ct += LastRays[j].Count;
                Thread.Sleep(500);
            } while (ct > 0);

            _ts = DateTime.Now - _st;
        }

        /// <summary>
        /// Called by each thread from the begin method.
        /// </summary>
        /// <param name="i">the object is type "Calc_Params" which holds all the necessary information to run a portion of the simulation.</param>
        public void Calculate_Conv(object i, object seed)
        {
            ///Rewrite to run all threads? Possibly keep a record for comparison in check, and run the next calc while checking?

            ///Homogeneous media only...
            LastRays = new BlockingCollection<OctaveRay>[_processorCt];
            Random rnd = new Random((int)seed);
            _st = DateTime.Now;
            //Thread[] lastraytracers = new Thread[_processorCt];
            bool conv = true;

            for (int P_I = 0; P_I < _processorCt; P_I++)
            {
                LastRays[P_I] = new BlockingCollection<OctaveRay>();
            }

            //for (int P_I = 0; P_I < Math.Max(1,_processorCt); P_I++)
            //{
            //    //for (int j = 0; j < 2; j++)
            //    //{
            //        System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { TraceWorker(rnd.Next(), P_I); });
            //        lastraytracers[P_I] = new System.Threading.Thread(TS);
            //        lastraytracers[P_I].Name = "Last Trace " + P_I;
            //        lastraytracers[P_I].Start();
            //    //}
            //}

            int rev = 0;

            do
            {
                rev++;
                if (conclude) 
                    break;
                for (int P_I = 0; P_I < _processorCt; P_I++)
                {
                    int start = (int)Math.Floor((double)P_I * 80 / _processorCt);
                    int end;
                    if (P_I == _processorCt - 1) end = 80;
                    else end = (P_I + 1) * (80 / _processorCt);
                    Calc_Params T = new Calc_Params(Room, start, end, P_I, rnd.Next());
                    System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate
                    {
                        Random Rnd = new Random(T.RandomSeed);

                        for (int ray = 0; ray < T.EndIndex - T.StartIndex; ray++)
                        {
                            //cast_ray(Params.ThreadID, ref Rnd);
                            cast_ray_distributed(T.Room, T.ThreadID, ref Rnd);
                            _currentRay[T.ThreadID]++;
                        }
                    });
                    _tlist[P_I + 1] = new System.Threading.Thread(TS);
                    _tlist[P_I + 1].Name = "Prime Trace " + P_I;
                    _tlist[P_I + 1].Start();
                }

                do
                {
                    Thread.Sleep(100);
                    int completed = 0;
                    for (int r = 0; r < _currentRay.Length; r++) completed +=_currentRay[r];
                    if (completed < 80 * rev) continue;
                    int remainder = 0;
                    int j;
                    for (j = 0; j < 3; j++)
                    {
                        for (int r = 0; r < LastRays.Length; r++) remainder += LastRays[r].Count;
                        if (remainder > 0) break;
                    }
                    if (j >= 2) break;
                } while (true);

                GC.Collect();
                conv = true;
                foreach (Convergence_Check c in check) if (c != null) conv &= c.Check();
            } while (!conv);

            Conclude_Simulation();
            foreach (Thread t in _tlist) t.Join();

            while (LastRays[0].Count > 0)
            {
                Thread.Sleep(200);
            }

            //foreach (Thread t in lastraytracers) t.Abort();
            _ts = DateTime.Now - _st;
        }

        //RayQueue<OctaveRay> LastRays;
        BlockingCollection<OctaveRay>[] LastRays;
        double thr2_template = Math.Pow(10, -4.5);

        private void cast_ray_distributed(Scene Room, int ThreadID, ref Random Rnd)
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
            //for (int oct = 0; oct < 8; oct++) R.Energy[oct] /= Raycount;
            for (int j = 0; j < 8; j++) _eInit[ThreadID] += R.Energy[j];
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
                        RecMain.CheckBroadbandRay(R, R.origin + R.direction * 1000000);
                    for (int j = 0; j < 8; j++) _lost[ThreadID] += R.Energy[j];
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
                R.origin = Start[0];

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
                            RecMain.CheckRay(OR, OR.origin + OR.direction * 1000000);
                            _lost[ThreadID] += OR.Intensity;
                            goto Do_Scattered;
                        }

                        RecMain.CheckRay(OR, Start[0]);

                        OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                        OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                        OR.origin = Start[0];

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

                    //OR.Decimation = Threshold_Power_3[OR.Octave];
                    //LastRays.Enqueue(OR);
                    LastRays[ThreadID].TryAdd(OR);
                }
                while (Rays.Count > 0);
            }
            
            if (LastRays[ThreadID].Count > 0)
            {
                OctaveRay OR;
                do
                {
                    if (LastRays[ThreadID].TryTake(out OR))
                    {
                        End_Trace(OR, ref Rnd, ref u, ref v, ref Start, ref leg, ref code);
                    }
                    if (LastRays[ThreadID].Count == 0)
                    {
                        TraceWorker(Rnd.Next(), ThreadID);
                        return;
                    }
                } while (true);
            }
        }

        private void TraceWorker(int randomseed, int threadid)
        {
            List<Point> Start = new List<Point>();
            List<double> leg = new List<double>();
            List<int> code = new List<int>();
            double u = 0, v = 0;
            Random Rnd = new Random(randomseed);
            int procmod = 0;

            while (true)
            {
                OctaveRay OR;
                while (true)
                {
                    if (LastRays[(int)Math.Floor((double)threadid + procmod) % _processorCt].TryTake(out OR))
                    {
                        _rayTotal[OR.ThreadID]++;
                        break;
                    }
                    procmod++;
                    if (procmod > _processorCt)
                    {
                        if (Raycount > 0)
                        {
                            //if (_currentRay.Sum() >= Raycount) return;
                            int remain = 0;
                            for (int i = 0; i < LastRays.Length; i++) remain += LastRays[i].Count;
                            if (remain == 0) return;
                        }
                        else return;
                        Thread.Sleep(100);
                    }
                }
                procmod = 0;
                End_Trace(OR, ref Rnd, ref u, ref v, ref Start, ref leg, ref code);
            }
        }

        private void End_Trace(OctaveRay OR, ref Random Rnd, ref double u, ref double v, ref List<Point> Start, ref List<double> leg, ref List<int> code)
        {
            do
            {
                OR.Ray_ID = Rnd.Next();

                if (!Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                {
                    RecMain.CheckRay(OR, OR.origin + OR.direction * 1000000);
                    _lost[OR.ThreadID] += OR.Intensity;
                    break;
                }

                RecMain.CheckRay(OR, Start[0]);
                OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                OR.origin = Start[0];

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
        }

        //private void cast_ray(int ThreadID, ref Random Rnd)   
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
            //for (int oct = 0; oct < 8; oct++) R.Energy[oct] /= Raycount;
            for (int j = 0; j < 8; j++) _eInit[ThreadID] += R.Energy[j];
            Threshold_Power_1 = R.Energy[h_oct] * 1E-2; //1E-3
            foreach (int o in _octaves)
            {
                Threshold_Power_2[o] = R.Energy[o] * 1E-4;//Cease splitting at 0.0001 sound intensity, or 40 dB down.
                Threshold_Power_3[o] = R.Energy[o] * 1E-10;//Finish 4 digits (more than generous) under 60 dB of decay. Full double Precision would be 15 digits.
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
                        RecMain.CheckBroadbandRay(R, R.origin + R.direction * 1000000);
                    for (int j = 0; j < 8; j++) _lost[ThreadID] += R.Energy[j];
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
                R.origin = Start[0];

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
            Do_Scattered:

            if (Rays.Count > 0)
            {
                do
                {
                    OctaveRay OR = Rays.Dequeue();
                    _rayTotal[ThreadID]++;
                    do
                    {
                        OR.Ray_ID = Rnd.Next();
                        if (!Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                        {
                            RecMain.CheckRay(OR, OR.origin + OR.direction * 1000000);
                            _lost[ThreadID] += OR.Intensity;
                            goto Do_Scattered;
                        }

                        RecMain.CheckRay(OR, Start[0]);

                        OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                        OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                        OR.origin = Start[0];

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
                            RecMain.CheckRay(OR, OR.origin + OR.direction * 1000000);
                            _lost[OR.ThreadID] += OR.Intensity;
                            goto Do_Scattered;
                        }

                        RecMain.CheckRay(OR, Start[0]);
                        //OctChecks.Enqueue(new object[] { OR.Clone(), new Hare.Geometry.Point(Start[0].x, Start[0].y, Start[0].z) });
                        OR.Intensity *= Math.Pow(10, -.1 * Room.Attenuation(code[0])[OR.Octave] * leg[0]);
                        OR.AddLeg(leg[0] / Room.Sound_speed(code[0]));
                        OR.origin = Start[0];

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
                int Ray_CT = 0;
                int subrays = 0;
                for (int i = 0; i < _processorCt; i++)
                {
                    Ray_CT += (int)_currentRay[i];
                    subrays += (int)_rayTotal[i];
                }
                if (Ray_CT == 0) return "Preparing Threads";
                _ts = new TimeSpan((long)(_ts.Ticks * (((double)(Raycount - Ray_CT)) / (Ray_CT + 1))));
                return string.Format("Calculating Ray {0} of {1} (sub-ray {2}). ({3} hours,{4} minutes,{5} seconds Remaining.) Press 'Esc' to Cancel...", Ray_CT, Raycount, subrays, _ts.Hours, _ts.Minutes, _ts.Seconds);
            }
            else
            {
                foreach (Convergence_Check c in check)
                {
                    if (c == null) break;
                    c.Update();
                    //if (!c.Check(this._currentRay.Sum())||conclude) continue;
                    if (!conclude) break;
                    this.Abort_Calculation();
                    //Convergence_Progress_WinForms.Instance.Hide();
                    return "Concluding Simulation...";
                }

                //Convergence_Progress_WinForms.Instance.Refresh();

                int Ray_CT = 0;
                int subrays = 0;
                for (int i = 0; i < _processorCt; i++)
                {
                    Ray_CT += (int)_currentRay[i];
                    subrays += (int)_rayTotal[i];
                }
                return string.Format("Revolution {0} - Ray {1} : sub {2} - Convergence not yet reached...", (int)Math.Ceiling((double)Ray_CT / 80), Ray_CT, subrays);
            }
        }

        /// <summary>
        /// Consolidates output from all threads into a single set of output.
        /// </summary>
        public override void Combine_ThreadLocal_Results()
        {
            RecMain.Scale(_currentRay.Sum());

            this._ts = DateTime.Now - _st;

            double LTotal = 0;
            double RTotal = 0;
            foreach (double L in _lost)
            {
                LTotal += L;
            }
            foreach (double T in _eInit)
            {
                RTotal += T;
            }
            _lost_percent = 100 * LTotal / RTotal;
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

            public Convergence_Check(Receiver_Bank R, int id, int _oct, int check_id)
            {
                oct = _oct;
                RunningSim = R.Rec_List[id].Recs.Energy[oct];
                check_no = check_id;
            }

            protected void Plot_Feedback(double[] diff, double conv, double ConvInf, int ID, int count, double corr)
            {
                On_Convergence_Check(diff, conv, ConvInf, ID, count, corr);
            }

            public abstract bool Check();
            public abstract bool Check(int No_of_Rays);

            public abstract bool Update();
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

            public Detailed_Convergence_Check(Receiver_Bank R, int id, int _oct, int check_id)
                : base(R, id, _oct, check_id)
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

            public override bool Update()
            {
                Plot_Feedback(diff, conv, 1000, 0, count, check_no);
                return true;
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

            public Minimum_Convergence_Check(Source S, Scene Sc, Receiver_Bank R, int id, int _oct, int check_id)
                : base(R, id, _oct, check_id)
            {
                SampleStart = (int)Math.Floor((S.H_Origin() - R.Origin(id)).Length() / Sc.Sound_speed(R.Origin(id)) * R.SampleRate);
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
                
                r = Schr_old == null? 0 : MathNet.Numerics.Statistics.Correlation.Spearman(Schr_new, Schr_old);
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

            public override bool Update()
            {
                Plot_Feedback(new double[1] { conv1 }, conv2, convinf, 0, check_no, r);
                return true;
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