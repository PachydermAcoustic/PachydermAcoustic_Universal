//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2018, Arthur van der Harten 
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
using System.Runtime.InteropServices;
using System.Linq;

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
                Raycount = int.MaxValue;
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
                check = (rayCountIn < 0) ? new Convergence_Check[2] { T_id < 0 ? null : new Minimum_Convergence_Check(this.Source, this.Room, receiverIn, T_id, h_oct), F_id < 0 ? null : new Minimum_Convergence_Check(this.Source, this.Room, receiverIn, F_id, h_oct) }
                : new Convergence_Check[2] { T_id < 0 ? null : new Detailed_Convergence_Check(receiverIn, T_id, h_oct), F_id < 0 ? null : new Detailed_Convergence_Check(receiverIn, F_id, h_oct) };
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
            Random Rnd = new Random();
            _processorCt = Pach_Properties.Instance.ProcessorCount();
            _currentRay = new int[_processorCt];
            _rayTotal = new int[_processorCt];
            _lost = new double[_processorCt];
            _eInit = new double[_processorCt];
            _u = new double[_processorCt];
            _v = new double[_processorCt];
            _tlist = new System.Threading.Thread[_processorCt];

            System.Runtime.GCSettings.LatencyMode = System.Runtime.GCLatencyMode.SustainedLowLatency;

            for (int P_I = 0; P_I < _processorCt; P_I++)
            {
                int start = (int)Math.Floor((double)P_I * Raycount / _processorCt);
                int end;
                if (P_I == _processorCt - 1) end = Raycount;
                else end = (P_I + 1) * (Raycount / _processorCt);
                Calc_Params T = new Calc_Params(Room, start, end, P_I, Rnd.Next());
                System.Threading.ParameterizedThreadStart TS = new System.Threading.ParameterizedThreadStart(delegate { Calculate(T); });
                _tlist[P_I] = new System.Threading.Thread(TS);
                _tlist[P_I].Start();
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
                if (T.ThreadState == System.Threading.ThreadState.Running) return System.Threading.ThreadState.Running;
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
            //foreach (System.Threading.Thread T in _tlist)
            //{
            //    T.Join();
            //}
        }

        /// <summary>
        /// Called by each thread from the begin method.
        /// </summary>
        /// <param name="i">the object is type "Calc_Params" which holds all the necessary information to run a portion of the simulation.</param>
        public void Calculate(object i)
        {
            ///Homogeneous media only...
            Calc_Params Params = (Calc_Params)i;
            Random Rnd = new Random(Params.RandomSeed);
            _st = DateTime.Now;
            Point Origin = Source.H_Origin();
            int Dir_ID = Params.StartIndex;
            List<Point> Start;
            List<double> leg;
            List<int> code;
            Queue<OctaveRay> Rays = new Queue<OctaveRay>();
            BroadRay R;
            int order;

            double Threshold_Power_1 = 0;
            double[] Threshold_Power_2 = new double[8];
            double[] Threshold_Power_3 = new double[8];
            for (_currentRay[Params.ThreadID] = 0; _currentRay[Params.ThreadID] < Params.EndIndex - Params.StartIndex; _currentRay[Params.ThreadID]++)
            {
                R = Source.Directions(_currentRay[Params.ThreadID] + Params.StartIndex, Params.ThreadID, ref Rnd, _octaves);
                //for (int oct = 0; oct < 8; oct++) R.Energy[oct] /= Raycount;
                for (int j = 0; j < 8; j++) _eInit[Params.ThreadID] += R.Energy[j];
                Threshold_Power_1 = R.Energy[h_oct] * 1E-3;
                foreach (int o in _octaves)
                {
                    Threshold_Power_2[o] = R.Energy[o] * 1E-6;//Cease splitting at 0.0001 sound intensity, or 40 dB down.
                    Threshold_Power_3[o] = R.Energy[o] * 1E-14;//Finish 8 digits (more than generous) under 60 dB of decay. Full double Precision would be 15 digits.
                }
                order = 0;
                double u = 0, v = 0;
                do
                {
                    R.Ray_ID = Rnd.Next();
                    if (!Params.Room.shoot(R, out u, out v, out R.Surf_ID, out Start, out leg, out code))
                    {
                        //Ray is lost... move on...
                        for (int j = 0; j < 8; j++) _lost[Params.ThreadID] += R.Energy[j];
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
                    foreach (int oct in _octaves)
                    {
                        if (Params.Room.TransmissionValue[R.Surf_ID][oct] > 0.0) Rays.Enqueue(R.SplitRay(oct, Params.Room.TransmissionValue[R.Surf_ID][oct]));
                    }
                    //3. Apply Scattering
                    Room.Scatter_Early(ref R, ref Rays, ref Rnd, cos_theta, u, v);
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
                        _rayTotal[Params.ThreadID]++;
                        do
                        {
                            OR.Ray_ID = Rnd.Next();
                            if (!Params.Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                            {
                                _lost[Params.ThreadID] += OR.Intensity;
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
                            if (Params.Room.TransmissionValue[OR.Surf_ID][OR.Octave] > 0.0)
                            {
                                //_rayTotal[Params.ThreadID]++;
                                //// a. Create new source for transmitted energy (E * Transmission).
                                //// b. Modify E (E * 1 - Transmission).
                                if (Rnd.NextDouble() < Params.Room.TransmissionValue[OR.Surf_ID][OR.Octave])
                                {
                                    OR.direction *= -1;
                                    Room.ReflectRay(ref OR, ref _u[OR.ThreadID], ref _v[OR.ThreadID], ref OR.Surf_ID);
                                }
                            }

                            // 3. Apply Scattering.
                            Room.Scatter_Late(ref OR, ref Rays, ref Rnd, cos_theta, u, v);
                        }
                        while (OR.t_sum < COTime && OR.Intensity > Threshold_Power_2[OR.Octave]);

                        do
                        {
                            OR.Ray_ID = Rnd.Next();
                            if (!Params.Room.shoot(OR, out u, out v, out OR.Surf_ID, out Start, out leg, out code))
                            {
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

                            // 2. Apply Transmission (if any).
                            if (Params.Room.TransmissionValue[OR.Surf_ID][OR.Octave] > 0.0)
                            {
                                //// a. Create new source for transmitted energy (E * Transmission).
                                //// b. Modify E (E * 1 - Transmission).
                                if (Rnd.NextDouble() < Params.Room.TransmissionValue[OR.Surf_ID][OR.Octave])
                                {
                                    OR.direction *= -1;
                                    Room.ReflectRay(ref OR, ref _u[OR.ThreadID], ref _v[OR.ThreadID], ref OR.Surf_ID);
                                }
                            }

                            // 3. Apply Scattering.
                            Room.Scatter_Simple(ref OR, ref Rnd, cos_theta, u, v);
                        }
                        while (OR.t_sum < COTime && OR.Intensity > Threshold_Power_3[OR.Octave]);
                    }
                    while (Rays.Count > 0);
                }
            }
            _ts = DateTime.Now - _st;
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
                for (int i = 0; i < _processorCt; i++)
                {
                    Ray_CT += (int)_currentRay[i];
                }
                if (Ray_CT == 0) return "Preparing Threads";
                _ts = new TimeSpan((long)(_ts.Ticks * (((double)(Raycount - Ray_CT)) / (Ray_CT + 1))));
                return string.Format("Calculating Ray {0} of {1}. ({2} hours,{3} minutes,{4} seconds Remaining.) Press 'Esc' to Cancel...", Ray_CT, Raycount, _ts.Hours, _ts.Minutes, _ts.Seconds);
            }
            else
            {
                foreach (Convergence_Check c in check)
                {                   
                    if (c != null && c.Check())
                    {
                        this.Abort_Calculation();
                        return "Concluding Simulation...";
                    }
                }
                int Ray_CT = 0;
                for (int i = 0; i < _processorCt; i++)
                {
                    Ray_CT += (int)_currentRay[i];
                }
                return string.Format("Ray {0} - Convergence not yet reached...", Ray_CT);
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
            
            public Convergence_Check(Receiver_Bank R, int id, int _oct)
            {
                oct = _oct;
                RunningSim = R.Rec_List[id].Recs.Energy[oct];
            }

            public abstract bool Check();
        }

        public class Detailed_Convergence_Check: Convergence_Check
        {
            //int MaxCheck;
            double[] Snapshot;
            int step;
            int binct;
            int count = 0;

            public Detailed_Convergence_Check(Receiver_Bank R, int id, int _oct)
                : base(R, id, _oct)
            {
                step = (R.SampleRate / 1000);
                binct = RunningSim.Length / step;
                Snapshot = new double[binct];
            }

            public override bool Check()
            {
                double[] sn = new double[Snapshot.Length];
                for (int i = 0; i < Snapshot.Length; i++)
                {
                    for (int j = 0; j < step; j++)
                    {
                        sn[i] += RunningSim[i*step + j];
                    }
                }

                for (int i = 0; i < sn.Length; i++)
                {
                    if (sn[i] == 0) continue;
                    double d = (sn[i] - Snapshot[i]) / Snapshot[i];
                    if (d > 0.1)
                    {
                        Snapshot = sn;
                        count = 0;
                        return false;
                    }
                }
                Snapshot = sn;
                count++;
                if (count > 10) return true;
                return false;
            }
        }

        public class Minimum_Convergence_Check: Convergence_Check
        {
            double snapshot50, snapshot80, snapshotinf;
            int SampleStart, Sample50, Sample80, SampleInf;
            int count = 0;

            public Minimum_Convergence_Check(Source S, Scene Sc, Receiver_Bank R, int id, int _oct)
                :base(R, id, _oct)
            {
                SampleStart =(int)Math.Floor((S.H_Origin() - R.Origin(id)).Length() / Sc.Sound_speed(R.Origin(id)) * R.SampleRate);
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
                    sn80 += RunningSim[i];
                }
                for (int i = Sample80; i < SampleInf; i++)
                {
                    sninf += RunningSim[i];
                }

                if (sn50 == 0
                    || sn50 / snapshot50 > 1.02
                    || sn80 == 0
                    || sn80 / snapshot80 > 1.02
                    || sninf == 0
                    || sninf / snapshotinf > 1.1)
                {
                    snapshot50 = sn50;
                    snapshot80 = sn80;
                    snapshotinf = sninf;
                    count = 0;
                    return false;
                }

                snapshot50 = sn50;
                snapshot80 = sn80;
                snapshotinf = sninf;
                count++;
                if (count > 10) return true;
                return false;
            }
        }
    }
}