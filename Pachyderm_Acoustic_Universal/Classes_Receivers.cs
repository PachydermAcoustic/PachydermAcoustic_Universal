//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2025, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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
using MathNet.Numerics;
using Pachyderm_Acoustic.Pach_Graphics;
using Pachyderm_Acoustic.Utilities;
using System.CodeDom;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using Eto.Forms;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        /// <summary>
        /// A bank for storing all receiver objects.
        /// </summary>
        [Serializable]
        public class Receiver_Bank
        {
            public Spherical_Receiver[] Rec_List;
            public Type Rec_Type;
            public double[] SWL;
            public string SrcType = "";
            public int SampleCT;
            public int SampleRate;
            public int RayCount;
            public double CutOffTime;
            public Hare.Geometry.Point Min;
            public Hare.Geometry.Point Max;
            public int[] Oct_choice;
            public double delay_ms;
            public bool Third_Oct = false;

            /// <summary>
            /// Specifies the kinds of receivers that can be placed in this receiver bank.
            /// </summary>
            public enum Type
            {
                Stationary,
                Variable
            }

            public Receiver_Bank()
            {
            }

            public Receiver_Bank(IList<Point> Pt, Source SrcPT, Scene Sc, int SampleRate_in, double COTime_in, Type Type, bool third_Oct)
            {
                SrcPT = SrcPT.ForBandResolution(third_Oct);
                SWL = SrcPT.SWL();
                SrcType = SrcPT.Type();
                Third_Oct = third_Oct;
                Oct_choice = Enumerable.Range(0, third_Oct ? 24 : 8).ToArray();
                SampleRate = SampleRate_in;
                SampleCT = (int)Math.Floor(COTime_in * SampleRate_in / 1000);
                this.CutOffTime = COTime_in;
                Rec_Type = Type;
                Point[] arrPts = Pt.ToArray<Point>();
                Rec_List = new Spherical_Receiver[arrPts.Length];
                Min = new Hare.Geometry.Point(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
                Max = new Hare.Geometry.Point(double.NegativeInfinity, double.NegativeInfinity, double.NegativeInfinity);

                for (int i = 0; i < arrPts.Length; i++)
                {
                    if (Type == Type.Stationary) Rec_List[i] = new Spherical_Receiver(arrPts[i], SrcPT.Origin, Sc.Attenuation(arrPts[i], third_Oct), Sc.Sound_speed(arrPts[i]), Sc.Rho(arrPts[i]), SampleRate_in, COTime_in);
                    if (Type == Type.Variable) Rec_List[i] = new Expanding_Receiver(arrPts[i], SrcPT.Origin, RayCount, Sc.Attenuation(arrPts[i], third_Oct), Sc.Sound_speed(arrPts[i]), Sc.Rho(arrPts[i]), SampleRate_in, COTime_in);

                    if (arrPts[i].x > Max.x) Max.x = arrPts[i].x;
                    if (arrPts[i].y > Max.y) Max.y = arrPts[i].y;
                    if (arrPts[i].z > Max.z) Max.z = arrPts[i].z;
                    if (arrPts[i].x < Min.x) Min.x = arrPts[i].x;
                    if (arrPts[i].y < Min.y) Min.y = arrPts[i].y;
                    if (arrPts[i].z < Min.z) Min.z = arrPts[i].z;
                }
                Third_Oct = third_Oct;
                Oct_choice = Enumerable.Range(0, third_Oct ? 24 : 8).ToArray();
            }

            public virtual Receiver_Bank Duplicate(Source Src, Scene Room)
            {
                SWL = Src.SWL();
                SrcType = Src.Type();
                List<Point> pts = new List<Point>();
                foreach (Spherical_Receiver s in Rec_List) pts.Add(s.Origin);
                return new Receiver_Bank(pts, Src, Room, this.SampleRate, this.CO_Time, this.Rec_Type, this.Third_Oct);
            }

            /// <summary>
            /// Private receiver bank constructer for file read-in.
            /// </summary>
            /// <param name="Rec_Count">The number of receivers specified by the user</param>
            /// <param name="SampleRate_in">the simulation histogram sampling frequency</param>
            /// <param name="CSound">the speed of sound in m/s</param>
            /// <param name="CO_Time_in">The Cut Off Time of the simulation in ms.</param>
            /// <param name="Type">The type of receiver to be used</param>
            private Receiver_Bank(int Rec_Count, int SampleRate_in, double CO_Time_in, double[] rho_c, int[] Octaves, Receiver_Bank.Type Type, bool third_Oct)
            {
                Third_Oct = third_Oct;
                Oct_choice = Enumerable.Range(0, third_Oct ? 24 : 8).ToArray();
                SampleRate = SampleRate_in;
                SampleCT = (int)Math.Ceiling(CO_Time_in * SampleRate_in / 1000);
                this.CutOffTime = CO_Time_in;
                Rec_List = new Spherical_Receiver[Rec_Count];
                Rec_Type = Type;
                for (int i = 0; i < Rec_Count; i++)
                {
                    Rec_List[i] = new Spherical_Receiver(SampleRate_in, CO_Time_in, rho_c[i], Oct_choice.Length);
                }
            }

            /// <summary>
            /// Receiver bank constructor.
            /// </summary>
            /// <param name="Pt">array of receiver Origin points</param>
            /// <param name="SrcPT">sound source point</param>
            /// <param name="Room">the acoustical scene</param>
            /// <param name="RCT">the number or rays emanating from the source</param>
            /// <param name="CSound">the speed of sound in m/s</param>
            /// <param name="SampleRate_in">the simulation histogram sampling frequency</param>
            /// <param name="COTime_in">the Cut Off Time in ms.</param>
            /// <param name="Type">the type of receivers contained in this receiver bank</param>
            public Receiver_Bank(IEnumerable<Point> Pt, Source SrcPT, Scene Sc, int SampleRate_in, double COTime_in, Type Type, bool third_Oct)
            {
                Third_Oct = third_Oct;
                Oct_choice = Enumerable.Range(0, third_Oct ? 24 : 8).ToArray();
                SrcPT = SrcPT.ForBandResolution(third_Oct);
                SWL = SrcPT.SWL();
                SampleRate = SampleRate_in;
                SampleCT = (int)Math.Floor(COTime_in * SampleRate_in / 1000);
                this.CutOffTime = COTime_in;
                Rec_Type = Type;
                Point[] arrPts = Pt.ToArray<Point>();
                Rec_List = new Spherical_Receiver[arrPts.Length];
                Min = new Hare.Geometry.Point(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
                Max = new Hare.Geometry.Point(double.NegativeInfinity, double.NegativeInfinity, double.NegativeInfinity);

                for (int i = 0; i < arrPts.Length; i++)
                {
                    if (Type == Type.Stationary) Rec_List[i] = new Spherical_Receiver(arrPts[i], SrcPT.Origin, Sc.Attenuation(arrPts[i], third_Oct), Sc.Sound_speed(arrPts[i]), Sc.Rho(arrPts[i]), SampleRate_in, COTime_in);
                    if (Type == Type.Variable) Rec_List[i] = new Expanding_Receiver(arrPts[i], SrcPT.Origin, RayCount, Sc.Attenuation(arrPts[i], third_Oct), Sc.Sound_speed(arrPts[i]), Sc.Rho(arrPts[i]), SampleRate_in, COTime_in);

                    if (arrPts[i].x > Max.x) Max.x = arrPts[i].x;
                    if (arrPts[i].y > Max.y) Max.y = arrPts[i].y;
                    if (arrPts[i].z > Max.z) Max.z = arrPts[i].z;
                    if (arrPts[i].x < Min.x) Min.x = arrPts[i].x;
                    if (arrPts[i].y < Min.y) Min.y = arrPts[i].y;
                    if (arrPts[i].z < Min.z) Min.z = arrPts[i].z;
                }
            }

            public double CO_Time
            {
                get
                {
                    return CutOffTime;
                }
                set
                {
                    CutOffTime = value;
                    if (Rec_List == null) return;
                    foreach (Spherical_Receiver S in Rec_List) S.CO_Time = value;
                }
            }

            /// <summary>
            /// this method checks all receivers in the bank for intersections with a broadband ray.
            /// </summary>
            /// <param name="Length">the distance traveled by the ray up to the reflection prior to potential receiver intersection</param>
            /// <param name="R">the broadband ray</param>
            /// <param name="EndPt">the point of intersection with the room after the potential receiver intersection</param>
            public virtual void CheckBroadbandRay(BroadRay R, Hare.Geometry.Point EndPt)
            {
                foreach (Spherical_Receiver S in Rec_List) { S.CheckBroadbandRay(R, EndPt); };
            }

            /// <summary>
            /// this method checks all receivers in the bank for intersections with a single octave ray.
            /// </summary>
            /// <param name="Length">the distance traveled by the ray up to the reflection prior to potential receiver intersection</param>
            /// <param name="R">the single octave band ray</param>
            /// <param name="EndPt">the point of intersection with the room after the potential receiver intersection</param>
            public virtual void CheckRay(OctaveRay R, Hare.Geometry.Point EndPt)
            {
                foreach (Spherical_Receiver S in Rec_List) { S.CheckRay(R, EndPt); };
            }

            /// <summary>
            /// Checks all receivers for intersections with a ray. This method does not record energy.
            /// </summary>
            /// <param name="R">the single octave band ray</param>
            /// <param name="EndPt">the point of intersection with the room after the potential receiver intersection</param>
            /// <param name="Length">the distance traveled by the ray up to the reflection prior to potential receiver intersection</param>
            /// <returns> an array of boolean values, indexed by receiver number</returns>
            public virtual bool[] SimpleCheck(BroadRay R, Hare.Geometry.Point EndPt, double Length)
            {
                bool[] B = new bool[Rec_List.Length];
                for (int i = 0; i < Rec_List.Length; i++) { B[i] = Rec_List[i].SimpleCheck(R, EndPt); };
                return B;
            }

            /// <summary>
            /// this method saves a completed receiver bank to memory.
            /// </summary>
            /// <param name="BW">the binary writer used to copy data to disk.</param>
            public void Write_Data(ref System.IO.BinaryWriter BW)
            {
                //1. Write an indicator to signify that there is ray traced data: BOOL
                BW.Write("Ray-Traced_Data");

                //2. Write the type of receivers used
                BW.Write((int)Rec_Type);

                //3. Write the sample rate:int
                BW.Write(SampleRate);

                //4. Write the number of samples:int
                BW.Write(SampleCT);

                //5. Write the Speed of Sound:double (deprecated as of v.1.6)
                //BW.Write(343.0d);

                //6. Write the cut off time:double
                BW.Write(CutOffTime);

                for (int q = 0; q < Rec_List.Length; q++)
                {
                    //7a. Write the Speed of Sound for this receiver: double
                    BW.Write(Rec_List[q].Sound_Speed);
                    double[][] Hist = new double[SWL.Length][];
                    double[][] HistP = new double[SWL.Length][];
                    for (int band = 0; band < SWL.Length; band++)
                    {
                        Hist[band] = Rec_List[q].GetEnergyHistogram(band);
                        HistP[band] = Rec_List[q].GetRMSPressureHistogram(band);
                    }

                    //7b. Write the echogram sample and directional sample:(double + 3double) * 8
                    for (int j = 0; j < SampleCT; j++)
                    {
                        for (int band = 0; band < SWL.Length; band++)
                        {

                            BW.Write(Hist[band][j]);
                            BW.Write(HistP[band][j]);

                            BW.Write(Rec_List[q].Directions_Pos(band, j, 0));
                            BW.Write(Rec_List[q].Directions_Pos(band, j, 1));
                            BW.Write(Rec_List[q].Directions_Pos(band, j, 2));
                            BW.Write(Rec_List[q].Directions_Neg(band, j, 0));
                            BW.Write(Rec_List[q].Directions_Neg(band, j, 1));
                            BW.Write(Rec_List[q].Directions_Neg(band, j, 2));


                        } 
                    }
                }
            }

            /// <summary>
            /// this method creates a complete Receiver bank from saved data.
            /// </summary>
            /// <param name="BR">the binary reader from which the data will come.</param>
            /// <param name="Rec_CT">The number of receivers</param>
            /// <param name="SampleRate">the histogram sampling frequency</param>
            /// <returns>a complete receiver bank</returns>
            public static Receiver_Bank Read_Data(ref System.IO.BinaryReader BR, double[] SWL, int Rec_CT, IEnumerable<Hare.Geometry.Point> RecPts, double[] rho_c, double delayms, ref int SampleRate, string version, IProgressFeedback VB)
            {
                //2. Write the type of receivers used
                Type Rec_typ = (Type)Enum.ToObject(typeof(Type), BR.ReadUInt32());

                //3. Write the sample rate:int
                SampleRate = BR.ReadInt32();

                //4. Write the number of samples:int
                int SampleCT = BR.ReadInt32();

                //5. Write the Speed of Sound:double (deprecated as of v1.6...)
                //double C_Sound = BR.ReadDouble();

                //6. Write the cut off time:double
                double CutoffTime = BR.ReadDouble();

                Receiver_Bank Rec = new Receiver_Bank(Rec_CT, SampleRate, CutoffTime, rho_c, new int[] { 0, 1, 2, 3, 4, 5, 6, 7 }, Rec_typ, SWL.Length == 24);

                double v = double.Parse(version.Substring(0, 3), System.Globalization.CultureInfo.InvariantCulture);

                for (int q = 0; q < Rec_CT; q++)
                {
                    Rec.Rec_List[q].Origin = RecPts.ElementAt<Hare.Geometry.Point>(q);
                    Rec.Rec_List[q].Recs = new Spherical_Receiver.Directional_Histogram(SampleRate, SampleCT, SWL.Length);
                    Rec.Rec_List[q].Reflection_Analysis_Data = new double[3][][];
                    for (int type = 0; type < 3; type++)
                    {
                        Rec.Rec_List[q].Reflection_Analysis_Data[type] = new double[SWL.Length][];
                        for (int band = 0; band < SWL.Length; band++) Rec.Rec_List[q].Reflection_Analysis_Data[type][band] = new double[SampleCT];
                    }

                    //if (v == 1.7)
                    //{
                    //    for (int i = 0; i < 16; i++) BR.ReadSingle();
                    //}
                    ////////Initiates the directional Histogram////////
                    //Rec.Rec_List[q].Recs = new Spherical_Receiver.Directional_Histogram(SampleRate, SampleCT);
                    //7a. Reads in the speed of sound.
                    Rec.Rec_List[q].Sound_Speed = BR.ReadDouble();

                    for (int j = 0; j < SampleCT; j++)
                    {
                        //7b. Write the echogram sample and directional sample:(double + 3double) * 8
                        for (int oct = 0; oct < SWL.Length; oct++)
                        {
                            Rec.Rec_List[q].Recs.Energy[oct][j] = BR.ReadDouble();
                            if (v == 1.7)
                            {
                                BR.ReadSingle(); BR.ReadSingle();
                                Rec.Rec_List[q].Recs.Add_direction(oct, j, new Vector(BR.ReadSingle(), BR.ReadSingle(), BR.ReadSingle()), new Vector(BR.ReadSingle(), BR.ReadSingle(), BR.ReadSingle()));
                            }
                            else if (v <= 2.1)
                            {
                                Rec.Rec_List[q].Recs.Pressure[oct][j] = Math.Sqrt(Rec.Rec_List[q].Recs.Energy[oct][j] * rho_c[q]);
                                Rec.Rec_List[q].Recs.Add_direction(oct, j, new Vector(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble()), new Vector(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble()));
                            }
                            else
                            {
                                //v2.5.x.x
                                Rec.Rec_List[q].Recs.Pressure[oct][j] = BR.ReadDouble();
                                Rec.Rec_List[q].Recs.Add_direction(oct, j, new Vector(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble()), new Vector(BR.ReadDouble(), BR.ReadDouble(), BR.ReadDouble()));
                            }
                        }
                    }
                }

                //Task t = new Task(() => { foreach (Receiver_Bank R in Receiver) R.Create_Filter(VB); });
                //t.Start();

                //while (! t.IsCompleted) await Task.Delay(1000);

                Rec.SampleCT = SampleCT;
                Rec.SWL = (double[])SWL.Clone();
                Rec.delay_ms = delayms;
                Rec.Create_Filter(VB);
                //while (! t.IsCompleted) await Task.Delay(1000);
                return Rec;
            }

            /// <summary>
            /// returns the raytraced energy response of a single receiver at an indexed octave band.
            /// </summary>
            /// <param name="Octave">octave band index</param>
            /// <param name="Rec_Index">the index of the receiver to be indexed</param>
            /// <returns>the energy response</returns>
            public virtual double[] GetEnergyHistogram(int Octave, double Delay_ms, int Rec_Index)
            {
                    double[] hist = (double[])Rec_List[Rec_Index].GetEnergyHistogram(Octave).Clone();
                    double[] pad = new double[(int)Math.Floor(Delay_ms / 1000.0 * SampleRate) + hist.Length];
                    Array.Copy(hist, 0, pad, (int)Math.Floor(Delay_ms / 1000.0 * SampleRate), hist.Length);
                    return pad;
            }

            public void reset_filter()
            {
                for (int i = 0; i < Rec_List.Length; i++)
                {
                    Rec_List[i].Recs.F = null;
                }
            }

            public virtual bool HasFilter()
            {
                return HasFilter(0);
            }

            public virtual bool HasFilter(int Rec_Index)
            {
                while (CreatingFilter) System.Threading.Thread.Sleep(500);
                return Rec_List[Rec_Index].Recs.HasFilter();
            }

            public virtual void GetFilter(int Rec_Index, out double[] Filter)
            {
                if (Rec_List[Rec_Index].Recs.F == null) Filter = new double[SampleCT];
                else
                {
                    Filter = (double[])Rec_List[Rec_Index].Recs.GetFilter().Clone();
                }
            }

            /// <summary>
            /// get a histogram of vectors indicating the integrated direction of sound at any given time.
            /// </summary>
            /// <param name="Octave">the indexed octave band</param>
            /// <param name="t">time index</param>
            /// <param name="Rec_Index">the index of the receiver</param>
            /// <returns>the integrated direction of sound at the chosen frequency and time</returns>
            public virtual Vector Directions_Pos(int Octave, int t, int Rec_Index)
            {
                Vector Dir = new Vector(Rec_List[Rec_Index].Directions_Pos(Octave, t, 0), Rec_List[Rec_Index].Directions_Pos(Octave, t, 1), Rec_List[Rec_Index].Directions_Pos(Octave, t, 2));
                return Dir;
            }

            /// <summary>
            /// get a histogram of vectors indicating the integrated direction of sound at any given time.
            /// </summary>
            /// <param name="Octave">the indexed octave band</param>
            /// <param name="t">time index</param>
            /// <param name="Rec_Index">the index of the receiver</param>
            /// <returns>the integrated direction of sound at the chosen frequency and time</returns>
            public virtual Vector Directions_Neg(int Octave, int t, int Rec_Index)
            {
                Vector Dir = new Vector(Rec_List[Rec_Index].Directions_Neg(Octave, t, 0), Rec_List[Rec_Index].Directions_Neg(Octave, t, 1), Rec_List[Rec_Index].Directions_Neg(Octave, t, 2));
                return Dir;
            }

            public virtual Vector Directions_Pos(int Octave, int t, int Rec_Index, double alt, double azi, bool degrees)
            {
                Vector V = new Vector(Rec_List[Rec_Index].Directions_Pos(Octave, t).dx, Rec_List[Rec_Index].Directions_Pos(Octave, t).dy, Rec_List[Rec_Index].Directions_Pos(Octave, t).dz);
                return Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V, azi, 0, degrees), 0, alt, degrees);
            }

            public virtual Vector Directions_Neg(int Octave, int t, int Rec_Index, double alt, double azi, bool degrees)
            {
                Vector V = new Vector(Rec_List[Rec_Index].Directions_Neg(Octave, t).dx, Rec_List[Rec_Index].Directions_Neg(Octave, t).dy, Rec_List[Rec_Index].Directions_Neg(Octave, t).dz);
                return Utilities.PachTools.Rotate_Vector(Utilities.PachTools.Rotate_Vector(V, azi, 0, degrees), 0, alt, degrees);
            }

            public virtual Vector Directions_Pos(int Octave, int t, int Rec_Index, Vector V)
            {
                double l = Math.Sqrt(V.dz * V.dz + V.dx * V.dx);
                double azi = Math.Asin(V.dy / l);
                double alt = Math.Atan2(V.dx, V.dz);
                return Directions_Pos(Octave, t, Rec_Index, alt, azi, false);
            }

            public virtual Vector Directions_Neg(int Octave, int t, int Rec_Index, Vector V)
            {
                double l = Math.Sqrt(V.dz * V.dz + V.dx * V.dx);
                double azi = Math.Asin(V.dy / l);
                double alt = Math.Atan2(V.dx, V.dz);
                return Directions_Neg(Octave, t, Rec_Index, alt, azi, false);
            }

            /// <summary>
            /// combines all threadlocal results
            /// </summary>
            /// <param name="Rec">a receiver bank clone to combine with this receiver bank instance</param>
            public virtual void Combine_Clones(Receiver_Bank Rec)
            {
                if (Third_Oct != Rec.Third_Oct || Rec.Count != Count || Rec.Duration() != Duration()) throw new ArgumentException("Receiver banks must have matching bands, receivers and duration.");
                for (int R = 0; R < Rec.Count; R++)
                {
                    for (int oct = 0; oct < SWL.Length; oct++)
                    {
                        for (int T = 0; T < Rec.Duration(); T++)
                        {
                            Rec_List[R].Combine_Sample(T, Rec.Rec_List[R].Energy(T, oct), Rec.Rec_List[R].Pressure_rms(T, oct), Rec.Rec_List[R].Directions_Pos(oct, T), Rec.Rec_List[R].Directions_Neg(oct, T), oct);
                            if (Rec_List[R].Reflection_Analysis_Data != null && Rec.Rec_List[R].Reflection_Analysis_Data != null)
                                for (int type = 0; type < 3; type++) Rec_List[R].Reflection_Analysis_Data[type][oct][T] += Rec.Rec_List[R].Reflection_Analysis_Data[type][oct][T];
                        }
                    }
                }
            }

            /// <summary>
            /// Gets the Origin of a receiver in the bank
            /// </summary>
            /// <param name="ID">the index or 'ID' of a receiver in the bank</param>
            /// <returns>the Origin point in Hare format</returns>
            public Hare.Geometry.Point Origin(int ID)
            {
                return Rec_List[ID].Origin;
            }

            public Point[] Origins()
            {
                Point[] P = new Point[Rec_List.Length];
                for (int i = 0; i < Rec_List.Length; i++)
                {
                    P[i] = Rec_List[i].Origin;
                }
                return P;
            }

            /// <summary>
            /// obtains the number of samples in histograms owned by this instance.
            /// </summary>
            /// <returns>Number of histogram samples</returns>
            public int Duration()
            {
                return SampleCT;
            }

            /// <summary>
            /// Gets the number of receivers in this bank.
            /// </summary>
            public virtual int Count
            {
                get
                {
                    return Rec_List.Length;
                }
            }

            /// <summary>
            /// Combines two Map Receivers with data in differing octave bands.
            /// </summary>
            /// <param name="A"></param>
            /// <param name="B"></param>
            /// <returns></returns>
            public static Receiver_Bank operator *(Receiver_Bank A, Receiver_Bank B)
            {
                Point[] a_pts = A.Origins();
                Point[] b_pts = B.Origins();
                if (A.Rec_Type != B.Rec_Type || a_pts.Length != b_pts.Length)
                {
                    Eto.Forms.MessageBox.Show("Data is for two different calculations. Simulations not Combined.");
                    return null;
                }

                for (int i = 0; i < a_pts.Length; i++)
                {
                    if (a_pts[i].GetHashCode() != b_pts[i].GetHashCode())
                    {
                        Eto.Forms.MessageBox.Show("Data is for two different calculations. Simulations not Combined.");
                        return null;
                    }
                }

                foreach (int a in A.Oct_choice)
                {
                    foreach (int b in B.Oct_choice)
                    {
                        if (a == b)
                        {
                            Eto.Forms.MessageBox.Show("Data Conflicts. Simulations not Combined.");
                            return null;
                        }
                    }
                }

                foreach (int oct in B.Oct_choice)
                {
                    for (int rec = 0; rec < a_pts.Length; rec++)
                    {
                        for (int t = 0; t < A.Rec_List[rec].Recs.SampleCT; t++) A.Rec_List[rec].Combine_Sample(t, B.Rec_List[rec].Energy(t, oct), B.Rec_List[rec].Pressure_rms(t, oct), B.Rec_List[rec].Directions_Pos(oct, t), B.Rec_List[rec].Directions_Neg(oct, t), oct);
                    }
                }

                return A;
            }

            private bool CreatingFilter = false;

            public async virtual void Create_Filter(IProgressFeedback VB = null)
            {
                CreatingFilter = true;
                //if (VB != null) VB.Report((int)(100 * (1f - ((float)CDE.CurrentCount / (float)IR.Length))));
                //if (VB != null) VB.change_title(string.Format("Filter Progress: Creating Positive {0}.", new string[3] { "X", "Y", "Z" }[dir]));
                //if (VB != null) VB.change_title(string.Format("Filter Progress: Creating Negative {0}.", new string[3] { "X", "Y", "Z" }[dir]));

                //VB.change_title(string.Format("Filter Progress: Creating Positive {0}.", new string[3] { "X", "Y", "Z" }[dir]));
                //VB.change_title(string.Format("Filter Progress: Creating Negative {0}.", new string[3] { "X", "Y", "Z" }[dir]));

                Task[] t = new Task[Rec_List.Length];

                for (int rec = 0; rec < Rec_List.Length; rec++)
                {
                    if (VB != null)
                    {
                        VB.change_title(string.Format("Creating Pressure for receiver {0}", rec));
                        VB.Report((int)(((double)(rec + 1) * 100) / ((double)Rec_List.Length)));
                    }
                    //t[rec] = Task.Run(() => { Rec_List[rec].Create_Filter(); });
                    //while (t[rec].Status != TaskStatus.RanToCompletion) System.Threading.Thread.Sleep(300);//await Task.Delay(500);
                    Rec_List[rec].Create_Filter();
                }
                //await Task.WhenAll(t);
                CreatingFilter = false;
            }

            public double[] Create_Filter(double[] SWL, int Rec_ID, int dim, IProgressFeedback VB = null)
            {
                CreatingFilter = true;
                return Rec_List[Rec_ID].Create_Filter(SWL, dim);
            }

            public double[][] Create_Filter(double[] SWL, int Rec_ID, IProgressFeedback VB = null)
            {
                CreatingFilter = true;
                return Rec_List[Rec_ID].Create_Filter(SWL);
            }

            public double[][] Filter_3Axis(int rec_id)
            {
                return Rec_List[rec_id].Filter3Axis();
            }

            public virtual void Set_Power(double[] factor)
            {
                if (factor == null || factor.Length != SWL.Length) throw new ArgumentException("Power factors must match the receiver band count.");
                for(int i = 0; i < SWL.Length; i++) SWL[i] = AcousticalMath.SPL_Intensity(AcousticalMath.Intensity_SPL(SWL[i]) * factor[i]);
                for (int i = 0; i < Rec_List.Length; i++) Rec_List[i].Set_Power(factor);
            }

            public void Scale(double Ray_CT)
            {
                for (int i = 0; i < Rec_List.Length; i++) Rec_List[i].Scale(Ray_CT);
            }

            /// <summary>
            /// The direct sound outputs a power mod factor. This takes that function's place, since there is no direct sound for this simulation type.
            /// </summary>
            /// <param name="new_SWL"></param>
            /// <returns></returns>
            public double[] PowerModFactor(double[] new_SWL)
            {
                double[] factor = new double[SWL.Length];
                if (this.SWL == null) SWL = new double[8] { 120, 120, 120, 120, 120, 120, 120, 120 };
                for (int i = 0; i < SWL.Length; i++) factor[i] = Utilities.AcousticalMath.Intensity_SPL(new_SWL[i]) / Utilities.AcousticalMath.Intensity_SPL(this.SWL[i]);
                return factor;
            }
        }

        /// <summary>
        /// 1m. stationary receiver
        /// </summary>
        [Serializable]
        public class Spherical_Receiver
        {
            public int[] Ray_ID;
            public Hare.Geometry.Point Origin;
            protected double C_Sound;
            public double Rho_C;
            protected double Inv_C_Sound;
            protected double[] Atten = new double[8];
            protected internal Histogram Recs;
            protected double D_Length;
            public int SampleRate;
            public double CO_Time;
            public double Radius;
            protected internal double Radius2;
            protected internal double SizeMod;
            public double[][][] Reflection_Analysis_Data;

        //    public int BandCount { get { return Recs.Energy.Length; } }
        //public void ConfigureBands(bool thirdOctave, double[] attenuation = null)
        //{
        //    if (attenuation != null) Atten = (double[])attenuation.Clone();
        //    if (Recs != null && Recs.Energy.Length == (thirdOctave ? 24 : 8)) return;
        //    Recs = thirdOctave ? (Histogram)new Directional_Histogram_ThirdOctave(SampleRate, CO_Time) : new Directional_Histogram(SampleRate, CO_Time);
        //    ConfigureReflectionAnalysis();
        //}

        //protected internal void ConfigureReflectionAnalysis()
        //{
        //    Reflection_Analysis_Data = new double[3][][];
        //    for (int type = 0; type < 3; type++)
        //    {
        //        Reflection_Analysis_Data[type] = new double[Recs.Energy.Length][];
        //        for (int b = 0; b < Recs.Energy.Length; b++) Reflection_Analysis_Data[type][b] = new double[Recs.SampleCT];
        //    }
        //}

        public Spherical_Receiver()
            { }

            /// <summary>
            /// Use this constructor when reading in a file...
            /// </summary>
            /// <param name="SampleRate_in">The sampling frequency at which detections will be recorded. Enter 0 for raw detection recording. Watch out for memory overflow...</param>
            /// <param name="CSound_in"> the speed of sound</param>
            /// <param name="CO_Time_in">the Cut Off Time in ms.</param>
            public Spherical_Receiver(int SampleRate_in, double CO_Time_in, double rho_c)
                : this(SampleRate_in, CO_Time_in, rho_c, 8)
            { }

            public Spherical_Receiver(int SampleRate_in, double CO_Time_in, double rho_c, int bandCount)
            {
                if (bandCount != 8 && bandCount != 24) throw new ArgumentOutOfRangeException(nameof(bandCount));
                Atten = new double[bandCount];
                Rho_C = rho_c;
                CO_Time = CO_Time_in;
                SampleRate = SampleRate_in;
                Recs = new Directional_Histogram(SampleRate, CO_Time, bandCount);
                Reflection_Analysis_Data = new double[3][][];
                Reflection_Analysis_Data[0] = new double[bandCount][];
                Reflection_Analysis_Data[1] = new double[bandCount][];
                Reflection_Analysis_Data[2] = new double[bandCount][];
                for(int oct = 0; oct < bandCount; oct++) { Reflection_Analysis_Data[0][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; Reflection_Analysis_Data[1][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; Reflection_Analysis_Data[2][oct] = new double[(int)(SampleRate * (CO_Time / 1000))]; }
            }

            /// <summary>
            /// Constructor which takes Rhino point input.
            /// </summary>
            /// <param name="Center"></param>
            /// <param name="SrcCenter"></param>
            /// <param name="room"></param>
            /// <param name="RCT"></param>
            /// <param name="C_Sound_in"></param>
            /// <param name="SampleRate_in"></param>
            /// <param name="COTime_in"></param>
            public Spherical_Receiver(Point Center, Point SrcCenter, double[] Attenuation, double C_Sound_in, double rho, int SampleRate_in, double COTime_in)
            {
                if (Attenuation == null || (Attenuation.Length != 8 && Attenuation.Length != 24)) throw new ArgumentException("Provide an octave or third-octave attenuation spectrum.", nameof(Attenuation));
                D_Length = (Center - SrcCenter).Length();
                Radius = 1;
                Radius2 = Radius * Radius;
                SampleRate = SampleRate_in;
                CO_Time = COTime_in;
                Origin = new Hare.Geometry.Point(Center.x, Center.y, Center.z);
                Sound_Speed = C_Sound_in;
                Rho_C = C_Sound_in * rho;
                Atten = (double[])Attenuation.Clone();
                SizeMod = 1 / Math.PI;
                Recs = new Directional_Histogram(SampleRate, CO_Time, Atten.Length);
                Reflection_Analysis_Data = new double[3][][];
                Reflection_Analysis_Data[0] = new double[Atten.Length][];
                Reflection_Analysis_Data[1] = new double[Atten.Length][];
                Reflection_Analysis_Data[2] = new double[Atten.Length][];
                for (int oct = 0; oct < Atten.Length; oct++) { Reflection_Analysis_Data[0][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; Reflection_Analysis_Data[1][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; Reflection_Analysis_Data[2][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; }
            }

            /// <summary>
            /// Constructor which takes Hare Point input.
            /// </summary>
            /// <param name="Center"></param>
            /// <param name="SrcCenter"></param>
            /// <param name="room"></param>
            /// <param name="RCT"></param>
            /// <param name="C_Sound_in"></param>
            /// <param name="SampleRate_in"></param>
            /// <param name="COTime_in"></param>
            public Spherical_Receiver(Hare.Geometry.Point Center, Hare.Geometry.Point SrcCenter, Scene room, double C_Sound_in, int SampleRate_in, double COTime_in)
            {
                D_Length = Math.Sqrt(Center.x * SrcCenter.x + Center.y * SrcCenter.y * Center.z * SrcCenter.z);
                Radius = 1;
                Radius2 = Radius * Radius;
                CO_Time = COTime_in;
                SampleRate = SampleRate_in;
                Origin = Center;
                Sound_Speed = C_Sound_in;
                Atten = room.Attenuation(Center);
                SizeMod = 1 / Math.PI;
                Recs = new Directional_Histogram(SampleRate, CO_Time);
                Reflection_Analysis_Data = new double[3][][];
                Reflection_Analysis_Data[0] = new double[8][];
                Reflection_Analysis_Data[1] = new double[8][];
                Reflection_Analysis_Data[2] = new double[8][];
                for (int oct = 0; oct < 8; oct++) { Reflection_Analysis_Data[0][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; Reflection_Analysis_Data[1][oct] = new double[(int)(SampleRate * (CO_Time / 1000.0))]; Reflection_Analysis_Data[2][oct] = new double[(int)(SampleRate * (CO_Time / 1000))]; }
            }

            /// <summary>
            /// This method checks a receiver for a ray using a Broadband Ray (these typically occur before rays are split in the Raytracing simulation.
            /// </summary>
            /// <param name="R">the ray.</param>
            /// <param name="EndPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            public virtual void CheckBroadbandRay(BroadRay R, Hare.Geometry.Point EndPt)
            {
                double mx = R.x - Origin.x;
                double my = R.y - Origin.y;
                double mz = R.z - Origin.z;
                double b = Hare_math.Dot(mx, my, mz, R.dx, R.dy, R.dz);
                double c = Hare_math.Dot(mx, my, mz, mx, my, mz)  - Radius2;
                if (c > 0 && b > 0) return;
                double discr = b * b - c;
                if (discr < 0) return;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);
                double t = (t1 + t2) * .5;
                if (t > 0 && t * t < SqDistance(EndPt.x,EndPt.y, EndPt.z, R.x, R.y, R.z))
                {
                    double tsphere = (t2 - t1) * 0.5;
                    double RayTime = t * Inv_C_Sound + R.t_sum;
                    double raydist = RayTime * C_Sound;
                    foreach (int group in R.Freq_Bands)
                    for (int component = 0; component < (R.ThirdOctave ? 3 : 1); component++)
                    {
                        int oct = R.ThirdOctave ? 3 * group + component : group;
                        double energy = R.Energy[oct] * Math.Pow(10, -.1 * Atten[oct] * raydist) * SizeMod * tsphere; //* t
                        Recs.Add(RayTime, energy, -R.dx, -R.dy, -R.dz, Rho_C, oct);
                        int sample = (int)(RayTime * SampleRate);
                        if (sample >= 0 && sample < Reflection_Analysis_Data[0][oct].Length)
                            lock (Reflection_Analysis_Data[0][oct]) Reflection_Analysis_Data[0][oct][sample] += energy;
                    }
                }
            }

            /// <summary>
            /// Checks receiver for an intersection with a ray
            /// </summary>
            /// <param name="Length">The length of the ray at the reflection point before potential intersection</param>
            /// <param name="R">the ray.</param>
            /// <param name="EndPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            public virtual void CheckRay(OctaveRay R, Hare.Geometry.Point EndPt)
            {
                double mx = R.x - Origin.x;
                double my = R.y - Origin.y;
                double mz = R.z - Origin.z;
                double b = Hare_math.Dot(mx, my, mz, R.dx, R.dy, R.dz);
                double c = Hare_math.Dot(mx, my, mz, mx, my, mz) - Radius2;
                if (c > 0 && b > 0) return;
                double discr = b * b - c;
                if (discr < 0) return;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);

                double t = (t1 + t2) *.5;
                if (t > 0 && t * t < SqDistance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z))
                {
                    double tsphere = (t2 - t1) * 0.5;
                    double Raytime = t * Inv_C_Sound + R.t_sum;
                    if (R.BandCount > 1)
                    {
                        int thirdSample = (int)(Raytime * SampleRate);
                        for (int cband = 0; cband < R.BandCount; cband++)
                        {
                            int band = R.BandStart + cband;
                            double e = (R.Intensity * R.Fractions[cband]) * Math.Pow(10, -.1 * Atten[band] * Raytime * C_Sound) * SizeMod * tsphere;
                            Recs.Add(Raytime, e, -R.dx, -R.dy, -R.dz, Rho_C, band);
                            int kind = R.IsSpecular ? (R.WasScattered ? 2 : 0) : 1;
                            if (thirdSample >= 0 && thirdSample < Reflection_Analysis_Data[kind][band].Length)
                                lock (Reflection_Analysis_Data[kind][band]) Reflection_Analysis_Data[kind][band][thirdSample] += e;
                        }
                        return;
                    }
                    double energy = R.Intensity * Math.Pow(10, -.1 * Atten[R.Octave] * Raytime * C_Sound) * SizeMod * tsphere; // t
                    Recs.Add(Raytime, energy, -R.dx, -R.dy, -R.dz, Rho_C, R.Octave);
                    int sample = (int)(Raytime * SampleRate);
                    if (sample < Reflection_Analysis_Data[0][0].Length)
                    {
                        if (R.IsSpecular)
                        {
                            if (R.WasScattered)
                            {
                                Reflection_Analysis_Data[2][R.Octave][sample] += energy;
                            }
                            else
                            {
                                Reflection_Analysis_Data[0][R.Octave][sample] += energy;
                            }
                        }
                        else
                        {
                            Reflection_Analysis_Data[1][R.Octave][sample] += energy;
                        }
                    }
                }
            }

            /// <summary>
            /// This method checks a receiver for an intersection with a ray.
            /// </summary>
            /// <param name="R">the ray.</param>
            /// <param name="EndPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            /// <param name="Length">The length of the ray at the reflection point before potential intersection</param>
            /// <returns>true if intersects, false if not.</returns>
            public virtual bool SimpleCheck(BroadRay R, Hare.Geometry.Point EndPt)
            {
                double mx = R.x - Origin.x, my = R.y - Origin.y, mz = R.z - Origin.z;
                double b = Hare_math.Dot(mx,my, mz, R.dx, R.dy, R.dz);
                double c = Hare_math.Dot(mx, my, mz, mx, my, mz) - Radius2;
                if (c > 0 && b > 0) return false;
                double discr = b * b - c;
                if (discr < 0) return false;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);
                double t = (t1 + t2) * .5;
                if (t > 0 && t * t < SqDistance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z)) return true;
                return false;
            }

            /// <summary>
            /// A quick way of getting a relative distance between two points. The equivalent of the pythagorean theorem without the last square root.
            /// </summary>
            /// <param name="StartPt">Point 1</param>
            /// <param name="EndPt">Point 2</param>
            /// <returns>the square of the distance</returns>
            protected double SqDistance(Hare.Geometry.Point StartPt, Hare.Geometry.Point EndPt)
            {
                double dx = EndPt.x - StartPt.x, dy = EndPt.y - StartPt.y, dz = EndPt.z - StartPt.z;
                return dx * dx + dy * dy + dz * dz;
            }

            protected double SqDistance(double x1, double y1, double z1, double x2, double y2, double z2)
            {
                double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
                return dx * dx + dy * dy + dz * dz;
            }

            /// <summary>
            /// returns the energy response of the simulation.
            /// </summary>
            /// <param name="Octave">an octave band index</param>
            /// <returns></returns>
            public double[] GetEnergyHistogram(int Octave)
            {
                return Recs.GetEnergyHistogram(Octave);
            }

            public double[][][] Analysis_Histogram
            {
                get
                {
                    return Reflection_Analysis_Data;
                }
            }

            /// <summary>
            /// returns the root mean square pressure histogram of the simulation.
            /// </summary>
            /// <param name="Octave">an octave band index</param>
            /// <returns></returns>
            public double[] GetRMSPressureHistogram(int Octave)
            {
                return Recs.GetRMSPressureHistogram(Octave);
            }

            /// <summary>
            /// returns the direction of sound at a given time.
            /// </summary>
            /// <param name="Octave">an octave band index</param>
            /// <param name="t">a point in time</param>
            /// <returns></returns>
            public virtual Vector Directions_Pos(int Octave, int t)
            {
                return Recs.Directions_Pos(Octave, t);
            }

            /// <summary>
            /// returns the direction of sound at a given time.
            /// </summary>
            /// <param name="Octave">an octave band index</param>
            /// <param name="t">a point in time</param>
            /// <returns></returns>
            public virtual Vector Directions_Neg(int Octave, int t)
            {
                return Recs.Directions_Neg(Octave, t);
            }

            ///// <summary>
            ///// returns the direction of sound at a given time.
            ///// </summary>
            ///// <param name="Octave">an octave band index</param>
            ///// <param name="t">a point in time</param>
            ///// <param name="c">0 for x, 1 for y, 2 for z</param>
            /// <returns></returns>
            public virtual double Directions_Pos(int Octave, int t, int c)
            {
                return Recs.Directions_Pos(Octave, t, c);
            }

            ///// <summary>
            ///// returns the direction of sound at a given time.
            ///// </summary>
            ///// <param name="Octave">an octave band index</param>
            ///// <param name="t">a point in time</param>
            ///// <param name="c">0 for x, 1 for y, 2 for z</param>
            /// <returns></returns>
            public virtual double Directions_Neg(int Octave, int t, int c)
            {
                return Recs.Directions_Neg(Octave, t, c);
            }

            /// <summary>
            /// An energy value at a given time and frequency.
            /// </summary>
            /// <param name="t">the index of the sample.</param>
            /// <param name="oct">the index of the octave band</param>
            /// <returns></returns>
            public virtual double Energy(int t, int oct)
            {
                if (oct < Recs.Energy.Length) return Recs.Energy[oct][t];
                double sum = 0;
                for (int band = 0; band < Recs.Energy.Length; band++) sum += Recs.Energy[band][t];
                return sum;
            }

            /// <summary>
            /// A root mean square pressure value at a given time and frequency.
            /// </summary>
            /// <param name="t">the index of the sample.</param>
            /// <param name="oct">the index of the octave band</param>
            /// <returns></returns>
            public virtual double Pressure_rms(int t, int oct)
            {
                if (oct < Recs.Pressure.Length) return Recs.Pressure[oct][t];
                double sum = 0;
                for (int band = 0; band < Recs.Pressure.Length; band++)
                    sum += Recs.Pressure[band][t];
                return sum;
            }

            /// <summary>
            /// Add an energy entry to this receiver's histograms.
            /// </summary>
            /// <param name="Distance">the distance traveled from the source to this receiver detection</param>
            /// <param name="Energy_in">the energy to be recorded</param>
            /// <param name="direction">The incoming direction of the ray</param>
            /// <param name="Octave">the octave band index</param>
            public void Add(double Distance, double Energy_in, Vector direction, int Octave)
            {
                Recs.Add(Distance, Energy_in, direction.dx, direction.dy, direction.dz, Rho_C, Octave);
            }

            /// <summary>
            /// Add an energy entry to this receiver's histograms using the actual sample point.
            /// </summary>
            /// <param name="Sample"></param>
            /// <param name="Energy_in"></param>
            /// <param name="direction"></param>
            /// <param name="Octave"></param>
            public void Add(int Sample, double Energy_in, Vector direction, int Octave)
            {
                Vector positive = new Vector(Math.Max(0, direction.dx) * Energy_in, Math.Max(0, direction.dy) * Energy_in, Math.Max(0, direction.dz) * Energy_in);
                Vector negative = new Vector(Math.Min(0, direction.dx) * Energy_in, Math.Min(0, direction.dy) * Energy_in, Math.Min(0, direction.dz) * Energy_in);
                Recs.Combine_Sample(Sample, Energy_in, Math.Sqrt(Math.Max(0, Energy_in) * Rho_C), positive, negative, Octave);
            }

            public void Combine_Sample(int sample, double Energy_in, double Pressure_in, Vector direction_pos, Vector direction_neg, int Octave)// float p_real, float p_imag,
            {
                Recs.Combine_Sample(sample, Energy_in, Pressure_in, direction_pos, direction_neg, Octave);// p_real, p_imag,
            }

            public double Sound_Speed
            {
                get
                {
                    return C_Sound;
                }
                set
                {
                    C_Sound = value;
                    Inv_C_Sound = 1 / value;
                }
            }

            public virtual void Create_Filter()
            {
                Recs.Create_Filter(Rho_C);
            }

            public virtual double[] Create_Filter(double[] SWL, int dim)
            {
                return Recs.Create_Filter(SWL, Rho_C, dim);
            }

            public virtual double[][] Create_Filter(double[] SWL)
            {
                return Recs.Create_Filter(SWL, Rho_C);
            }

            public double[][] Filter3Axis()
            {
                if (Recs is Directional_Histogram)
                {
                    return (Recs as Directional_Histogram).Fdir;
                }
                else return null;
            }

            public virtual void Set_Power(double[] factor)
            {
                Recs.Set_Power(factor);
                if (Reflection_Analysis_Data != null)
                    for (int type = 0; type < 3; type++) for (int band = 0; band < Reflection_Analysis_Data[type].Length; band++)
                        for (int t = 0; t < Reflection_Analysis_Data[type][band].Length; t++) Reflection_Analysis_Data[type][band][t] *= factor[band];
            }

            public void Scale(double ray_ct)
            {
                Recs.Scale(ray_ct);
                if (Reflection_Analysis_Data == null) return;
                for (int i = 0; i < Reflection_Analysis_Data.Length; i++)
                {
                    for (int j = 0; j < Reflection_Analysis_Data[i].Length; j++)
                    {
                        for (int k = 0; k < Reflection_Analysis_Data[i][j].Length; k++)
                        {
                            Reflection_Analysis_Data[i][j][k] /= ray_ct;
                        }
                    }
                }
            }

            /// <summary>
            /// Class used to store the energy histogram of the receiver.
            /// </summary>
            [Serializable]
            protected internal class Histogram
            {
                protected internal double[][] Energy;
                protected internal double[][] Pressure;
                protected internal double[] F;
                protected internal double CO_Time;
                protected internal int SampleRate;
                protected internal int SampleCT;

                public Histogram(int SampleRate_in, int SampleCT)
                    : this(SampleRate_in, SampleCT, 8)
                { }

                public Histogram(int SampleRate_in, int SampleCT, int bandCount)
                {
                    SampleRate = SampleRate_in;
                    this.SampleCT = SampleCT;
                    CO_Time = (double)SampleCT / SampleRate;
                    Energy = new double[bandCount][];
                    Pressure = new double[Energy.Length][];

                    for (int o = 0; o < Energy.Length; o++)
                    {
                        Energy[o] = new double[SampleCT];
                        Pressure[o] = new double[SampleCT];
                    }
                }

                public Histogram(int SampleRate_in, double COTime_in)
                    : this(SampleRate_in, COTime_in, 8)
                { }

                public Histogram(int SampleRate_in, double COTime_in, int bandCount)
                {
                    SampleRate = SampleRate_in;
                    CO_Time = COTime_in / 1000;
                    SampleCT = (int)(CO_Time * SampleRate);
                    Energy = new double[bandCount][];
                    Pressure = new double[Energy.Length][];

                    for (int o = 0; o < Energy.Length; o++)
                    {
                        Energy[o] = new double[SampleCT];
                        Pressure[o] = new double[SampleCT];
                    }
                }

                public virtual string Type()
                {
                    return "Type;Map_Data_NoDir";
                }

                /// <summary>
                /// Adds energy detections to this histogram.
                /// </summary>
                /// <param name="Sample">the time point of arrival</param>
                /// <param name="Energy_in">the energy to be recorded</param>
                /// <param name="direction">the incoming direction of the ray</param>
                /// <param name="Octave">the index of the octave band.</param>
                public virtual void Add(double time, double Energy_in, double dx, double dy, double dz, double Rho_C, int Octave)
                {
                    lock (Energy[Octave])
                    {
                    int sample = (int)(time * SampleRate);
                    if (sample >= Energy[Octave].Length) return;
                    if (sample < 0) return;
                    Energy[Octave][sample] += Energy_in;
                    Pressure[Octave][sample] += Math.Sqrt(Math.Max(0, Energy_in) * Rho_C);

                    }
                }

                public virtual void Combine_Sample(int Sample, double Energy_in, double Pressure_in, Vector direction_pos, Vector direction_neg, int Octave)
                {
                    lock (Energy[Octave])
                    {
                    if (Sample >= Energy[Octave].Length) return;
                    Energy[Octave][Sample] += Energy_in;
                    Pressure[Octave][Sample] += Pressure_in;

                    }
                }

                public virtual void Create_Filter(double Rho_C)
                {
                    if (Energy[0].Length == 1)
                    {
                        return;
                    }

                    F = Audio.Pach_SP.ETCToFilter(this.Pressure, Enumerable.Repeat(120.0, Energy.Length).ToArray(), this.SampleRate, 44100);
                    Array.Resize(ref F, F.Length - 4096);
                }

                public virtual double[][] Create_Filter(double[] SWL, double Rho_C)
                {
                    if (Energy[0].Length == 1)
                    {
                        return null;
                    }

                    double[][] F = new double[1][] { Audio.Pach_SP.ETCToFilter(this.Pressure, SWL, this.SampleRate, 44100) };
                    Array.Resize(ref F[0], F[0].Length - 4096);
                    return F;
                }

                public virtual double[] Create_Filter(double[] SWL, double Rho_C, int dim)
                {
                    if (Energy[0].Length == 1)
                    {
                        return null;
                    }

                    double[] F = Audio.Pach_SP.ETCToFilter(this.Pressure, SWL, this.SampleRate, 44100);
                    Array.Resize(ref F, F.Length - 4096);
                    return F;
                }

                public virtual void Scale(double ray_ct)
                {
                    for (int oct = 0; oct < Energy.Length; oct++) for (int i = 0; i < Energy[oct].Length; i++)
                        {
                            Pressure[oct][i] /= ray_ct;
                            Energy[oct][i] /= ray_ct;
                        }
                }

                public virtual void Set_Power(double[] factor)
                {
                    double[] p_factor = new double[Energy.Length];
                    for(int oct = 0; oct < Energy.Length; oct++) p_factor[oct] = Math.Sqrt(factor[oct]);
                    for (int i = 0; i < Energy.Length; i++) for (int j = 0; j < Energy[i].Length; j++)
                        {
                            Energy[i][j] *= factor[i];
                            Pressure[i][j] *= p_factor[i];
                        }
                }

                /// <summary>
                /// Returns the root mean square pressure response of the simulation.
                /// </summary>
                /// <param name="Octave">the octave band index</param>
                /// <returns></returns>
                public virtual double[] GetRMSPressureHistogram(int Octave)
                {
                    if (Octave < Energy.Length)
                    {
                        return Pressure[Octave];
                    }
                    else
                    {
                        double[] Sum_Pressure = new double[(Energy[0].Length)];
                        for (int Oct = 0; Oct < Energy.Length; Oct++)
                        {
                            for (int idx = 0; idx < Pressure[Oct].Length; idx++)
                            {
                                Sum_Pressure[idx] += Pressure[Oct][idx];
                            }
                        }
                        return Sum_Pressure;
                    }
                }

                /// <summary>
                /// Returns the energy response of the simulation.
                /// </summary>
                /// <param name="Octave">the octave band index</param>
                /// <returns></returns>
                public virtual double[] GetEnergyHistogram(int Octave)
                {
                    if (Octave < Energy.Length)
                    {
                        return Energy[Octave];
                    }
                    else
                    {
                        double[] Sum_Energy = new double[(Energy[0].Length)];
                        for (int Oct = 0; Oct < Energy.Length; Oct++)
                        {
                            for (int idx = 0; idx < Energy[Oct].Length; idx++)
                            {
                                Sum_Energy[idx] += Energy[Oct][idx];
                            }
                        }
                        return Sum_Energy;
                    }
                }

                public virtual double[] GetFilter()
                {
                    return F;
                }

                public virtual double[][] GetFilter3Axis(int rec_id)
                {
                    return null;
                }

                public virtual void Divide(double divisor)
                {
                    for (int i = 0; i < Energy.Length; i++)
                    {
                        for (int j = 0; j < Energy[i].Length; j++)
                        {
                            Energy[i][j] /= divisor;
                        }
                    }
                }

                /// <summary>
                /// Direction addiion empty for this operation.
                /// </summary>
                /// <param name="octave"></param>
                /// <param name="sample"></param>
                /// <param name="V"></param>
                public virtual void Add_direction(int octave, int sample, Vector VPos, Vector VNeg)
                {
                    return;
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <returns></returns>
                public virtual Vector Directions_Pos(int oct, int t)
                {
                    return new Vector(0,0,-1);
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <returns></returns>
                public virtual Vector Directions_Neg(int oct, int t)
                {
                    return new Vector(0, 0, -1);
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <param name="c">0 for x, 1 for y, 2 for z</param>
                /// <returns></returns>
                public virtual double Directions_Pos(int oct, int t, int c)
                {
                    return 0;
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <param name="c">0 for x, 1 for y, 2 for z</param>
                /// <returns></returns>
                public virtual double Directions_Neg(int oct, int t, int c)
                {
                    return 0;
                }

                public virtual bool HasFilter()
                {
                    return F != null;
                }
            }

            /// <summary>
            /// Class used ot store the energy histogram of the receiver.
            /// </summary>
            [Serializable]
            protected internal class Histogram_1Pt:Histogram
            {
                public Histogram_1Pt(int SampleRate_in, int SampleCT)
                :base(1, 1)
                {
                }

                public Histogram_1Pt(int SampleRate_in, double COTime_in)
                :base(1, 1)
                {
                }

                public virtual string Receiver_Type()
                {
                    return "Type;Map_Data_NoDir";
                }

                /// <summary>
                /// Adds energy detections to this histogram.
                /// </summary>
                /// <param name="Sample">the time point of arrival</param>
                /// <param name="Energy_in">the energy to be recorded</param>
                /// <param name="direction">the incoming direction of the ray</param>
                /// <param name="Octave">the index of the octave band.</param>
                public override void Add(double time, double Energy_in, double dx, double dy, double dz, double Rho_C, int Octave)
                {
                    lock (Energy[Octave])
                    {
                    Energy[Octave][0] += Energy_in;

                    }
                }

                public override bool HasFilter()
                {
                    return false;
                }
            }

            /// <summary>
            /// This class stores all energy detections in a compact, reduced sampling rate series of histograms.
            /// </summary>
            [Serializable]
            protected internal class Directional_Histogram_ThirdOctave : Directional_Histogram
            {
                public Directional_Histogram_ThirdOctave(int rate, int samples) : base(rate, samples, 24) { }
                public Directional_Histogram_ThirdOctave(int rate, double cutoff) : base(rate, cutoff, 24) { }
            }

            [Serializable]
            protected internal class Directional_Histogram : Histogram
            {
                public double[][][] Dir_Rec_Pos;
                public double[][][] Dir_Rec_Neg;
                public double[][] Fdir;
                public Directional_Histogram(int SampleRate_in, int SampleCT)
                    : this(SampleRate_in, SampleCT, 8)
                { }

                public Directional_Histogram(int SampleRate_in, int SampleCT, int bandCount)
                : base(SampleRate_in, SampleCT, bandCount)
                {
                    SampleRate = SampleRate_in;
                    this.SampleCT = SampleCT;
                    Energy = new double[bandCount][];

                    for (int o = 0; o < Energy.Length; o++)
                    {
                        Energy[o] = new double[SampleCT];
                    }
                    Dir_Rec_Pos = new double[3][][];
                    Dir_Rec_Neg = new double[3][][];
                    for (int i = 0; i < 3; i++)
                    {
                        Dir_Rec_Pos[i] = new double[Energy.Length][];
                        Dir_Rec_Neg[i] = new double[Energy.Length][];
                        for (int oct = 0; oct < Energy.Length; oct++)
                        {
                            Dir_Rec_Pos[i][oct] = new double[SampleCT];
                            Dir_Rec_Neg[i][oct] = new double[SampleCT];
                        }
                    }
                }

                /// <summary>
                /// gets the histogram of aggregated vector directions.
                /// </summary>
                /// <param name="SampleRate_in"></param>
                /// <param name="C_Sound_in"></param>
                /// <param name="COTime_in"></param>
                public Directional_Histogram(int SampleRate_in, double COTime_in)
                    : this(SampleRate_in, COTime_in, 8)
                { }

                public Directional_Histogram(int SampleRate_in, double COTime_in, int bandCount)
                : base(SampleRate_in, COTime_in, bandCount)
                {
                    SampleRate = SampleRate_in;
                    CO_Time = COTime_in / 1000;
                    SampleCT = (int)(CO_Time * SampleRate);
                    Energy = new double[bandCount][];
                    for (int o = 0; o < Energy.Length; o++)
                    {
                        Energy[o] = new double[(int)(SampleRate * CO_Time)];
                    }
                    Dir_Rec_Pos = new double[3][][];
                    Dir_Rec_Neg = new double[3][][];
                    for (int i = 0; i < 3; i++)
                    {
                        Dir_Rec_Pos[i] = new double[Energy.Length][];
                        Dir_Rec_Neg[i] = new double[Energy.Length][];
                        for (int oct = 0; oct < Energy.Length; oct++)
                        {
                            Dir_Rec_Pos[i][oct] = new double[SampleCT];
                            Dir_Rec_Neg[i][oct] = new double[SampleCT];
                        }
                    }
                }

                public override string Type()
                {
                    return "Type;Map_Data";
                }

                public override void Combine_Sample(int Sample, double Energy_in, double Pressure_in, Vector direction_pos, Vector direction_neg, int Octave)// float p_real, float p_imag,
                {
                    lock (Energy[Octave])
                    {
                    if (Sample >= Energy[Octave].Length) return;
                    if (Sample < 0) Sample = 0;
                    Energy[Octave][Sample] += Energy_in;
                    Pressure[Octave][Sample] += Pressure_in;
                    Dir_Rec_Pos[0][Octave][Sample] += (direction_pos.dx);
                    Dir_Rec_Pos[1][Octave][Sample] += (direction_pos.dy);
                    Dir_Rec_Pos[2][Octave][Sample] += (direction_pos.dz);
                    Dir_Rec_Neg[0][Octave][Sample] += (direction_neg.dx);
                    Dir_Rec_Neg[1][Octave][Sample] += (direction_neg.dy);
                    Dir_Rec_Neg[2][Octave][Sample] += (direction_neg.dz);

                    }
                }

                public override void Scale(double ray_ct)
                {
                    for (int oct = 0; oct < Energy.Length; oct++) for (int i = 0; i < Energy[oct].Length; i++)
                        {
                            Energy[oct][i] /= ray_ct;
                            Pressure[oct][i] /= Math.Sqrt(ray_ct);
                        }
                    for (int dir = 0; dir < 3; dir++) for (int oct = 0; oct < Dir_Rec_Pos[dir].Length; oct++) for (int i = 0; i < Dir_Rec_Pos[dir][oct].Length; i++)
                            {
                                Dir_Rec_Pos[dir][oct][i] /= ray_ct;
                                Dir_Rec_Neg[dir][oct][i] /= ray_ct;
                            }
                }

                /// <summary>
                /// Returns the energy response of the simulation.
                /// </summary>
                /// <param name="Octave">the octave band index</param>
                /// <returns></returns>
                public override double[] GetEnergyHistogram(int Octave)
                {
                    if (Octave < Energy.Length)
                    {
                        return Energy[Octave];
                    }
                    else
                    {
                        double[] Sum_Energy = new double[(Energy[0].Length)];
                        for (int Oct = 0; Oct < Energy.Length; Oct++)
                        {
                            for (int idx = 0; idx < Energy[Oct].Length; idx++)
                            {
                                Sum_Energy[idx] += Energy[Oct][idx];
                            }
                        }
                        return Sum_Energy;
                    }
                }

                public override void Add(double time, double Energy_in, double dx, double dy, double dz, double Rho_C, int Octave)
                {
                    lock (Energy[Octave])
                    {
                    int sample = (int)(time * SampleRate);
                    if (sample < 0 || sample >= Energy[Octave].Length) return;
                    Energy[Octave][sample] += Energy_in;
                    Pressure[Octave][sample] += Math.Sqrt(Energy_in * Rho_C);
                    if (dx > 0) Dir_Rec_Pos[0][Octave][sample] += (float)(dx * Energy_in);
                    else Dir_Rec_Neg[0][Octave][sample] += (float)(dx * Energy_in);
                    if (dy > 0) Dir_Rec_Pos[1][Octave][sample] += (float)(dy * Energy_in);
                    else Dir_Rec_Neg[1][Octave][sample] += (float)(dy * Energy_in);
                    if (dz > 0) Dir_Rec_Pos[2][Octave][sample] += (float)(dz * Energy_in);
                    else Dir_Rec_Neg[2][Octave][sample] += (float)(dz * Energy_in);

                    }
                }

                /// <summary>
                /// Adds the direction of the arrival.
                /// </summary>
                /// <param name="octave">the octave of the arriving ray.</param>
                /// <param name="sample">the time point of arrival</param>
                /// <param name="V">the direction of the arrival.</param>
                public override void Add_direction(int octave, int sample, Vector VPos, Vector VNeg)
                {
                    this.Dir_Rec_Pos[0][octave][sample] = (float)VPos.dx;
                    this.Dir_Rec_Pos[1][octave][sample] = (float)VPos.dy;
                    this.Dir_Rec_Pos[2][octave][sample] = (float)VPos.dz;
                    this.Dir_Rec_Neg[0][octave][sample] = (float)VNeg.dx;
                    this.Dir_Rec_Neg[1][octave][sample] = (float)VNeg.dy;
                    this.Dir_Rec_Neg[2][octave][sample] = (float)VNeg.dz;
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <returns></returns>
                public override Vector Directions_Pos(int oct, int t)
                {
                    if (oct < Energy.Length) return new Vector(Dir_Rec_Pos[0][oct][t], Dir_Rec_Pos[1][oct][t], Dir_Rec_Pos[2][oct][t]);
                    Vector result = new Vector(0, 0, 0);
                    for (int b = 0; b < Energy.Length; b++) result += Directions_Pos(b, t);
                    return result;
            }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <returns></returns>
                public override Vector Directions_Neg(int oct, int t)
                {
                    if (oct < Energy.Length) return new Vector(Dir_Rec_Neg[0][oct][t], Dir_Rec_Neg[1][oct][t], Dir_Rec_Neg[2][oct][t]);
                    Vector result = new Vector(0, 0, 0);
                    for (int b = 0; b < Energy.Length; b++) result += Directions_Neg(b, t);
                    return result;
            }

                public override void Divide(double divisor)
                {
                    for (int i = 0; i < Energy.Length; i++)
                    {
                        for (int j = 0; j < Energy[i].Length; j++)
                        {
                            this.Energy[i][j] /= divisor;
                            this.Pressure[i][j] /= divisor;
                            for (int k = 0; k < 3; k++) this.Dir_Rec_Pos[k][i][j] /= (float)divisor;
                            for (int k = 0; k < 3; k++) this.Dir_Rec_Neg[k][i][j] /= (float)divisor;
                        }
                    }
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <param name="c">0 for x, 1 for y, 2 for z</param>
                /// <returns></returns>
                public override double Directions_Pos(int oct, int t, int c)
                {
                    return Dir_Rec_Pos[c][oct][t];
                }

                /// <summary>
                /// Creates a filter with a custom leve/EQ. This is useful for auralization using a source with a colored spectrum.
                /// </summary>
                /// <param name="SWL"></param>
                /// <param name="Rho_C"></param>
                /// <returns></returns>
                public override double[][] Create_Filter(double[] SWL, double Rho_C)
                {
                    //DirPos & DirNeg refer to arrival directions of incoming sound. Filter refers to direction sound comes from from reference of the receiver.
                    double[][] F_ = new double[7][];
                    F_[0] = Audio.Pach_SP.ETCToFilter(this.Pressure, SWL, this.SampleRate, 44100);
                    Array.Resize(ref F_[0], F_[0].Length - 4096);

                    for (int dir = 0; dir < 3; dir++)
                    {
                        double[][] temp_ptcP = new double[Energy.Length][];
                        double[][] temp_ptcN = new double[Energy.Length][];
                        for (int oct = 0; oct < Energy.Length; oct++)
                        {
                            temp_ptcP[oct] = new double[this.Pressure[oct].Length];
                            temp_ptcN[oct] = new double[this.Pressure[oct].Length];
                            for (int t = 0; t < temp_ptcP[oct].Length; t++)
                            {
                                temp_ptcP[oct][t] = this.Energy[oct][t] == 0 ? 0 : this.Pressure[oct][t] * this.Dir_Rec_Neg[dir][oct][t] / this.Energy[oct][t];
                                temp_ptcN[oct][t] = this.Energy[oct][t] == 0 ? 0 : this.Pressure[oct][t] * this.Dir_Rec_Pos[dir][oct][t] / this.Energy[oct][t];
                            }
                        }
                        F_[2 * dir + 1] = Audio.Pach_SP.ETCToFilter(temp_ptcP, SWL, this.SampleRate, 44100);
                        Array.Resize(ref F_[2 * dir + 1], F_[2 * dir + 1].Length - 4096);
                        F_[2 * dir + 2] = Audio.Pach_SP.ETCToFilter(temp_ptcN, SWL, this.SampleRate, 44100);
                        Array.Resize(ref F_[2 * dir + 2], F_[2 * dir + 2].Length - 4096);
                    }
                    return F_;
                }

                public override double[] Create_Filter(double[] SWL, double Rho_C, int dim)
                {
                    //DirPos & DirNeg refer to arrival directions of incoming sound. Filter refers to direction sound comes from from reference of the receiver.
                    double[][] temp = new double[Energy.Length][];
                    string msg;
                    if (dim < 1 || dim > 6)
                    {
                        temp = this.Pressure;
                        msg = "Filter Progress: Creating Omnidirectional..";
                    }
                    else
                    {
                        int d = (int)Math.Ceiling((double)dim / 2) - 1;
                        double[][] tref = dim.IsOdd() ? Dir_Rec_Neg[d] : Dir_Rec_Pos[d];
                        msg = string.Format(dim.IsOdd() ? "Filter Progress: Creating Positive {0}." : "Filter Progress: Creating Negative {0}.", new string[3] { "X", "Y", "Z" }[d]);

                        for (int oct = 0; oct < Energy.Length; oct++)
                        {
                            temp[oct] = new double[Pressure[oct].Length];
                            for (int t = 0; t < temp[oct].Length; t++)
                            {
                                temp[oct][t] = this.Energy[oct][t] == 0 ? 0 : this.Pressure[oct][t] * tref[oct][t] / this.Energy[oct][t];
                            }
                        }
                    }
                    double[] F = Audio.Pach_SP.ETCToFilter(temp, SWL, this.SampleRate, 44100);
                    Array.Resize(ref F, F.Length - 4096);
                    return F;
                }

                public override void Create_Filter(double Rho_C)
                {
                    //DirPos & DirNeg refer to arrival directions of incoming sound. Filter refers to direction sound comes from from reference of the receiver.
                    F = Audio.Pach_SP.ETCToFilter(this.Pressure, Enumerable.Repeat(120.0, Energy.Length).ToArray(), this.SampleRate, 44100);
                    Array.Resize(ref F, F.Length - 4096);
                    Fdir = new double[6][];

                    for (int dir = 0; dir < 3; dir++)
                    {
                        double[][] temp_ptcP = new double[Energy.Length][];
                        double[][] temp_ptcN = new double[Energy.Length][];
                        for (int oct = 0; oct < Energy.Length; oct++)
                        {
                            temp_ptcP[oct] = new double[this.Pressure[oct].Length];
                            temp_ptcN[oct] = new double[this.Pressure[oct].Length];
                            for (int t = 0; t < temp_ptcP[oct].Length; t++)
                            {
                                temp_ptcP[oct][t] = this.Energy[oct][t] == 0 ? 0 : this.Pressure[oct][t] * this.Dir_Rec_Neg[dir][oct][t] / this.Energy[oct][t];
                                temp_ptcN[oct][t] = this.Energy[oct][t] == 0 ? 0 : this.Pressure[oct][t] * this.Dir_Rec_Pos[dir][oct][t] / this.Energy[oct][t];
                            }
                        }

                        Fdir[2 * dir] = Audio.Pach_SP.ETCToFilter(temp_ptcP, Enumerable.Repeat(120.0, Energy.Length).ToArray(), this.SampleRate, 44100);
                        Array.Resize(ref Fdir[2*dir], Fdir[2 * dir].Length - 4096);
                        Fdir[2 * dir + 1] = Audio.Pach_SP.ETCToFilter(temp_ptcN, Enumerable.Repeat(120.0, Energy.Length).ToArray(), this.SampleRate, 44100);
                        Array.Resize(ref Fdir[2 * dir + 1], Fdir[2 * dir + 1].Length - 4096);
                    }
                }

                /// <summary>
                /// the direction of the sound at a given instant in the energy response.
                /// </summary>
                /// <param name="oct">octave band index</param>
                /// <param name="t">the sample in time of the energy response to reference</param>
                /// <param name="c">0 for x, 1 for y, 2 for z</param>
                /// <returns></returns>
                public override double Directions_Neg(int oct, int t, int c)
                {
                    return Dir_Rec_Neg[c][oct][t];
                }

                public override double[][] GetFilter3Axis(int rec_id)
                {
                    return Fdir;
                }

                public override void Set_Power(double[] factor)
                {
                    base.Set_Power(factor);
                    for(int i = 0; i < Dir_Rec_Pos.Length; i++) for(int j = 0; j < Dir_Rec_Pos[i].Length; j++) for(int k = 0; k < Dir_Rec_Pos[i][j].Length; k++)
                            {
                                Dir_Rec_Pos[i][j][k] *= factor[j];
                                Dir_Rec_Neg[i][j][k] *= factor[j];
                            }
                }
                public override bool HasFilter()
                {
                    return Fdir != null && Fdir[5] != null && Fdir[5].Length > 0;
                }
            }
        }

        [Serializable]
        public class Expanding_Receiver : Spherical_Receiver
        {
            int RayCount;

            public void SetRayCount(int count) { RayCount = count; }

            public Expanding_Receiver(Point Point, Point SrcPoint, int RCT, double[] Attenuation, double Speed_of_Sound, double rho, int SampleRate, double COTime)
                : base(Point, SrcPoint, Attenuation, Speed_of_Sound, rho, SampleRate, COTime)
            {
                RayCount = RCT;
            }

            public Expanding_Receiver(Hare.Geometry.Point Point, Hare.Geometry.Point SrcPoint, Scene room, int RCT, double Speed_of_Sound, int SampleRate_in, double COTime_in)
                :base(Point, SrcPoint, room, Speed_of_Sound,SampleRate_in,COTime_in)
            {
                RayCount = RCT;
            }

            /// <summary>
            /// Checks receiver for an intersection with a ray
            /// </summary>
            /// <param name="R">the ray.</param>
            /// <param name="EndPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            public override void CheckRay(OctaveRay R, Hare.Geometry.Point EndPt)
            {
                double Radius = (R.t_sum * C_Sound) * Math.Sqrt(6.28 / RayCount);
                double Radius2 = Radius * Radius;
                double mx = R.x - Origin.x;
                double my = R.y - Origin.y;
                double mz = R.z - Origin.z;
                double RayLength = Hare_math.distance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z);
                double dx = (EndPt.x - R.x) / RayLength, dy = (EndPt.y - R.y)/RayLength, dz = (EndPt.z - R.z) / RayLength;
                double b = Hare_math.Dot(mx, my, mz, dx, dy, dz);
                double c = Hare_math.Dot(mx, my, mz, mx, my, mz) - (Radius2);
                if (c > 0 && b > 0) return;
                double discr = b * b - c;
                if (discr < 0) return;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);
                double t = (t1 + t2) * .5;

                if (t > 0 && t * t < SqDistance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z))
                {
                    double Raydist = t * Inv_C_Sound + R.t_sum;
                    if (R.BandCount > 1)
                    {
                        for (int cband = 0; cband < R.BandCount; cband++)
                        {
                            int band = R.BandStart + cband;
                            double energy = (R.Intensity * R.Fractions[cband]) * Math.Pow(10, -.1 * Atten[band] * t) / (Math.PI * Radius2);
                            Recs.Add(Raydist, energy, -R.dx, -R.dy, -R.dz, Rho_C, band);
                            int kind = R.IsSpecular ? (R.WasScattered ? 2 : 0) : 1;
                            int sample = (int)(Raydist * Recs.SampleRate);
                            if (sample >= 0 && sample < Reflection_Analysis_Data[kind][band].Length)
                                lock (Reflection_Analysis_Data[kind][band]) Reflection_Analysis_Data[kind][band][sample] += energy;
                        }
                        return;
                    }
                    double Ei = R.Intensity * Math.Pow(10,-.1 * Atten[R.Octave] * t) / (Math.PI * Radius2);
                    Recs.Add(Raydist, Ei, -R.dx, -R.dy, -R.dz, Rho_C, R.Octave);
                }
                return;
            }

            /// <summary>
            /// This method checks a receiver for a ray using a Broadband Ray (these typically occur before rays are split in the Raytracing simulation.
            /// </summary>
            /// <param name="Length">The length of the ray at the reflection point before potential intersection</param>
            /// <param name="R">the ray.</param>
            /// <param name="EndPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            public override void CheckBroadbandRay(BroadRay R, Hare.Geometry.Point EndPt)
            {
                double Radius = (R.t_sum * C_Sound) * Math.Sqrt(6.28 / RayCount);
                double Radius2 = Radius * Radius;
                double mx = R.x - Origin.x, my = R.y - Origin.y, mz = R.z - Origin.z;
                double RayLength = Hare_math.distance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z);
                double dx = (EndPt.x - R.x) / RayLength;
                double dy = (EndPt.y - R.y) / RayLength;
                double dz = (EndPt.z - R.z) / RayLength;
                double b = Hare_math.Dot(mx, my, mz, dx, dy, dz);
                double c = Hare_math.Dot(mx, my, mz, mx, my, mz) - (Radius2);
                if (c > 0 && b > 0) return;
                double discr = b * b - c;
                if (discr < 0) return;
                double t1 = -b - Math.Sqrt(discr);
                double t2 = -b + Math.Sqrt(discr);
                double t = (t1 + t2) * .5;
                if (t > 0 && t * t < SqDistance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z))
                {
                    double t_ACC = Math.Abs(t1 - t2);
                    double RayTime = t * Inv_C_Sound + R.t_sum;
                    double Area = Math.PI * Radius2;

                    foreach (int group in R.Freq_Bands)
                        for (int component = 0; component < (R.ThirdOctave ? 3 : 1); component++)
                        {
                            int band = R.ThirdOctave ? 3 * group + component : group;
                            double energy = R.Energy[band] * Math.Pow(10, -.1 * Atten[band] * t) / Area;
                            Recs.Add(RayTime, energy, -R.dx, -R.dy, -R.dz, Rho_C, band);
                            if (R.ThirdOctave)
                            {
                                int sample = (int)(RayTime * Recs.SampleRate);
                                if (sample >= 0 && sample < Reflection_Analysis_Data[0][band].Length)
                                    lock (Reflection_Analysis_Data[0][band]) Reflection_Analysis_Data[0][band][sample] += energy;
                            }
                        }
                }
                return;
            }

            /// <summary>
            /// This method checks a receiver for an intersection with a ray.
            /// </summary>
            /// <param name="R">the ray.</param>
            /// <param name="EndPt">The point at which the ray intersects the model after potential receiver intersection.</param>
            /// <param name="Length">The length of the ray at the reflection point before potential intersection</param>
            /// <returns>true if intersects, false if not.</returns>
            public override bool SimpleCheck(BroadRay R, Hare.Geometry.Point EndPt)
            {
                double Radius = R.t_sum * Math.Sqrt(6.28 / RayCount);
                double mx = R.x - Origin.x, my = R.y - Origin.y, mz = R.z - Origin.z;
                double b = Hare_math.Dot(mx, my, mz, R.dx, R.dy, R.dz);
                double c = Hare_math.Dot(mx, my, mz, mx, my, mz) - (Radius * Radius);
                if (c > 0 && b > 0) return false;
                double discr = b * b - c;
                if (discr < 0) return false;
                double t = (-b - Math.Sqrt(discr)) + (-b + Math.Sqrt(discr)) / 2;
                if (t * t < SqDistance(EndPt.x, EndPt.y, EndPt.z, R.x, R.y, R.z)) return true;
                return false;
            }
        }
    }
}
