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

using System.Collections.Generic;
using System;
using Pachyderm_Acoustic.Environment;
using System.Linq;
using Pachyderm_Acoustic.Pach_Graphics;

namespace Pachyderm_Acoustic
{
    namespace Utilities
    {
        public static partial class Numerics
        {
            public static double rt2 = Math.Sqrt(2);
            public static double PiX2 = 2 * Math.PI;
            public static double Pi_180 = Math.PI / 180;
            public static double[] angularFrequency_Octave = new double[] { 62.5 * PiX2, 125 * PiX2, 250 * PiX2, 500 * PiX2, 1000 * PiX2, 2000 * PiX2, 4000 * PiX2, 8000 * PiX2 };
            public static double[] angularFrequency_ThirdOctave = new double[] { 50 * PiX2, 62.5 * PiX2, 80 * PiX2, 100 * PiX2, 125 * PiX2, 160 * PiX2, 200 * PiX2, 250 * PiX2, 315 * PiX2, 400 * PiX2, 500 * PiX2, 630 * PiX2, 800 * PiX2, 1000 * PiX2, 1250 * PiX2, 1600 * PiX2, 2000 * PiX2, 2500 * PiX2, 3150 * PiX2, 4000 * PiX2, 5000 * PiX2, 6300 * PiX2, 8000 * PiX2, 10000 * PiX2 };

            public enum ComplexComponent
            {
                Real,
                Imaginary,
                Magnitude
            }

            public static void ExpComplex(float re, float im, out float re_out, out float im_out)
            {
                float exp = (float)Math.Exp(re);
                re_out = exp * (float)Math.Cos(im);
                im_out = exp * (float)Math.Sin(im);
            }

            public static void ExpComplex(double re, double im, out double re_out, out double im_out)
            {
                float exp = (float)Math.Exp(re);
                re_out = exp * (float)Math.Cos(im);
                im_out = exp * (float)Math.Sin(im);
            }

            public static float Abs(float re, float im)
            {
                return (float)Math.Sqrt((re * re) + (im * im));
            }

            public static double Abs(double re, double im)
            {
                return Math.Sqrt((re * re) + (im * im));
            }

            public static System.Numerics.Complex jBessel(int order, System.Numerics.Complex X)
            {
                //Asymptotic Solution
                //Get Angle
                double Arg = X.Phase;
                int asy_sign = (Arg >= 0) ? -1 : 1;

                return (1 / System.Numerics.Complex.Sqrt(2 * Math.PI * X)) * System.Numerics.Complex.Exp(asy_sign * System.Numerics.Complex.ImaginaryOne * (X - order * Math.PI / 2 - Math.PI / 4));
            }
        }

        ///<summary>
        /// Simple Mathematical tools to estimate certain acoustical parameters...
        ///</summary>
        public class AcousticalMath
        {
            /// <summary>
            /// Calculate Sound Intensity from sound pressure level (dB).
            /// </summary>
            /// <param name="Level"></param>
            /// <returns></returns>
            public static double Intensity_SPL(double Level)
            {
                return Math.Pow(10, Level / 10) * 1E-12;
            }

            /// <summary>
            /// Calculate Sound pressure from sound pressure level (dB).
            /// </summary>
            /// <param name="Level"></param>
            /// <returns></returns>
            public static double Pressure_SPL(double Level)
            {
                return Math.Pow(10, (Level-93) / 20);// * 20E-6;
            }

            /// <summary>
            /// Calculate Sound Pressure Level (in dB) from intensity/energy.
            /// </summary>
            /// <param name="Intensity"></param>
            /// <returns></returns>
            public static double SPL_Intensity(double Intensity)
                {
                    double SPL = 10 * Math.Log10(Intensity / 1E-12);
                    if (SPL < 0) return 0;
                    return SPL;
                }

            public static double SPL_Pressure(double Pressure)
                {
                    return 20 * Math.Log10(Math.Abs(Pressure) / 2E-5);
                }

            public static double Intensity_Pressure(double Pressure, double Rho_C = 405)
            {
                return Pressure * Pressure / Rho_C;
            }

            public static double Pressure_Intensity(double Intensity, double Rho_C = 405)
            {
                return Math.Sqrt(Intensity * Rho_C);
            }

            /// <summary>
            /// Calculate Sound Pressure Level (in dB) from intensity/energy for entire signal.
            /// </summary>
            /// <param name="Intensity"></param>
            /// <returns></returns>
            public static double[] SPL_Intensity_Signal(double[] Intensity)
                {
                    for (int i = 0; i < Intensity.Length; i++)
                    {
                        Intensity[i] = SPL_Intensity(Intensity[i]);
                    }
                    return Intensity;
                }

            public static double[] SPL_Pressure_Signal(double[] P)
                {
                    double[] SPL = new double[P.Length];

                    for (int i = 0; i < P.Length; i++)
                    {
                        SPL[i] = P[i] == 0 ? 0 : SPL_Pressure(P[i]);
                    }
                    return SPL;
                }
                
            /// <summary>
                /// Calculates SPL-time curve from simulation output.
                /// </summary>
                /// <param name="Direct"></param>
                /// <param name="ISData"></param>
                /// <param name="RTData"></param>
                /// <param name="CO_Time"></param>
                /// <param name="samplerate"></param>
                /// <param name="Octave">the chosen octave band.</param>
                /// <param name="Rec_ID">the id of the selected receiver.</param>
                /// <returns></returns>
            public static double[] SPLTCurve_Intensity(Direct_Sound[] Direct, ImageSourceData[] ISData, Environment.Receiver_Bank[] RTData, double CO_Time, int samplerate, int Octave, int Rec_ID, int SrcID, bool Start_at_Zero)
                {
                    double[] Energy = Utilities.IR_Construction.ETCurve(Direct, ISData, RTData, CO_Time, samplerate, Octave, Rec_ID, SrcID, Start_at_Zero);
                    double[] SPL = new double[Energy.Length];

                    for (int i = 0; i < Energy.Length; i++)
                    {
                        SPL[i] = SPL_Intensity(Energy[i]);
                    }
                    return SPL;
                }

            public static double[] SPLTCurve_Intensity(Direct_Sound Direct, ImageSourceData ISData, Environment.Receiver_Bank RTData, double CO_Time, int samplerate, int Octave, int Rec_ID, bool StartAtZero)
                {
                    Direct_Sound[] ArrDirect = new Direct_Sound[1];
                    ArrDirect[0] = Direct;
                    ImageSourceData[] ArrIS = new ImageSourceData[1];
                    ArrIS[0] = ISData;
                    Environment.Receiver_Bank[] ArrRT = new Environment.Receiver_Bank[1];
                    ArrRT[0] = RTData;
                    return SPLTCurve_Intensity(ArrDirect, ArrIS, ArrRT, CO_Time, samplerate, Octave, Rec_ID, 0, StartAtZero);
                }

            //Sound Pressure Level at 1 m. of several kinds of people.
            public static double[][] Females = new double[][] { new double[] { 20, 36, 48, 49, 42, 39, 35, 35 }, new double[] { 20, 37, 51, 53, 49, 44, 42, 37 }, new double[] { 20, 35, 56, 59, 57, 53, 48, 43 }, new double[] { 20, 34, 58, 64, 66, 63, 56, 51 }, new double[] { 20, 30, 56, 69, 76, 75, 69, 58 } };
            public static double[][] Males = new double[][] { new double[] { 20, 45, 49, 50, 42, 41, 38, 35 }, new double[] { 20, 51, 56, 57, 50, 47, 43, 36 }, new double[] { 20, 53, 59, 64, 58, 54, 49, 43 }, new double[] { 20, 56, 64, 72, 70, 66, 60, 50 }, new double[] { 20, 45, 70, 80, 84, 80, 72, 63 } };
            public static double[][] Children = new double[][] { new double[] { 20, 27, 48, 52, 44, 41, 38, 38 }, new double[] { 20, 30, 53, 56, 50, 45, 43, 42 }, new double[] { 20, 31, 56, 60, 60, 55, 51, 46 }, new double[] { 20, 30, 56, 63, 66, 65, 57, 51 }, new double[] { 20, 45, 55, 69, 75, 72, 70, 58 } };

            public static double Sound_Pressure_Level_A (double[] unweighted_SPL)
            {
                double[] AFactors = new double[8] { Math.Pow(10, (-26.2 / 10)), Math.Pow(10, (-16.1 / 10)), Math.Pow(10, (-8.6 / 10)), Math.Pow(10, (-3.2 / 10)), 1, Math.Pow(10, (1.2 / 10)), Math.Pow(10, (1 / 10)), Math.Pow(10, (-1.1 / 10)) };
                double SW = 0;
                for (int f = 0; f < unweighted_SPL.Length; f++)
                {
                    double s = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL(unweighted_SPL[f]);
                    SW += s * AFactors[f];
                }
                return Pachyderm_Acoustic.Utilities.AcousticalMath.SPL_Intensity(SW);
            }

            ///<summary>
            /// Sabine Reverberation Time : T60 = 0.161V/A
            ///</summary>
            public static void Sabine(Environment.Scene Room, double volume, ref double[] T60)
            {
                T60 = new double[8];
                if (Room == null) return;
                double TA;
                for (int t = 0; t <= 7; t++)
                {
                    TA = 0;
                    for (int q = 0; q <= Room.Count() - 1; q++)
                    {
                        TA += (Room.SurfaceArea(q) * Room.AbsorptionValue[q].Coefficient_A_Broad(t));
                    }
                    T60[t] = (0.161 * (volume / (TA + (4 * Room.Attenuation(0)[t] * volume))));
                }
            }

            ///<summary>
            /// Total Absorption A for a given scene object.
            ///</summary>
            public static void Absorption_Total(Environment.Scene Room, out double[] Alpha)
            {
                Alpha = new double[8];
                if (Room == null) return;
                for (int t = 0; t <= 7; t++)
                {
                    for (int q = 0; q <= Room.Count() - 1; q++)
                    {
                        Alpha[t] += (Room.SurfaceArea(q) * Room.AbsorptionValue[q].Coefficient_A_Broad(t));
                    }
                }
            }
            
            ///<summary>
            /// Sabine Reverberation Time : T60 = 0.161V/ln(1-a)S
            ///</summary>
            public static void Eyring(Environment.Scene Room, double volume, ref double[] T60)
            {
                T60 = new double[8];
                if (Room == null) return;
                double TA;
                for (int t = 0; t <= 7; t++)
                {
                    TA = 0;
                    for (int q = 0; q <= Room.Count() - 1; q++)
                    {
                        TA += (Room.SurfaceArea(q) * (-System.Math.Log(1 - Room.AbsorptionValue[q].Coefficient_A_Broad(t),System.Math.E)));
                    }
                    T60[t] = 0.161 * volume / (TA + 4 * Room.Attenuation(0)[t] * volume);
                }
            }

            /// <summary>
            /// The speed of sound in m/s.
            /// </summary>
            /// <param name="Air_Temp"></param>
            /// <returns></returns>
            public static double SoundSpeed(double Air_Temp)
            {
                return 331 + 0.6 * Air_Temp;
            }

            /// <summary>
            /// Calculates the schroeder integral from the energy time curve.
            /// </summary>
            /// <param name="etc">energy time curve.</param>
            /// <returns></returns>
            public static double[] Schroeder_Integral(double[] etc)
            {
                double sum_energy = 0;
                double[] build_up = new double[etc.Length];
                for (int index = 0; index < etc.Length; index++)
                {
                    build_up[index] = sum_energy += etc[index];
                }

                double[] schroed = new double[etc.Length];
                if (sum_energy == 0) return schroed;
                for (int index = 0; index < etc.Length; index++)
                {
                    schroed[index] = 1 - (build_up[index] / sum_energy);
                }

                return schroed;
            }

            public static double[] Schroeder_Integral(double[] etc, double T_Limit)
            {
                double sum_energy = 0;

                if (T_Limit > etc.Length) T_Limit = etc.Length;

                double[] build_up = new double[etc.Length];
                for (int index = 0; index < T_Limit; index++)
                {
                    build_up[index] = sum_energy += etc[index];
                }

                double[] schroed = new double[etc.Length];
                if (sum_energy == 0) return schroed;
                for (int index = 0; index < T_Limit; index++)
                {
                    schroed[index] = 1 - (build_up[index] / sum_energy);
                }

                return schroed;
            }


            /// <summary>
            /// Gets the logarithmic value of the signal.
            /// </summary>
            /// <param name="etc">input data</param>
            /// <param name="limit_in_db">the minimum value of the result.</param>
            /// <returns></returns>
            public static double[] Log10Data(double[] etc, double limit_in_db)
            {
                double[] log_data = new double[etc.Length];
                double limit_val = Math.Pow(10, (limit_in_db / 10));
                for (int index = 1; index <= etc.Length - 1; index++)
                {
                    if (etc[index] >= limit_val)
                    {
                        log_data[index] = 10 * Math.Log10(etc[index]);
                    }
                    else
                    {
                        log_data[index] = limit_in_db;
                    }
                }
                return log_data;
            }

            /// <summary>
            /// Find the start and end of the linear least square fit line, based on the decay interfval chosen.
            /// </summary>
            /// <param name="log_sch">logarithmic schroeder integral</param>
            /// <param name="first">the start of the linear least square fit line</param>
            /// <param name="second">the end of the linear least square fit line</param>
            /// <returns>the array of start and end values.</returns>
            public static int[] Schroeder_Limits(double[] log_sch, double first, double second)
            {
                int[] Limits = new int[2];
                for (int index = 1; index <= log_sch.Length - 1; index++)
                {
                    if (log_sch[index] < first)
                    {
                        Limits[0] = index;
                        break;
                    }
                }

                for (int index = Limits[0]; index <= log_sch.Length - 1; index++)
                {
                    if (log_sch[index] < second)
                    {
                        Limits[1] = index;
                        break;
                    }
                }
                return Limits;
            }

            /// <summary>
            /// Calculates linear least square fit curve.
            /// </summary>
            /// <param name="Magnitude">Section schroeder integral within limits established for lsf operation.</param>
            /// <param name="sample_f">sample frequency.</param>
            /// <returns></returns>
            public static double[] Least_Square_Fit(double[] Magnitude, double sample_f)
            {
                double[] Time = new double[Magnitude.Length];

                int length = Time.Length;

                double Time_Sum = 0;
                double Mag_Sum = 0;
                double TT_Sum = 0;
                double TM_Sum = 0;
                int intIndex = 0;

                for (intIndex = 0; intIndex <= length - 1; intIndex++)
                {
                    double q = intIndex / sample_f;
                    Time_Sum += (q);
                    Mag_Sum += Magnitude[intIndex];

                    TT_Sum += q * q;
                    TM_Sum += q * Magnitude[intIndex];
                }

                double Slope = ((length * TM_Sum) - (Time_Sum * Mag_Sum)) / ((length * TT_Sum) - (Time_Sum * Time_Sum));

                double Intercept = (Mag_Sum - (Slope * Time_Sum)) / length;
                double[] dblResult = { Slope, Intercept };
                return dblResult;
            }

            /// <summary>
            /// normalizes the SPL-Time curve to relative SPL.
            /// </summary>
            /// <param name="SPL_Curve">Spl-Time curve</param>
            /// <returns>relative spl-time curve.</returns>
            public static double[] Normalize_Function(double[] SPL_Curve, List<double> secondary = null)
            {
                double max = 0;
                foreach (double x in SPL_Curve)
                {
                    if (max < x) max = x;
                }

                for (int i = 0; i < SPL_Curve.Length; i++)
                {
                    SPL_Curve[i] -= max;
                }

                if (secondary != null)
                {
                    for (int i = 0; i < secondary.Count; i++)
                    {
                        secondary[i] -= max;
                    }
                }

                return SPL_Curve;
            }

            /// <summary>
            /// Calculate reverberation time. (T-X)
            /// </summary>
            /// <param name="Schroeder_Integral"></param>
            /// <param name="Decay_Depth">For T-10, 10, For T-20, 20, etc...</param>
            /// <param name="sample_f">The number of bins/second</param>
            /// <returns></returns>
            public static double T_X(double[] Schroeder_Integral, int Decay_Depth, int sample_f)
            {
                double[] log_sch = Log10Data(Schroeder_Integral, -60);
                int[] Limits = Schroeder_Limits(log_sch, -5, -5 - Math.Abs(Decay_Depth));

                double[] snippet = new double[Limits[1] - Limits[0] + 1];
                for (int index = Limits[0]; index <= Limits[1]; index++)
                {
                    snippet[index - Limits[0]] = log_sch[index];
                }

                double[] Coefficients = Least_Square_Fit(snippet, sample_f);
                double Value = (-60 / Coefficients[0]);

                return Value;
            }

            /// <summary>
            /// Calculates early decay time
            /// </summary>
            /// <param name="Schroeder_Integral"></param>
            /// <param name="sample_f">sample frequency.</param>
            /// <returns></returns>
            public static double EarlyDecayTime(double[] Schroeder_Integral, int sample_f)
            {
                double[] log_sch = Log10Data(Schroeder_Integral, -60);
                int[] Limits = Schroeder_Limits(log_sch, -1, -10.0001);

                double[] snippet = new double[Limits[1] - Limits[0] + 1];
                for (int index = Limits[0]; index <= Limits[1]; index++)
                {
                    snippet[index - Limits[0]] = log_sch[index];
                }

                double[] Coefficients = Least_Square_Fit(snippet, sample_f);
                double edt = (-60 / Coefficients[0]);

                return edt;
            }

            /// <summary>
            /// Calculate clarity from energy time curve.
            /// </summary>
            /// <param name="etc">energy time curve.</param>
            /// <param name="sample_f"></param>
            /// <param name="seconds"></param>
            /// <param name="startTime">direct arrival time.</param>
            /// <returns></returns>
            public static double Clarity(double[] etc, double sample_f, double seconds, double startTime, bool pressure)
            {
                if (pressure) seconds += (double)1024 / sample_f;

                double Binwidth = 1 / sample_f;

                int StartIndex = (int)Math.Floor((startTime) * sample_f);
                int EndIndex = (int)Math.Floor((seconds + startTime) * sample_f);
                double Sum_Early = 0;
                double Sum_Late = 0;

                for (int q = StartIndex; q < EndIndex; q++)
                {
                    Sum_Early += etc[q];
                }

                double bin_split_pt = (seconds + startTime - (EndIndex * Binwidth)) / Binwidth;
                //Yields percentage of last bin to give to early energy calculation.
                Sum_Early += bin_split_pt * etc[EndIndex];

                for (int q = EndIndex + 1; q < etc.Length; q++)
                {
                    Sum_Late += etc[q];
                }

                Sum_Late += (1 - bin_split_pt) * etc[EndIndex];
                return 10 * Math.Log10(Sum_Early / Sum_Late);
            }

            public static double[] abs_rt = new double[7]{46, 27, 12, 6.5, 7.5, 8, 12};

            public static double[] Modulation_Transfer_Index(double[][] ETC, double rhoC, double[] Noise, int samplefreq)
            {
                for (int i = 0; i < 7; i++) if (ETC[i].Length < samplefreq*1.6) Array.Resize<double>(ref ETC[i], (int)(samplefreq * 1.6));

                double[] MTI = new double[7];
                double[] fm = new double[14] { .63, .8, 1.0, 1.25, 1.6, 2.0, 2.5, 3.15, 4.0, 5.0, 6.3, 8.0, 10.0, 12.5 };
                double I_LowerBand = 0;
                double[] etc = ETC[0];

                for (int s = 0; s < ETC[0].Length; s++)
                {
                    I_LowerBand += etc[s];
                }

                I_LowerBand *= rhoC;

                for (int oct = 1; oct < 8; oct++)
                {
                    etc = ETC[oct];// AcousticalMath.Add_Noise_I(ETC[oct], Noise[oct]);
                    double[] mk = new double[14];
                    double[] ISin = new double[14], ICos = new double[14];
                    double sumI = 0;

                    for (int s = 0; s < etc.Length; s++)
                    {
                        double t = s/(double)samplefreq;
                        double e = etc[s];
                        sumI += e;
                        for (int mct = 0; mct < 14; mct++)
                        {
                            double p = (Utilities.Numerics.PiX2 * fm[mct] * t);
                            ISin[mct] += e * Math.Sin(p);
                            ICos[mct] += e * Math.Cos(p);
                        }
                    }

                    sumI *= rhoC;

                    for (int mct = 0; mct < 14; mct++)
                    {
                        ISin[mct] *= rhoC;
                        ICos[mct] *= rhoC;
                    }

                    double Irtk = Math.Pow(10, abs_rt[oct-1]/10) * 1e-12 * rhoC;
                    double amf;
                    double LowerbandL = AcousticalMath.SPL_Intensity(I_LowerBand);
                    if (LowerbandL < 63)
                    {
                        amf = Math.Pow(10,(0.5 * LowerbandL - 65)/10);
                    }
                    else if(LowerbandL < 67)
                    {
                        amf = Math.Pow(10,(1.8 * LowerbandL - 146.9)/10);
                    }
                    else if(LowerbandL <100)
                    {
                        amf = Math.Pow(10,(0.5 * LowerbandL - 59.8)/10);
                    }
                    else
                    {
                        amf = Math.Pow(10,-10/10);
                    }

                    double I_noise = Math.Pow(10, Noise[oct] / 10) * 1E-12 * rhoC;
                    double Iamk = (I_LowerBand + I_noise) * amf;

                    double Msnr = 1d / (1 + I_noise / sumI);

                    for (int mct = 0; mct < 14; mct++)
                    {
                        mk[mct] = Math.Sqrt(ISin[mct] * ISin[mct] + ICos[mct] * ICos[mct]) / sumI;
                        mk[mct] *= sumI / (sumI + Iamk +Irtk) * Msnr;
                        double TI = ((10 * Math.Log10(mk[mct] / (1 - mk[mct])) + 15.0) / 30.0);
                        MTI[oct - 1] += TI < 0 ? 0 : TI > 1 ? 1 : TI;
                    }

                    I_LowerBand = sumI;
                }

                for (int i = 0; i < MTI.Length; i++) MTI[i] /= 14;
                return MTI;
            }

            public static double[] Speech_Transmission_Index(double[][] ETC, double rhoC, double[] Noise, int samplefreq)
            {
                double[] MTI = Modulation_Transfer_Index(ETC, rhoC, Noise, samplefreq);

                double[] alpha2003 = new double[7] { 0.13, 0.14, 0.11, 0.12, 0.19, 0.17, 0.14 };
                double[] alphaMale = new double[7] { 0.085, 0.127, 0.23, 0.233, 0.309, 0.224, 0.173 };
                double[] alphaFemale = new double[7] { 0, 0.117, 0.223, 0.216, 0.328, 0.25, 0.194 };
                double[] BetaMale = new double[6] { 0.085, 0.078, 0.065, 0.011, 0.047, 0.095 };
                double[] BetaFemale = new double[6] { 0, 0.099, 0.066, 0.062, 0.025, 0.076};

                double[] STI = new double[3];
                for (int oct = 0; oct < 7; oct++) 
                {
                    STI[0] += MTI[oct] * alpha2003[oct];
                    STI[1] += MTI[oct] * alphaMale[oct];
                    STI[2] += MTI[oct] * alphaFemale[oct];
                }
                for (int oct = 0; oct < 6; oct++)
                {
                    double R = Math.Sqrt(MTI[oct] * MTI[oct + 1]);
                    STI[1] -= BetaMale[oct] * R;
                    STI[2] -= BetaFemale[oct] * R;
                }
                //STI:
                //1 = general (2003)
                //2 = Male
                //3 = Female
                return STI;
            }

            public static double Lateral_Parameter(double[][] Dir_ETC, double[] Total_ETC, double LowerBound_s, double UpperBound_s, double sample_f, double startTime, bool pressure)
            {
                if (pressure) startTime += (double)1024 / sample_f;

                double sum_Lateral = 0, sum_Total = 0;
                int i = (int)Math.Floor(startTime * sample_f);
                while (i < sample_f * (LowerBound_s + startTime))
                {
                    sum_Total += Total_ETC[i];
                    i++;
                }
                while (i <= Math.Floor(sample_f * (UpperBound_s + startTime)))
                {
                    sum_Lateral += Dir_ETC[1][i];
                    sum_Total += Total_ETC[i];
                    i++;
                }

                return sum_Lateral / sum_Total;
            }

            public static double Lateral_Parameter(double[] Lateral_ETC, double[] Total_ETC, double LowerBound_s, double UpperBound_s, double sample_f, double startTime, bool pressure)
            {
                double sum_Lateral = 0, sum_Total = 0;
                int i = (int)Math.Floor(startTime * sample_f);
                while(i < sample_f * (LowerBound_s + startTime))
                {
                    sum_Total += Total_ETC[i];
                    i++;
                }
                while(i <= Math.Floor(sample_f * (UpperBound_s + startTime)))
                {
                    sum_Lateral += Math.Abs(Lateral_ETC[i]);
                    sum_Total += Total_ETC[i];
                    i++;
                }

                return sum_Lateral / sum_Total;
            }

            public static double Lateral_Fraction(double[][] ETC, double sample_f, double startTime, bool pressure)
            {
                double[] Total_ETC = new double[ETC[0].Length];
                for (int i = 0; i < ETC[0].Length; i++) Total_ETC[i] = Math.Sqrt(ETC[0][i] * ETC[0][i] + ETC[1][i] * ETC[1][i] + ETC[2][i] * ETC[2][i]);

                return Lateral_Parameter(ETC[1], Total_ETC, 0.005, 0.08, sample_f, startTime, pressure);
            }

            public static double Lateral_Fraction(double[] Total_ETC, double[] Lateral_ETC, double sample_f, double startTime, bool pressure)
            {
                return Lateral_Parameter(Lateral_ETC, Total_ETC, 0.005, 0.08, sample_f, startTime, pressure);
            }

            public static double Lateral_Fraction(double[] Total_ETC, double[][] Lateral_ETC, double sample_f, double startTime, bool pressure)
            {
                return Lateral_Parameter(Lateral_ETC, Total_ETC, 0.005, 0.08, sample_f, startTime, pressure);
            }

            public static double Lateral_Efficiency(double[] Total_ETC, double[] Lateral_ETC, double sample_f, double startTime, bool pressure)
            {
                return Lateral_Parameter(Lateral_ETC, Total_ETC, 0.025, 0.08, sample_f, startTime, pressure);
            }

            /// <summary>
            /// Calculates the Echo Kriterion (Dietsch & Kraak) of the impulse response.
            /// </summary>
            /// <param name="PTC">Pressure time curve</param>
            /// <param name="sample_freq">usually 44100 Hz.</param>
            /// <param name="Direct_Time">Time of arrival of the direct sound.</param>
            /// <param name="Speech">Is it speech (true) or music (false)?</param>
            /// <param name="EKGrad">The actual EKGrad values...</param>
            /// <param name="PercentEcho">The linearly interpolated percentage of people that would perceive an echo at that location.</param>
            /// <param name="Echo10">Do 10% of people perceive an echo?</param>
            /// <param name="Echo50">Do 50% of people perceive an echo?</param>
            public static void EchoCriterion(double[] PTC, int sample_freq, double Direct_Time, bool Speech, out double[] EKGrad, out double[] PercentEcho, out bool Echo10, out bool Echo50)
            {
                EKGrad = new double[PTC.Length];
                double[] EK = new double[2];
                double dte;
                double n;

                if (Speech)
                {
                    //Speech
                    EK = new double[]{.9, 1};
                    dte = 0.009;
                    n = 2f/3f;
                }
                else
                {           
                    //Music
                    EK = new double[2]{1.5, 1.8};
                    dte = 0.014;
                    n = 1;
                }

                double[] time = new double[PTC.Length];
                double[] num = new double[PTC.Length];
                double[] denom = new double[PTC.Length];
                double[] ts = new double[PTC.Length];
                PercentEcho = new double[PTC.Length];

                num[0] = 0;
                denom[0] = 0;
                ts[0] = 0;

                for (int t = (int)(Direct_Time * sample_freq); t < PTC.Length; t++)
                {
                    time[t] = time[t-1] + 1 / (double)sample_freq;
                    double Pn = Math.Pow(Math.Abs(PTC[t]), n);
                    denom[t] = denom[t-1] + Pn;
                    num[t] = num[t-1] + Pn * time[t];                    
                    ts[t] = num[t] / denom[t];
                }
                
                double dEK = (EK[1] - EK[0]) / 40;

                for (int t = (int)(dte * sample_freq); t < PTC.Length; t++)
                {
                    EKGrad[t] = (ts[t] - ts[t - (int)(dte*sample_freq)] )/ dte;
                    //Linear interpolation of Echo percentage... (for prudence, not to exceed 50).
                    PercentEcho[t] = (EKGrad[t] - EK[0]) * dEK + 10;
                }

                double max = PercentEcho.Max();

                Echo10 = false; Echo50 = false;

                if (max > 10) Echo10 = true;
                if (max > 50) Echo50 = true;
            }

            /// <summary>
            /// Generic logarithmic energy ratio calculation
            /// </summary>
            /// <param name="etc">energy time curve.</param>
            /// <param name="sample_f"></param>
            /// <param name="num_start">numerator early bound.</param>
            /// <param name="num_end">numberator late bound.</param>
            /// <param name="denom_start">denominator early bound.</param>
            /// <param name="denom_end">denominator late bound.</param>
            /// <param name="startTime">direct arrival time.</param>
            /// <returns></returns>
            public static double Energy_Ratio(double[] etc, double sample_f, double num_start, double num_end, double denom_start, double denom_end, double startTime)
            {
                double Binwidth = 1 / sample_f;

                int EarlyStartindex;
                double ESbinsplt;
                if (num_start == 0)
                {
                    EarlyStartindex = 0;
                    ESbinsplt = 1;
                }
                else
                {
                    EarlyStartindex = (int)Math.Floor((num_start + startTime) * sample_f);
                    ESbinsplt = 1 - (num_start + startTime - (EarlyStartindex * Binwidth)) * sample_f;
                }

                int LateStartindex;
                double LSbinsplt;
                if (denom_start == 0)
                {
                    LateStartindex = 0;
                    LSbinsplt = 1;
                }
                else
                {
                    LateStartindex = (int)Math.Floor((denom_start + startTime) * sample_f);
                    LSbinsplt = 1 - (denom_start + startTime - (LateStartindex * Binwidth)) * sample_f;
                }

                int EarlyEndindex;
                double EEbinsplt;
                if (num_end == 0)
                {
                    EarlyEndindex = 0;
                    EEbinsplt = 1;
                }
                else
                {
                    EarlyEndindex = (int)Math.Floor((num_end + startTime) * sample_f);
                    EEbinsplt = (num_end + startTime - (EarlyEndindex * Binwidth)) * sample_f;
                }

                int LateEndindex;
                double LEbinsplt;
                if (denom_end == 0)
                {
                    LateEndindex = 0;
                    LEbinsplt = 0;
                }
                else
                {
                    LateEndindex = (int)Math.Floor((denom_end + startTime) * sample_f);
                    LEbinsplt = (denom_end + startTime - (LateEndindex * Binwidth)) * sample_f;
                }

                double Sum_Early = 0;
                double Sum_Late = 0;

                for (int q = EarlyStartindex + 1; q < EarlyEndindex; q++)
                {
                    Sum_Early += etc[q];
                }

                //Yields percentage of first and last bins to give to numerator calculation.
                Sum_Early += etc[EarlyStartindex] * ESbinsplt + etc[EarlyEndindex] * EEbinsplt;

                for (int q = LateStartindex + 1; q < LateEndindex; q++)
                {
                    Sum_Late += etc[q];
                }

                //Yields percentage of first and last bins to give to numerator calculation.
                Sum_Late += etc[LateStartindex] * LSbinsplt + etc[LateEndindex] * LEbinsplt;

                return 10 * Math.Log10(Sum_Early / Sum_Late);
            }

            /// <summary>
            /// Calculate definition from energy time curve.
            /// </summary>
            /// <param name="etc">energy time curve.</param>
            /// <param name="sample_f"></param>
            /// <param name="seconds"></param>
            /// <param name="StartTime"></param>
            /// <returns></returns>
            public static double Definition(double[] etc, double sample_f, double seconds, double StartTime, bool pressure)
            {
                if (pressure) StartTime += (double)1024 / sample_f;

                double Binwidth = 1 / sample_f;

                int StartIndex = 0;
                StartIndex = (int)Math.Floor(StartTime * sample_f);

                int EndIndex = (int)Math.Floor(seconds * sample_f) + StartIndex;
                double Sum_Early = 0;
                double Sum_All = 0;

                for (int q = 0; q < EndIndex; q++)
                {
                    Sum_Early += etc[q];
                }
                double bin_split_pt = (seconds + StartTime - (EndIndex * Binwidth)) / Binwidth;

                //Yields percentage of last bin to give to early energy calculation.
                Sum_Early += bin_split_pt * etc[EndIndex];

                for (int q = 0; q < etc.Length; q++)
                {
                    Sum_All += etc[q];
                }

                return Sum_Early / Sum_All * 100;
            }

            /// <summary>
            /// Calculate center time - center of gravity of impulse response.
            /// </summary>
            /// <param name="etc">energy time curve.</param>
            /// <param name="sample_f"></param>
            /// <param name="Direct_time"></param>
            /// <returns></returns>
            public static double Center_Time(double[] etc, int sample_f, double Direct_time)
            {
                double sumPT=0, sumT=0;
                double BW = 1 / (double)sample_f;

                int start = (int)Math.Floor(Direct_time * sample_f);

                for(int i = start; i < etc.Length; i++)
                {
                    sumPT += etc[i] * (i*BW - Direct_time);
                    sumT += etc[i]; 
                }
                return sumPT / sumT;
            }

            /// <summary>
            /// Calculate sound strength.
            /// </summary>
            /// <param name="etc"></param>
            /// <returns></returns>
            public static double Strength(double[] etc, double SWL, bool pressure)
            {
                double Sum_All = 0;
                for (int q = 0; q < etc.Length; q++)
                {
                    Sum_All += etc[q];
                }

                if (pressure) return AcousticalMath.SPL_Pressure(Sum_All) - SWL + 31;
                else return AcousticalMath.SPL_Intensity(Sum_All) - SWL + 31;
            }

            /// <summary>
            /// Calculate Initial Time Delay Gap (ITDG). **Caution: Unbenchmarked**
            /// </summary>
            /// <param name="etc"></param>
            /// <returns>ITDG in milliseconds</returns>
            public static double InitialTimeDelayGap(double[] etc, int Sample_Frequency)
            {
                int t1 = 0, t2 = 0;
                double maxdiff = 0, nextdiff = 0;
                for (int t = 0; t < 100; t++)
                {
                    double diff = etc[t + 1] - etc[t];

                    if (diff > maxdiff)
                    {
                        maxdiff = diff;
                        t1 = t;
                        continue;
                    }

                    if (diff > nextdiff)
                    {
                        nextdiff = diff;
                        t2 = t;
                    }
                }

                return 1000 * (t2 - t1) / Sample_Frequency;
            }
        }

        public class IR_Construction
        {
            public static double[] ETCurve(Direct_Sound[] Direct, ImageSourceData[] ISData, Environment.Receiver_Bank[] RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, List<int> SrcIDs, bool StartAtZero)
            {
                //TODO: Make provisions for specifying source delays...
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;

                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData[0] is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }

                maxdelay *= (double)Sampling_Frequency / 1000.0;

                double[] Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[] IR = ETCurve(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Octave, Rec_ID, s, StartAtZero);
                    for (int i = 0; i < IR.Length; i++)
                    {
                        Histogram[i] += IR[i];
                    }
                }
                return Histogram;
            }

            public static double[] ETCurve(Direct_Sound Direct, ImageSourceData ISData, Environment.Receiver_Bank RTData, double CO_Time_ms, int samplerate, int Octave, int Rec_ID, bool StartAtZero)
            {
                Direct_Sound[] ArrDirect = new Direct_Sound[1];
                ArrDirect[0] = Direct;
                ImageSourceData[] ArrIS = new ImageSourceData[1];
                ArrIS[0] = ISData;
                Receiver_Bank[] ArrRT = new Environment.Receiver_Bank[1];
                ArrRT[0] = RTData;
                return ETCurve(ArrDirect, ArrIS, ArrRT, CO_Time_ms, samplerate, Octave, Rec_ID, 0, StartAtZero);
            }

            /// <summary>
            /// Calculates Energy-time curve from simulation output.
            /// </summary>
            /// <param name="Direct"></param>
            /// <param name="ISData"></param>
            /// <param name="RTData"></param>
            /// <param name="CO_Time"></param>
            /// <param name="Sampling_Frequency"></param>
            /// <param name="Octave">the chosen octave band.</param>
            /// <param name="Rec_ID">the id of the selected receiver.</param>
            /// <returns></returns>
            public static double[] ETCurve(Direct_Sound[] Direct, ImageSourceData[] ISData, Environment.Receiver_Bank[] RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, int Src_ID, bool Start_at_Zero)
            {
                double[] Histogram = null;

                double delay = 0;
                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    delay = Direct[Src_ID].Delay_ms;
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData[0] is PachMapReceiver)
                {
                    delay = (RTData[Src_ID] as PachMapReceiver).delay_ms;
                }

                if (RTData[Src_ID] != null)
                {
                    Histogram = RTData[Src_ID].GetEnergyHistogram(Octave, delay, Rec_ID);
                }
                else
                {
                    Histogram = new double[(int)((CO_Time_ms + (int)delay) * 0.001 * Sampling_Frequency)];
                }

                if (Direct[Src_ID] != null && Direct[Src_ID].IsOccluded(Rec_ID))
                {
                    int D_Start = 0;
                    if (!Start_at_Zero) D_Start = (int)Math.Ceiling(((double)delay * 0.001 + Direct[Src_ID].Time(Rec_ID)) * Sampling_Frequency);
                    for (int i = 0; i < Direct[Src_ID].Io[Rec_ID][0].Length; i++)
                    {
                        double DirectValue = 0;
                        switch (Octave)
                        {
                            case 8:
                                DirectValue = Direct[Src_ID].EnergySum(Rec_ID, i);
                                break;
                            default:
                                DirectValue = Direct[Src_ID].EnergyValue(Octave, Rec_ID)[i];
                                break;
                        }
                        Histogram[D_Start + i] += DirectValue;
                    }
                }

                if (ISData[Src_ID] != null)
                {
                    foreach (Deterministic_Reflection value in ISData[Src_ID].Paths[Rec_ID])
                    {
                        int place = (int)Math.Ceiling(Sampling_Frequency * (value.TravelTime + (double)Direct[Src_ID].Delay_ms * 0.001));
                        if (place < Histogram.Length - 1 && place > 0)
                        {
                            double[] e = value.Energy(Octave, Sampling_Frequency);
                            for (int t = 0; t < e.Length; t++) if (place + t < Histogram.Length - 1) Histogram[place + t] += e[t];
                        }
                    }
                }

                return Histogram;
            }

            public static double[][] ETCurve_1d_Tight(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees)
            {
                double[][] Histogram = new double[3][];

                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;

                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }

                maxdelay *= Sampling_Frequency / 1000;

                Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[][] IR = ETCurve_1d_Tight(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Octave, Rec_ID, StartAtZero, alt, azi, degrees);
                    for (int d = 0; d < 3; d++)
                    {
                        for (int i = 0; i < IR[0].Length; i++)
                        {
                            Histogram[d][i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[d][i];
                        }
                    }
                }

                return Histogram;
            }

            public static double[][] ETCurve_1d_Tight(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, bool Start_at_Zero, double alt, double azi, bool degrees)
            {
                double[][] Histogram = new double[3][];

                if (Octave < 8)
                {
                    //Get Power...
                    double power = 0;
                    if (Direct != null)
                    {
                        power = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL(Direct.SWL[Octave]);
                    }
                    else if (RTData != null && RTData is PachMapReceiver)
                    {
                        power = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL((RTData as PachMapReceiver).SWL[Octave]);
                    }

                    Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                    Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                    Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                    if (RTData != null)
                    {
                        for (int i = 0; i < Histogram[0].Length; i++)
                        {
                            Hare.Geometry.Vector Vpos = RTData.Directions_Pos(Octave, i, Rec_ID, alt, azi, degrees);
                            Hare.Geometry.Vector Vneg = RTData.Directions_Neg(Octave, i, Rec_ID, alt, azi, degrees);

                            Hare.Geometry.Vector VTot = new Hare.Geometry.Vector(-Math.Abs(Vpos.dx) + Math.Abs(Vneg.dx), -Math.Abs(Vpos.dy) + Math.Abs(Vneg.dy), -Math.Abs(Vpos.dz) + Math.Abs(Vneg.dz));
                            Histogram[0][i] = Math.Abs(VTot.dx);
                            Histogram[1][i] = Math.Abs(VTot.dy);
                            Histogram[2][i] = Math.Abs(VTot.dz);
                        }
                    }

                    if (Direct != null && Direct.IsOccluded(Rec_ID))
                    {
                        int D_Start = 0;
                        if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                        Hare.Geometry.Vector[] DirectValue;
                        switch (Octave)
                        {
                            case 8:
                                DirectValue = Direct.Dir_Energy_Sum(Rec_ID, alt, azi, degrees);
                                break;
                            default:
                                DirectValue = Direct.Dir_Energy(Octave, Rec_ID, alt, azi, degrees);
                                break;
                        }

                        for (int i = 0; i < DirectValue.Length; i++)
                        {
                            Histogram[0][D_Start + i] += DirectValue[i].dx;// * DirectValue[i].x / E);
                            Histogram[1][D_Start + i] += DirectValue[i].dy;// * DirectValue[i].y / E);
                            Histogram[2][D_Start + i] += DirectValue[i].dz;// * DirectValue[i].z / E);
                        }
                    }

                    if (ISData != null)
                    {
                        switch (Octave)
                        {
                            case 8:
                                foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                                {
                                    if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                                    {
                                        Hare.Geometry.Vector[] E_Sum = value.Dir_EnergySum(alt, azi, degrees);
                                        for (int i = 0; i < E_Sum.Length; i++)
                                        {
                                            Histogram[0][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += E_Sum[i].dx;// E_Sum[i].x;// * E_Sum[i].x / E);
                                            Histogram[1][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += E_Sum[i].dy;// E_Sum[i].y;// * E_Sum[i].y / E);
                                            Histogram[2][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += E_Sum[i].dz;// E_Sum[i].z;// * E_Sum[i].z / E);
                                        }
                                    }
                                }
                                break;
                            default:
                                foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                                {
                                    if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                                    {
                                        Hare.Geometry.Vector[] E_Dir = value.Dir_Energy(Octave, alt, azi, degrees);
                                        for (int i = 0; i < E_Dir.Length; i++)
                                        {
                                            Histogram[0][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += E_Dir[i].dx;// * E_Dir[i].x / E);
                                            Histogram[1][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += E_Dir[i].dy;// * E_Dir[i].y / E);
                                            Histogram[2][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += E_Dir[i].dz;// * E_Dir[i].z / E);
                                        }
                                    }
                                }
                                break;
                        }
                    }

                    for (int i = 0; i < Histogram.Length; i++)
                    {
                        for (int j = 0; j < Histogram[i].Length; j++)
                        {
                            Histogram[i][j] = Math.Abs(Histogram[i][j]);
                        }
                    }
                    for (int d = 0; d < 3; d++) for (int i = 0; i < Histogram[d].Length; i++) Histogram[d][i] *= power;

                }
                else
                {
                    //Take Sum
                    for (int oct = 0; oct < 8; oct++)
                    {
                        double[][] Hist = ETCurve_1d_Tight(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, oct, Rec_ID, Start_at_Zero, alt, azi, degrees);
                        if (Histogram[0] == null) for (int d = 0; d < 3; d++) Histogram[d] = new double[Hist[d].Length];
                        for (int d = 0; d < 3; d++) for (int i = 0; i < Histogram[d].Length; i++) Histogram[d][i] += Hist[d][i];
                    }
                }

                return Histogram;
            }

            public static double[][] ETCurve_1d(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees)
            {
                double[][] Histogram = new double[3][];

                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;

                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[][] IR = ETCurve_1d(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Octave, Rec_ID, StartAtZero, alt, azi, degrees);
                    for (int d = 0; d < 3; d++)
                    {
                        for (int i = 0; i < IR[0].Length; i++)
                        {
                            Histogram[d][i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[d][i];
                        }
                    }
                }
                return Histogram;
            }

            public static double[][] ETCurve_1d(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, bool Start_at_Zero, double alt, double azi, bool degrees)
            {

                double[][] Histogram = new double[3][];

                if (Octave < 8)
                {
                    //Get Power...
                    double power = 0;
                    if (Direct != null)
                    {
                        power = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL(Direct.SWL[Octave]);
                    }
                    else if (RTData != null && RTData is PachMapReceiver)
                    {
                        power = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL((RTData as PachMapReceiver).SWL[Octave]);
                    }
                    Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                    Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                    Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                    if (RTData != null)
                    {
                        for (int i = 0; i < Histogram[0].Length; i++)
                        {
                            Hare.Geometry.Vector Vpos = RTData.Directions_Pos(Octave, i, Rec_ID, alt, azi, degrees);
                            Hare.Geometry.Vector Vneg = RTData.Directions_Neg(Octave, i, Rec_ID, alt, azi, degrees);

                            double E = RTData.Rec_List[Rec_ID].Energy(i, Octave);
                            Hare.Geometry.Vector VTot = new Hare.Geometry.Vector(Math.Abs(Vpos.dx) - Math.Abs(Vneg.dx), Math.Abs(Vpos.dy) - Math.Abs(Vneg.dy), Math.Abs(Vpos.dz) - Math.Abs(Vneg.dz));
                            VTot.Normalize();
                            VTot *= E;

                            Histogram[0][i] = VTot.dx;
                            Histogram[1][i] = VTot.dy;
                            Histogram[2][i] = VTot.dz;
                        }
                    }

                    if (Direct != null && Direct.IsOccluded(Rec_ID))
                    {
                        int D_Start = 0;
                        if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                        Hare.Geometry.Vector[] DirectValue;
                        switch (Octave)
                        {
                            case 8:
                                DirectValue = Direct.Dir_Energy_Sum(Rec_ID, alt, azi, degrees);
                                break;
                            default:
                                DirectValue = Direct.Dir_Energy(Octave, Rec_ID, alt, azi, degrees);
                                break;
                        }

                        for (int i = 0; i < DirectValue.Length; i++)
                        {
                            Histogram[0][D_Start + i] += Math.Abs(DirectValue[i].dx);
                            Histogram[1][D_Start + i] += Math.Abs(DirectValue[i].dy);
                            Histogram[2][D_Start + i] += Math.Abs(DirectValue[i].dz);
                        }
                    }

                    if (ISData != null)
                    {
                        switch (Octave)
                        {
                            case 8:
                                foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                                {
                                    if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                                    {
                                        Hare.Geometry.Vector[] E_Sum = value.Dir_EnergySum(alt, azi, degrees);
                                        for (int i = 0; i < E_Sum.Length; i++)
                                        {
                                            Histogram[0][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += Math.Abs(E_Sum[i].dx);
                                            Histogram[1][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += Math.Abs(E_Sum[i].dy);
                                            Histogram[2][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i / Sampling_Frequency)] += Math.Abs(E_Sum[i].dz);
                                        }
                                    }
                                }
                                break;
                            default:
                                foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                                {
                                    if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                                    {
                                        Hare.Geometry.Vector[] E_Dir = value.Dir_Energy(Octave, alt, azi, degrees, true);
                                        for (int i = 0; i < E_Dir.Length; i++)
                                        {
                                            Histogram[0][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i)] += Math.Abs(E_Dir[i].dx);
                                            Histogram[1][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i)] += Math.Abs(E_Dir[i].dy);
                                            Histogram[2][(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i)] += Math.Abs(E_Dir[i].dz);
                                        }
                                    }
                                }
                                break;
                        }
                    }
                    for (int d = 0; d < 3; d++) for (int i = 0; i < Histogram[d].Length; i++) Histogram[d][i] *= power;
                }
                else
                {
                    //Take Sum
                    for (int oct = 0; oct < 8; oct++)
                    {
                        double[][] Hist = ETCurve_1d(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, oct, Rec_ID, Start_at_Zero, alt, azi, degrees);
                        if (Histogram[0] == null) for (int d = 0; d < 3; d++) Histogram[d] = new double[Hist[d].Length];
                        for (int d = 0; d < 3; d++) for (int i = 0; i < Histogram[d].Length; i++) Histogram[d][i] += Hist[d][i];
                    }
                }

                return Histogram;
            }

            public static double[] ETCurve_Directional(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees)
            {
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;
                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                double[] Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[] IR = ETCurve_Directional(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Octave, Rec_ID, StartAtZero, alt, azi, degrees);
                    for (int i = 0; i < IR.Length; i++)
                    {
                        Histogram[i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[i];
                    }
                }
                return Histogram;
            }

            public static double[] ETCurve_Directional(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, bool Start_at_Zero, double alt, double azi, bool degrees)
            {
                double[] Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                if (Octave < 8)
                {
                    //Get Power...
                    double power = 0;
                    if (Direct != null)
                    {
                        power = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL(Direct.SWL[Octave]);
                    }
                    else if (RTData != null && RTData is PachMapReceiver)
                    {
                        power = Pachyderm_Acoustic.Utilities.AcousticalMath.Intensity_SPL((RTData as PachMapReceiver).SWL[Octave]);
                    }

                    if (RTData != null)
                    {
                        for (int i = 0; i < Histogram.Length; i++)
                        {
                            Hare.Geometry.Vector Vpos = RTData.Directions_Pos(Octave, i, Rec_ID, alt, azi, degrees);
                            Hare.Geometry.Vector Vneg = RTData.Directions_Neg(Octave, i, Rec_ID, alt, azi, degrees);

                            double E = RTData.Rec_List[Rec_ID].Energy(i, Octave);
                            Hare.Geometry.Vector VTot = new Hare.Geometry.Vector(Math.Abs(Vpos.dx) - Math.Abs(Vneg.dx), Math.Abs(Vpos.dy) - Math.Abs(Vneg.dy), Math.Abs(Vpos.dz) - Math.Abs(Vneg.dz));

                            if (Vpos.dx > 0)
                            {
                                Histogram[i] += Math.Abs(Vpos.dx);
                            }
                            if (Vneg.dx > 0)
                            {
                                Histogram[i] += Math.Abs(Vneg.dx);
                            }

                            double L = VTot.Length();
                            if (L > 0) Histogram[i] *= E / L;
                        }
                    }

                    if (Direct != null && Direct.IsOccluded(Rec_ID))
                    {
                        int D_Start = 0;
                        if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                        Hare.Geometry.Vector[] DirectValue;
                        switch (Octave)
                        {
                            case 8:
                                DirectValue = Direct.Dir_Energy_Sum(Rec_ID, alt, azi, degrees);
                                break;
                            default:
                                DirectValue = Direct.Dir_Energy(Octave, Rec_ID, alt, azi, degrees);
                                break;
                        }

                        for (int i = 0; i < DirectValue.Length; i++)
                        {
                            if (DirectValue[i].dx > 0) Histogram[D_Start + i] += Math.Abs(DirectValue[i].dx);
                        }
                    }

                    if (ISData != null)
                    {
                        switch (Octave)
                        {
                            case 8:
                                foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                                {
                                    if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram.Length - 1)
                                    {
                                        Hare.Geometry.Vector[] E_Sum = value.Dir_EnergySum(alt, azi, degrees);
                                        for (int i = 0; i < E_Sum.Length; i++)
                                        {
                                            if (E_Sum[i].dx > 0) Histogram[(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i)] += Math.Abs(E_Sum[i].dx);
                                        }
                                    }
                                }
                                break;
                            default:
                                foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                                {
                                    if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram.Length - 1)
                                    {
                                        Hare.Geometry.Vector[] E_Dir = value.Dir_Energy(Octave, alt, azi, degrees);
                                        for (int i = 0; i < E_Dir.Length; i++)
                                        {
                                            if (E_Dir[i].dx > 0) Histogram[(int)Math.Ceiling(Sampling_Frequency * value.TravelTime + i)] += Math.Abs(E_Dir[i].dx);
                                        }
                                    }
                                }
                                break;
                        }
                    }

                    for (int i = 0; i < Histogram.Length; i++) Histogram[i] *= power;

                }
                else
                {
                    //Take Sum
                    for (int oct = 0; oct < 8; oct++)
                    {
                        double[] Hist = ETCurve_Directional(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, oct, Rec_ID, Start_at_Zero, alt, azi, degrees);
                        if (Histogram == null) Histogram = new double[Hist.Length];
                        for (int i = 0; i < Histogram.Length; i++) Histogram[i] += Hist[i];
                    }
                }

                return Histogram;
            }
            
            /// <summary>
            /// Calculates pressure impulse response from simulation output.
            /// </summary>
            /// <param name="Direct"></param>
            /// <param name="ISData"></param>
            /// <param name="RTData"></param>
            /// <param name="CO_Time"></param>
            /// <param name="Sampling_Frequency"></param>
            /// <param name="Octave">the chosen octave band.</param>
            /// <param name="Rec_ID">the id of the selected receiver.</param>
            /// <returns></returns>
            public static double[] Auralization_Filter(Direct_Sound[] Direct, ImageSourceData[] ISData, Environment.Receiver_Bank[] RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, int Src_ID, int dim, bool Start_at_Zero, bool flat, IProgressFeedback VB)
            {
                double[] SWL = Direct[Src_ID].SWL;
                double[] F;

                if (RTData[Src_ID] != null)
                {
                    if (flat) RTData[Src_ID].GetFilter(Rec_ID, out F);
                    else F = RTData[Src_ID].Create_Filter(Direct[Src_ID].SWL, Rec_ID, 0, VB);
                }
                else
                {
                    F = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency)];
                }

                if (Direct[Src_ID] != null && Direct[Src_ID].IsOccluded(Rec_ID))
                {
                    int D_Start = 0;
                    if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct[Src_ID].Time(Rec_ID) * Sampling_Frequency);

                    double[] Filter = flat ? Direct[Src_ID].Get_Filter(Rec_ID, Sampling_Frequency)[0] : Direct[Src_ID].Create_Filter(Direct[Rec_ID].SWL, Rec_ID, dim, Sampling_Frequency);

                    for (int i = 0; i < Filter.Length; i++)
                    {
                        F[D_Start + i] += Filter[i];
                    }
                }

                if (ISData[Src_ID] != null)
                {
                    foreach (Deterministic_Reflection value in ISData[Src_ID].Paths[Rec_ID])
                    {
                        if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < F.Length - 1)
                        {
                            double[] f_value = flat ? value.Filter: value.Create_Filter(SWL, Sampling_Frequency, 0, 4096, 0);

                            int end = value.Filter.Length < F.Length - (int)Math.Ceiling(Sampling_Frequency * value.TravelTime) ? value.Filter.Length : F.Length - (int)Math.Ceiling(Sampling_Frequency * value.TravelTime);
                            int R_start = (int)Math.Ceiling(Sampling_Frequency * value.TravelTime);
                            for (int t = 0; t < end; t++)
                            {
                                int t_s = R_start + t;
                                if (t_s >= 0) F[t_s] += (float)value.Filter[t];
                            }
                        }
                    }
                }
                return F;
            }

            public static double[] Aurfilter_Directional(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, bool Start_at_Zero, double alt, double azi, bool degrees, bool flat, IProgressFeedback VB)
            {
                double[] Histogram;
                if (RTData != null)
                {
                    int[] ids = new int[3];
                    ids[0] = (azi > 90 && azi < 270) ? 1 : 0;
                    ids[1] = (azi <= 180) ? 2 : 3;
                    ids[2] = (alt < 0) ? 5 : 4;
                    int SIGN = 1;
                    for (int i = 1; i < 2; i++) SIGN *= (ids[i] % 2 == 1) ? -1 : 1;

                    double[][] hist_temp = flat ? RTData.Filter_3Axis(Rec_ID) : RTData.Create_Filter(Direct.SWL, Rec_ID, VB);
                    Histogram = new double[hist_temp[0].Length];
                    for (int i = 0; i < hist_temp[0].Length; i++)
                    {
                        Hare.Geometry.Vector V = PachTools.Rotate_Vector(PachTools.Rotate_Vector(new Hare.Geometry.Vector(hist_temp[ids[0]][i], hist_temp[ids[1]][i], hist_temp[ids[2]][i]), azi, 0, true), 0, alt, true);
                        Histogram[i] = V.dx * SIGN;
                    }
                }
                else
                {
                    Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                }

                if (Direct != null && Direct.IsOccluded(Rec_ID))
                {
                    int D_Start = 0;
                    if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                    double[][] V = Direct.Dir_Filter(Rec_ID, alt, azi, degrees,Sampling_Frequency, false, flat);
                    for (int i = 0; i < V.Length; i++)
                    {
                        Histogram[D_Start + i] += V[i][0];
                    }
                }

                if (ISData != null)
                {
                    foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                    {
                        if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram.Length - 1)
                        {
                            int R_Start = (int)Math.Ceiling(Sampling_Frequency * value.TravelTime);

                            double[] V = value.Dir_Filter(Direct.SWL, alt, azi, degrees, false, Sampling_Frequency, flat);
                            for (int i = 0; i < value.Filter.Length; i++)
                            {
                                Histogram[R_Start + i] += V[i];
                            }
                        }
                    }
                }
                return Histogram;
            }

            public static double[][] AurFilter_Fig8_3Axis(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, bool Start_at_Zero, double xpos_alt, double xpos_azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[3][];
                if (RTData != null)
                {
                    double[][] hist_temp = RTData.Filter_3Axis(Rec_ID);
                    Histogram[0] = new double[hist_temp[0].Length];
                    Histogram[1] = new double[hist_temp[0].Length];
                    Histogram[2] = new double[hist_temp[0].Length];
                    for (int i = 0; i < hist_temp[0].Length; i++)
                    {
                        Hare.Geometry.Vector V = PachTools.Rotate_Vector(PachTools.Rotate_Vector(new Hare.Geometry.Vector(hist_temp[0][i] - hist_temp[1][i], hist_temp[2][i] - hist_temp[3][i], hist_temp[4][i] - hist_temp[5][i]), xpos_azi, 0, true), 0, xpos_alt, true);
                        Histogram[0][i] = V.dx;
                        Histogram[1][i] = V.dy;
                        Histogram[2][i] = V.dz;
                    }
                }
                else
                {
                    Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                }

                if (Direct != null && Direct.IsOccluded(Rec_ID))
                {
                    int D_Start = 0;
                    if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                    double[][] V = Direct.Dir_Filter(Rec_ID, xpos_alt, xpos_azi, degrees, Sampling_Frequency, true, flat);
                    for (int i = 0; i < V.Length; i++)
                    {
                        Histogram[0][D_Start + i] += V[i][0];
                        Histogram[1][D_Start + i] += V[i][1];
                        Histogram[2][D_Start + i] += V[i][2];
                    }
                }

                if (ISData != null)
                {
                    foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                    {
                        if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                        {
                            int R_Start = (int)Math.Ceiling(Sampling_Frequency * value.TravelTime);

                            double[][] V = value.Dir_Filter(Direct.SWL, xpos_alt, xpos_azi, degrees, Sampling_Frequency, flat);

                            //Hare.Geometry.Vector dir = value.Path[0][value.Path[0].Length - 1] - value.Path[0][value.Path[0].Length - 2];
                            //dir.Normalize();
                            //for (int i = 0; i < value.Filter.Length; i++)
                            for (int i = 0; i < V.Length; i++)
                            {
                                Histogram[0][R_Start + i] += V[i][0];
                                Histogram[1][R_Start + i] += V[i][1];
                                Histogram[2][R_Start + i] += V[i][2];
                            }
                        }
                    }
                }
                return Histogram; //XYZ - Furse Malham (FUMA)
            }

            /// <summary>
            /// Approximation of the 5 second order ambisonics channels. (Fig8 3Axis and omni are the first four).
            /// </summary>
            /// <param name="Direct"></param>
            /// <param name="ISData"></param>
            /// <param name="RTData"></param>
            /// <param name="CO_Time"></param>
            /// <param name="Sampling_Frequency"></param>
            /// <param name="Rec_ID"></param>
            /// <param name="Start_at_Zero"></param>
            /// <param name="xpos_alt"></param>
            /// <param name="xpos_azi"></param>
            /// <param name="degrees"></param>
            /// <returns></returns>
            public static double[][] AurFilter_Ambisonics2(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, bool Start_at_Zero, double xpos_alt, double xpos_azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[5][];
                if (RTData != null)
                {
                    double[][] hist_temp = RTData.Filter_3Axis(Rec_ID);
                    Histogram[0] = new double[hist_temp[0].Length];
                    Histogram[1] = new double[hist_temp[0].Length];
                    Histogram[2] = new double[hist_temp[0].Length];
                    Histogram[3] = new double[hist_temp[0].Length];
                    Histogram[4] = new double[hist_temp[0].Length];
                    for (int i = 0; i < hist_temp[0].Length; i++)
                    {
                        Hare.Geometry.Vector Vpos = PachTools.Rotate_Vector(PachTools.Rotate_Vector(new Hare.Geometry.Vector(hist_temp[0][i], hist_temp[2][i], hist_temp[4][i]), xpos_azi, 0, true), 0, xpos_alt, true);
                        Hare.Geometry.Vector Vneg = PachTools.Rotate_Vector(PachTools.Rotate_Vector(new Hare.Geometry.Vector(-hist_temp[1][i], -hist_temp[3][i], -hist_temp[5][i]), xpos_azi, 0, true), 0, xpos_alt, true);
                        double magpos = Math.Sqrt(Vpos.dx * Vpos.dx + Vpos.dy * Vpos.dy + Vpos.dz * Vpos.dz);
                        double magneg = Math.Sqrt(Vneg.dx * Vneg.dx + Vneg.dy * Vneg.dy + Vneg.dz * Vneg.dz);
                        double magxyp = Math.Sqrt(Vpos.dx * Vpos.dx + Vpos.dy * Vpos.dy);
                        double magxyn = Math.Sqrt(Vneg.dx * Vneg.dx + Vneg.dy * Vneg.dy);
                        double phipos = Math.Atan(Vpos.dz / (magxyp == 0? 1:magxyp));
                        double phineg = Math.Atan(Vneg.dz / (magxyn == 0? 1:magxyn));
                        double thetapos = Math.Asin(Vpos.dy / (magxyp == 0 ? 1 : magxyp));
                        double thetaneg = Math.Asin(Vneg.dy / (magxyn == 0 ? 1 : magxyn));
                        double rt3_2 = Math.Sqrt(3) / 2;

                        double sin2phpos = Math.Sin(2 * phipos);
                        double sin2phneg = Math.Sin(2 * phineg);
                        double cossqphpos = Math.Cos(phipos) * Math.Cos(phipos);
                        double cossqphneg = Math.Cos(phineg) * Math.Cos(phineg);

                        Histogram[0][i] = magpos * (3 * (Math.Sin(phipos) * Math.Sin(phipos) - 1) / 2 + magneg * 3 * Math.Sin(phineg) * Math.Sin(phineg) - 1) / 2; //R
                        Histogram[1][i] = rt3_2 * (Math.Cos(thetapos) * sin2phpos * magpos + Math.Cos(thetaneg) * sin2phneg * magneg);  //S
                        Histogram[2][i] = rt3_2 * (Math.Sin(thetapos) * sin2phpos * magpos + Math.Sin(thetaneg) * sin2phneg * magneg);  //T
                        Histogram[3][i] = rt3_2 * (Math.Cos(2 * thetapos) * cossqphpos * magpos + Math.Cos(2 * thetaneg) * cossqphneg * magneg);  //U
                        Histogram[4][i] = rt3_2 * (Math.Sin(2 * thetapos) * cossqphpos * magpos + Math.Sin(2 * thetaneg) * cossqphneg * magneg);  //V
                    }
                }
                else
                {
                    Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[4] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[5] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                }

                if (Direct != null && Direct.IsOccluded(Rec_ID))
                {
                    int D_Start = 0;
                    if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                    double[][] V = Direct.Dir_Filter(Rec_ID, xpos_alt, xpos_azi, degrees, Sampling_Frequency, true, flat);
                    for (int i = 0; i < V.Length; i++)
                    {
                        double mag = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1] + V[i][2] * V[i][2]);
                        double magxy = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1]);
                        double phi = Math.Atan(V[i][2] / (magxy == 0 ? 1 : magxy));// Math.Asin(Vpos.z / (magpos == 0? 1 : magpos));
                        double theta = Math.Asin(V[i][1] / (magxy == 0 ? 1 : magxy));//phipos / Math.Cos(phipos));
                        double rt3_2 = Math.Sqrt(3) / 2;

                        double sin2phi = Math.Sin(2 * phi);
                        double cossqphi = Math.Cos(phi) * Math.Cos(phi);

                        Histogram[0][i + D_Start] += mag * 3 * (Math.Sin(phi) * Math.Sin(phi) - 1) / 2;
                        Histogram[1][i + D_Start] += rt3_2 * Math.Cos(theta) * sin2phi * mag;
                        Histogram[2][i + D_Start] += rt3_2 * Math.Sin(theta) * sin2phi * mag;
                        Histogram[3][i + D_Start] += rt3_2 * Math.Cos(2 * theta) * cossqphi * mag;
                        Histogram[4][i + D_Start] += rt3_2 * Math.Sin(2 * theta) * cossqphi * mag;
                    }
                }

                if (ISData != null)
                {
                    foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                    {
                        if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                        {
                            int R_Start = (int)Math.Ceiling(Sampling_Frequency * value.TravelTime);

                            double[][] V = value.Dir_Filter(Direct.SWL, xpos_alt, xpos_azi, degrees, Sampling_Frequency, flat);

                            Hare.Geometry.Vector dir = value.Path[0][value.Path[0].Length - 1] - value.Path[0][value.Path[0].Length - 2];
                            dir.Normalize();
                            for (int i = 0; i < V.Length; i++)
                            {
                                double mag = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1] + V[i][2] * V[i][2]);
                                double magxy = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1]);
                                double phi = Math.Atan(V[i][2] / (magxy == 0 ? 1 : magxy));
                                double theta = Math.Asin(V[i][1] / (magxy == 0 ? 1 : magxy));
                                double rt3_2 = Math.Sqrt(3) / 2;

                                double sin2phi = Math.Sin(2 * phi);
                                double cossqphi = Math.Cos(phi) * Math.Cos(phi);

                                Histogram[0][i + R_Start] += mag * 3 * (Math.Sin(phi) * Math.Sin(phi) - 1) / 2;
                                Histogram[1][i + R_Start] += rt3_2 * Math.Cos(theta) * sin2phi * mag;
                                Histogram[2][i + R_Start] += rt3_2 * Math.Sin(theta) * sin2phi * mag;
                                Histogram[3][i + R_Start] += rt3_2 * Math.Cos(2 * theta) * cossqphi * mag;
                                Histogram[4][i + R_Start] += rt3_2 * Math.Sin(2 * theta) * cossqphi * mag;
                            }
                        }
                    }
                }
                return Histogram;
            }

            public static double[][] AurFilter_Ambisonics3(Direct_Sound Direct, ImageSourceData ISData, Receiver_Bank RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, bool Start_at_Zero, double xpos_alt, double xpos_azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[7][];
                if (RTData != null)
                {
                    double[][] hist_temp = RTData.Filter_3Axis(Rec_ID);
                    Histogram[0] = new double[hist_temp[0].Length];
                    Histogram[1] = new double[hist_temp[0].Length];
                    Histogram[2] = new double[hist_temp[0].Length];
                    Histogram[3] = new double[hist_temp[0].Length];
                    Histogram[4] = new double[hist_temp[0].Length];
                    Histogram[5] = new double[hist_temp[0].Length];
                    Histogram[6] = new double[hist_temp[0].Length];

                    for (int i = 0; i < hist_temp[0].Length; i++)
                    {
                        Hare.Geometry.Vector Vpos = PachTools.Rotate_Vector(PachTools.Rotate_Vector(new Hare.Geometry.Vector(hist_temp[0][i], hist_temp[2][i], hist_temp[4][i]), xpos_azi, 0, true), 0, xpos_alt, true);
                        Hare.Geometry.Vector Vneg = PachTools.Rotate_Vector(PachTools.Rotate_Vector(new Hare.Geometry.Vector(-hist_temp[1][i], -hist_temp[3][i], -hist_temp[5][i]), xpos_azi, 0, true), 0, xpos_alt, true);
                        double magpos = Math.Sqrt(Vpos.dx * Vpos.dx + Vpos.dy * Vpos.dy + Vpos.dz * Vpos.dz);
                        double magneg = Math.Sqrt(Vneg.dx * Vneg.dx + Vneg.dy * Vneg.dy + Vneg.dz * Vneg.dz);
                        double magxyp = Math.Sqrt(Vpos.dx * Vpos.dx + Vpos.dy * Vpos.dy);
                        double magxyn = Math.Sqrt(Vneg.dx * Vneg.dx + Vneg.dy * Vneg.dy);
                        double phipos = Math.Atan(Vpos.dz / (magxyp == 0 ? 1 : magxyp));
                        double phineg = Math.Atan(Vneg.dz / (magxyn == 0 ? 1 : magxyn));
                        double thetapos = Math.Asin(Vpos.dy / (magxyp == 0 ? 1 : magxyp));
                        double thetaneg = Math.Asin(Vneg.dy / (magxyn == 0 ? 1 : magxyn));
                        double rt3_8 = Math.Sqrt(3.0 / 8.0);
                        double rt15_2 = Math.Sqrt(15.0) / 2.0;
                        double rt5_8 = Math.Sqrt(5.0 / 8.0);

                        double LM_compos = Math.Cos(phipos) * (5 * Math.Pow(Math.Sin(phipos), 2) - 1);
                        double LM_comneg = Math.Cos(phineg) * (5 * Math.Pow(Math.Sin(phineg), 2) - 1);
                        double NO_compos = Math.Sin(phipos) * Math.Pow(Math.Cos(phipos), 2);
                        double NO_comneg = Math.Sin(phineg) * Math.Pow(Math.Cos(phineg), 2);
                        double PQ_compos = Math.Pow(Math.Cos(phipos), 3);
                        double PQ_comneg = Math.Pow(Math.Cos(phipos), 3);

                        Histogram[0][i] = magpos * Math.Sin(phipos) * (5 * Math.Sin(phipos) * Math.Sin(phipos) - 3) / 2 + magneg * Math.Sin(phineg) * (5 * Math.Sin(phineg) * Math.Sin(phineg) - 3) / 2; //K
                        Histogram[1][i] = rt3_8 * (magpos * Math.Cos(thetapos) * LM_compos + magneg * Math.Cos(thetaneg) * LM_comneg); //L
                        Histogram[2][i] = rt3_8 * (magpos * Math.Sin(thetapos) * LM_compos + magneg * Math.Sin(thetaneg) * LM_comneg); //M
                        Histogram[3][i] = rt15_2 * (magpos * Math.Cos(2 * thetapos) * NO_compos + magneg * Math.Cos(2 * thetaneg) * NO_compos); //N
                        Histogram[4][i] = rt15_2 * (magpos * Math.Sin(2 * thetapos) * NO_compos + magneg * Math.Sin(2 * thetaneg) * NO_compos); //O
                        Histogram[5][i] = rt5_8 * (magpos * Math.Cos(3 * thetapos) * PQ_compos + magneg * Math.Cos(3 * thetaneg) * PQ_compos); //P
                        Histogram[6][i] = rt5_8 * (magpos * Math.Sin(3 * thetapos) * PQ_compos + magneg * Math.Sin(3 * thetaneg) * PQ_compos); //Q
                    }
                }
                else
                {
                    Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[3] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[4] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[5] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                    Histogram[6] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096];
                }

                if (Direct != null && Direct.IsOccluded(Rec_ID))
                {
                    int D_Start = 0;
                    if (!Start_at_Zero) D_Start = (int)Math.Ceiling(Direct.Time(Rec_ID) * Sampling_Frequency);

                    double[][] V = Direct.Dir_Filter(Rec_ID, xpos_alt, xpos_azi, degrees, Sampling_Frequency, true, flat);
                    for (int i = 0; i < V.Length; i++)
                    {
                        double mag = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1] + V[i][2] * V[i][2]);
                        double magxy = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1]);
                        double phi = Math.Atan(V[i][2] / (magxy == 0 ? 1 : magxy));
                        double theta = Math.Asin(V[i][1] / (magxy == 0 ? 1 : magxy));
                        double rt3_8 = Math.Sqrt(3.0 / 8.0);
                        double rt15_2 = Math.Sqrt(15.0) / 2.0;
                        double rt5_8 = Math.Sqrt(5.0 / 8.0);

                        double LM_com = Math.Cos(phi) * (5 * Math.Pow(Math.Sin(phi), 2) - 1);
                        double NO_com = Math.Sin(phi) * Math.Pow(Math.Cos(phi), 2);
                        double PQ_com = Math.Pow(Math.Cos(phi), 3);

                        Histogram[0][i + D_Start] += mag * Math.Sin(phi) * (5 * Math.Sin(phi) * Math.Sin(phi) - 3) / 2;
                        Histogram[1][i + D_Start] += rt3_8 * mag * Math.Cos(theta) * LM_com;
                        Histogram[2][i + D_Start] += rt3_8 * mag * Math.Sin(theta) * LM_com;
                        Histogram[3][i + D_Start] += rt15_2 * mag * Math.Cos(2 * theta) * NO_com;
                        Histogram[4][i + D_Start] += rt15_2 * mag * Math.Sin(2 * theta) * NO_com;
                        Histogram[5][i + D_Start] += rt5_8 * mag * Math.Cos(3 * theta) * PQ_com;
                        Histogram[6][i + D_Start] += rt5_8 * mag * Math.Sin(3 * theta) * PQ_com;
                    }
                }

                if (ISData != null)
                {
                    foreach (Deterministic_Reflection value in ISData.Paths[Rec_ID])
                    {
                        if (Math.Ceiling(Sampling_Frequency * value.TravelTime) < Histogram[0].Length - 1)
                        {
                            int R_Start = (int)Math.Ceiling(Sampling_Frequency * value.TravelTime);

                            double[][] V = value.Dir_Filter(Direct.SWL, xpos_alt, xpos_azi, degrees, Sampling_Frequency, flat);

                            Hare.Geometry.Vector dir = value.Path[0][value.Path[0].Length - 1] - value.Path[0][value.Path[0].Length - 2];
                            dir.Normalize();
                            for (int i = 0; i < V.Length; i++)
                            {
                                double mag = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1] + V[i][2] * V[i][2]);
                                double magxy = Math.Sqrt(V[i][0] * V[i][0] + V[i][1] * V[i][1]);
                                double phi = Math.Atan(V[i][2] / (magxy == 0 ? 1 : magxy));// Math.Asin(Vpos.z / (magpos == 0? 1 : magpos));
                                double theta = Math.Asin(V[i][1] / (magxy == 0 ? 1 : magxy));//phipos / Math.Cos(phipos));
                                double rt3_8 = Math.Sqrt(3.0 / 8.0);
                                double rt15_2 = Math.Sqrt(15.0) / 2.0;
                                double rt5_8 = Math.Sqrt(5.0 / 8.0);

                                double LM_com = Math.Cos(phi) * (5 * Math.Pow(Math.Sin(phi), 2) - 1);
                                double NO_com = Math.Sin(phi) * Math.Pow(Math.Cos(phi), 2);
                                double PQ_com = Math.Pow(Math.Cos(phi), 3);

                                Histogram[0][i + R_Start] += mag * Math.Sin(phi) * (5 * Math.Sin(phi) * Math.Sin(phi) - 3) / 2;
                                Histogram[1][i + R_Start] += rt3_8 * mag * Math.Cos(theta) * LM_com;
                                Histogram[2][i + R_Start] += rt3_8 * mag * Math.Sin(theta) * LM_com;
                                Histogram[3][i + R_Start] += rt15_2 * mag * Math.Cos(2 * theta) * NO_com;
                                Histogram[4][i + R_Start] += rt15_2 * mag * Math.Sin(2 * theta) * NO_com;
                                Histogram[5][i + R_Start] += rt5_8 * mag * Math.Cos(3 * theta) * PQ_com;
                                Histogram[6][i + R_Start] += rt5_8 * mag * Math.Sin(3 * theta) * PQ_com;
                            }
                        }
                    }
                }
                return Histogram;
            }
            
            public static double[] Auralization_Filter(Direct_Sound[] Direct, ImageSourceData[] ISData, Environment.Receiver_Bank[] RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, bool flat, IProgressFeedback VB)
            {
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;

                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData[0] is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }

                maxdelay *= Sampling_Frequency / 1000;

                double[] Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 * 2 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[] P = Auralization_Filter(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, s, 0, StartAtZero, flat, VB);
                    for (int i = 0; i < P.Length; i++)
                    {
                        Histogram[i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += P[i];
                    }
                }
                return Histogram;
            }

            public static double[] PressureTimeCurve(Direct_Sound[] Direct, ImageSourceData[] ISData, Environment.Receiver_Bank[] RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, bool flat, IProgressFeedback VB)
            {
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;

                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData[0] is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }

                if (RTData.ElementAt<Receiver_Bank>(0) != null)
                {
                    if (CO_Time_ms == 0) CO_Time_ms = RTData[0].CO_Time;
                }
                else
                {
                    if (CO_Time_ms == 0) CO_Time_ms = 1000;
                }

                maxdelay *= Sampling_Frequency / 1000;

                double[] Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 * 2 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[] P = Auralization_Filter(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, s, 0, StartAtZero, flat, VB);
                    P = Audio.Pach_SP.Filter2Signal(P, Direct[s].SWL, Sampling_Frequency, 0);
                    for (int i = 0; i < P.Length; i++)
                    {
                        Histogram[i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += P[i];
                    }
                }
                return Histogram;
            }

            public static double[][] AurFilter_Ambisonics2(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[5][];

                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;
                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[3] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[4] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[][] IR = AurFilter_Ambisonics2(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Rec_ID, StartAtZero, alt, azi, degrees, flat);
                    for (int d = 0; d < IR.Length; d++)
                    {
                        for (int i = 0; i < IR[0].Length; i++)
                        {
                            Histogram[d][i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[d][i];
                        }
                    }
                }
                return Histogram;
            }

            public static double[][] AurFilter_Ambisonics3(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[7][];

                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;
                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[3] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[4] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[5] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[6] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[][] IR = AurFilter_Ambisonics3(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Rec_ID, StartAtZero, alt, azi, degrees, flat);
                    for (int d = 0; d < IR.Length; d++)
                    {
                        for (int i = 0; i < IR[0].Length; i++)
                        {
                            Histogram[d][i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[d][i];
                        }
                    }
                }

                return Histogram;
            }

            public enum Ambisonics_Component_Order { FuMa, SID, ACN }

            public static double[][] AurFilter_Fig8_3Axis(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat, Ambisonics_Component_Order order)
            {
                double[][] filter = AurFilter_Fig8_3Axis(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, SrcIDs, StartAtZero, alt, azi, degrees, flat);
                if (order == Ambisonics_Component_Order.ACN)
                {
                    double[][] final = new double[filter.Length][];
                    final[0] = filter[2];
                    final[1] = filter[0];
                    final[2] = filter[1];
                    return final;
                }
                else return filter; //FuMa,SID
            }

            public static double[][] AurFilter_Ambisonics2(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat, Ambisonics_Component_Order order)
            {
                double[][] filter = AurFilter_Ambisonics2(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, SrcIDs, StartAtZero, alt, azi, degrees, flat);
                if (order == Ambisonics_Component_Order.ACN)
                {
                    double[][] final = new double[filter.Length][];
                    final[0] = filter[4];//4
                    final[1] = filter[0];//5
                    final[2] = filter[3];//6
                    final[3] = filter[1];//7
                    final[4] = filter[2];//8
                    return final;
                }
                else return filter; //FuMa,SID
            }

            public static double[][] AurFilter_Ambisonics3(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat, Ambisonics_Component_Order order)
            {
                double[][] filter = AurFilter_Ambisonics3(Direct, ISData, RTData, CO_Time_ms, Sampling_Frequency, Rec_ID, SrcIDs, StartAtZero, alt, azi, degrees, flat);
                if (order == Ambisonics_Component_Order.ACN)
                {
                    double[][] final = new double[filter.Length][];
                    final[0] = filter[6];//9
                    final[1] = filter[0];//10
                    final[2] = filter[5];//11
                    final[3] = filter[1];//12
                    final[4] = filter[4];//13
                    final[5] = filter[2];//14
                    final[6] = filter[3];//15
                    return final;
                }
                else return filter; //FuMa,SID
            }

            public static double[][] AurFilter_Fig8_3Axis(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[3][];
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];
                double maxdelay = 0;
                List<double> delays = new List<double>();
                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[][] IR = AurFilter_Fig8_3Axis(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Rec_ID, StartAtZero, alt, azi, degrees, flat);
                    for (int d = 0; d < 3; d++)
                    {
                        for (int i = 0; i < IR[0].Length; i++)
                        {
                            Histogram[d][i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[d][i];
                        }
                    }
                }
                return Histogram;
            }

            public static double[][] PTC_Fig8_3Axis(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat)
            {
                double[][] Histogram = new double[3][];
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];
                double maxdelay = 0;
                List<double> delays = new List<double>();
                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                Histogram[0] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[1] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];
                Histogram[2] = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[][] IR = AurFilter_Fig8_3Axis(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Rec_ID, StartAtZero, alt, azi, degrees, flat);
                    for (int d = 0; d < 3; d++)
                    {
                        IR[d] = Audio.Pach_SP.Filter2Signal(IR[d], Direct.ElementAt<Direct_Sound>(s).SWL, Sampling_Frequency, 0);
                        for (int i = 0; i < IR[0].Length; i++)
                        {
                            Histogram[d][i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[d][i];
                        }
                    }
                }
                return Histogram;
            }


            public static double[] AurFilter_Directional(IEnumerable<Direct_Sound> Direct, IEnumerable<ImageSourceData> ISData, IEnumerable<Environment.Receiver_Bank> RTData, double CO_Time_ms, int Sampling_Frequency, int Octave, int Rec_ID, List<int> SrcIDs, bool StartAtZero, double alt, double azi, bool degrees, bool flat, IProgressFeedback VB)
            {
                if (Direct == null) Direct = new Direct_Sound[SrcIDs[SrcIDs.Count - 1] + 1];
                if (ISData == null) ISData = new ImageSourceData[SrcIDs[SrcIDs.Count - 1] + 1];
                if (RTData == null) RTData = new Environment.Receiver_Bank[SrcIDs[SrcIDs.Count - 1] + 1];

                double maxdelay = 0;
                List<double> delays = new List<double>();

                if (Direct.ElementAt<Direct_Sound>(0) != null)
                {
                    foreach (Direct_Sound d in Direct)
                    {
                        delays.Add(d.Delay_ms);
                        maxdelay = Math.Max(maxdelay, d.Delay_ms);
                    }
                }
                else if (RTData.ElementAt<Receiver_Bank>(0) != null && RTData.ElementAt(0) is PachMapReceiver)
                {
                    foreach (Receiver_Bank r in RTData)
                    {
                        delays.Add((r as PachMapReceiver).delay_ms);
                        maxdelay = Math.Max(maxdelay, (r as PachMapReceiver).delay_ms);
                    }
                }
                maxdelay *= Sampling_Frequency / 1000;

                double[] Histogram = new double[(int)(CO_Time_ms * 0.001 * Sampling_Frequency) + 4096 + (int)Math.Ceiling(maxdelay)];

                foreach (int s in SrcIDs)
                {
                    double[] IR = Aurfilter_Directional(Direct.ElementAt<Direct_Sound>(s), ISData.ElementAt<ImageSourceData>(s), RTData.ElementAt<Receiver_Bank>(s), CO_Time_ms, Sampling_Frequency, Rec_ID, StartAtZero, alt, azi, degrees, flat, VB);
                    for (int i = 0; i < IR.Length; i++)
                    {
                        Histogram[i + (int)Math.Ceiling(delays[s] / 1000 * Sampling_Frequency)] += IR[i];
                    }
                }
                return Histogram;
            }
        }

        public class PachTools
        {
            public sealed class RandomNumberGenerator : Random
            {
                private static Random _global = new Random();
                [ThreadStatic]
                private static Random _localInstance;

                public RandomNumberGenerator()
                : base()
                {

                }

                public static Random Instance
                {
                    get
                    {
                        Random inst = _localInstance;
                        if (inst == null)
                        {
                            int seed;
                            lock (_global) seed = _global.Next();
                            _localInstance = inst = new Random(seed);
                        }
                        return _localInstance;
                    }
                }
            }

            /// <summary>
            /// Gets the version of this library...
            /// </summary>
            /// <returns></returns>
            public static string Version()
            {
                System.Reflection.Assembly me = System.Reflection.Assembly.GetExecutingAssembly();
                return me.GetName().Version.ToString();
            }

            /// <summary>
            /// Searches for the minimum value in an array of numbers.
            /// </summary>
            /// <param name="ArrNumbers"></param>
            /// <returns></returns>
            public static double Minimum(double[] ArrNumbers)
            {
                double number = ArrNumbers[0];
                foreach (double value in ArrNumbers)
                {
                    if (value < number) number = value;
                }
                return number;
            }

            public static void World_Angles(Hare.Geometry.Point Src, Hare.Geometry.Point Rec, bool degrees, out double alt, out double azi)
            {
                Hare.Geometry.Vector V = Src - Rec;
                V.Normalize();
                azi = Math.Atan2(V.dy, V.dx);
                alt = Math.Asin(V.dz);

                while (azi < 0) azi += 2 * Math.PI;

                if (degrees)
                {
                    azi *= (180 / Math.PI);
                    alt *= (180 / Math.PI);
                }
            }

            public static bool check_circumsphere(Hare.Geometry.Point ctr, Hare.Geometry.Point Vertex, double r2)
            {
                Hare.Geometry.Vector vd = Vertex - ctr;
                if (vd.dx * vd.dx + vd.dy * vd.dy + vd.dz * vd.dz < r2) return false;
                return true;
            }

            public static void CircumCircleRadius(Hare.Geometry.Point a, Hare.Geometry.Point b, Hare.Geometry.Point c, out double radius, out Hare.Geometry.Point ctr)
            {
                Hare.Geometry.Point acmid = (a + c) / 2;
                Hare.Geometry.Point abmid = (a + b) / 2;
                Hare.Geometry.Vector ac = c - a;
                Hare.Geometry.Vector ab = a - b;
                Hare.Geometry.Vector abXac = Hare.Geometry.Hare_math.Cross(ab, ac);
                Hare.Geometry.Vector Pab = Hare.Geometry.Hare_math.Cross(abXac, ab);
                Hare.Geometry.Vector Pac = Hare.Geometry.Hare_math.Cross(abXac, ac);

                // this is the vector from a TO the circumsphere center
                radius = (Pab.dx * acmid.y - Pab.dx * abmid.y + Pab.dy * abmid.x - Pab.dy * acmid.x) / (Pac.dx * Pab.dy - Pab.dx * Pac.dy);
                ctr = acmid + Pac * radius;
                radius = (ctr - a).Length();
            }

            public static Hare.Geometry.Vector Rotate_Vector(Hare.Geometry.Vector V, double azi, double alt, bool degrees)
            {
                double yaw, pitch;
                if (degrees)
                {
                    yaw = Math.PI * alt / 180.0;
                    pitch = Math.PI * azi / 180.0;
                }
                else
                {
                    yaw = alt;
                    pitch = azi;
                }

                ///Implicit Sparse Rotation Matrix
                double[] r1 = new double[3] { Math.Cos(pitch) * Math.Cos(yaw), Math.Sin(pitch) * Math.Cos(yaw), Math.Sin(yaw) };
                Hare.Geometry.Vector fwd = new Hare.Geometry.Vector(r1[0], r1[1], r1[2]);
                Hare.Geometry.Vector up = new Hare.Geometry.Vector(0, 0, 1) - Hare.Geometry.Hare_math.Dot(new Hare.Geometry.Vector(0, 0, 1), fwd) * fwd;
                up.Normalize();
                double[] r3 = new double[3] { up.dx, up.dy, up.dz };
                Hare.Geometry.Vector right = Hare.Geometry.Hare_math.Cross(up, fwd);
                double[] r2 = new double[3] { right.dx, right.dy, right.dz };

                return (new Hare.Geometry.Vector(r1[0] * V.dx + r1[1] * V.dy + r1[2] * V.dz, r2[0] * V.dx + r2[1] * V.dy + r2[2] * V.dz, r3[0] * V.dx + r3[1] * V.dy + r3[2] * V.dz));
            }

            public static double Polygon_Closest_Distance(Hare.Geometry.Point p, Hare.Geometry.Point a, Hare.Geometry.Point b, Hare.Geometry.Point c )
            {
                Hare.Geometry.Vector ab = b - a;
                Hare.Geometry.Vector ac = c - a;
                Hare.Geometry.Vector bc = c - b;
                // Compute parametric position s for projection P’ of P on AB,
                // P’ = A + s*AB, s = snom/(snom+sdenom)
                double snom = Hare.Geometry.Hare_math.Dot(p - a, ab), sdenom = Hare.Geometry.Hare_math.Dot(p - b, a - b);
                // Compute parametric position t for projection P’ of P on AC,
                // P’ = A + t*AC, s = tnom/(tnom+tdenom)
                double tnom = Hare.Geometry.Hare_math.Dot(p - a, ac), tdenom = Hare.Geometry.Hare_math.Dot(p - c, a - c);
                if (snom <= 0.0f && tnom <= 0.0f) return (p - a).Length(); // Vertex region early out
                // Compute parametric position u for projection P’ of P on BC,
                // P’ = B + u*BC, u = unom/(unom+udenom)
                double unom = Hare.Geometry.Hare_math.Dot(p - b, bc), udenom = Hare.Geometry.Hare_math.Dot(p - c, b - c);
                if (sdenom <= 0.0f && unom <= 0.0f) return (p - b).Length(); // Vertex region early out
                if (tdenom <= 0.0f && udenom <= 0.0f) return (p - c).Length(); // Vertex region early out
                // P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
                Hare.Geometry.Vector n = Hare.Geometry.Hare_math.Cross(b - a, c - a);
                double vc = Hare.Geometry.Hare_math.Dot(n, Hare.Geometry.Hare_math.Cross(a - p, b - p));
                // If P outside AB and within feature region of AB,
                // return projection of P onto AB
                Hare.Geometry.Point temp;
                if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
                {
                    temp = a + (snom / (snom + sdenom) * ab);
                    return Math.Sqrt(temp.x*temp.x + temp.y * temp.y + temp.z * temp.z);
                }
                // P is outside (or on) BC if the triple scalar product [N PB PC] <= 0
                double va = Hare.Geometry.Hare_math.Dot(n, Hare.Geometry.Hare_math.Cross(b - p, c - p));
                // If P outside BC and within feature region of BC,
                // return projection of P onto BC
                if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
                {
                    temp = b + unom / (unom + udenom) * bc;
                    return Math.Sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
                }
                // P is outside (or on) CA if the triple scalar product [N PC PA] <= 0
                double vb = Hare.Geometry.Hare_math.Dot(n, Hare.Geometry.Hare_math.Cross(c - p, a - p));
                // If P outside CA and within feature region of CA,
                // return projection of P onto CA
                if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
                {
                    temp = a + tnom / (tnom + tdenom) * ac;
                    return Math.Sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
                }
                // P must project inside face region. Compute Q using barycentric coordinates
                double u = va / (va + vb + vc);
                double v = vb / (va + vb + vc);
                double w = 1.0f - u - v; // = vc / (va + vb + vc)
                temp = (u* a + v* b + w* c);
                return Math.Sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
            }

            public static double[] NormalDistribution(int samplect, double sum)
            {
                double[] function = new double[samplect];
                double sigma2 = samplect / 2;
                sigma2 *= sigma2;

                double k = 1 / Math.Sqrt(Utilities.Numerics.PiX2 * sigma2);

                for(int i = 0; i < samplect; i++)
                {
                    double num = i - (double)samplect / 2;
                    function[i] = k * Math.Exp(- num * num / (2 * sigma2));
                }

                //double max = function.Max();
                //for (int i = 0; i < samplect; i++) function[i] *= sum / max;
                return function;
            }

            public static double Polygon_Farthest_Distance(Hare.Geometry.Point p, Hare.Geometry.Point a, Hare.Geometry.Point b, Hare.Geometry.Point c)
            {
                Hare.Geometry.Vector pa = p - a, pb = p - b, pc = p - c;
                double dpa = pa.dx * pa.dx + pa.dy * pa.dy + pa.dz * pa.dz;
                double dpb = pb.dx * pb.dx + pb.dy * pb.dy + pb.dz * pb.dz;
                if (dpa < dpb)
                {
                    double dpc = pc.dx * pc.dx + pc.dy * pc.dy + pc.dz * pc.dz;
                    return dpb < dpc ? Math.Sqrt(dpc): Math.Sqrt(dpb);
                }
                else
                {
                    double dpc = pc.dx * pc.dx + pc.dy * pc.dy + pc.dz * pc.dz;
                    return dpa < dpc ? Math.Sqrt(dpc) : Math.Sqrt(dpa);
                }
            }

            public static void Euler_Pitch_Yaw(Hare.Geometry.Vector V1, Hare.Geometry.Vector V2, out double yaw, out double pitch)
            {
                Hare.Geometry.Vector V = V2 - V1;
                yaw = Math.Atan2(V.dx, V.dz);
                double padj = Math.Sqrt(V.dx * V.dx + V.dz * V.dz);
                pitch = Math.Atan2(padj, V.dy);
            }

            public static void RotateMatrix(double phi_radians, ref MathNet.Numerics.LinearAlgebra.Matrix<double> matrix)
            {
                if (matrix.ColumnCount != 2 || matrix.RowCount != 2) throw new Exception("Matrix Rotation was designed for 2x2 gaussian curvature matrices...");
                double sinphi = Math.Sin(phi_radians), cosphi = Math.Cos(phi_radians);
                double m00 = matrix[0, 0], m01 = matrix[0, 1], m10 = matrix[1, 0], m11 = matrix[1, 1];
                double n00 = m00 * cosphi - m10 * sinphi, n01 = m01 * cosphi - m11 * sinphi, n10 = m00 * sinphi + m10 * cosphi, n11 = m01 * sinphi + m11 * cosphi;
                matrix[0, 0] = n00 * cosphi - n01 * sinphi;
                matrix[0, 1] = n00 * sinphi + n01 * cosphi;
                matrix[1, 0] = n10 * cosphi - n11 * sinphi;
                matrix[1, 1] = n10 * sinphi + n11 * cosphi;
            }

            /// <summary>
            /// obtains the octave band or third octave band index by string input.
            /// </summary>
            /// <param name="choice"></param>
            /// <returns></returns>
            public static int OctaveStr2Int(string choice, bool thirdoct)
            {
                if (thirdoct) return ThirdOctaveStr2Int(choice); else return OctaveStr2Int(choice); 
            }

                /// <summary>
                /// obtains the octave band index by string input.
                /// </summary>
                /// <param name="choice"></param>
                /// <returns></returns>
                public static int OctaveStr2Int(string choice)
                {
                switch (choice)
                {
                    case "Summation: All Octaves":
                        return 8;
                    case "62.5 Hz.":
                    case "63":
                        return 0;
                    case "125 Hz.":
                    case "125":
                        return 1;
                    case "250 Hz.":
                    case "250":
                        return 2;
                    case "500 Hz.":
                    case "500":
                        return 3;
                    case "1 kHz.":
                    case "1k":
                    case "1000":
                        return 4;
                    case "2 kHz.":
                    case "2k":
                    case "2000":
                        return 5;
                    case "4 kHz.":
                    case "4k":
                    case "4000":
                        return 6;
                    case "8 kHz.":
                    case "8k":
                    case "8000":
                        return 7;
                }
                return -1;
            }

            /// <summary>
            /// obtains the third-octave band index by string input.
            /// </summary>
            /// <param name="choice"></param>
            /// <returns></returns>
            public static int ThirdOctaveStr2Int(string choice)
            {
                switch (choice)
                {
                    case "Summation: All Octaves":
                        return 8;
                    case "62.5 Hz.":
                    case "63":
                        return 0;
                    case "80 Hz.":
                    case "80":
                        return 1;
                    case "100 Hz.":
                    case "100":
                        return 2;
                    case "125 Hz.":
                    case "125":
                        return 3;
                    case "160 Hz.":
                    case "160":
                        return 4;
                    case "200 Hz.":
                    case "200":
                        return 5;
                    case "250 Hz.":
                    case "250":
                        return 6;
                    case "315 Hz.":
                    case "315":
                        return 7;
                    case "400 Hz.":
                    case "400":
                        return 8;
                    case "500 Hz.":
                    case "500":
                        return 9;
                    case "630 Hz.":
                    case "630":
                        return 10;
                    case "800 Hz.":
                    case "800":
                        return 11;
                    case "1 kHz.":
                    case "1k":
                    case "1000":
                        return 12;
                    case "1.25 kHz.":
                    case "1.25k":
                    case "1250":
                        return 13;
                    case "1.6 kHz.":
                    case "1.6k":
                    case "1600":
                        return 14;
                    case "2 kHz.":
                    case "2k":
                    case "2000":
                        return 15;
                    case "2.5 kHz.":
                    case "2.5k":
                    case "2500":
                        return 16;
                    case "3.15 kHz.":
                    case "3.15k":
                    case "3150":
                        return 17;
                    case "4 kHz.":
                    case "4k":
                    case "4000":
                        return 18;
                    case "5 kHz.":
                    case "5k":
                    case "5000":
                        return 19;
                    case "6.13 kHz.":
                    case "6.13k":
                    case "6130":
                        return 20;
                    case "8 kHz.":
                    case "8k":
                    case "8000":
                        return 21;
                    case "10 kHz.":
                    case "10k":
                    case "10000":
                        return 22;
                }
                return -1;
            }

            public static Eto.Drawing.Color HsvColor(double h, double S, double V)
            {
                // ######################################################################
                // T. Nathan Mundhenk
                // mundhenk@usc.edu
                // C/C++ Macro HSV to RGB
                double H = h * 180 / Math.PI;
                while (H < 0) { H += 360; };
                while (H >= 360) { H -= 360; };
                double R, G, B;
                if (V <= 0)
                { R = G = B = 0; }
                else if (S <= 0)
                {
                    R = G = B = V;
                }
                else
                {
                    double hf = H / 60.0;
                    int i = (int)Math.Floor(hf);
                    double f = hf - i;
                    double pv = V * (1 - S);
                    double qv = V * (1 - S * f);
                    double tv = V * (1 - S * (1 - f));
                    switch (i)
                    {
                        // Red is the dominant color
                        case 0:
                            R = V;
                            G = tv;
                            B = pv;
                            break;
                        // Green is the dominant color
                        case 1:
                            R = qv;
                            G = V;
                            B = pv;
                            break;
                        case 2:
                            R = pv;
                            G = V;
                            B = tv;
                            break;
                        // Blue is the dominant color
                        case 3:
                            R = pv;
                            G = qv;
                            B = V;
                            break;
                        case 4:
                            R = tv;
                            G = pv;
                            B = V;
                            break;
                        // Red is the dominant color
                        case 5:
                            R = V;
                            G = pv;
                            B = qv;
                            break;
                        // Just in case we overshoot on our math by a little, we put these here. Since it's a switch it won't slow us down at all to put these here.
                        case 6:
                            R = V;
                            G = tv;
                            B = pv;
                            break;
                        case -1:
                            R = V;
                            G = pv;
                            B = qv;
                            break;

                        // The color is not defined, we should throw an error.

                        default:
                            //LFATAL("i Value error in Pixel conversion, Value is %d", i);
                            R = G = B = V; // Just pretend its black/white
                            break;
                    }
                }
                int r = (int)(R * 255.0);
                int g = (int)(G * 255.0);
                int b = (int)(B * 255.0);
                return Eto.Drawing.Color.FromArgb(255, r<0?0:r>255?255:r, g<0?0:g>255?255:g, b<0?0:b>255?255:b);
            }

            /// <summary>
            /// Codes a string for the user key "SourcePower"
            /// </summary>
            /// <param name="SWL">8 number double array of sound power levels.</param>
            public static string EncodeSourcePower(double[] SWL)
            {
                string code = "";
                for (int i = 0; i < SWL.Length; i++)
                {
                    code += SWL[i].ToString() + ";";
                }
                return code;
            }

            public static string EncodeSourcePower(float[] SWL)
            {
                string code = "";
                for (int i = 0; i < SWL.Length; i++)
                {
                    code += SWL[i].ToString() + ";";
                }
                return code;
            }

            public static string EncodeAcoustics(int[] Absorption, int[] Scattering, int[] Transparency)
            {
                string Code = "";
                for (int q = 0; q < Absorption.Length; q++)
                {
                    string Temp = null;
                    if (Absorption[q] > 991)
                    {
                        Temp = "xxx";
                    }
                    else if (Absorption[q] < 10)
                    {
                        Temp = string.Concat("00", Absorption[q]);
                    }
                    else if (Absorption[q] < 100)
                    {
                        Temp = string.Concat("0", Absorption[q]);
                    }
                    else
                    {
                        Temp = Absorption[q].ToString();
                    }
                    Code = string.Concat(Code, Temp);
                }
                for (int q = 0; q < 8; q++)
                {
                    string Temp = null;
                    if (Scattering[q] > 99)
                    {
                        Temp = "xx";
                    }
                    else if (Scattering[q] < 10)
                    {
                        Temp = string.Concat("0", Scattering[q].ToString());
                    }
                    else
                    {
                        Temp = Scattering[q].ToString();
                    }
                    Code = string.Concat(Code, Temp);
                }

                if (Transparency != null && Transparency.Length == 8)

                    for (int q = 0; q < 8; q++)
                    {
                        string Temp = null;
                        if (Transparency[q] > 99)
                        {
                            Temp = "xx";
                        }
                        else if (Transparency[q] < 10)
                        {
                            Temp = string.Concat("0", Transparency[q].ToString());
                        }
                        else
                        {
                            Temp = Transparency[q].ToString();
                        }
                        Code = string.Concat(Code, Temp);
                    }
                return Code;
            }

            public static string EncodeAcoustics(double[] Absorption, double[] Scattering, double[] Transparency)
            {
                string Code = "";
                for (int q = 0; q < 8; q++)
                {
                    int abs = (int)(Absorption[q] * 10);
                    string Temp = null;
                    if (abs > 991)
                    {
                        Temp = "xxx";
                    }
                    else if (abs < 10)
                    {
                        Temp = string.Concat("00", abs);
                    }
                    else if (abs < 100)
                    {
                        Temp = string.Concat("0", abs);
                    }
                    else
                    {
                        Temp = abs.ToString();
                    }
                    Code = string.Concat(Code, Temp);
                }
                for (int q = 0; q < 8; q++)
                {
                    string Temp = null;
                    if (Scattering[q] > 99)
                    {
                        Temp = "xx";
                    }
                    else if (Scattering[q] < 10)
                    {
                        Temp = string.Concat("0", ((int)Scattering[q]).ToString());
                    }
                    else
                    {
                        Temp = ((int)Scattering[q]).ToString();
                    }
                    Code = string.Concat(Code, Temp);
                }

                if (Transparency != null && Transparency.Length == 8)

                    for (int q = 0; q < 8; q++)
                    {
                        string Temp = null;
                        if (Transparency[q] > 99)
                        {
                            Temp = "xx";
                        }
                        else if (Transparency[q] < 10)
                        {
                            Temp = string.Concat("0", ((int)Transparency[q]).ToString());
                        }
                        else
                        {
                            Temp = ((int)Transparency[q]).ToString();
                        }
                        Code = string.Concat(Code, Temp);
                    }
                return Code;
            }


            public static bool DecodeAcoustics(string Code, ref double[] Absorption, ref double[] Scattering, ref double[] Transparency)
            {
                if (Code == null) return false;
                int abslength = Code.Length > 48 ? 24 : 8;  
                Absorption = new double[abslength];
                Scattering = new double[8];
                Transparency = new double[8];
                int mod = 0;

                for (int q = 0; q < abslength; q++)
                {
                    if (Code.Length == 48 || Code.Length == 32)
                    {
                        string Temp = string.Concat(Code[2 * q], Code[2 * q + 1]);
                        if (Temp == "xx")
                        { 
                        }
                        else
                        {
                            Absorption[q] = Double.Parse(Temp) / 100;
                        }
                    }
                    else if (Code.Length == 56 || Code.Length == 40)
                    {
                        string Temp = string.Concat(Code[3 * q], Code[3 * q + 1], Code[3 * q + 2]);
                        if (Temp == "xxx")
                        {
                            Absorption[q] = 1;
                        }
                        else
                        {
                            Absorption[q] = Double.Parse(Temp) / 1000;
                        }
                        mod = 8;
                    }
                }
                for (int q = 0; q < 8; q++)
                {
                    string Temp = string.Concat(Code[2 * q + 16 + mod], Code[2 * q + 17 + mod]);
                    if (Temp == "xx")
                    {
                        Scattering[q] = 1;
                    }
                    else
                    {
                        Scattering[q] = Double.Parse(Temp)/100;
                    }
                }
                if (Code.Length == 48 || Code.Length == 56)
                {
                    for (int q = 0; q < 8; q++)
                    {
                        string Temp = string.Concat(Code[2 * q + 32 + mod], Code[2 * q + 33 + mod]);
                        if (Temp == "xx")
                        {
                            Transparency[q] = 1;
                        }
                        else
                        {
                            Transparency[q] = Double.Parse(Temp) / 100;
                        }
                    }
                }
                return true;
            }

            /// <summary>
            /// Takes a SourcePower user string and extracts sound power levels.
            /// </summary>
            /// <param name="code">Rhino source power user string.</param>
            /// <returns></returns>
            public static double[] DecodeSourcePower(string code)
            {
                if (code == null) return new double[] { 120, 120, 120, 120, 120, 120, 120, 120 };
                string[] SWLCodes = code.Split(";".ToCharArray());
                if (SWLCodes.Length < 8) return new double[] { 120, 120, 120, 120, 120, 120, 120, 120 };
                double[] SWL = new double[SWLCodes.Length-1];
                for (int i = 0; i < SWL.Length; i++)
                {                                         
                    SWL[i] = double.Parse(SWLCodes[i]);
                }
                return SWL;
            }

            public static double[] DecodeTransmissionLoss(string code)
            {
                if (string.IsNullOrEmpty(code)) return new double[] { 0, 0, 0, 0, 0, 0, 0, 0 };
                string[] TLCodes = code.Split(";".ToCharArray());
                if (TLCodes.Length < 8) return new double[] { 0, 0, 0, 0, 0, 0, 0, 0 };
                double[] TL = new double[TLCodes.Length];
                for (int i = 0; i < TL.Length; i++)
                {
                    TL[i] = double.Parse(TLCodes[i]);
                }
                return TL;
            }

            public static string EncodeTransmissionLoss(double[] TL)
            {
                string code = "";
                for (int i = 0; i < TL.Length - 1; i++)
                {
                    code += TL[i].ToString() + ";";
                }
                code += TL[7].ToString();
                return code;
            }
        }

        public class StandardConstructions
        {
            public static double[] WelshTraffic = new double[] { -32, -16, -4, -3, -10, -13, -17, -17 };
            /// FHWA traffic coefficients. [Vehicle Type i][Pavement Type p][Full Throttle][A	B	C	D1	D2	E1	E2	F1	F2	G1	G2	H1	H2	I1	I2	J1	J2]
            public static double[][][][] FHWATraffic = new double[][][][] {
                new double[][][]{new double[][]{new double[]{41.740807,1.148546, 67,-7516.580054,-9.7623,16460.1,11.65932,-14823.9,-1.233347,7009.474786,-4.327918,-1835.189815,2.579086,252.418543,-0.573822,-14.268316,0.045682},
                new double[]{41.740807,1.148546,50.128316,-7516.580054,-9.7623,16460.1,11.65932,-14823.9,-1.23334,7009.474786,-4.327918,-1835.189815,2.579086,252.418543,-0.573822,-14.268316,0.045682}},
                new double[][]{new double[]{41.740807,0.494698,67,-7313.985627,-19.697019,16009.5,34.363901,-14414.4,-22.462943,6814.317463,6.093141,-1783.723974,-0.252834,245.299562,-0.170266,-13.86487,0.022131},
                new double[]{41.740807,0.494698,50.128316,-7313.985627,-19.697019,16009.5,34.363901,-14414.4,-22.462943,6814.317463,6.093141,-1783.723974,-0.252834,245.299562,-0.170266,-13.86487,0.022131}},
                new double[][]{new double[]{41.740807,-1.065026,67,-9549.987851,-146.173482,21064,340.622686,-19060.8,-324.802942,9032.990872,161.886578,-2363.810485,-44.454426,324.077238,6.378783,-18.21167,-0.373971},
                new double[]{41.740807,-1.065026,50.128316,-9549.987851,-146.173482,21064,340.622686,-19060.8,-324.802942,9032.990872,161.886578,-2363.810485,-44.454426,324.077238,6.378783,-18.21167,-0.373971}},
                new double[][]{new double[]{41.740807,3.520004,67,-2027.8376,-70.674562,3728.329033,155.109567,-2768.001364,-138.780925,1030.541403,64.525774,-195.32456,-16.430316,16.418899,2.17435,-0.339616,-0.117021},
                new double[]{41.740807,3.520004,50.128316,-2027.8376,-70.674562,3728.329033,155.109567,-2768.001364,-138.780925,1030.541403,64.525774,-195.32456,-16.430316,16.418899,2.17435,-0.339616,-0.117021}}},
                new double[][][]{new double[][]{new double[]{33.918713,20.591046,74,-8997.974274,96.301703,19015.4,-196.241744,-16587,162.56952,7627.874332,-70.394575,-1950.412341,16.876826,263.093464,-2.132793,-14.645109,0.111404},
                new double[]{33.918713,20.591046,68.002978,-1238.353632,-68.218944,2532.436947,151.781493,-2124.165806,-140.388413,919.784302,68.545463,-215.745405,-18.551234,25.909788,2.634001,-1.244253,-0.153272}},
                new double[][]{new double[]{33.918713,19.903775,74,-8997.974274,96.301703,19015.4,-196.241744,-16587,162.56952,7627.874332,-70.394575,-1950.412341,16.876826,263.093464,-2.132793,-14.645109,0.111404},
                new double[]{33.918713,19.903775,68.002978,-230.440015,-82.783198,172.725033,186.80143,131.655819,-174.718246,-207.664798,86.12481,95.139145,-23.513441,-18.96669,3.366475,1.407549,-0.197472}},
                new double[][]{new double[]{33.918713,19.345214,74,-8997.974274,96.301703,19015.4,-196.241744,-16587,162.56952,7627.874332,-70.394575,-1950.412341,16.876826,263.093464,-2.132793,-14.645109,0.111404},
                new double[]{33.918713,19.345214,68.002978,-234.711357,-103.147894,162.036132,244.033651,133.970948,-237.867685,-196.613672,121.527971,87.517298,-34.222359,-17.12562,5.031804,1.253128,-0.301914}},
                new double[][]{new double[]{33.918713,22.141611,74,-8997.974274,96.301703,19015.4,-196.241744,-16587,162.56952,7627.874332,-70.394575,-1950.412341,16.876826,263.093464,-2.132793,-14.645109,0.111404},
                new double[]{33.918713,22.141611,68.002978,-139.27717,-132.207111,97.357937,296.574807,65.350117,-273.981431,-104.555273,132.85439,47.637332,-35.600554,-9.424641,4.997542,0.689877,-0.287335}}},
                new double[][][]{new double[][]{new double[]{35.87985,21.019665,80,-6864.586846,-94.379848,14368.7,226.701375,-12459.2,-220.015419,5710.525999,110.518825,-1458.340416,-30.365892,196.811136,4.33716,-10.977676,-0.252197},
                new double[]{35.87985,21.019665,74.298135,1468.440649,-235.319117,-3852.393214,537.981518,3886.430673,-502.160068,-1986.858782,244.714955,549.002247,-65.686556,-78.239429,9.217734,4.509121,-0.529106}},
                new double[][]{new double[]{35.87985,20.358498,80,-6864.586846,-94.379848,14368.7,226.701375,-12459.2,-220.015419,5710.525999,110.518825,-1458.340416,-30.365892,196.811136,4.337165,-10.977676,-0.252197},
                new double[]{35.87985,20.358498,74.298135,-290.277032,-196.828915,156.854882,450.144699,151.082001,-420.250062,-168.033708,204.806845,60.772941,-54.968455,-9.681901,7.711617,0.570105,-0.442469}},
                new double[][]{new double[]{35.87985,19.107151,80,-6864.586846,-94.379848,14368.7,226.701375,-12459.2,-220.015419,5710.525999,110.518825,-1458.340416,-30.365892,196.811136,4.337165,-10.977676,-0.252197},
                new double[]{35.87985,19.107151,74.298135,-258.941348,-255.205946,135.514216,587.489921,132.973712,-552.824216,-151.366531,272.102657,57.66924,-73.912732,-9.928293,10.514055,0.649271,-0.612569}},
                new double[][]{new double[]{35.87985,21.822818,80,-6864.586846,-94.379848,14368.7,226.701375,-12459.2,-220.015419,5710.525999,110.518825,-1458.340416,-30.365892,196.811136,4.337165,-10.977676,-0.252197},
                new double[]{35.87985,21.822818,74.298135,87.378338,-224.132311,-497.410428,509.705253,579.584033,-473.326603,-298.5689955,229.5809,78.021585,-61.374037,-10.058424,8.58403,0.498685,-0.49149}}},
                new double[][][]{new double[][]{new double[]{23.47953,38.006238,74,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557},
                new double[]{23.47953,38.006238,68.002978,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557}},
                new double[][]{new double[]{23.47953,37.318967,74,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557},
                new double[]{23.47953,37.318967,68.002978,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557}},
                new double[][]{new double[]{23.47953,36.760406,74,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557},
                new double[]{23.47953,36.760406,68.002978,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557}},
                new double[][]{new double[]{23.47953,39.556803,74,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557},
                new double[]{23.47953,39.556803,68.002978,4621.365424,-123.140566,-11601.5,284.796174,11535.3,-267.623062,-5896.461017,130.822488,1645.797051,-35.139019,-238.929963,4.927783,14.139828,-0.282557}}},
                new double[][][]{new double[][]{new double[]{41.022542,10.013879,67,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005},
                new double[]{41.022542,10.013879,56,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005}},
                new double[][]{new double[]{41.022542,10.013879,67,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005},
                new double[]{41.022542,10.013879,56,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005}},
                new double[][]{new double[]{41.022542,10.013879,67,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005},
                new double[]{41.022542,10.013879,56,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005}},
                new double[][]{new double[]{41.022542,10.013879,67,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005},
                new double[]{41.022542,10.013879,56,7546.65902,-8.870177,-17396,7.899209,16181.8,2.526152,-7828.632535,-5.314462,2085.468458,2.344913,-290.816544,-0.435913,16.614043,0.03005}}}};

            public enum Pavement { Average_DGAC_PCC = 0,
                DGAC_Asphalt = 1,
                PCC_Concrete = 2,
                OGAC_OpenGradedAsphalt = 3}

            public static double[] FHWA_Welsh_SoundPower(double SPLW)
            {
                double[] SWL = new double[8];
                for (int oct = 0; oct < 8; oct++) SWL[oct] = SPLW + WelshTraffic[oct];
                return SWL;
            }

            public static double[] FHWA_TNM10_SoundPower(double speed_kph, int pavement, int automobile, int med_truck, int heavy_truck, int buses, int motorcycles, bool full_throttle)
            {
                ///Described at:
                ///http://www.fhwa.dot.gov/environment/noise/traffic_noise_model/old_versions/tnm_version_10/tech_manual/tnm03.cfm#tnma2

                double s = speed_kph;
                //int i = 0;

                int t = (full_throttle) ? 1 : 0;
                double root2 = Math.Sqrt(2);
                double vtot = 0;
                double[] Es = new double[8] { 0, 0, 0, 0, 0, 0, 0, 0 };

                double[] Veh = new double[5] { automobile, med_truck, heavy_truck, buses, motorcycles };

                for (int v = 0; v < 5; v++)
                {
                    double A = FHWATraffic[v][pavement][t][0];
                    double B = FHWATraffic[v][pavement][t][1];
                    double C = FHWATraffic[v][pavement][t][2];
                    double D1 = FHWATraffic[v][pavement][t][3];
                    double D2 = FHWATraffic[v][pavement][t][4];
                    double E1 = FHWATraffic[v][pavement][t][5];
                    double E2 = FHWATraffic[v][pavement][t][6];
                    double F1 = FHWATraffic[v][pavement][t][7];
                    double F2 = FHWATraffic[v][pavement][t][8];
                    double G1 = FHWATraffic[v][pavement][t][9];
                    double G2 = FHWATraffic[v][pavement][t][10];
                    double H1 = FHWATraffic[v][pavement][t][11];
                    double H2 = FHWATraffic[v][pavement][t][12];
                    double I1 = FHWATraffic[v][pavement][t][13];
                    double I2 = FHWATraffic[v][pavement][t][14];
                    double J1 = FHWATraffic[v][pavement][t][15];
                    double J2 = FHWATraffic[v][pavement][t][16];

                    
                    vtot += automobile + med_truck + heavy_truck + buses + motorcycles;

                    for (int oct = 0; oct < 8; oct++)
                    {
                        double f = 62.5 * Math.Pow(2, oct);
                        double[] freq = new double[3] { f / root2, f, f * root2 };

                        for (int oct3 = 0; oct3 < 3; oct3++)
                        {
                            double Ea = Math.Pow(0.6214 * s, A / 10) * Math.Pow(10, B / 10) + Math.Pow(10, C / 10);
                            double logf = Math.Log10(freq[oct3]);
                            double Ls = 10 * Math.Log10(Ea) + (D1 + 0.6214 * D2 * s) + (E1 + 0.6214 * E2 * s) * logf
                                + (F1 + 0.6214 * F2 * s) * logf * logf + (G1 + 0.6214 * G2 * s) * logf * logf * logf
                                + (H1 + 0.6214 * H2 * s) * logf * logf * logf * logf + (I1 + 0.6214 * I2 * s) * logf * logf * logf * logf * logf
                                + (J1 + 0.6214 * J2 * s) * logf * logf * logf * logf * logf * logf;
                            Es[oct] += 0.0476 * Math.Pow(10, Ls / 10) * Veh[v] / s;
                        }
                    }
                }

                double[] Awt = new double[8] { -26, -16, -9, -3, 0, 1.2, 1, -1 };
                double dmod = 10 * Math.Log10(1 / (Utilities.Numerics.PiX2 * 15));
                double[] SWL = new double[8];

                for (int oct = 0; oct < 8; oct++)
                {
                    SWL[oct] = 10 * Math.Log10(Es[oct]) - Awt[oct] - dmod;//
                }
                return SWL;
            }

            public enum Ancon_runway_use {Takeoff = 0, Landing = 1, Both = 2 }

            public static double[] Ancon_SoundPower(double SWLA, double velocity, double slant_delta, Ancon_runway_use TLChoice)
            {
                double[][] Aircraft_Normalization = new double[3][] {
                                //new double[8]{ -12, -10.5, -12, -15, -20, -27, -40, -44},
                                //new double[8]{-11, -13, -12, -13.5, -18, -21, -25, -35},
                                //new double[8]{-11, -10.5, -12, -13.5, -18, -21, -25, -35}
                            new double[8] { 6.191472203, 7.691472203, 6.191472203, 3.191472203, -1.808527797, -8.808527797,-21.8085278, -25.8085278},
                            new double[8] { 5.6783811710, 3.6783811710, 4.678381171, 3.178381171, -1.321618829, -4.321618829, -8.321618829, -18.32161883},
                            new double[8] { 5.678381171, 6.178381171, 4.678381171, 3.178381171, -1.321618829, -4.321618829, -8.321618829, -18.32161883}
                            };

                double[] SWL = new double[8];

                for (int oct = 0;                                                                                                                        oct < 8; oct++)
                {
                    SWL[oct] = SWLA + Aircraft_Normalization[(int)TLChoice][oct];//
                }
                return SWL;
            }
        }
    }
}