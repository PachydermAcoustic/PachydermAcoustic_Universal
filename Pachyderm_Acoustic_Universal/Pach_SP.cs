﻿//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2023, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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
using System.Numerics;
using FFTWSharp;
using System.Linq;
using MathNet.Numerics.Statistics;
using System.ComponentModel.Design.Serialization;
using Pachyderm_Acoustic.Utilities;
using Pachyderm_Acoustic.Pach_Graphics;
using System.Drawing.Printing;
using System.Threading.Tasks;
using System.Runtime.CompilerServices;
using Eto.Forms;
using System.Threading;

namespace Pachyderm_Acoustic
{
    namespace Audio
    {
        public partial class Pach_SP
        {
            static System.Threading.Semaphore S = new System.Threading.Semaphore(1, 1);
            //For all standard 4096 sample FFTs...
            static fftw_complexarray[] FFT_ArrayIn4096;
            static fftw_complexarray[] FFT_ArrayOut4096;
            static fftw_plan[] FFT_Plan4096;
            static object[] fftlock4096;
            static fftw_complexarray[] IFFT_ArrayIn4096;
            static fftw_complexarray[] IFFT_ArrayOut4096;
            static fftw_plan[] IFFT_Plan4096;
            static object[] ifftlock4096;


            //For all others...
            static fftw_complexarray[] FFT_ArrayIn;
            static fftw_complexarray[] FFT_ArrayOut;
            static fftw_plan[] FFT_Plan;
            static object[] fftlock;
            static fftw_complexarray[] IFFT_ArrayIn;
            static fftw_complexarray[] IFFT_ArrayOut;
            static fftw_plan[] IFFT_Plan;
            static object[] ifftlock;


            //Persistent Filter settings
            public static DSP_System Filter;

            /// <summary>
            /// A function used to describe the relationship between the ends of a raised 
            /// cosine filter and the amount of spill into a neighboring octave band.
            /// </summary>
            static MathNet.Numerics.Interpolation.CubicSpline RCos_Integral;


            public static void Initialize_FFTW()
            {
                int proc = System.Environment.ProcessorCount;
                FFT_ArrayIn4096 = new fftw_complexarray[proc];
                FFT_ArrayOut4096 = new fftw_complexarray[proc];
                FFT_Plan4096 = new fftw_plan[proc];
                fftlock4096 = new object[proc];

                IFFT_ArrayIn4096 = new fftw_complexarray[proc];
                IFFT_ArrayOut4096 = new fftw_complexarray[proc];
                IFFT_Plan4096 = new fftw_plan[proc];
                ifftlock4096 = new object[proc];

                FFT_ArrayIn = new fftw_complexarray[proc];
                FFT_ArrayOut = new fftw_complexarray[proc];
                FFT_Plan = new fftw_plan[proc];
                fftlock = new object[proc];
                IFFT_ArrayIn = new fftw_complexarray[proc];
                IFFT_ArrayOut = new fftw_complexarray[proc];
                IFFT_Plan = new fftw_plan[proc];
                ifftlock = new object[proc];

                for (int i = 0; i < proc; i++)
                {
                    FFT_ArrayIn4096[i] = new fftw_complexarray(4096);
                    FFT_ArrayOut4096[i] = new fftw_complexarray(4096);
                    FFT_Plan4096[i] = fftw_plan.dft_1d(4096, FFT_ArrayIn4096[i], FFT_ArrayOut4096[i], fftw_direction.Forward, fftw_flags.Estimate);
                    fftlock4096[i] = new object();
                    IFFT_ArrayIn4096[i] = new fftw_complexarray(4096);
                    IFFT_ArrayOut4096[i] = new fftw_complexarray(4096);
                    IFFT_Plan4096[i] = fftw_plan.dft_1d(4096, IFFT_ArrayIn4096[i], IFFT_ArrayOut4096[i], fftw_direction.Forward, fftw_flags.Estimate);
                    ifftlock4096[i] = new object();
                    fftlock[i] = new object();
                    ifftlock[i] = new object();
                }

                Initialize_filter_functions();
            }

            public abstract class DSP_System
            {
                public abstract double[] Signal(double[] OctavePressure, int SampleFrequency, int LengthStartToFinish, int Threadid);
                public abstract double[] Transfer_Function(double[] OctavePressure, int SampleFrequency, int LengthStartToFinish, int Threadid);
                public abstract double[] Response(double[] Spectrum, int SampleFrequency, int Threadid);
                public abstract Complex[] Spectrum(double[] Octave_Pessure, int SampleFrequency, int LengthStartToFinish, int Threadid);
            }

            public class Minimum_Phase_System : DSP_System
            {
                public override double[] Signal(double[] OctavePressure, int SampleFrequency, int LengthStartToFinish, int Threadid)
                {
                    return Pach_SP.Minimum_Phase_Signal(OctavePressure, SampleFrequency, LengthStartToFinish, Threadid);
                }

                public override double[] Transfer_Function(double[] OctavePressure, int SampleFrequency, int LengthStartToFinish, int Threadid)
                {
                    double[] m_spec = Pach_SP.Magnitude_Filter(OctavePressure, SampleFrequency, LengthStartToFinish, Threadid);
                    return Pach_SP.Minimum_Phase_Response(m_spec, SampleFrequency, Threadid);
                }

                public override double[] Response(double[] Spectrum, int SampleFrequency, int Threadid)
                {
                    return Pach_SP.Minimum_Phase_Response(Spectrum, SampleFrequency, Threadid);
                }

                public override Complex[] Spectrum(double[] Octave_Pessure, int SampleFrequency, int LengthStartToFinish, int Threadid)
                {
                    return Minimum_Phase_Spectrum(Octave_Pessure, SampleFrequency, LengthStartToFinish, Threadid);
                }
            }

            public class Linear_Phase_System : DSP_System
            {
                public override double[] Signal(double[] OctavePressure, int SampleFrequency, int LengthStartToFinish, int Threadid)
                {
                    return Pach_SP.Linear_Phase_Signal(OctavePressure, SampleFrequency, LengthStartToFinish, Threadid);
                }

                public override double[] Transfer_Function(double[] OctavePressure, int SampleFrequency, int LengthStartToFinish, int Threadid)
                {
                    double[] m_spec = Pach_SP.Magnitude_Filter(OctavePressure, SampleFrequency, LengthStartToFinish, Threadid);
                    double[] response = Pach_SP.Linear_Phase_Response(m_spec, SampleFrequency, Threadid);
                    Raised_Cosine_Window(ref response);
                    return response;
                }

                public override double[] Response(double[] Spectrum, int SampleFrequency, int Threadid)
                {
                    return Pach_SP.Linear_Phase_Response(Spectrum, SampleFrequency, Threadid);
                }

                public override Complex[] Spectrum(double[] Octave_Pessure, int SampleFrequency, int LengthStartToFinish, int Threadid)
                {
                    return Linear_Phase_Spectrum(Octave_Pessure, SampleFrequency, LengthStartToFinish, Threadid);
                }
            }

            public static void Initialize_filter_functions()
            {
                int n = 100;
                double[] A = new double[n];
                double[] ph = new double[n];

                for (int i = 0; i < 100; i++)
                {
                    ph[i] = Math.PI / 2 * i;
                    A[i] = 0.25 * (ph[i] - Math.Cos(ph[i])) / (Math.PI * 2);
                }

                RCos_Integral = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(A, ph);
            }

            public static void Raised_Cosine_Window(ref double[] h)
            {
                double sum_before = 0;
                double sum_after = 0;

                int div = 4;

                int idl = (div / 2 - 1) * h.Length / div;
                int idu = (div / 2 + 1) * h.Length / div;

                for (int i = 0; i < idl; i++)
                {
                    h[i] = 0;
                }

                for (int i = idu; i < h.Length; i++)
                {
                    h[i] = 0;
                }

                for (int i = idl; i < idu; i++)
                {
                    sum_before += h[i] * h[i];
                    double weight = (-Math.Cos(((double)(i - idl) / (2 * h.Length / div)) * 2 * Math.PI) + 1) * .5;
                    h[i] *= weight;
                    sum_after += h[i] * h[i];
                }

                double factor = Math.Sqrt(sum_before / sum_after);
                for (int i = 0; i < h.Length; i++) h[i] *= factor;
            }

            public static void Raised_HCosine_Window(ref double[] h)
            {
                double sum_before = 0;
                double sum_after = 0;

                for (int i = 0; i < h.Length; i++)
                {
                    sum_before += h[i] * h[i];
                    double weight = i < h.Length / 2 ? (Math.Cos((1 - Math.Pow(1 - ((double)i / h.Length / 2), 1)) * Math.PI) + 1) * .5 : 0;
                    h[i] *= weight;
                    sum_after += h[i] * h[i];
                }

                double factor = Math.Sqrt(sum_before / sum_after);
                for (int i = 0; i < h.Length; i++) h[i] *= factor;
            }

            public static double[] MaximumLengthSequence(int length)
            {
                bool[] MLS = new bool[length];
                double[] output = new double[length];
                MLS[0] = true; MLS[1] = true; MLS[2] = false;
                for (int n = 3; n < length; n++)
                {
                    MLS[n] = MLS[n - 3] ^ MLS[n - 1];
                }
                for (int n = 0; n < length; n++)
                {
                    output[n] = MLS[n] ? -1 : 1;
                }
                return output;
            }
            public static double[] MLS_Reverb(double length_s, double[] RT, int sampling_frequency, double[] magnitude)
            {
                double[] MLS = MaximumLengthSequence((int)(length_s * sampling_frequency));
                double[] signal = new double[MLS.Length];
                double[] output = new double[MLS.Length];
                for (int oct = 0; oct < 8; oct++)
                {
                    signal = FIR_Bandpass(MLS, oct, sampling_frequency, 0);
                    for (int i = 0; i < signal.Length; i++)
                    {
                        double slope = Math.Pow(10, ((-60d / RT[oct]) * ((double)i / sampling_frequency)) / 10);
                        output[i] += magnitude[oct] * signal[i] * slope;
                    }
                }
                return output;
            }

            public static double[] FIR_Third_Bandpass(double[] h, int octave_index, int Sample_Freq, int thread)
            {
                int length = h.Length;

                if (length != 4096)
                {
                    Array.Resize(ref h, (int)Math.Pow(2, Math.Ceiling(Math.Log(h.Length, 2)) + 1));
                }

                double ctr = 50 * Math.Pow(2, (double)octave_index / 3.0);
                double freq_l = 1.09 * ctr / Utilities.Numerics.rt2;
                double freq_u = 0.92 * ctr * Utilities.Numerics.rt2;
                //int idl = (int)Math.Round((h.Length * freq_l) / (Sample_Freq));
                //int idu = (int)Math.Round((h.Length * freq_u) / (Sample_Freq));

                //(double[] A, double[] B) = MathNet.Filtering.Butterworth.IirCoefficients.BandPass(freq_l * .8/22050d, freq_l/22050d, freq_u/22050d, freq_u / 0.8/22050d, .1, 100);
                //Complex[] Ac = new Complex[A.Length], Bc = new Complex[B.Length];
                double[] freq = new double[h.Length / 2];
                //for (int i = 0; i < A.Length; i++) Ac[i] = new Complex(A[i],0);
                //for (int i = 0; i < B.Length; i++) Bc[i] = new Complex(B[i], 0);
                for (int i = 0; i < h.Length / 2; i++) freq[i] = i * (double)Sample_Freq / h.Length;
                //Complex[] magspec = Pach_SP.IIR_Design.AB_FreqResponse(Bc.ToList(), Ac.ToList(), freq);
                //Complex[] magspec = Pach_SP.IIR_Design.Butter_FreqResponse(freq, 10, freq_l, freq_u);
                //double[] magspec = Pach_SP.IIR_Design.Butter_FreqResponse(105, freq_l, freq_u, Sample_Freq,h.Length);
                double[] magspec = Pach_SP.IIR_Design.Butter_FreqResponse(105, freq_l, freq_u, Sample_Freq, h.Length / 2);

                /////////////Use Zero Phase Bandpass//////////////////////
                //Array.Resize(ref magspec, h.Length / 2);
                System.Numerics.Complex[] filter = Mirror_Spectrum(magspec);

                //double[] f_r = new double[filter.Length];
                //double[] f_i = new double[filter.Length];
                //double[] f_m = new double[filter.Length];

                //for (int i = 0; i < filter.Length; i++) { f_r[i] = filter[i].Real; f_i[i] = filter[i].Imaginary; f_m[i] = filter[i].Magnitude; }

                /////////////Minimum Phase Bandpass///////////////////////
                //double[] response = Minimum_Phase_Response(magspec, 44100, 0);
                //filter = FFT_General(response, 0);

                //Convolve signal with Bandpass Filter.
                Complex[] freq_h = FFT_General(h, thread);

                for (int i = 0; i < freq_h.Length; i++)
                {
                    freq_h[i] *= filter[i];
                }

                double[] h_oct = IFFT_Real_General(freq_h, thread);
                Scale(ref h_oct);

                Array.Resize(ref h_oct, length);

                return h_oct;
            }

            public static double[] FIR_Bandpass(double[] h, int octave_index, int Sample_Freq, int thread)
            {
                int length = h.Length;

                if (length != 4096)
                {
                    Array.Resize(ref h, (int)Math.Pow(2, Math.Ceiling(Math.Log(h.Length, 2)) + 1));
                }

                double ctr = 62.5 * Math.Pow(2, octave_index);
                double freq_l = 1.25 * ctr / Utilities.Numerics.rt2;
                double freq_u = 0.75 * ctr * Utilities.Numerics.rt2;
                //int idl = (int)Math.Round((h.Length * freq_l) / (Sample_Freq));
                //int idu = (int)Math.Round((h.Length * freq_u) / (Sample_Freq));

                //(double[] A, double[] B) = MathNet.Filtering.Butterworth.IirCoefficients.BandPass(freq_l * .8/22050d, freq_l/22050d, freq_u/22050d, freq_u / 0.8/22050d, .1, 100);
                //Complex[] Ac = new Complex[A.Length], Bc = new Complex[B.Length];
                double[] freq = new double[h.Length / 2];
                //for (int i = 0; i < A.Length; i++) Ac[i] = new Complex(A[i],0);
                //for (int i = 0; i < B.Length; i++) Bc[i] = new Complex(B[i], 0);
                for (int i = 0; i < h.Length / 2; i++) freq[i] = i * (double)Sample_Freq / h.Length;
                //Complex[] magspec = Pach_SP.IIR_Design.AB_FreqResponse(Bc.ToList(), Ac.ToList(), freq);
                //Complex[] magspec = Pach_SP.IIR_Design.Butter_FreqResponse(freq, 10, freq_l, freq_u);
                //double[] magspec = Pach_SP.IIR_Design.Butter_FreqResponse(105, freq_l, freq_u, Sample_Freq,h.Length);
                double[] magspec = Pach_SP.IIR_Design.Butter_FreqResponse(105, freq_l, freq_u, Sample_Freq, h.Length / 2);

                /////////////Use Zero Phase Bandpass//////////////////////
                //Array.Resize(ref magspec, h.Length / 2);
                System.Numerics.Complex[] filter = Mirror_Spectrum(magspec);

                //double[] f_r = new double[filter.Length];
                //double[] f_i = new double[filter.Length];
                //double[] f_m = new double[filter.Length];

                //for (int i = 0; i < filter.Length; i++) { f_r[i] = filter[i].Real; f_i[i] = filter[i].Imaginary; f_m[i] = filter[i].Magnitude; }

                /////////////Minimum Phase Bandpass///////////////////////
                //double[] response = Minimum_Phase_Response(magspec, 44100, 0);
                //filter = FFT_General(response, 0);

                //Convolve signal with Bandpass Filter.
                Complex[] freq_h = FFT_General(h, thread);

                for (int i = 0; i < freq_h.Length; i++)
                {
                    freq_h[i] *= filter[i];
                }

                double[] h_oct = IFFT_Real_General(freq_h, thread);
                Scale(ref h_oct);

                Array.Resize(ref h_oct, length);

                return h_oct;
            }

            public static double[] Filter2Signal(double[] Filter_in, double[] OctaveSWL, int SampleFrequency, int threadid)
            {
                double[] octave_Pressure = new double[8];
                for (int i = 0; i < 8; i++) octave_Pressure[i] = Math.Pow(10, (OctaveSWL[i]) / 20); //Math.Pow(10,(SWL[5] - SWL[i])/20);

                int length = 8192;
                double[] pulse = (Audio.Pach_SP.Filter is Audio.Pach_SP.Linear_Phase_System) ? Audio.Pach_SP.Minimum_Phase_Signal(octave_Pressure, SampleFrequency, length, threadid) : Audio.Pach_SP.Linear_Phase_Signal_Centered(octave_Pressure, SampleFrequency, length, threadid);//
                //for (int i = 0; i < pulse.Length; i++) pulse[i] /= Math.Sqrt(2);
                //pulse = pulse.Reverse().ToArray();
                Raised_HCosine_Window(ref pulse);// new double[pulse.Length / 2];
                //Array.Resize(ref pulse2, 4096);
                //for (int i = 0; i < pulse.Length / 2; i++) pulse2[i] = pulse[pulse.Length / 2 - 1 + i];
                //pulse2 = pulse2.Reverse().ToArray();
                double[] Signal = FFT_Convolution_double(Filter_in, pulse, threadid);
                int qw = pulse.Length / 4;
                if (Audio.Pach_SP.Filter is Audio.Pach_SP.Minimum_Phase_System) for (int i = qw; i < Signal.Length; i++) Signal[i - qw] = Signal[i];
                Array.Resize(ref Signal, Filter_in.Length);
                return Signal;
            }

            public static System.Numerics.Complex[] FFT4096(double[] Signal, int threadid)
            {
                lock (fftlock4096[threadid])
                {

                    //FFTW.Net Setup//
                    FFT_ArrayIn4096[threadid].SetData(Signal);
                    FFT_Plan4096[threadid].Execute();

                    System.Numerics.Complex[] Out = FFT_ArrayOut4096[threadid].GetData_Complex();
                    return Out;
                }
            }

            public static double[] IFFT_Real4096(System.Numerics.Complex[] spectrum, int threadid)
            {
                lock (ifftlock4096[threadid])
                {

                    //FFTW.Net Setup//
                    IFFT_ArrayIn4096[threadid].SetData(spectrum);
                    IFFT_Plan4096[threadid].Execute();

                    double[] Out = IFFT_ArrayOut4096[threadid].GetData_Real();

                    return Out;
                }
            }

            public static Complex[] IFFT4096(System.Numerics.Complex[] spectrum, int threadid)
            {
                lock (ifftlock[threadid])
                {
                    //FFTW.Net Setup//
                    IFFT_ArrayIn4096[threadid].SetData(spectrum);
                    IFFT_Plan4096[threadid].Execute();

                    System.Numerics.Complex[] Out = IFFT_ArrayOut4096[threadid].GetData_Complex();

                    return Out;
                }
            }

            public static System.Numerics.Complex[] FFT_General(Complex[] Signal, int threadid)
            {
                //FFTW.Net Setup//
                lock (fftlock[threadid])
                {
                    FFT_ArrayIn[threadid] = new fftw_complexarray(Signal);
                    FFT_ArrayOut[threadid] = new fftw_complexarray(Signal.Length);
                    FFT_Plan[threadid] = fftw_plan.dft_1d(Signal.Length, FFT_ArrayIn[threadid], FFT_ArrayOut[threadid], fftw_direction.Forward, fftw_flags.Estimate);
                    FFT_Plan[threadid].Execute();

                    System.Numerics.Complex[] Out = FFT_ArrayOut[threadid].GetData_Complex();
                    return Out;
                }
            }

            public static System.Numerics.Complex[] FFT_General(double[] Signal, int threadid)
            {
                lock (fftlock[threadid])
                {
                    double[] Sig_complex = new double[Signal.Length * 2];
                    for (int i = 0; i < Signal.Length; i++) Sig_complex[i * 2] = Signal[i];

                    //FFTW.Net Setup//
                    FFT_ArrayIn[threadid] = new fftw_complexarray(Sig_complex);
                    FFT_ArrayOut[threadid] = new fftw_complexarray(Signal.Length);
                    FFT_Plan[threadid] = fftw_plan.dft_1d(Signal.Length, FFT_ArrayIn[threadid], FFT_ArrayOut[threadid], fftw_direction.Forward, fftw_flags.Estimate);

                    FFT_Plan[threadid].Execute();

                    System.Numerics.Complex[] Out = FFT_ArrayOut[threadid].GetData_Complex();
                    return Out;
                }
            }

            public static double[] IFFT_Real_General(System.Numerics.Complex[] spectrum, int threadid)
            {
                lock (ifftlock[threadid])
                {
                    //FFTW.Net Setup//
                    IFFT_ArrayIn[threadid] = new fftw_complexarray(spectrum);
                    IFFT_ArrayOut[threadid] = new fftw_complexarray(spectrum.Length);
                    IFFT_Plan[threadid] = fftw_plan.dft_1d(spectrum.Length, IFFT_ArrayIn[threadid], IFFT_ArrayOut[threadid], fftw_direction.Backward, fftw_flags.Estimate);

                    IFFT_Plan[threadid].Execute();

                    double[] Out = IFFT_ArrayOut[threadid].GetData_Real();

                    return Out;
                }
            }

            public static Complex[] IFFT_General(System.Numerics.Complex[] spectrum, int threadid)
            {
                lock (fftlock[threadid])
                {
                    //FFTW.Net Setup//
                    IFFT_ArrayIn[threadid] = new fftw_complexarray(spectrum);
                    IFFT_ArrayOut[threadid] = new fftw_complexarray(spectrum.Length);
                    IFFT_Plan[threadid] = fftw_plan.dft_1d(spectrum.Length, IFFT_ArrayIn[threadid], IFFT_ArrayOut[threadid], fftw_direction.Backward, fftw_flags.Estimate);

                    IFFT_Plan[threadid].Execute();

                    System.Numerics.Complex[] Out = IFFT_ArrayOut[threadid].GetData_Complex();

                    return Out;
                }
            }

            public static Complex[] Mirror_Spectrum(double[] spectrum)
            {
                //Mirror Spectrum
                Complex[] samplep = new Complex[spectrum.Length * 2];
                samplep[0] = spectrum[0];
                samplep[spectrum.Length] = spectrum[spectrum.Length - 1];
                for (int i = 1; i < spectrum.Length; i++)
                {
                    samplep[i] = spectrum[i];
                    samplep[samplep.Length - i] = spectrum[i];
                }
                return samplep;
            }

            public static Complex[] Mirror_Spectrum(Complex[] spectrum)
            {
                //Mirror Spectrum
                Complex[] samplep = new Complex[spectrum.Length * 2];
                samplep[0] = spectrum[0];
                samplep[spectrum.Length] = spectrum[spectrum.Length - 1];
                for (int i = 1; i < spectrum.Length; i++)
                {
                    samplep[i] = spectrum[i];
                    samplep[samplep.Length - i] = Complex.Conjugate(spectrum[i]);
                }
                return samplep;
            }

            /// <summary>
            /// Standard resampling algorithm...
            /// </summary>
            /// <param name="Input">input signal</param>
            /// <param name="inputFS">input sampling frequency</param>
            /// <param name="outputFS">target sampling frequency</param>
            /// <param name="threadid">thread id if applicable</param>
            /// <param name="normalized">is this signal normalized? If so, we will not scale it to match the energy of the original.</param>
            /// <returns>resampled signal</returns>
            public static double[] Resample(double[] Input, int inputFS, int outputFS, int threadid, bool normalized)
            {
                double[] IN = Input.Clone() as double[];
                Array.Resize(ref IN, Input.Length * 2 - 1);
                Complex[] Spec = Pach_SP.FFT_General(IN, threadid);
                Spec.Reverse();
                Array.Resize(ref Spec, (int)(((double)outputFS / (double)inputFS) * IN.Length / 2));
                Spec.Reverse();
                Spec = Pach_SP.Mirror_Spectrum(Spec);
                double[] s_out = Pach_SP.IFFT_Real_General(Spec, threadid);
                double factor = normalized ? 1 : Math.Sqrt(s_out.Length) * Math.Sqrt(IN.Length);
                for (int i = 0; i < s_out.Length; i++) s_out[i] /= factor;
                return s_out;
            }

            /// <summary>
            /// Resampling algorithm using cubic spline interpolation...
            /// </summary>
            /// <param name="Input">input signal</param>
            /// <param name="inputFS">input sampling frequency</param>
            /// <param name="outputFS">target sampling frequency</param>
            /// <param name="threadid">thread id if applicable</param>
            /// <param name="normalized">is this signal normalized? If so, we will not scale it to match the energy of the original.</param>
            /// <returns>resampled signal</returns>
            public static double[] Resample_Cubic(double[] Input, int inputFS, int outputFS, int threadid, bool normalized = false)
            {
                double[] Time = new double[Input.Length];
                double maxtime = (double)Input.Length / inputFS;
                int max_sample = (int)Math.Floor(outputFS * maxtime);
                for (int i = 0; i < Input.Length; i++) Time[i] = (double)i / inputFS;
                MathNet.Numerics.Interpolation.CubicSpline res = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(Time, Input);
                double Factor = normalized ? 1 : (double)outputFS / inputFS;
                double[] output = new double[max_sample];
                for (int i = 0; i < max_sample; i++)
                {
                    output[i] = res.Interpolate((double)i / outputFS) / Factor;
                }
                return output;
            }

            public static void ScaleRoot(ref Complex[] IN)
            {
                double factor = Math.Sqrt(IN.Length);
                for (int i = 0; i < IN.Length; i++) IN[i] /= factor;
            }

            public static void ScaleRoot(ref double[] IN)
            {
                double factor = Math.Sqrt(IN.Length);
                for (int i = 0; i < IN.Length; i++) IN[i] /= factor;
            }

            public static void Scale(ref Complex[] IN)
            {
                for (int i = 0; i < IN.Length; i++) IN[i] /= IN.Length;
            }

            public static void Scale(ref double[] IN)
            {
                for (int i = 0; i < IN.Length; i++) IN[i] /= IN.Length;
            }

            public static double[] Magnitude_Spectrum(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                double[] p_i = new double[length_starttofinish / 2];
                double df = (double)(sample_frequency) / (double)length_starttofinish;

                int oct = 0;
                ///////////////////////////////////////////////////////////////////
                ///Do the bottom frequencies
                double ctr = 62.5 * Math.Pow(2, oct);
                double freq_l = ctr / Utilities.Numerics.rt2;
                double freq_u = ctr * Utilities.Numerics.rt2;
                int idl = 0;
                int idu = (int)Math.Round(freq_u / df);
                int ctr_id = ((idu + (int)Math.Round(freq_l / df)) / 2);
                double l = Math.Sqrt(idu - idl) * Math.Sqrt(2);
                int length;
                if (Octave_pressure[1] == 0)
                {
                    //do nothing. boundaries stay as they are...
                    length = idu;
                }
                else if (Octave_pressure[1] * .25 / Octave_pressure[0] > 0.002503912)
                {
                    ///Use default overflow...
                    //TODO: Default overflow needs work...
                    length = (int)(idu / (1 - 1d / 16));
                }
                else
                {
                    //calculate how much overflow can be tolerated...
                    double ph = RCos_Integral.Integrate(Octave_pressure[1] * .25 / Octave_pressure[0]);
                    //Idu is at ph. end of total bin is 0 to idu + ph.
                    length = (int)(idu / (1 - ph));
                }

                double[] magspec = new double[length];

                double mod = 0;

                int t_length = ctr_id - 0;
                for (int i = 0; i < t_length; i++)
                {
                    double v = 0;
                    v = (-Math.Cos(Math.PI * (double)i / t_length) + 1);
                    magspec[i] = v;
                    mod += magspec[i];
                }

                t_length = magspec.Length - ctr_id;
                for (int i = 0; i < t_length; i++)
                {
                    double v = (-Math.Cos(Math.PI + Math.PI * (double)i / t_length) + 1);
                    magspec[i + ctr_id] = v * v;
                    mod += magspec[i];
                }

                //TODO: handle step discontinuity resuling at high end of this bin.

                mod = Octave_pressure[oct] * Octave_pressure[oct] / (mod * 2);
                for (int i = 0; i < magspec.Length; i++)
                {
                    if ((idl + i) <= p_i.Length)
                        p_i[idl + i] += Math.Sqrt(magspec[i] * mod);
                }

                ///Do the top frequencies...
                oct = 7;
                ctr = 62.5 * Math.Pow(2, oct);
                freq_l = ctr / Utilities.Numerics.rt2;
                freq_u = ctr * Utilities.Numerics.rt2;
                idl = (int)Math.Round(freq_l / df);
                idu = (int)Math.Round(freq_u / df);
                ctr_id = (idu + idl) / 2;
                l = Math.Sqrt(idu - idl) * Math.Sqrt(2);
                length = (int)(length_starttofinish / 2) - idl;
                if (Octave_pressure[6] == 0)
                {
                    //do nothing. boundaries stay as they are...
                    length = (idu - idl);
                }
                else if (Octave_pressure[6] * .25 / Octave_pressure[7] > 0.019572846)
                {
                    ///Use default overflow... 
                    /// TODO: Default Overflow is incorrect...
                    length = (int)((idu - idl) / (1 - 1d / 8)); //Check that this does not overload the neighboring bin...
                }
                else
                {
                    //calculate how much overflow  can be tolerated...
                    double ph = RCos_Integral.Integrate(Octave_pressure[6] * .25 / Octave_pressure[7]);
                    //Idu is at ph. end of total bin is 0 to idu + ph.
                    int add = (int)((idu - idl) / (1 - ph)) - (idu - idl);
                    length += add;
                    idl -= add;
                }

                magspec = new double[length];

                mod = 0;
                t_length = ctr_id - idl;
                for (int i = 0; i < t_length; i++)
                {
                    double v = 0;
                    v = (-Math.Cos(Math.PI * (double)i / t_length) + 1);
                    magspec[i] = v * v;
                    mod += magspec[i];
                }

                t_length = idu - ctr_id;
                for (int i = 0; i < t_length; i++)
                {
                    double v = 0;
                    v = (-Math.Cos(Math.PI + Math.PI * (double)i / t_length) + 1);
                    magspec[i + ctr_id - idl] = v * v;
                    if (i < (idu - idl)) mod += magspec[i];
                }
                mod = Octave_pressure[oct] * Octave_pressure[oct] / (mod * 2);
                for (int i = 0; i < magspec.Length; i++)
                {
                    if ((idl + i) > p_i.Length) break;
                    p_i[idl + i] += Math.Sqrt(magspec[i] * mod);
                }

                //////////////////////////////////////////////////////////////////
                for (oct = 1; oct < 7; oct++)
                {
                    if (Octave_pressure[oct] == 0) continue;

                    ctr = 62.5 * Math.Pow(2, oct);
                    freq_l = ctr / Utilities.Numerics.rt2;
                    freq_u = ctr * Utilities.Numerics.rt2;
                    idl = (int)Math.Round(freq_l / df);
                    idu = (int)Math.Round(freq_u / df);
                    l = Math.Sqrt(idu - idl) * Math.Sqrt(2);
                    length = idu - idl;
                    ctr_id = (idu + idl) / 2;
                    double ph_l = 0;
                    double ph_u = 0;
                    if (Octave_pressure[oct + 1] == 0)
                    {
                        //do nothing. boundaries stay as they are...
                    }
                    else if (Octave_pressure[oct + 1] * .25 / Octave_pressure[oct] > 0.019572846) //The area under the curve at phase = pi/4
                    {
                        ///Use default overflow...
                        //length += ((idu - (int)(ctr / df)) / 2);
                        ph_u = 1d / 8;
                    }
                    else
                    {
                        //calculate how much overflow  can be tolerated...
                        double ph = RCos_Integral.Integrate(Octave_pressure[oct + 1] * .25 / Octave_pressure[oct]);
                        //Idu is at ph. end of total bin is 0 to idu + ph.
                        ph_u = ph;
                    }

                    if (Octave_pressure[oct - 1] == 0)
                    {
                        //do nothing. boundaries stay as they are...
                    }
                    else if (Octave_pressure[oct - 1] * .25 / Octave_pressure[oct] > 0.002503912)//The area under the curve at phase = pi/8
                    {
                        ///Use default overflow...
                        //length += ((idu - (int)(ctr / df)) / 2);
                        ph_l += 1d / 16;
                    }
                    else
                    {
                        //calculate how much overflow  can be tolerated...
                        double ph = RCos_Integral.Integrate(Octave_pressure[oct - 1] * .25 / Octave_pressure[oct]);
                        //Idu is at ph. end of total bin is 0 to idu + ph.
                        ph_l += ph;
                    }

                    length = (int)(length / (1 - ph_l - ph_u));
                    idl -= (int)(length * ph_l);
                    magspec = new double[length];

                    mod = 0;
                    t_length = ctr_id - idl;
                    for (int i = 0; i < t_length; i++)
                    {
                        double v = 0;
                        v = (-Math.Cos(Math.PI * (double)i / t_length) + 1);
                        magspec[i] = v * v;
                        mod += magspec[i];
                    }

                    t_length = idu - ctr_id;
                    for (int i = 0; i < t_length; i++)
                    {
                        double v = 0;
                        v = (-Math.Cos(Math.PI + Math.PI * (double)i / t_length) + 1);
                        magspec[i + ctr_id - idl] = v * v;
                        mod += magspec[i];
                    }

                    mod = Octave_pressure[oct] * Octave_pressure[oct] / (mod * 2);
                    for (int i = 0; i < magspec.Length; i++)
                    {
                        if ((idl + i) >= p_i.Length) break;
                        p_i[idl + i] += Math.Sqrt(magspec[i] * mod);
                    }
                }
                return p_i;
            }

            public static double[] Magnitude_Filter(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                int spec_length = length_starttofinish / 2;
                List<double> f = new List<double>();
                List<double> pr = new List<double>();

                if (Octave_pressure.Length == 8)
                {
                    f.Add(0);

                    for (int oct = 0; oct < 9; oct++)
                    {
                        f.Add(62.5 * Math.Pow(2, oct) / Utilities.Numerics.rt2);
                        f.Add(62.5 * Math.Pow(2, oct) * .8);
                        f.Add(62.5 * Math.Pow(2, oct));
                        f.Add(62.5 * Math.Pow(2, oct) * 1.2);
                    }
                    f.Add(sample_frequency / 2);

                    pr.Add(Octave_pressure[0] * 0.0000000001);
                    pr.Add(Octave_pressure[0] * 0.00000001);
                    pr.Add(Octave_pressure[0] * 0.001);
                    pr.Add(Octave_pressure[0]);
                    pr.Add(Octave_pressure[0]);

                    for (int oct = 1; oct < 7; oct++)
                    {
                        double min = Math.Min(Octave_pressure[oct], pr.Last());
                        double max = Math.Max(Octave_pressure[oct], pr.Last());
                        pr.Add((max * .95 + min * .05));
                        pr.Add(Octave_pressure[oct]);
                        pr.Add(Octave_pressure[oct]);
                        pr.Add(Octave_pressure[oct]);
                    }
                    if (pr.Count < f.Count) pr.Add(Octave_pressure[7]);//8k
                    if (pr.Count < f.Count) pr.Add(Octave_pressure[7] * .9);//10k
                    if (pr.Count < f.Count) pr.Add(Octave_pressure[7] * .5);//12k
                    if (pr.Count < f.Count) pr.Add(Octave_pressure[7] * .25);//16k
                    int mod = 0;
                    while (pr.Count < f.Count)
                    {
                        mod++;
                        pr.Add(Octave_pressure[7] * .125 * Math.Pow(.5, mod));
                    }
                }
                else if (Octave_pressure.Length == 24)
                {
                    if (Octave_pressure.Length == 8)
                    {
                        f.Add(20);
                        pr.Add(Octave_pressure[0] * 0.0001);

                        for (int oct = 0; oct < 24; oct++)
                        {
                            f.Add(20 * Math.Pow(2, 1 + oct / 3.0));
                            pr.Add(Octave_pressure[oct]);
                        }
                        f.Add(sample_frequency / 2);
                        pr.Add(Octave_pressure[23] * 0.0001);

                    }
                }
                else throw new Exception();
                //MathNet.Numerics.Interpolation.CubicSpline prm = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkimaSorted(f.ToArray(), pr.ToArray(),);
                MathNet.Numerics.Interpolation.LinearSpline prm = MathNet.Numerics.Interpolation.LinearSpline.InterpolateSorted(f.ToArray(), pr.ToArray());

                double[] p_i = new double[spec_length];

                for (int j = 0; j < spec_length; j++)
                {
                    double fr = (j + 1) * (sample_frequency / 2) / spec_length;
                    p_i[j] = prm.Interpolate(fr);// / (fr);
                    if (p_i[j] < 0) p_i[j] = 0;
                }

                return p_i;
            }

            public static double[] Magnitude_Filter_BW(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                //double rt2 = Math.Sqrt(2);
                int spec_length = length_starttofinish / 2;

                //List<double> f = new List<double>();
                //List<double> fl = new List<double>();
                //List<double> fu = new List<double>();
                //f.Add(0);
                //f.Add(31.25 * rt2);

                //for (int oct = 0; oct < 9; oct++)
                //{
                //    f.Add(62.5 * Math.Pow(2, oct));
                //    f.Add(rt2 * 62.5 * Math.Pow(2, oct));
                //}
                //f.Add(sample_frequency / 2);

                double[] output = new double[spec_length];
                double[] samplep = new double[spec_length];
                double fs_l = sample_frequency / spec_length;

                for (int oct = 0; oct < 8; oct++)
                {
                    double f = 62.5 * Math.Pow(2, oct);
                    double fl = 1.1 * f / Utilities.Numerics.rt2;
                    double fu = .9 * f * Utilities.Numerics.rt2;
                    samplep = Audio.Pach_SP.IIR_Design.Butter_FreqResponse(10, fl, fu, sample_frequency, spec_length);
                    //double m = samplep.Max();
                    for (int i = 0; i < spec_length; i++)
                    {
                        //double w_f = i * fs_l;
                        //output[i] = Octave_pressure[oct] * 1025 / (Math.Sqrt(1 + Math.Pow(w_f / fl, 20)) * Math.Sqrt(1 + Math.Pow(fu / w_f, 20)));
                        output[i] = Math.Max(Octave_pressure[oct] * samplep[i], output[i]);///m;
                    }
                }

                return output;

                //List<double> pr = new List<double>();
                //pr.Add(Octave_pressure[0]);
                //pr.Add(Octave_pressure[0]);

                //for (int oct = 0; oct < 7; oct++)
                //{
                //    pr.Add(Octave_pressure[oct]);
                //    pr.Add((Octave_pressure[oct] + Octave_pressure[oct + 1]) / 2);
                //}
                //if (pr.Count < f.Count) pr.Add(Octave_pressure[7]);//8k
                //if (pr.Count < f.Count) pr.Add(Octave_pressure[7]*.5);//10k
                //if (pr.Count < f.Count) pr.Add(Octave_pressure[7]*.25);//12k
                //if (pr.Count < f.Count) pr.Add(Octave_pressure[7]*.125);//16k
                //int mod = 0;
                //while (pr.Count < f.Count)
                //{
                //    mod++;
                //    pr.Add(Octave_pressure[7] * .125 * Math.Pow(.5,mod));
                //}

                //MathNet.Numerics.Interpolation.CubicSpline prm = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkimaSorted(f.ToArray(), pr.ToArray());

                //double[] p_i = new double[spec_length];

                //for (int j = 0; j < spec_length; j++)
                //{
                //    double fr = (j + 1) * (sample_frequency / 2) / spec_length;
                //    p_i[j] = prm.Interpolate(fr);// / (fr);
                //    if (p_i[j] < 0) p_i[j] = 0;
                //}

                //return p_i;
            }

            public static System.Numerics.Complex[] Minimum_Phase_Spectrum(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                double[] M_spec = Magnitude_Spectrum(Octave_pressure, sample_frequency, length_starttofinish, threadid);

                System.Numerics.Complex[] logspec = new System.Numerics.Complex[M_spec.Length];
                for (int i = 0; i < M_spec.Length; i++)
                {
                    if (M_spec[i] == 0) logspec[i] = -100;
                    else logspec[i] = Math.Log(M_spec[i]);
                }

                double[] real_cepstrum = IFFT_Real_General(Mirror_Spectrum(logspec), threadid);
                Scale(ref real_cepstrum);

                double[] ym = new double[length_starttofinish / 2];
                ym[0] = real_cepstrum[0];

                for (int i = 1; i < length_starttofinish / 4; i++)
                {
                    ym[i] = 2 * real_cepstrum[i];
                }
                ym[length_starttofinish / 4] = real_cepstrum[length_starttofinish / 4];
                System.Numerics.Complex[] ymspec = FFT_General(ym, threadid);

                for (int i = 0; i < ymspec.Length; i++)
                {
                    ymspec[i] = Complex.Exp(ymspec[i]);
                }

                return ymspec;
            }

            public static double[] Minimum_Phase_Response(double[] M_spec, int sample_frequency, int threadid)
            {
                System.Numerics.Complex[] logspec = new System.Numerics.Complex[M_spec.Length];
                for (int i = 0; i < M_spec.Length; i++)
                {
                    if (M_spec[i] == 0) logspec[i] = -100;
                    else logspec[i] = Math.Log(M_spec[i]);
                }

                double[] real_cepstrum = IFFT_Real4096(Mirror_Spectrum(logspec), threadid);
                Scale(ref real_cepstrum);

                double[] ym = new double[real_cepstrum.Length];
                ym[0] = real_cepstrum[0];

                for (int i = 1; i < M_spec.Length; i++)
                {
                    ym[i] = 2 * real_cepstrum[i];
                }
                ym[M_spec.Length] = real_cepstrum[M_spec.Length];
                System.Numerics.Complex[] ymspec = FFT4096(ym, threadid);

                for (int i = 0; i < ymspec.Length; i++)
                {
                    ymspec[i] = Complex.Exp(ymspec[i]);
                }

                double[] Signal = IFFT_Real4096(ymspec, threadid);
                Scale(ref Signal);

                double[] S2 = new double[Signal.Length];
                for (int i = 0; i < Signal.Length; i++) S2[i] = Signal[(Signal.Length + i + 1) % Signal.Length];
                S2 = S2.Reverse().ToArray();

                return S2;
            }

            public static double[] Minimum_Phase_Signal(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                double[] M_spec = Magnitude_Spectrum(Octave_pressure, sample_frequency, length_starttofinish, threadid);

                double sum_start = 0;
                for (int i = 0; i < M_spec.Length; i++) sum_start += M_spec[i] * M_spec[i];

                System.Numerics.Complex[] logspec = new System.Numerics.Complex[M_spec.Length];

                double min = double.MaxValue;

                for (int i = 0; i < M_spec.Length; i++)
                {
                    if (M_spec[i] != 0) min = Math.Min(M_spec[i], min);
                }

                for (int i = 0; i < M_spec.Length; i++)
                {
                    if (M_spec[i] == 0) logspec[i] = Math.Log(min * 0.01);
                    else logspec[i] = Math.Log(M_spec[i]);
                }

                double[] real_cepstrum = IFFT_Real_General(Mirror_Spectrum(logspec), threadid);
                Scale(ref real_cepstrum);

                double[] ym = new double[length_starttofinish];
                ym[0] = real_cepstrum[0];

                for (int i = 1; i < length_starttofinish / 2; i++)
                {
                    ym[i] = 2 * real_cepstrum[i];
                }
                ym[length_starttofinish / 2] = real_cepstrum[length_starttofinish / 2];
                System.Numerics.Complex[] ymspec = FFT_General(ym, threadid);

                for (int i = 0; i < ymspec.Length; i++)
                {
                    ymspec[i] = Complex.Exp(ymspec[i]);
                }

                double[] Signal = IFFT_Real_General(ymspec, threadid);
                Scale(ref Signal);
                double[] S2 = new double[Signal.Length];
                for (int i = 0; i < Signal.Length; i++) S2[i] = Signal[i];

                return S2;
            }

            public static double[] Linear_Phase_Signal(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                double[] M_s = Magnitude_Spectrum(Octave_pressure, sample_frequency, length_starttofinish, threadid);
                Complex[] M_spec = new Complex[M_s.Length];

                for (int i = 0; i < M_s.Length; i++)
                {
                    M_spec[i] = M_s[i];
                }

                ///////////////Create Zero Phase Time Domain Filter/////////////
                double[] prefilter = IFFT_Real_General(Mirror_Spectrum(M_spec), threadid);
                double scale = Math.Pow(prefilter.Length, 2.0);
                //Rotate filter in time domain to center.
                int hw = prefilter.Length / 2;
                double[] filter = new double[prefilter.Length];

                for (int i = 0; i < prefilter.Length; i++)
                {
                    filter[i] = prefilter[i] / scale;
                }

                return filter;
            }

            public static double[] Linear_Phase_Signal_Centered(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                double[] M_s = Magnitude_Spectrum(Octave_pressure, sample_frequency, length_starttofinish, threadid);
                Complex[] M_spec = new Complex[M_s.Length];

                for (int i = 0; i < M_s.Length; i++)
                {
                    M_spec[i] = M_s[i];
                }

                ///////////////Create Zero Phase Time Domain Filter/////////////
                double[] prefilter = IFFT_Real_General(Mirror_Spectrum(M_spec), threadid);
                double scale = Math.Pow(prefilter.Length, 2.0);
                //Rotate filter in time domain to center.
                int hw = prefilter.Length / 2;
                double[] filter = new double[prefilter.Length];

                for (int i = 0; i < prefilter.Length; i++)
                {
                    filter[i] = prefilter[(i + hw) % prefilter.Length] / scale;
                }

                return filter;
            }

            public static Complex[] Linear_Phase_Spectrum(double[] Octave_pressure, int sample_frequency, int length_starttofinish, int threadid)
            {
                double[] M_s = Magnitude_Spectrum(Octave_pressure, sample_frequency, length_starttofinish, threadid);
                Complex[] M_spec = new Complex[M_s.Length];

                for (int i = 0; i < M_s.Length; i++)
                {
                    M_spec[i] = M_s[i];
                }

                ///////////////Create Zero Phase Time Domain Filter/////////////
                double[] prefilter = IFFT_Real4096(Mirror_Spectrum(M_spec), threadid);
                double scale = 2 * Math.Pow(prefilter.Length, 3.0 / 2.0);
                //Rotate filter in time domain to center.
                int hw = prefilter.Length / 2;
                double[] filter = new double[prefilter.Length];

                for (int i = 0; i < prefilter.Length; i++)
                {
                    filter[i] = prefilter[(i + hw) % prefilter.Length] / scale;
                }

                Complex[] spectrum = FFT4096(filter, threadid);
                Array.Resize<Complex>(ref spectrum, hw);
                return spectrum;
            }

            public static double[] Linear_Phase_Response(double[] M_spec, int sample_frequency, int threadid)
            {
                ///////////////Create Zero Phase Time Domain Filter/////////////
                double[] prefilter = IFFT_Real4096(Mirror_Spectrum(M_spec), threadid);
                double scale = 2 * Math.Pow(prefilter.Length, 2.0);
                //Rotate filter in time domain to center.
                int hw = prefilter.Length / 2;
                double[] filter = new double[prefilter.Length];
                //double[] SPL = new double[prefilter.Length];

                for (int i = 0; i < prefilter.Length; i++)
                {
                    filter[i] = prefilter[(i + hw) % prefilter.Length] / scale;
                    //SPL[i] = 10 * Math.Log10(Math.Pow(filter[i] / 20E-6, 2));
                    //if (double.IsInfinity(SPL[i])) SPL[i] = -170;
                }

                //for (int i = 0; i < filter.Length; i++)
                //{
                //    SPL[i] = 10 * Math.Log10(Math.Pow(filter[i] / 20E-6, 2));
                //    if (double.IsInfinity(SPL[i])) SPL[i] = -170;
                //}

                return filter;
            }

            public static double[] Arbitrary_Phase_Response(System.Numerics.Complex[] M_spec, int sample_frequency, int threadid)
            {
                ///////////////Create Zero Phase Time Domain Filter/////////////
                double[] prefilter = (M_spec.Length == 2048) ? IFFT_Real4096(Mirror_Spectrum(M_spec), threadid) : IFFT_Real_General(Mirror_Spectrum(M_spec), threadid);
                double scale = 2 * Math.Pow(prefilter.Length, 2.0);
                //Rotate filter in time domain to center.
                int hw = prefilter.Length / 2;
                double[] filter = new double[prefilter.Length];

                for (int i = 0; i < prefilter.Length; i++)
                {
                    filter[i] = prefilter[(i + hw) % prefilter.Length] / scale;
                }

                return filter;
            }

            public static double[] ETCToFilter(double[][] Octave_PRMS, double[] SWL, int sample_frequency_in = 44100, int sample_frequency_out = 44100)
            {
                int length = 4096;
                double BW = (double)sample_frequency_out / (double)sample_frequency_in;
                double[] IR = new double[(int)Math.Floor((double)(Octave_PRMS[0].Length * BW)) + (int)length];

                double[] p_mod = new double[8];
                for (int i = 0; i < 8; i++) p_mod[i] = Math.Pow(10, (120 - SWL[i]) / 20);

                //Convert to Pressure & Interpolate full resolution IR
                int proc = System.Environment.ProcessorCount;
                double[][] output = new double[proc][];
                double[][] samplep = new double[proc][];
                System.Threading.Thread[] T = new System.Threading.Thread[proc];
                int[] to = new int[proc];
                int[] from = new int[proc];

                Task[] tsk = new Task[proc];

                //Action<int> processchunk = (p) =>
                Parallel.For(0,proc, p =>
                {
                    output[p] = new double[length];
                    samplep[p] = new double[length * 2];
                    to[p] = p * Octave_PRMS[0].Length / proc;
                    from[p] = (p + 1) * Octave_PRMS[0].Length / proc;

                    for (int t = to[p]; t < from[p]; t++)
                    {
                        double[] pr = new double[8];
                        for (int oct = 0; oct < 8; oct++) pr[oct] = Math.Abs(Octave_PRMS[oct][t]) * p_mod[oct];
                        double sum = 0;
                        foreach (double d in pr) sum += d;
                        if (sum > 0)
                        {
                            output[p] = Filter.Transfer_Function(pr, sample_frequency_out, length, p);
                            for (int k = 0; k < length; k++)
                            {
                                IR[(int)Math.Floor(t * BW) + k] += output[p][k];
                            }
                        }
                    }
                });

                    //for (int p = 0; p < proc; p++)
                    //{
                    //    int id = (int)p;
                    //    //tsk[p] = Task.Run(() => processchunk(id));
                    //    tsk[p] = Task.Factory.StartNew(() => processchunk(id), CancellationToken.None, TaskCreationOptions.LongRunning, TaskScheduler.Default);
                    //}

                //Task.WaitAll(tsk);

                //System.Threading.CountdownEvent CDE = new System.Threading.CountdownEvent(Octave_PRMS[0].Length);

                //for (int p = 0; p < proc; p++)
                //{
                //    output[p] = new double[length];
                //    samplep[p] = new double[length * 2];
                //    to[p] = p * Octave_PRMS[0].Length / proc;
                //    from[p] = (p + 1) * Octave_PRMS[0].Length / proc;

                //    T[p] = new System.Threading.Thread((thread) =>
                //    {
                //        int thr = (int)thread;
                //        for (int t = to[thr]; t < from[thr]; t++)
                //        {
                //            double[] pr = new double[8];
                //            for (int oct = 0; oct < 8; oct++) pr[oct] = Math.Abs(Octave_PRMS[oct][t]) * p_mod[oct];
                //            double sum = 0;
                //            foreach (double d in pr) sum += d;
                //            if (sum > 0)
                //            {
                //                output[thr] = Filter.Transfer_Function(pr, sample_frequency_out, length, thr);
                //                for (int k = 0; k < length; k++)
                //                {
                //                    IR[(int)Math.Floor(t * BW) + k] += output[thr][k];
                //                }
                //            }
                //            CDE.Signal();
                //        }
                //    });
                //    T[p].Start(p);
                //}

                //do
                //{
                //    if (CDE.IsSet)
                //    {
                //        break;
                //    }
                //    System.Threading.Thread.Sleep(500);
                //} while (true);

                return IR;
            }

            /// <summary>
            /// Frequency domain convolution.
            /// </summary>
            /// <param name="SignalBuffer">the dry signal.</param>
            /// <param name="Filter">the pressure domain impulse response.</param>
            /// <returns>the convolved signal.</returns>
            public static float[] FFT_Convolution(double[] SignalBuffer, double[] Filter, int threadid)
            {
                if (SignalBuffer == null) return null;
                int minlength = SignalBuffer.Length > Filter.Length ? SignalBuffer.Length : Filter.Length;

                int W = (int)Math.Pow(2, Math.Ceiling(Math.Log(minlength, 2)));

                if (SignalBuffer.Length < W) Array.Resize(ref SignalBuffer, W);
                if (Filter.Length < W) Array.Resize(ref Filter, W);

                System.Numerics.Complex[] freq1 = FFT_General(SignalBuffer, threadid);
                System.Numerics.Complex[] freq2 = FFT_General(Filter, threadid);
                System.Numerics.Complex[] freq3 = new System.Numerics.Complex[W];

                for (int i = 0; i < freq1.Length; i++) freq3[i] = freq1[i] * freq2[i];

                double[] conv = IFFT_Real_General(freq3, threadid);

                float[] output = new float[conv.Length];
                //double mod = 1d / Math.Sqrt(conv.Length);
                double mod = 1d / conv.Length;
                for (int i = 0; i < conv.Length; i++) output[i] = (float)(conv[i] * mod);// * mod);

                double maxfilt = Filter.Max();
                double maxsig = SignalBuffer.Max();
                double max = conv.Max();
                double outmax = output.Max();
                max++;
                maxsig++;
                maxfilt++;
                return output;
            }

            /// <summary>
            /// Frequency domain convolution.
            /// </summary>
            /// <param name="SignalBuffer">the dry signal.</param>
            /// <param name="Filter">the pressure domain impulse response.</param>
            /// <returns>the convolved signal.</returns>
            public static double[] FFT_Convolution_double(double[] SignalBuffer, double[] Filter, int threadid)
            {
                if (SignalBuffer == null) return null;
                int minlength = SignalBuffer.Length > Filter.Length ? SignalBuffer.Length : Filter.Length;

                int W = (int)Math.Pow(2, Math.Ceiling(Math.Log(minlength, 2)));

                if (SignalBuffer.Length < W) Array.Resize(ref SignalBuffer, W);
                if (Filter.Length < W) Array.Resize(ref Filter, W);

                System.Numerics.Complex[] freq1 = FFT_General(SignalBuffer, threadid);
                System.Numerics.Complex[] freq2 = FFT_General(Filter, threadid);

                System.Numerics.Complex[] freq3 = new System.Numerics.Complex[W];

                for (int i = 0; i < freq1.Length; i++) freq3[i] = freq1[i] * freq2[i];

                double[] conv = IFFT_Real_General(freq3, threadid);

                double[] output = new double[conv.Length];
                double mod = 1d / conv.Length;
                for (int i = 0; i < conv.Length; i++) output[i] = conv[i] * mod;// * mod;

                double maxfilt = Filter.Max();
                double maxsig = SignalBuffer.Max();
                double max = conv.Max();
                double outmax = output.Max();
                max++;
                maxsig++;
                maxfilt++;
                return output;
            }

            /// <summary>
            /// Frequency domain deconvolution.
            /// </summary>
            /// <param name="Recording">the recorded signal</param>
            /// <param name="Signal">the signal used to excite the space.</param>
            /// <returns>the deconvolved IR (hopefully).</returns>
            public static double[] FFT_Deconvolution(double[] Recording, double[] Signal)
            {
                if (Recording == null) return null;
                int minlength = Recording.Length > Signal.Length ? Recording.Length : Signal.Length;

                int W = (int)Math.Pow(2, Math.Ceiling(Math.Log(minlength, 2)));

                if (Recording.Length < W) Array.Resize(ref Recording, W);
                if (Signal.Length < W) Array.Resize(ref Signal, W);

                System.Numerics.Complex[] Recording_FD = FFT_General(Recording, 0);
                System.Numerics.Complex[] Signal_FD = FFT_General(Signal, 0);

                Complex[] Quotient = new Complex[W];

                for (int f = 0; f < W; f++)
                {
                    Quotient[f] = Recording_FD[f] / Signal_FD[f];
                }

                double[] Complete = IFFT_Real_General(Quotient, 0);

                double[] signalout = new double[W];

                double maxsignal = Complete.Max();

                for (int q = 0; q < Recording.Length; q++) signalout[q] = (float)(Complete[q] * ((0.75) / maxsignal));

                //crop zeros from final signal 
                for (int q = signalout.Length - 1; q > 0; q--)
                {
                    if (signalout[q] != 0)
                    {
                        Array.Resize(ref signalout, q + 2);
                        break;
                    }
                }
                return signalout;
            }

            public static double[] Filter_Interpolation(double[] SWL, double[][] ETC, int SampleRate_IN, int SampleRate_Out, double Rho_C, IProgressFeedback VB = null)
            {
                double[][] SPLetc = new double[8][];
                double[] Total_E = new double[8];

                double[] time = new double[ETC[0].Length];
                double dtIn = 1.0f / SampleRate_IN;
                double dtOut = 1.0f / SampleRate_Out;
                for (int i = 1; i < time.Length; i++) time[i] = (double)i / SampleRate_IN;

                for (int oct = 0; oct < 8; oct++)
                {
                    SPLetc[oct] = new double[ETC[0].Length];
                    for (int i = 0; i < time.Length; i++)
                    {
                        Total_E[oct] += ETC[oct][i];
                        if (ETC[oct][i] > 0) SPLetc[oct][i] = Math.Log10(ETC[oct][i]);
                        else SPLetc[oct][i] = Math.Log10(1E-12);
                    }
                }

                int NewLength = (int)Math.Ceiling((double)(SampleRate_Out * ETC[0].Length) / SampleRate_IN);
                double[][] NewETC = new double[8][];
                double[] NewTE = new double[8];
                for (int oct = 0; oct < 8; oct++)
                {
                    MathNet.Numerics.Interpolation.CubicSpline CS = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(time, SPLetc[oct]);
                    NewETC[oct] = new double[NewLength];
                    for (int i = 0; i < NewLength; i++)
                    {
                        NewETC[oct][i] = Math.Pow(10, CS.Interpolate(i * dtOut));
                        NewTE[oct] += NewETC[oct][i];
                    }

                    for (int i = 0; i < NewLength; i++)
                    {
                        NewETC[oct][i] *= Total_E[oct] / NewTE[oct];
                    }
                }

                if (VB != null)
                {
                    VB.change_title("Interpolating to Pressure...");
                } 
                return ETCToFilter(NewETC, SWL, SampleRate_Out, SampleRate_Out);
            }

            public static class Wave
            {
                public static double[][] ReadtoDouble(string Path, bool Normalize, out int Sample_Frequency, bool tryunsupported = false)
                {
                    int[][] data;
                    try
                    {
                        data = ReadtoInt(Path, out Sample_Frequency, tryunsupported);
                    }
                    catch (Exception x)
                    {
                        throw x;
                    }

                    double[][] signal = new double[data.Length][];

                    for (int c = 0; c < data.Length; c++)
                    {
                        signal[c] = new double[data[c].Length];
                        if (Normalize)
                        {   
                            int max = Math.Max(data[c].Max(), Math.Abs(data[c].Min()));
                            for (int i = 0; i < data[c].Length; i++)
                            {
                                signal[c][i] = (double)data[c][i] / max;
                            }
                        }
                        else
                        {
                            for (int i = 0; i < data[c].Length; i++)
                            {
                                signal[c][i] = (double)data[c][i] / int.MaxValue;
                            }
                        }
                    }
                    return signal;
                }

                public static int[][] ReadtoInt(string Path, bool Normalize, out int Sample_Frequency, bool tryunsupported = false)
                {

                    int[][] data;

                    try
                    {
                        data = ReadtoInt(Path, out Sample_Frequency, tryunsupported);
                    }
                    catch (Exception x)
                    {
                        MessageBox.Show(x.Message);
                        Sample_Frequency = 44100;
                        return new int[1][];
                    }

                    if (Normalize)
                    {
                        for (int c = 0; c < data.Length; c++)
                        {
                            double mod = (double)int.MaxValue / data[c].Max();
                            for (int i = 0; i < data[c].Length; i++)
                            {
                                data[c][i] = (int)((double)data[c][i] * mod);
                            }
                        }
                    }

                    return data;
                }

                public static int[][] ReadtoInt(string Path, out int Sample_Frequency, bool try_unsupported_data = false)
                {
                    //Open a stream:
                    bool RF64 = false;
                    System.IO.BinaryReader wav = new System.IO.BinaryReader(new System.IO.FileStream(Path, System.IO.FileMode.Open));
                    //Read the header
                    string Riffhdr = new string(wav.ReadChars(4));
                    long filesize = wav.ReadUInt32();
                    long dataChunkLength;
                    if (Riffhdr == "RF64")
                    {
                        RF64 = true;
                        if (new string(wav.ReadChars(4)) != "ds64") throw new Exception("RF64 files must have ds64 chunks. This file lacks one.");
                        int chunklength = wav.ReadInt32();
                        filesize = wav.ReadInt64();
                        dataChunkLength = wav.ReadInt64();
                        long sampleCount = wav.ReadInt64();
                        wav.ReadBytes(chunklength - 24);
                    }
                    else if (Riffhdr != "RIFF") throw new Exception("Wave files other than Riff or RF64 type not supported.");

                    long end = Math.Min(filesize + 8, wav.BaseStream.Length);
                    string wave = new string(wav.ReadChars(4));

                    //Read WAVE header
                    if (wave != "WAVE") throw new FormatException("WAVE header not present.");

                    List<int> chunklengths = new List<int>();

                    UInt16 FormatTag;
                    Int16 no_of_channels;
                    Int32 BytesPerSecond;
                    Int16 blockAlign;
                    Int16 bitsPerSample;

                    //Read in Format information
                    if (new string(wav.ReadChars(4)) == "fmt ")
                    {
                        UInt32 fmt_len = wav.ReadUInt32();
                        if (fmt_len > Int32.MaxValue || fmt_len < 16) throw new Exception("Invalid format chunk...");
                        FormatTag = wav.ReadUInt16();
                        no_of_channels = wav.ReadInt16();
                        Sample_Frequency = wav.ReadInt32();
                        BytesPerSecond = wav.ReadInt32();
                        blockAlign = wav.ReadInt16();
                        bitsPerSample = wav.ReadInt16();
                        if (fmt_len > 16)
                        {
                            Int16 extraSize = wav.ReadInt16();
                            if (extraSize != fmt_len - 18)
                            {
                                extraSize = (short)(fmt_len - 18);
                            }
                            if (extraSize > 0)
                            {
                                byte[] additionaldata = new byte[extraSize];
                                wav.Read(additionaldata, 0, extraSize);
                            }
                        }
                    }
                    else
                    {
                        throw new Exception("Invalid wave format information...");
                    }

                    int[][] data = new int[no_of_channels][];
                    string str = new string(wav.ReadChars(4));

                    //Sometimes a LIST section is included...
                    if (str == "LIST")
                    {
                        UInt32 lst_len = wav.ReadUInt32();
                        byte[] additionaldata = new byte[lst_len];
                        wav.Read(additionaldata, 0, (int)lst_len);
                        str = new string(wav.ReadChars(4));
                    }

                    //The data chunk should be headed with "data" but occuasionally not. Read in as normal if so. Otherwise, ask the user.
                    if (str != "data" && !try_unsupported_data)
                    {
                        throw new Exception("This wavefile has an unsupported structure, which may lead to misread data. Try anyway?");
                    }
                    long bytelength = wav.BaseStream.Length - wav.BaseStream.Position;
                    UInt32 ch_len = wav.ReadUInt32();
                    for (int c = 0; c < no_of_channels; c++) data[c] = new int[(bytelength / (bitsPerSample / 8)) / no_of_channels];
                    int i = 0;
                    while (wav.BaseStream.Position < end - 8)
                    {
                        i++;
                        for (int c = 0; c < no_of_channels; c++)
                        {
                            byte[] temp = new byte[4];
                            int st = 4 - bitsPerSample / 8;
                            for (int j = st; j < 4; j++) temp[j] = (byte)wav.BaseStream.ReadByte();
                            data[c][i] = BitConverter.ToInt32(temp, 0);
                        }
                    }
                    wav.Close();
                    return data;
                }

                public static bool Write(float[][] Unit_signal, int sample_frequency, string Path, int bitrate = 32)
                {
                    int datalength = (bitrate / 8) * Unit_signal.Length * Unit_signal[0].Length;

                    //Open a stream.
                    System.IO.BinaryWriter wav = new System.IO.BinaryWriter(new System.IO.FileStream(Path, System.IO.FileMode.Create));
                    wav.Write("RIFF".ToCharArray());//4
                    wav.Write((UInt32)datalength + 32);//Filesize//4
                    wav.Write("WAVE".ToCharArray());//4
                    wav.Write("fmt ".ToCharArray());//4
                    wav.Write(bitrate==32?(UInt32)18:(UInt32)16);//Chunk Length//4
                    wav.Write(bitrate==32?(UInt16)3:(UInt16)1);//Format Tag//2
                    wav.Write((Int16)Unit_signal.Length);//Number of channels//2
                    wav.Write((Int32)sample_frequency);//Sample Frequency//4
                    wav.Write((Int32)(sample_frequency * (bitrate / 8) * Unit_signal.Length));// Bytes per second//4
                    wav.Write((Int16)((bitrate / 8) * Unit_signal.Length));//Block Align//2
                    wav.Write((Int16)bitrate);//Bits per sample//2
                    if (bitrate ==32) wav.Write((Int16)0); //cbSize//2
                    //wav.Write("fact");
                    //wav.Write((Int32)4);//ChunkSize//4
                    //wav.Write((Int32)Unit_signal.Length * Unit_signal[0].Length);//No_of_channles * No_of_Samples//4
                    wav.Write("data".ToCharArray());//4
                    wav.Write((UInt32) datalength);//Length of data portion of file.//4
                    bool clipped = false;

                    switch (bitrate)
                    {
                        case 32:
                            for (int i = 0; i < Unit_signal[0].Length; i++)
                            {
                                for (int c = 0; c < Unit_signal.Length; c++)
                                {
                                    clipped = clipped || Math.Abs(Unit_signal[c][i]) >= 1;
                                    wav.Write(Unit_signal[c][i]);
                                }
                            }
                            break;
                        case 24:
                            for (int i = 0; i < Unit_signal[0].Length; i++)
                            {
                                for (int c = 0; c < Unit_signal.Length; c++)
                                {
                                    var value = BitConverter.GetBytes((Int32)(Unit_signal[c][i] * (Math.Pow(2,22)-1)));
                                    byte[] sig_bts = new byte[4];
                                    sig_bts[0] = value[0];
                                    sig_bts[1] = value[1];
                                    sig_bts[2] = value[2];
                                    Int32 v = Math.Abs(BitConverter.ToInt32(new byte[4] { value[0], value[1], value[2], 0 }, 0)) / 2;
                                    clipped = clipped || v > 8388607 * 1.5;
                                    for (int s = 0; s < 3; s++) wav.Write(sig_bts[s]);
                                }
                            }
                            break;
                        case 16:
                            for (int i = 0; i < Unit_signal[0].Length; i++)
                            {
                                for (int c = 0; c < Unit_signal.Length; c++)
                                {
                                    Int16 v = (Int16)(Unit_signal[c][i] * Math.Pow(2, 14));
                                    clipped = clipped || Math.Abs(v) > Int16.MaxValue * 0.9;
                                    wav.Write(v);
                                }
                            }
                            break;
                        default:
                            wav.Close();
                            throw new Exception("invalid bitrate");
                    }
                    wav.Close();

                    return !clipped;
                }
            }
        }
    }
}