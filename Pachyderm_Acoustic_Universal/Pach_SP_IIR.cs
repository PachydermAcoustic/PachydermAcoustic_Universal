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
using System.Numerics;
using System.Linq;

namespace Pachyderm_Acoustic
{
    namespace Audio
    {
        public partial class Pach_SP
        {
            public class IIR_Design
            {
                public static Complex[] AutoCorrelation_Coef(Complex[] X, int maxlag)
                {
                    Complex[] r_l = new Complex[maxlag];

                    for (int lag = 1; lag <= maxlag; lag++)
                    {
                        Complex mean = 0;
                        for (int i = 0; i < X.Length; i++) mean += X[i];
                        mean /= X.Length;
                        Complex denom = 0;
                        Complex num = 0;

                        int N_k = X.Length - lag;

                        for (int i = 0; i < N_k; i++)
                        {
                            Complex x_X = X[i] - mean;
                            denom += x_X * x_X;
                            num += x_X * Complex.Conjugate(X[i + lag] - mean);
                        }

                        for (int i = 0; i < lag; i++)
                        {
                            int idx = N_k + i;
                            Complex x_X = X[idx] - mean;
                            denom += x_X * x_X;
                        }

                        r_l[lag - 1] = num / denom;
                    }

                    return r_l;
                }

                public static Complex[] RLC_FreqResponse(List<double> R, List<double> L, List<double> C, double[] frequencies)
                {
                    Complex[] response = new Complex[frequencies.Length];

                    if (R.Count != L.Count && L.Count != C.Count) throw new Exception("Poles and Zeros must be arranged in pairs...");

                    for (int f = 0; f < frequencies.Length; f++)
                    {
                        Complex w = Complex.ImaginaryOne * Utilities.Numerics.PiX2 * frequencies[f];
                        response[f] = 1;
                        for (int i = 0; i < R.Count; i++) response[f] *= L[i] * w + R[i] + C[i] / w;
                    }
                    return response;
                }

                public static Complex[] PZ_FreqResponse(List<Complex> poles, List<Complex> zeros, double gain, double[] frequencies)
                {
                    Complex[] response = new Complex[frequencies.Length];

                    for (int f = 0; f < frequencies.Length; f++)
                    {
                        Complex s = Complex.ImaginaryOne * Utilities.Numerics.PiX2 * frequencies[f];
                        response[f] = gain;
                        for (int i = 0; i < poles.Count; i++) response[f] /= (s - poles[i]);
                        for (int i = 0; i < zeros.Count; i++) response[f] *= (s - zeros[i]);
                    }

                    return response;
                }

                public static Complex[] AB_FreqResponse(List<Complex> b, List<Complex> a, double[] frequencies)
                {
                    Complex[] response = new Complex[frequencies.Length];

                    for (int f = 0; f < frequencies.Length; f++)
                    {
                        Complex s = Complex.ImaginaryOne * Utilities.Numerics.PiX2 * frequencies[f];
                        Complex num = 0, den = 0;
                        for (int i = 0; i < b.Count(); i++) num += b[i] * Complex.Pow(s, b.Count - i - 1);
                        for (int i = 0; i < a.Count(); i++) den += a[i] * Complex.Pow(s, a.Count - i - 1);
                        response[f] = num / den;
                    }

                    return response;
                }

                public static Complex[] Butter_FreqResponse(double[] frequencies, int order, double freq_l, double freq_u)
                {
                    Complex[] H = new Complex[frequencies.Length];
                    double w_l = Utilities.Numerics.PiX2 * freq_l;
                    double w_u = Utilities.Numerics.PiX2 * freq_u;
                    double max = 0;

                    for(int i = 0; i < frequencies.Length; i++)
                    {
                        double w_f = Utilities.Numerics.PiX2 * frequencies[i];
                        H[i] = 1025 / (Math.Sqrt(1 + Math.Pow(w_f / w_l, 2 * order)) * Math.Sqrt(1 + Math.Pow(w_u / w_f, 2 * order)));
                        max = Math.Max(max, H[i].Real);
                    }

                    return H;
                }

                public static double[] Butter_FreqResponse(int order, double freq_l, double freq_u, double sampling_freq, int length)
                {
                    double[] H = new double[length];
                    double w_l = Utilities.Numerics.PiX2 * freq_l;
                    double w_u = Utilities.Numerics.PiX2 * freq_u;
                    double max = 0;

                    for (int i = 0; i < length; i++)
                    {
                        double w_f = Utilities.Numerics.PiX2 * i * sampling_freq / length;
                        H[i] = 1 / (Math.Sqrt(1 + Math.Pow(w_f / w_l, 2 * order)) * Math.Sqrt(1 + Math.Pow(w_u / w_f, 2 * order)));
                        max = max > H[i]? max : H[i];
                    }
                    for (int i = 0; i < length; i++) H[i] /= max;

                    return H;
                }

                public static Complex[] AB_FreqResponse(List<double> b, List<double> a, double[] frequencies, double fs)
                {
                    int N = frequencies.Length;
                    Complex[] H_f = new Complex[N];

                    for (int i = 0; i < N; i++)
                    {
                        double f = frequencies[i];
                        double omega = 2 * Math.PI * f / fs;
                        Complex z = Complex.Exp(new Complex(0, -omega));

                        // Compute B(z) and A(z)
                        Complex B_z = 0;
                        Complex A_z = 0;

                        for (int k = 0; k < b.Count; k++)
                        {
                            B_z += b[k] * Complex.Pow(z, -k);
                        }

                        for (int k = 0; k < a.Count; k++)
                        {
                            A_z += a[k] * Complex.Pow(z, -k);
                        }

                        // Compute the frequency response H(f) = B(z) / A(z)
                        H_f[i] = B_z / A_z;
                    }

                    return H_f;
                }

                //public class IIR_Fitting_Spectrum_Objective : LibOptimization.Optimization.absObjectiveFunction
                //{
                //    int fs;
                //    int n;
                //    int[] sampled;
                //    Complex[] Spectrum;
                //    Complex[] sel_Spectrum;
                //    Complex[] abs;
                //    double[] freq;
                //    double[] sel_freq;
                //    int p_conj;

                //    public IIR_Fitting_Spectrum_Objective(int order, Complex[] BasisSpectrum, int samplingFrequency)
                //    {
                //        //Find all local maxima (within 10 points)
                //        Spectrum = BasisSpectrum;
                //        fs = samplingFrequency;
                //        n = order;
                //        p_conj = (int)(2 * Math.Floor((double)n / 2));
                //        for (int i = 0; i < freq.Length; i++) freq[i] = (double)(i * samplingFrequency) / freq.Length;
                //    }

                //    public override double F(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //        Complex[] p = new Complex[n];
                //        int p_conj = (int)(2 * Math.Floor((double)n / 2));
                //        for (int i = 0; i < p_conj; i++)
                //        {
                //            p[2 * i] = Complex.FromPolarCoordinates(x[2 * i], x[2 * i + 1] * Math.PI * 2);
                //            p[2 * i + 1] = Complex.Conjugate(p[i]);
                //        }
                //        for (int i = 0; i < n - p_conj; i++)
                //        {
                //            p[i + p_conj] = x[p_conj + i];
                //        }
                //        Complex[] a;
                //        P2poly(p, out a
                //            );
                //        Complex[] Spectrum_Mod = new Complex[Spectrum.Length];
                //        for (int i = 0; i < Spectrum.Length; i++)
                //        {
                //            Complex denom = 0;
                //            Complex s = Complex.ImaginaryOne * freq[i] * Utilities.Numerics.PiX2;
                //            Complex Spow = s;

                //            for (int j = n-1; j >= 0; i++)
                //            {
                //                denom += a[j] * Spow;
                //                Spow *= s;
                //            }
                //            Spectrum_Mod[i] = Spectrum[i] * denom;
                //        }

                //        //TODO: find a polynomial fit operation that works on complex numbers...
                //        //MathNet.Numerics.Fit.Polynomial(freq, Spectrum_Mod, 3);
                //    }

                //    public override List<double> Gradient(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //    }

                //    public override List<List<double>> Hessian(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //    }

                //    public override int NumberOfVariable
                //    {
                //        get
                //        {
                //            return n + p_conj;
                //        }
                //    }
                //}

                //public class IIR_fd_Spectrum_Objective : LibOptimization.Optimization.absObjectiveFunction
                //{
                //    int fs;
                //    int n;
                //    int[] sampled;
                //    Complex[] Spectrum;
                //    Complex[] sel_Spectrum;
                //    Complex[] abs;
                //    double[] freq;
                //    double[] sel_freq;
                //    public double[] pole_freqs;

                //    public IIR_fd_Spectrum_Objective(Complex[] BasisSpectrum, int samplingFrequency)
                //    {
                //        //Find all local maxima (within 10 points)
                //        List<int> maxima = new List<int>();
                //        freq = new double[BasisSpectrum.Length];

                //        fs = samplingFrequency;

                //        for (int i = 0; i < freq.Length; i++) freq[i] = (double)(i * samplingFrequency) / freq.Length;

                //        Spectrum = BasisSpectrum;
                //        n = BasisSpectrum.Length;
                //        fs = samplingFrequency;
                //        sampled = new int[24];
                //        sel_Spectrum = new Complex[24];
                //        sel_freq = new double[24];
                //        abs = new Complex[24];

                //        for (int i = 0; i < 24; i++)
                //        {
                //            double f = 24.803 * Math.Pow(2, (double)i / 3);
                //            sampled[i] = (int)(f * 2 * BasisSpectrum.Length / samplingFrequency);
                //            sel_Spectrum[i] = Spectrum[sampled[i]];
                //            abs[i] = AbsorptionModels.Operations.Absorption_Coef(sel_Spectrum[i]);
                //            sel_freq[i] = freq[sampled[i]];
                //        }

                //        double[] alpha = AbsorptionModels.Operations.Absorption_Coef(BasisSpectrum);
                //        double[] Derivative_Spectrum = new double[BasisSpectrum.Length - 1];
                //        for (int i = 0; i < Derivative_Spectrum.Length; i++) Derivative_Spectrum[i] = alpha[i + 1] - alpha[i];

                //        //for (int i = 10; i < Derivative_Spectrum.Length - 10; i++)
                //        //{
                //        //    bool max = true;
                //        //    for (int s = -10; s < 10; s++)
                //        //    {
                //        //        if (Derivative_Spectrum[i] < Derivative_Spectrum[i + s])
                //        //        {
                //        //            max = false;
                //        //            continue;
                //        //        }
                //        //    }
                //        //    if (max) maxima.Add(i);
                //        //}

                //        //Check for a zero pole...
                //        if (alpha[0] - alpha[3] > 0)
                //        {
                //            maxima.Add(0);
                //        }

                //        //Check for any other poles...
                //        for (int i = 3; i < BasisSpectrum.Length - 3; i++)
                //        {
                //            bool max = true;
                //            for (int s = -3; s < 3; s++)
                //            {
                //                if (alpha[i] > alpha[i + s])
                //                {
                //                    max = false;
                //                    continue;
                //                }
                //            }
                //            if (max) maxima.Add(i);
                //        }

                //        pole_freqs = new double[maxima.Count()];
                //        for (int i = 0; i < pole_freqs.Length; i++) pole_freqs[i] = Utilities.Numerics.PiX2 * (double)(maxima[i] * samplingFrequency) / BasisSpectrum.Length;
                //    }

                //    public override double F(List<double> x)
                //    {
                //        double[] a, b;
                //        DW2AB(x.ToArray(), pole_freqs, out a, out b);
                //        Complex[] spec = AB_FreqResponse(new List<double>(b), new List<double>(a), sel_freq);

                //        Complex[] alpha = new Complex[spec.Length];
                //        double[] alpha_d = new double[spec.Length];

                //        for (int i = 0; i < spec.Length; i++)
                //        {
                //            alpha[i] = AbsorptionModels.Operations.Absorption_Coef(spec[i]);
                //            alpha_d[i] = alpha[i].Real;
                //            if (alpha[i].Real > 1 || alpha[i].Real < 0) return double.PositiveInfinity;
                //        }

                //        //double ans = -Math.Abs(CrossCorrelation_Coef(spec, sel_Spectrum).Real);

                //        double ans = CompareFunctions(abs, alpha, sampled);
                //        //if (ans > -0.9) return 1000;
                //        return ans;
                //    }

                //    public override List<double> Gradient(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //    }

                //    public override List<List<double>> Hessian(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //    }

                //    public override int NumberOfVariable
                //    {
                //        get
                //        {
                //            return pole_freqs.Length;
                //        }
                //    }
                //}

                //public class IIR_Spectrum_Objective:LibOptimization.Optimization.absObjectiveFunction
                //{
                //    int fs;
                //    int n;
                //    int[] sampled;
                //    Complex[] Spectrum;
                //    Complex[] sel_Spectrum;
                //    Complex[] abs;
                //    public int pzct;
                //    int pno_conj;
                //    int pno_conj_phase;
                //    int pno_real;
                //    int zno_real;
                //    double[] freq;
                //    double[] sel_freq;
                //    public double[] pole_phase;

                //    public IIR_Spectrum_Objective(int Filter_Order, Complex[] BasisSpectrum, int samplingFrequency)
                //    {
                //        //Find all local maxima (within 10 points)
                //        List<int> maxima = new List<int>();
                //        freq = new double[BasisSpectrum.Length];
                //        for (int i = 10; i < BasisSpectrum.Length-10; i++)
                //        {
                //            bool max = true;
                //            for (int s = -10; s < 10; s++)
                //            {
                //                if (BasisSpectrum[i].Real < BasisSpectrum[i + s].Real)
                //                {
                //                    max = false;
                //                    continue;
                //                }
                //            }
                //            if (max) maxima.Add(i);
                //        }

                //        fs = samplingFrequency;

                //        for (int i = 0; i < freq.Length; i++) freq[i] = (double)(i * samplingFrequency) / freq.Length;

                //        pole_phase = new double[maxima.Count];
                //        for (int i = 0; i < maxima.Count; i++) pole_phase[i] = 2 * Math.PI * maxima[i] / BasisSpectrum.Count();

                //        if (Filter_Order < maxima.Count * 2) Filter_Order = maxima.Count * 2;

                //        pzct = Filter_Order;
                //        pno_conj = (int)Math.Floor((double)pzct / 2);
                //        pno_conj_phase = (int)Math.Floor((double)pzct / 2 - maxima.Count());
                //        pno_real = pzct % 2;
                //        zno_real = pzct;
                //        Spectrum = BasisSpectrum;
                //        n = BasisSpectrum.Length;
                //        fs = samplingFrequency;
                //        sampled = new int[24];
                //        sel_Spectrum = new Complex[24];
                //        sel_freq = new double[24];
                //        abs = new Complex[24];

                //        for (int i = 0; i < 24; i++)
                //        {
                //            double f = 24.803 * Math.Pow(2, (double)i/3);
                //            sampled[i] = (int)(f * 2 * BasisSpectrum.Length / samplingFrequency);
                //            sel_Spectrum[i] = Spectrum[sampled[i]];
                //            abs[i] = AbsorptionModels.Operations.Absorption_Coef(sel_Spectrum[i]);
                //            sel_freq[i] = freq[sampled[i]];
                //        }
                //    }

                //    public override double F(List<double> x)
                //    {
                //        //foreach (int i in x) if (i > 1 || i < -1) return double.PositiveInfinity;

                //        //Complex[] p = new Complex[(int)Math.Floor((double)x.Count / 4)];
                //        //Complex[] z = new Complex[(int)Math.Floor((double)x.Count / 4)];
                //        Complex[] p = new Complex[2*pno_conj + pno_real];
                //        Complex[] z = new Complex[zno_real];

                //        //for (int i = 0; i < p.Length; i++) p[i] = Complex.FromPolarCoordinates(x[2*i], x[2*i + 1] * Math.PI * 2);
                //        //for (int i = p.Length; i < 2*p.Length; i++) z[i - p.Length] = Complex.FromPolarCoordinates(x[2*i], x[2*i+1] * Math.PI * 2);
                //        //for (int i = 0; i < p.Length; i++)
                //        //{
                //        //    p[i] = Complex.FromPolarCoordinates(x[2 * i], x[2 * i + 1] * Math.PI * 2);
                //        //    z[i] = Complex.Conjugate(p[i]);
                //        //}
                //        for (int i = 0; i < pole_phase.Length; i++)
                //        {
                //            p[2 * i] = Complex.FromPolarCoordinates(1, pole_phase[i]);//x[i]
                //            p[2 * i + 1] = Complex.Conjugate(p[2*i]);
                //        }
                //        for (int i = 0; i < pno_conj_phase; i++)
                //        {
                //            p[2 * i + pole_phase.Length] = Complex.FromPolarCoordinates(x[pole_phase.Length + 2 * i], x[pole_phase.Length + 2 * i + 1] * Math.PI * 2);
                //            p[2 * i + pole_phase.Length + 1] = Complex.Conjugate(p[i]);
                //        }
                //        for (int i = 0; i < pno_real; i++)
                //        {
                //            p[i + pno_conj] = x[pno_conj + pno_conj_phase + i];
                //        }

                //        for (int i = 0; i < zno_real; i++)
                //        {
                //            z[i] = x[pno_conj + pno_conj_phase + pno_real + i];
                //        }

                //        //for (int i = 0; i < pole_phase.Length; i++)
                //        //{
                //        //    p[2 * i] = Complex.FromPolarCoordinates(x[2 * i], pole_phase[i]);
                //        //    p[2 * i + 1] = Complex.Conjugate(p[i]);
                //        //}
                //        //for (int i = 0; i < pno_conj_phase; i++)
                //        //{
                //        //    p[2 * i + pole_phase.Length] = Complex.FromPolarCoordinates(x[2 * i], x[2 * i + 1] * Math.PI * 2);
                //        //    p[2 * i + pole_phase.Length + 1] = Complex.Conjugate(p[i]);
                //        //}
                //        //for (int i = 0; i < pno_real; i++)
                //        //{
                //        //    p[i + 2 * pno_conj] = x[2 * pno_conj + i];
                //        //}

                //        //for (int i = 0; i < zno_real; i++)
                //        //{
                //        //    z[i] = x[2 * pno_conj + pno_real + i];
                //        //}

                //        Complex[] a, b;
                //        PZ2AB(p, z, 1, out a, out b);
                //        //Complex[] spec = SpectrumFromIIR(a, b, fs, n);
                //        Complex[] spec = AB_FreqResponse(new List<Complex>(b), new List<Complex>(a), sel_freq);

                //        Complex[] alpha = new Complex[spec.Length];
                //        double[] alpha_d = new double[spec.Length];

                //        for(int i = 0; i < spec.Length; i++)
                //        { 
                //            alpha[i] = AbsorptionModels.Operations.Absorption_Coef(spec[i]);
                //            alpha_d[i] = alpha[i].Real;
                //            if (alpha[i].Real > 1 || alpha[i].Real < 0) return double.PositiveInfinity;
                //        }

                //        //double ans = -Math.Abs(CrossCorrelation_Coef(spec, sel_Spectrum).Real);

                //        double ans = CompareFunctions(abs, alpha, sampled) - 1;
                //        //if (ans > -0.9) return 1000;
                //        return ans;
                //    }

                //    public override List<double> Gradient(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //    }

                //    public override List<List<double>> Hessian(List<double> x)
                //    {
                //        throw new NotImplementedException();
                //    }

                //    public override int NumberOfVariable
                //    {
                //        get
                //        {
                //            //return 2*pno_conj + pno_real + zno_real;
                //            return pno_conj_phase + pno_conj + pno_real + zno_real;
                //        }
                //    }
                //}

                ////public static void PZtoAB(Complex[] p, Complex[] z, out double[] a, out double[] b)
                ////{
                ////    Complex[] a1, b1;
                ////    PZtoAB(p, z, out a1, out b1);
                ////    a = new double[a1.Length];
                ////    b = new double[b1.Length];
                ////    for (int i = 0; i < a.Length; i++)
                ////    {
                ////        a[i] = a1[i].Magnitude;
                ////        b[i] = b1[i].Magnitude;
                ////    }
                ////}

                ////public static void ABtpPZ(out Complex[] a, out Complex[] b, Complex[] p, Complex[] z)
                ////{

                ////}

                //public static void DW2AB(double[] d, double[] w, out double[] a, out double[] b)
                //{
                //    if (d.Length != w.Length) throw new Exception("Angular frequency and damping coefficients must be organized in pairs...");

                //    b = new double[2] { 2 * d[0] * w[0], 0 };
                //    a = new double[3] { 1, 2 * d[0] * w[0], w[0] * w[0] };

                //    for (int i = 1; i < d.Length; i++)
                //    {
                //        double coef = 2 * d[i] * w[i];
                //        double[] btemp = new double[2] { coef, 0 };
                //        double[] atemp = new double[3] { 1, coef, w[0] * w[0] };

                //        b = Polynomial_Add(Polynomial_Multiply(b, atemp), Polynomial_Multiply(btemp, a));
                //        a = Polynomial_Multiply(a, atemp);
                //    }
                //}

                //public static void DW2PZ(double[] d, double[] w, out Complex[] p, out double[] z, out double[] g)
                //{
                //    if (d.Length != w.Length) throw new Exception("Angular frequency and damping coefficients must be organized in pairs...");
                //    p = new Complex[w.Length * 4];
                //    z = new double[w.Length];
                //    g = new double[w.Length];

                //    for (int i = 0; i < d.Length; i++)
                //    {
                //        double _dw = -d[i] * w[i];
                //        double rtwd2_1 = Math.Sqrt(_dw * _dw - w[i]);
                //        p[2 * i] = _dw + rtwd2_1;
                //        p[2 * i + 1] = _dw - rtwd2_1;
                //        z[i] = 0;
                //        g[i] = -2 * _dw;
                //    }
                //}

                //public static void P2poly(Complex[] p, out Complex[] poly)
                //{
                //    poly = new Complex[p.Length + 1];

                //    poly[0] = 1;

                //    if (p.Length == 1) poly[1] = -p[0];
                //    else
                //    {
                //        poly[1] = -2 * p[0] * p[1];
                //        poly[2] = p[0] * p[1];
                //        for (int n = 2; n < p.Length; n++)
                //        {
                //            for (int m = 1; m < poly.Length; m++) poly[m] += p[n] * -poly[m - 1];
                //        }
                //    }
                //}


                //public static void PZ2AB(Complex[] p, Complex[] z, double k, out Complex[] a, out Complex[] b)
                //{
                //    a = new Complex[p.Length + 1];
                //    b = new Complex[z.Length + 1];

                //    a[0] = 1;
                //    b[0] = k;

                //    if (p.Length == 1)
                //    {
                //        a[1] = -p[0];
                //        b[1] = -z[0];
                //    }
                //    else
                //    {
                //        a[1] = -2 * p[0] * p[1];
                //        b[1] = -2 * z[0] * z[1];
                //        a[2] = p[0] * p[1];
                //        b[2] = z[0] * z[1];

                //        for (int n = 2; n < p.Length; n++)
                //        {
                //            for (int m = 1; m < a.Length; m++)
                //            {
                //                a[m] += p[n] * -a[m - 1];
                //                b[m] += z[n] * -b[m - 1];
                //            }
                //        }
                //    }
                //}

                //private static double[] Polynomial_Add(double[] A, double[] B)
                //{
                //    double[] sum = new double[Math.Max(A.Length, B.Length)];
                //    int min = Math.Min(A.Length, B.Length) - 1;
                //    for (int i = min; i >= 0; i--) sum[i] = A[i] + B[i];
                //    return sum;
                //}

                //private static Complex[] Polynomial_Add(Complex[] A, Complex[] B)
                //{
                //    Complex[] sum = new Complex[Math.Max(A.Length, B.Length)];
                //    int min = Math.Min(A.Length, B.Length)-1;
                //    for (int i = min; i >= 0; i--) sum[i] = A[i] + B[i];
                //    return sum;
                //}

                //private static double[] Polynomial_Multiply(double[] A, double[] B)
                //{
                //    double[] prod = new double[A.Length + B.Length - 1];
                //    // Multiply two polynomials term by term
                //    for (int i = 0; i < A.Length; i++)
                //    {
                //        for (int j = 0; j < B.Length; j++)
                //            prod[i + j] += A[i] * B[j];
                //    }

                //    return prod;
                //}

                //private static Complex[] Polynomial_Multiply(Complex[] A, Complex[] B)
                //{
                //    Complex[] prod = new Complex[A.Length + B.Length - 1];
                //    // Multiply two polynomials term by term
                //    for (int i = 0; i < A.Length; i++)
                //    {
                //        for (int j = 0; j < B.Length; j++)
                //            prod[i + j] += A[i] * B[j];
                //    }

                //    return prod;
                //}

                //public static void OptimizeIIR(Complex[] Spectrum, int Sample_Freq, int filter_order, out double[] a, out double[] b)
                //{
                //    IIR_Spectrum_Objective obj = new IIR_Spectrum_Objective(filter_order, Spectrum, Sample_Freq);
                //    //LibOptimization.Optimization.clsOptRealGAREX SA = new LibOptimization.Optimization.clsOptRealGAREX(obj);
                //    //LibOptimization.Optimization.clsOptRealGABLX SA = new LibOptimization.Optimization.clsOptRealGABLX(obj);
                //    LibOptimization.Optimization.clsOptSimulatedAnnealing SA = new LibOptimization.Optimization.clsOptSimulatedAnnealing(obj);
                //    //LibOptimization.Optimization.clsOptNelderMead SA = new LibOptimization.Optimization.clsOptNelderMead(obj);
                //    //SA.UseEliteStrategy(10);
                //    //SA.PARAM_PopulationSize = Filter_Order * 50;
                //    SA.InitialValueRange = 0.00001;
                //    //SA.PARAM_InitRange = 1;
                //    //SA.PARAM_ChildrenSize = 20;
                //    //SA.PARAM_MAX_ITERATION = 1000;
                //    SA.Init();
                //    SA.DoIteration();

                //    DW2AB( SA.Result.ToArray(), null, out a, out b);
                //}

                //public static void OptimizeIIR(Complex[] Spectrum, int Filter_Order, int Sample_Freq, out Complex[] a, out Complex[] b)
                //{
                //    IIR_Spectrum_Objective obj = new IIR_Spectrum_Objective(Filter_Order, Spectrum, Sample_Freq);
                //    //LibOptimization.Optimization.clsOptRealGAREX SA = new LibOptimization.Optimization.clsOptRealGAREX(obj);
                //    LibOptimization.Optimization.clsOptRealGABLX SA = new LibOptimization.Optimization.clsOptRealGABLX(obj);
                //    //SA.UseEliteStrategy(10);
                //    //SA.PARAM_PopulationSize = Filter_Order * 50;

                //    SA.InitialValueRange = 1;
                //    SA.ChildrenSize = 20;
                //    //SA.PARAM_MAX_ITERATION = 1000;
                //    SA.Init();
                //    SA.DoIteration();

                //    Filter_Order = obj.pzct;

                //    //Complex[] p = new Complex[(int)Math.Floor((double)SA.Result.Count / 4)];
                //    //Complex[] z = new Complex[(int)Math.Ceiling((double)SA.Result.Count / 4)];
                //    Complex[] p = new Complex[Filter_Order];
                //    Complex[] z = new Complex[Filter_Order];

                //    int pzct = Filter_Order;
                //    int pno_conj = (int)Math.Floor((double)pzct / 2);
                //    int pno_conj_phase = pno_conj - obj.pole_phase.Length;
                //    int pno_real = pzct % 2;
                //    int zno_real = pzct;

                //    for (int i = 0; i < obj.pole_phase.Length; i++)
                //    {
                //        p[2 * i] = Complex.FromPolarCoordinates(SA.Result[i], obj.pole_phase[i]);
                //        p[2 * i + 1] = Complex.Conjugate(p[2*i]);
                //    }
                //    for (int i = 0; i < pno_conj_phase; i++)
                //    {
                //        p[2 * i + obj.pole_phase.Length] = Complex.FromPolarCoordinates(SA.Result[obj.pole_phase.Length + 2 * i], SA.Result[obj.pole_phase.Length + 2 * i + 1] * Math.PI * 2);
                //        p[2 * i + obj.pole_phase.Length + 1] = Complex.Conjugate(p[i]);
                //    }
                //    for (int i = 0; i < pno_real; i++)
                //    {
                //        p[i + pno_conj] = SA.Result[pno_conj + pno_conj_phase + i];
                //    }

                //    for (int i = 0; i < zno_real; i++)
                //    {
                //        z[i] = SA.Result[pno_conj + pno_conj_phase + pno_real + i];
                //    }

                //    //for (int i = 0; i < p.Length; i++) p[i] = Complex.FromPolarCoordinates(SA.Result[2 * i], SA.Result[2 * i + 1] * Math.PI * 2);
                //    //for (int i = p.Length; i < 2 * p.Length; i++) z[i-p.Length] = Complex.FromPolarCoordinates(SA.Result[2 * i], SA.Result[2 * i + 1] * Math.PI * 2);
                //    PZ2AB(p, z, 1, out a, out b);

                //    //Normalize
                //    //Complex[] spec = SpectrumFromIIR(a, b, Sample_Freq, Spectrum.Length);
                //    //double ms = 0;
                //    //foreach (Complex s in spec) ms = ((ms > s.Magnitude) ? ms : s.Magnitude);
                //    //double mb = 0;
                //    //foreach (Complex s in Spectrum) mb = ((mb > s.Magnitude) ? mb : s.Magnitude);

                //    //b[0] *= mb / ms;
                //}

                ////                public static void invfreq(Complex[] Freq_Response, double[] Frequencies, int nB, int nA, out double[] b, out double[] a)
                ////                {
                ////                    //int Zb = 0;
                ////                    int n = Math.Max(nA, nB);
                ////                    int m = n + 1, mA = nA + 1, mB = nB + 1;
                ////                    int nF = Frequencies.Length;
                ////                    if (Frequencies.Length != Freq_Response.Length) throw new Exception("Length of Freq_Response and Frequenices must be the same");
                ////                    //if nargin< 5 || isempty(W), W = ones(1, nF); endif
                ////                    //if nargin< 6, iter = []; endif
                ////                    //if nargin< 7  tol = []; endif
                ////                    //if nargin< 8 || isempty(tr), tr = ''; endif
                ////                    //if nargin< 9, plane = 'z'; endif
                ////                    //if nargin< 10, varargin = {}; endif
                ////                    //if iter ~=[], disp('no implementation for iter yet'),endif
                ////                    //if tol ~=[], disp('no implementation for tol yet'),endif
                ////                    //if (plane ~= 'z' && plane ~= 's'), disp('invfreqz: Error in plane argument'), endif

                ////                    //[reg, prop] = parseparams(varargin);
                ////                    //## should we normalise freqs to avoid matrices with rank deficiency ?
                ////                    //bool norm = false;
                ////                    //## by default, use Ordinary Least Square to solve normal equations
                ////                    //string method = 'LS';

                ////                    MathNet.Numerics.LinearAlgebra.Matrix<Complex> Ruu = MathNet.Numerics.LinearAlgebra.Complex.Matrix.Build.Dense(mB, mB);
                ////                    MathNet.Numerics.LinearAlgebra.Matrix<Complex> Ryy = MathNet.Numerics.LinearAlgebra.Complex.Matrix.Build.Dense(nA, nA);
                ////                    MathNet.Numerics.LinearAlgebra.Matrix<Complex> Ryu = MathNet.Numerics.LinearAlgebra.Complex.Matrix.Build.Dense(nA, mB);
                ////                    MathNet.Numerics.LinearAlgebra.Vector<Complex> Pu = MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.Dense(mB);
                ////                    MathNet.Numerics.LinearAlgebra.Vector<Complex> Py = MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.Dense(nA);

                ////                    Complex[] s = new Complex[Freq_Response.Length];
                ////                    double[] F = new double[Freq_Response.Length];

                ////                    for (int i = 0; i < Freq_Response.Length; i++)
                ////                    {
                ////                        F[i] = (double)i * Math.PI / Freq_Response.Length;
                ////                        s[i] = Complex.Exp(-F[i] * Complex.ImaginaryOne);
                ////                    }

                ////                    for (int k = 0; k < nF; k++)
                ////                    {
                ////                        Complex[] Zk = new Complex[n];
                ////                        for (int i = 0; i < n; i++)
                ////                        {
                ////                            Zk[i] = Complex.Pow(Frequencies[k] * Complex.ImaginaryOne, i);
                ////                        }
                ////                        Complex aHks = Freq_Response[k] * Complex.Conjugate(Freq_Response[k]);
                ////                        MathNet.Numerics.LinearAlgebra.Vector<Complex> Zkv = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.DenseOfArray(Zk);

                ////                        MathNet.Numerics.LinearAlgebra.Matrix<Complex> Rk = Zkv.ToColumnMatrix() * Zkv.ToRowMatrix();
                ////                        MathNet.Numerics.LinearAlgebra.Matrix<Complex> rRk = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.DenseOfMatrix(Rk);
                ////                        for (int x = 0; x < Rk.ColumnCount; x++) for (int y = 0; y < Rk.RowCount; y++) rRk[x, y] = rRk[x, y].Real;

                ////                        Ruu = Ruu + rRk.SubMatrix(0, mB, 0, mB);
                ////                        Ryy += aHks * rRk.SubMatrix(2, mA, 2, mA);
                ////                        Ryu += (Freq_Response[k] * Rk.SubMatrix(2, mA, 1, mB));
                ////                        for (int i = 0; i < mB; i++) Pu[i] += (Complex.Conjugate(Freq_Response[k]) * Zk[i]).Real;
                ////                        for (int i = 1; i < mA; i++) Py[i] += aHks * Zk[i].Real;
                ////                    }

                ////                    MathNet.Numerics.LinearAlgebra.Matrix<Complex> Rr = MathNet.Numerics.LinearAlgebra.Complex.Matrix.Build.Dense(s.Length, mB+nA, 1);
                ////                    MathNet.Numerics.LinearAlgebra.Vector<Complex> Zk = MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.DenseOfArray(s);

                ////                    for (int k = 1; k < Math.Min(nA, nB); k++)
                ////                    {
                ////                        Rr.SetSubMatrix(0, 1 + k, Zk.ToRowMatrix());
                ////                        Rr.SetSubMatrix(0, mB + k, (-Zk.PointwiseMultiply(MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.DenseOfArray(Freq_Response)).ToColumnMatrix()));
                ////                        Zk = Zk.PointwiseMultiply(MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.DenseOfArray(s));
                ////                    }

                ////                    int k;
                ////                    for (k = 1 + Math.Min(nA, nB); k < Math.Max(nA, nB) - 1; k++)
                ////                    {
                ////                        if (k <= nB) Rr.SetSubMatrix(0, Rr.RowCount, 1 + k, 1, Zk.ToColumnMatrix());
                ////                        if(k <= nA) Rr.SetSubMatrix(0, Rr.RowCount, mB + k, 1, -Zk.PointwiseMultiply(MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.DenseOfArray(Freq_Response)).ToColumnMatrix());
                ////                        Zk = Zk.PointwiseMultiply(MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.DenseOfArray(s));
                ////                    }
                ////                    k++;

                ////                    if(k <= nB) Rr.SetSubMatrix(0, 1 + k, Zk.ToColumnMatrix());
                ////                    if(k <= nA) Rr.SetSubMatrix(0, mB + k, (-Zk.PointwiseMultiply(MathNet.Numerics.LinearAlgebra.Complex.Vector.Build.DenseOfArray(Freq_Response)).ToColumnMatrix()));
                ////                    //Rr(:, mB + k) = -Zk.* H;

                //////## complex to real equation system -- this ensures real solution
                ////                    Rr = Rr.SubMatrix(0, Rr.RowCount, 1, Rr.ColumnCount);

                ////                    Rr = [real(Rr); imag(Rr)]; Pr = [real(H(:)); imag(H(:))];
                //////## normal equations -- keep for ref
                //////## Rn= [Ruu(1+zB:mB, 1+zB:mB), -Ryu(:, 1+zB:mB)';  -Ryu(:, 1+zB:mB), Ryy];
                //////## Pn= [Pu(1+zB:mB); -Py];

                //////## avoid scaling errors with Theta = R\P;
                //////## [Q, R] = qr([Rn Pn]); Theta = R(1:end, 1:end-1)\R(1:end, end);
                ////[Q, R] = qr([Rr Pr], 0); Theta = R(1:end-1, 1:end-1)\R(1:end-1, end);
                //////## SigN = R(end, end-1);
                ////                //SigN = R(end, end);

                ////                B = [zeros(0, 1); Theta(1:mB)];
                ////                A = [1; Theta(mB + (1:nA))];
                ////            }

                public static double CompareFunctions(Complex[] f1, Complex[] f2, int[] sampled)
                {
                    if (f1.Length != f2.Length) throw new Exception("f1 and f2 must be equal in length");

                    double Corr = 0;

                    //foreach (int f in sampled)
                    for (int i = 0; i < f1.Length; i++)
                    {
                        //if (f > f1.Length) continue;
                        //Corr += Math.Abs((f1[f].Magnitude - f2[f].Magnitude))/Math.Max(f1[f].Magnitude, f2[f].Magnitude);
                        Corr += Math.Pow(f1[i].Real - f2[i].Real, 2) + Math.Pow(f1[i].Imaginary - f2[i].Imaginary, 2);
                        //Corr += Math.Abs((f1[i].Real - f2[i].Real)) / Math.Max(f1[i].Real, f2[i].Real);
                        //Corr = Math.Max(Corr, Math.Abs((f1[f].Imaginary - f2[f].Imaginary)) / Math.Max(f1[f].Imaginary, f2[f].Imaginary));
                    }

                    //Corr /= sampled.Length;
                    return Corr;
                }

                public static Complex CrossCorrelation_Coef(Complex[] f1, Complex[] f2, int[] sampled)
                {
                    //if (f1.Length != f2.Length) throw new Exception("f1 and f2 must be equal in length");
                    int min = Math.Min(f1.Length, f2.Length);
                    int lost = 0;

                    Complex mf1 = 0, mf2 = 0;
                    foreach (int i in sampled)
                    {
                        if (i > min)
                        {
                            lost++;
                            continue;
                        }
                        mf1 += f1[i];
                        mf2 += f2[i];
                    }

                    mf1 /= f1.Length - lost;
                    mf2 /= f2.Length - lost;

                    Complex Corr = 0;
                    Complex sx = 0, sy = 0;
                    foreach (int f in sampled)
                    {
                        if (f > min) continue;
                        Complex x1 = f1[f] - mf1;
                        Complex x2 = f2[f] - mf2;
                        Corr += Complex.Conjugate(x1) * x2;
                        sx += f1[f];
                        sy += f2[f];
                    }

                    Corr /= Complex.Sqrt(sx * sy);
                    return Corr;
                }

                public static Complex CrossCorrelation_Coef(Complex[] f1, Complex[] f2)
                {
                    if (f1.Length != f2.Length) throw new Exception("f1 and f2 must be equal in length");
                    Complex mf1 = 0, mf2 = 0;
                    for (int i = 0; i < f1.Length; i++)
                    {
                        mf1 += f1[i];
                        mf2 += f2[i];
                    }

                    mf1 /= f1.Length;
                    mf2 /= f2.Length;

                    Complex Corr = 0;
                    Complex sx = 0, sy = 0;
                    for (int f = 0; f < f1.Length; f++)
                    {
                        Complex x1 = f1[f] - mf1;
                        Complex x2 = f2[f] - mf2;
                        Corr += Complex.Conjugate(x1) * x2;
                        sx += f1[f];
                        sy += f2[f];
                    }

                    Corr /= Complex.Sqrt(sx * sy);
                    return Corr;
                }

                ///// <summary>
                ///// Designs a band-pass Butterworth filter according to the specification.
                ///// </summary>
                ///// <param name="lowStopbandFreq">Lower stopband corner frequency, normalized to Nyquist frequency.</param>
                ///// <param name="lowPassbandFreq">Lower passband corner frequency, normalized to Nyquist frequency.</param>
                ///// <param name="highPassbandFreq">Higher passband corner frequency, normalized to Nyquist frequency.</param>
                ///// <param name="highStopbandFreq">Higher stopband corner frequency, normalized to Nyquist frequency.</param>
                ///// <param name="passbandRipple">Maximum allowed passband ripple.</param>
                ///// <param name="stopbandAttenuation">Minimum required stopband attenuation.</param>
                ///// <returns>Minimum required filter order and computed cutoff frequency.</returns>
                ///// <exception cref="ArgumentException">Lower passband corner frequency is not in [0,1] range.</exception>
                ///// <exception cref="ArgumentException">Lower stopband corner frequency is not in [0,1] range.</exception>
                ///// <exception cref="ArgumentException">Higher passband corner frequency is not in [0,1] range.</exception>
                ///// <exception cref="ArgumentException">Higher stopband corner frequency is not in [0,1] range.</exception>
                ///// <exception cref="ArgumentException">Passband ripple must be a finite number.</exception>
                ///// <exception cref="ArgumentException">Stopband attenuation must be a finite number.</exception>
                ///// <exception cref="ArgumentException">Lower stopband corner frequency must be lesser than lower passband corner frequency.</exception>
                ///// <exception cref="ArgumentException">Lower passband corner frequency must be lesser than higher passband corner frequency.</exception>
                ///// <exception cref="ArgumentException">Higher stopband corner frequency must be greater than higher passband corner frequency.</exception>
                ///// <exception cref="ArgumentException">Stopband attenuation must be greater than passband ripple.</exception>
                //public static (byte n, double wc1, double wc2) BandPass(double lowStopbandFreq, double lowPassbandFreq, double highPassbandFreq, double highStopbandFreq, double passbandRipple, double stopbandAttenuation)
                //{
                //    Helpers.Validators.CheckFrequency(lowStopbandFreq, nameof(lowStopbandFreq));
                //    Helpers.Validators.CheckFrequency(lowPassbandFreq, nameof(lowPassbandFreq));
                //    Helpers.Validators.CheckFrequency(highPassbandFreq, nameof(highPassbandFreq));
                //    Helpers.Validators.CheckFrequency(highPassbandFreq, nameof(highPassbandFreq));
                //    Helpers.Validators.CheckDouble(passbandRipple, nameof(passbandRipple));
                //    Helpers.Validators.CheckDouble(stopbandAttenuation, nameof(stopbandAttenuation));

                //    if (lowStopbandFreq > lowPassbandFreq)
                //    {
                //        throw new ArgumentException("Lower stopband corner frequency must be lesser than lower passband corner frequency.", nameof(lowStopbandFreq));
                //    }

                //    if (lowPassbandFreq > highPassbandFreq)
                //    {
                //        throw new ArgumentException("Lower passband corner frequency must be lesser than higher passband corner frequency.", nameof(highPassbandFreq));
                //    }

                //    if (highStopbandFreq < highPassbandFreq)
                //    {
                //        throw new ArgumentException("Higher stopband corner frequency must be greater than higher passband corner frequency.", nameof(highStopbandFreq));
                //    }

                //    if (stopbandAttenuation < passbandRipple)
                //    {
                //        throw new ArgumentException("Stopband attenuation must be greater than passband ripple.", nameof(stopbandAttenuation));
                //    }

                //    var wwp1 = Math.Tan(Math.PI * lowPassbandFreq / 2);
                //    var wwp2 = Math.Tan(Math.PI * highPassbandFreq / 2);
                //    var wws1 = Math.Tan(Math.PI * lowStopbandFreq / 2);
                //    var wws2 = Math.Tan(Math.PI * highStopbandFreq / 2);

                //    var w02 = wwp1 * wwp2;
                //    if (w02 < wws1 * wws2)
                //    {
                //        wws2 = w02 / wws1;
                //    }
                //    else
                //    {
                //        wws1 = w02 / wws2;
                //    }

                //    const double wwp = 1d;
                //    var wws = (wws2 - wws1) / (wwp2 - wwp1);

                //    var qp = Math.Log(Math.Pow(10, passbandRipple / 10) - 1);
                //    var qs = Math.Log(Math.Pow(10, stopbandAttenuation / 10) - 1);
                //    var n = (byte)Math.Ceiling((qs - qp) / (2 * (Math.Log(wws) - Math.Log(wwp))));

                //    var wpp1 = Math.Exp(Math.Log(wwp1) - (qp / 2 / n));

                //    var (wb, wa) = CutoffFrequencies(wwp1, wwp2, wpp1);

                //    return (n, Math.Atan(wb) * 2 / Math.PI, Math.Atan(wa) * 2 / Math.PI);
                //}

                //public static (double[] numerator, double[] denominator) BandPass(double lowStopbandFreq, double lowPassbandFreq, double highPassbandFreq, double highStopbandFreq, double passbandRipple, double stopbandAttenuation)
                //{
                //    var (n, wc1, wc2) = BandPass(lowStopbandFreq, lowPassbandFreq, highPassbandFreq, highStopbandFreq, passbandRipple, stopbandAttenuation);

                //    const double T = 2;
                //    var (gain, zeros, poles) = TransferFunction(n);

                //    wc1 = WarpFrequency(wc1, T);
                //    wc2 = WarpFrequency(wc2, T);
                //    (gain, zeros, poles) = TFBandPass(gain, zeros, poles, wc1, wc2);

                //    return Coefficients(gain, zeros, poles, T);
                //}

                //internal static (double gain, Complex[] zeros, Complex[] poles) TFBandPass(double gain, Complex[] zeros, Complex[] poles, double wc1, double wc2)
                //{
                //    var z = zeros.Length;
                //    var p = poles.Length;

                //    gain *= Math.Pow(1 / (wc2 - wc1), z - p);

                //    { // scopes b
                //        var b = poles.Select(pole => pole * ((wc2 - wc1) / 2)).ToArray();

                //        poles = new Complex[2 * p];

                //        for (int i = p - 1; i >= 0; i--)
                //        {
                //            var sqrt = Complex.Sqrt((b[i] * b[i]) - (wc2 * wc1));
                //            poles[i] = b[i] + sqrt;
                //            poles[i + p] = b[i] - sqrt;
                //        }
                //    }

                //    if (z == 0)
                //    {
                //        zeros = Enumerable.Repeat(Complex.Zero, p).ToArray();
                //    }
                //    else
                //    {
                //        var a = zeros.Select(zero => zero * ((wc2 - wc1) / 2)).ToArray();

                //        zeros = new Complex[2 * z];

                //        for (int i = z - 1; i >= 0; i--)
                //        {
                //            var sqrt = Complex.Sqrt((a[i] * a[i]) - (wc2 * wc1));
                //            zeros[i] = a[i] + sqrt;
                //            zeros[i + z] = a[i] - sqrt;
                //        }

                //        if (poles.Length > zeros.Length)
                //        {
                //            var tmp = zeros;
                //            zeros = new Complex[poles.Length];
                //            Array.Copy(tmp, zeros, zeros.Length);
                //        }
                //    }

                //    return (gain, zeros, poles);
                //}

                //public static (double gain, Complex[] zeros, Complex[] poles) Bilinear(double gain, Complex[] zeros, Complex[] poles, double samplingTime)
                //{
                //    //Helpers.Validators.CheckDouble(gain, nameof(gain));
                //    //Helpers.Validators.CheckDouble(samplingTime, nameof(samplingTime));
                //    //Helpers.Validators.CheckNull(zeros, nameof(zeros));
                //    //Helpers.Validators.CheckNull(poles, nameof(poles));

                //    var z = zeros.Length;
                //    var p = poles.Length;

                //    if (z > p || p == 0)
                //    {
                //        throw new ArgumentException("The number of poles must at least be equal to the number of zeros.", nameof(poles));
                //    }

                //    var nProd = zeros.Aggregate(Complex.One, (acc, zero) => acc *= (2 - (zero * samplingTime)) / samplingTime);
                //    var dProd = poles.Aggregate(Complex.One, (acc, pole) => acc *= (2 - (pole * samplingTime)) / samplingTime);

                //    gain = (gain * nProd / dProd).Real;

                //    for (int i = p - 1; i >= 0; i--)
                //    {
                //        poles[i] = (2 + (poles[i] * samplingTime)) / (samplingTime - (poles[i] * samplingTime));
                //    }

                //    if (z == 0)
                //    {
                //        zeros = Enumerable.Repeat(-Complex.One, p).ToArray();
                //    }
                //    else
                //    {
                //        var tmp = zeros;
                //        zeros = new Complex[p];

                //        for (int i = z - 1; i >= 0; i--)
                //        {
                //            zeros[i] = (2 + (tmp[i] * samplingTime)) / (samplingTime - (tmp[i] * samplingTime));
                //        }

                //        Array.Copy(Enumerable.Repeat(-Complex.One, p - z).ToArray(), 0, zeros, z, p - z);
                //    }

                //    return (gain, zeros, poles);
                //}
                //private static (double[] numerator, double[] denominator) Coefficients(double gain, Complex[] zeros, Complex[] poles, double T)
                //{
                //    (gain, zeros, poles) = Bilinear(gain, zeros, poles, T);

                //    double[] numerator = Generate.Map(PolynomialCoefficients(zeros), num => (num * gain).Real);
                //    double[] denominator = Generate.Map(PolynomialCoefficients(poles), den => den.Real);

                //    return (numerator, denominator);
                //}

                //private static (double gain, Complex[] zeros, Complex[] poles) TransferFunction(uint n)
                //{
                //    var zeros = new Complex[0];
                //    var poles = new Complex[n];
                //    for (int i = 0; i < n; i++)
                //    {
                //        poles[i] = Complex.Exp(Complex.ImaginaryOne * Math.PI * ((2 * (i + 1)) + n - 1) / (2 * n));
                //    }

                //    if ((n & 1) == 1) // if n is odd
                //    {
                //        var target = (n + 1) / 2;
                //        poles[target - 1] = -1;
                //    }

                //    return (1, zeros, poles);
                //}

                ///// <summary>
                ///// Computes the cutoff frequencies for a band-pass or a band-stop filter.
                ///// </summary>
                //private static (double wb, double wa) CutoffFrequencies(double wwp1, double wwp2, double wpp1)
                //{
                //    var w0 = Math.Sqrt(wwp1 * wwp2);
                //    var Q = w0 / (wwp2 - wwp1);
                //    var wc1 = wwp1;

                //    var wp = wpp1 / wc1;
                //    var wa = Math.Abs(wp + Math.Sqrt((wp * wp) + (4 * Q * Q))) / (2 * Q / w0);
                //    var wb = Math.Abs(wp - Math.Sqrt((wp * wp) + (4 * Q * Q))) / (2 * Q / w0);

                //    return (wb, wa);
                //}

                //internal static double WarpFrequency(double cutoffFrequency, double samplingTime)
                //{
                //    return Math.Tan(Math.PI * cutoffFrequency / samplingTime);
                //}

                ///// <summary>
                ///// Computes the coefficients of the polynomial whose solutions are roots are given as input parameter.
                ///// </summary>
                ///// <param name="roots">Roots of the polynomial.</param>
                ///// <returns>Polynomial coefficients.</returns>
                //internal static Complex[] PolynomialCoefficients(Complex[] roots)
                //{
                //    Complex[] y = Generate.Repeat(roots.Length + 1, Complex.Zero);
                //    y[0] = Complex.One;

                //    for (int i = 0; i < roots.Length; i++)
                //    {
                //        for (int j = i; j >= 0; j--)
                //        {
                //            y[j + 1] -= roots[i] * y[j];
                //        }
                //    }

                //    //    return y;
                //    //}

                //    public static Butterworth_Coefficients(double f_u, double f_l, out double[] A, out double[] B)
                //    {

                //    }


                //    //Estimate the coeffients of a band-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
                //    public double[][] Lp2bp(double W_f1, double W_f2, int order_filt)
                //    {
                //        double[][] save_filt_coeff = new double[2][];
                //        for (int kk = 0; kk < 2; kk++) save_filt_coeff[kk] = new double[2 * order_filt + 1];
                //        int fs = 2;
                //        int type_filt = 0;
                //        double u_f1 = 0, u_f2 = 0;

                //        //Step 1: get analog, pre - warped frequencies
                //        //Freq_pre_wrapped(type_filt, W_f1, W_f2);
                //        switch (type_filt)
                //        {
                //            //Band-pass
                //            case 0:
                //                u_f1 = 2 * fs * Math.Tan(Math.PI * W_f1 / fs);
                //                u_f2 = 2 * fs * Math.Tan(Math.PI * W_f2 / fs);
                //                break;
                //            //Band-stop
                //            case 1:
                //                u_f1 = 2 * fs * Math.Tan(Math.PI * W_f1 / fs);
                //                u_f2 = 2 * fs * Math.Tan(Math.PI * W_f2 / fs);
                //                break;
                //            //Low-pass
                //            case 2:
                //                u_f2 = 2 * fs * Math.Tan(Math.PI * W_f2 / fs);
                //                break;
                //            //High-pass
                //            case 3:
                //                u_f1 = 2 * fs * Math.Tan(Math.PI * W_f1 / fs);
                //                break;
                //        }

                //        //Step 2: convert to low-pass prototype estimate
                //        //Wn_f1_Wn_f2(type_filt, u_f1, u_f2);
                //        double Bw = 0, Wn= 0;
                //        switch (type_filt)
                //        {
                //            //Band-pass
                //            case 0:
                //                Bw = u_f2 - u_f1;
                //                Wn = Math.Sqrt(u_f1 * u_f2);
                //                break;
                //            //Band-stop
                //            case 1:
                //                Bw = u_f2 - u_f1;
                //                Wn = Math.Sqrt(u_f1 * u_f2);
                //                break;
                //            //Low-pass
                //            case 2:
                //                Wn = u_f2;
                //                break;
                //            //High-pass
                //            case 3:
                //                Wn = u_f1;
                //                break;
                //        }

                //        //Step 3: Get N - th order Butterworth analog lowpass prototype
                //        //Buttap(order_filt, out p);
                //        double order_filt_exp = (double)order_filt;
                //        int temp_length_vec = 0;
                //        int kkk = 1;
                //        do
                //        {
                //            temp_length_vec++;
                //            kkk += 2;
                //        } while (kkk <= order_filt - 1);

                //        Complex[] ptemp = new Complex[temp_length_vec];

                //        int track_cell = 0;
                //        for (double kk = 0; kk < (double)order_filt - 1; kk += 2)
                //        {
                //            ptemp[track_cell] = Complex.Exp(Complex.ImaginaryOne * ((Math.PI * (kk + 1)) / (2 * order_filt_exp) + Math.PI / 2));
                //            track_cell++;
                //        }

                //        Complex[] p = new Complex[order_filt];
                //        Math.DivRem(order_filt, 2, out int temp_rem);

                //        if (temp_rem != 0)
                //        {
                //            for (int kk = 0; kk < order_filt; kk++)
                //                if (kk < order_filt - 1) p[kk] = Complex.ImaginaryOne; 
                //                else p[kk] = -Complex.One;
                //        }
                //        else
                //        {
                //            for (int kk = 0; kk < order_filt; kk++) p[kk] = Complex.ImaginaryOne;
                //        }

                //        if (order_filt > 1)
                //        {
                //            track_cell = 0;
                //            for (int kk = 0; kk < temp_length_vec * 2; kk += 2)
                //            {
                //                p[kk] = ptemp[track_cell];
                //                p[kk + 1] = Complex.Conjugate(ptemp[track_cell]);
                //                track_cell++;
                //            }
                //        }

                //        //Step 4: Transform to state-space
                //        int temp_dim_arr_matr;
                //        Matrix<Complex> a;
                //        Complex[] b, c;
                //        double d;
                //        Zp2ss(order_filt, out temp_dim_arr_matr, out a, out b, out c, out d);


                //        if (order_filt > 1)
                //        {
                //            temp_dim_arr_matr -= 2;
                //        }
                //        else
                //        {
                //            temp_dim_arr_matr = order_filt;
                //        }

                //        //Copy the values of the matrix/arrays "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
                //        Matrix<Complex> a_arma = Matrix<Complex>.Build.Dense(2 * temp_dim_arr_matr, 2 * temp_dim_arr_matr);
                //        Matrix<Complex> a_arma_p_eye = Matrix<Complex>.Build.DenseIdentity(temp_dim_arr_matr, temp_dim_arr_matr);
                //        Matrix<Complex> a_arma_n_eye = -Matrix<Complex>.Build.DenseIdentity(temp_dim_arr_matr, temp_dim_arr_matr);
                //        Matrix<Complex> b_arma = Matrix<Complex>.Build.Dense(2 * temp_dim_arr_matr, 1);
                //        Matrix<Complex> c_arma = Matrix<Complex>.Build.Dense(1, 2 * temp_dim_arr_matr);
                //        Matrix<Complex> d_arma = Matrix<Complex>.Build.Dense(1, 1);

                //        double q = Wn / Bw;
                //        for (int kk = 0; kk < 2 * temp_dim_arr_matr; kk++)
                //        {
                //            if (kk < temp_dim_arr_matr)
                //            {
                //                b_arma[kk, 0] = b[kk] * Wn / q;
                //                c_arma[0, kk] = c[kk];
                //            }

                //            for (int ll = 0; ll < 2 * temp_dim_arr_matr; ll++)
                //            {
                //                if (kk < temp_dim_arr_matr)
                //                {
                //                    if (ll < temp_dim_arr_matr)
                //                    {
                //                        a_arma[kk, ll] = Wn * a[kk, ll] / q;
                //                    }
                //                    else
                //                    {
                //                        a_arma[kk, ll] = Wn * a_arma_p_eye[kk, ll - temp_dim_arr_matr];
                //                    }
                //                }
                //                else
                //                {
                //                    if (ll < temp_dim_arr_matr)
                //                    {
                //                        a_arma[kk, ll] = Wn * a_arma_n_eye[kk - temp_dim_arr_matr, ll];
                //                    }
                //                }
                //            }
                //        }

                //        //Step 5: Use Bilinear transformation to find discrete equivalent
                //        //Bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt);
                //        double t_arma;
                //        double r_arma;

                //        Matrix<Complex> t1_arma_eye;
                //        Matrix<Complex> t2_arma_eye;
                //        Matrix<Complex> t1_arma, t2_arma, ad_arma, bd_arma, cd_arma, dd_arma;

                //        if (type_filt > 1)
                //        {
                //            t1_arma = Matrix<Complex>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
                //            t1_arma_eye = Matrix<Complex>.Build.DenseIdentity(temp_dim_arr_matr, temp_dim_arr_matr);
                //            t2_arma = Matrix<Complex>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
                //            t2_arma_eye = Matrix<Complex>.Build.DenseIdentity(temp_dim_arr_matr, temp_dim_arr_matr);
                //            ad_arma = Matrix<Complex>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
                //            bd_arma = Matrix<Complex>.Build.Dense(temp_dim_arr_matr, 1);
                //            cd_arma = Matrix<Complex>.Build.Dense(1, temp_dim_arr_matr);
                //            dd_arma = Matrix<Complex>.Build.Dense(1, 1);
                //        }

                //        else
                //        {

                //            t1_arma = Matrix<Complex>.Build.DenseIdentity(2 * temp_dim_arr_matr, 2 * temp_dim_arr_matr);
                //            t1_arma_eye = Matrix<Complex>.Build.DenseIdentity(2 * temp_dim_arr_matr, 2 * temp_dim_arr_matr);
                //            t2_arma = Matrix<Complex>.Build.DenseIdentity(2 * temp_dim_arr_matr, 2 * temp_dim_arr_matr);
                //            t2_arma_eye = Matrix<Complex>.Build.DenseIdentity(2 * temp_dim_arr_matr, 2 * temp_dim_arr_matr);
                //            ad_arma = Matrix<Complex>.Build.Dense(2 * temp_dim_arr_matr, 2 * temp_dim_arr_matr);
                //            bd_arma = Matrix<Complex>.Build.Dense(2 * temp_dim_arr_matr, 1);
                //            cd_arma = Matrix<Complex>.Build.Dense(1, 2 * temp_dim_arr_matr);
                //            dd_arma = Matrix<Complex>.Build.Dense(1, 1);

                //        }

                //        t_arma = (1 / fs);
                //        r_arma = Math.Sqrt(t_arma);
                //        t1_arma = t1_arma_eye + a_arma * t_arma * 0.5; //t1_arma.eye() 
                //        t2_arma = t2_arma_eye - a_arma * t_arma * 0.5;
                //        ad_arma = t1_arma * (t2_arma.PseudoInverse());
                //        bd_arma = (t_arma / r_arma) * t2_arma.Solve(b_arma);//arma::solve(t2_arma, b_arma_f);
                //        cd_arma = (r_arma * c_arma) * (t2_arma.PseudoInverse());
                //        dd_arma = (c_arma * t2_arma.PseudoInverse()) * b_arma * (t_arma / 2) + d_arma;



                //        //Step 6: Transform to zero-pole-gain and polynomial forms
                //        //Zero_pole_gain(ad_arma, type_filt, order_filt, Wn, Bw);
                //        int dim_array;
                //        double[] num_filt, den_filt;
                //        if (type_filt > 1)
                //        {
                //            //Initialize the vectors "num_filt" and "den_filt"
                //            num_filt = new double[order_filt + 1];
                //            den_filt = new double[order_filt + 1];
                //            dim_array = temp_dim_arr_matr;
                //        }
                //        else
                //        {
                //            //Initialize the vectors "num_filt" and "den_filt"
                //            num_filt = new double[2 * order_filt + 1];
                //            den_filt = new double[2 * order_filt + 1];
                //            dim_array = 2 * temp_dim_arr_matr;
                //        }

                //        //Extract the coefficients of the denumerator
                //        Complex[] coeff_pol = new Complex[temp_dim_arr_matr + 1];

                //        if (type_filt > 1) coeff_pol = Char_poly(a_arma_f, temp_dim_arr_matr);
                //        else coeff_pol = Char_poly(a_arma_f, 2 * temp_dim_arr_matr);

                //        for (int qq = 0; qq < dim_array + 1; qq++)
                //        {
                //            den_filt[qq] = coeff_pol[qq].Real;
                //            save_filt_coeff[1][qq] = den_filt[qq];
                //        }

                //        //Extract the coefficients of the denominator
                //        double w = 0;
                //        Wn = 2 * Math.Atan2(Wn, 4);
                //        Complex[] r;

                //        switch (type_filt)
                //        {
                //            case 0: // band-pass
                //                r = new Complex[dim_array + 1];
                //                for (int kk = 0; kk < dim_array; kk++)
                //                {
                //                    if (kk < temp_dim_arr_matr) r[kk] = 1;
                //                    else r[kk] = -1;
                //                }
                //                w = Wn;
                //                break;
                //            case 1: // band-stop
                //                r = new Complex[dim_array + 1];
                //                for (int kk = 0; kk < dim_array; kk++) r[kk] = Complex.Exp(Complex.ImaginaryOne * Wn * Math.Pow(-1, kk));
                //                w = 0;
                //                break;
                //            case 2: // low-pass
                //                r = new Complex[dim_array + 1];
                //                for (int kk = 0; kk < dim_array; kk++) r[kk] = -1;
                //                w = 0;
                //                break;
                //            case 3: //high-pass
                //                r = new Complex[dim_array + 1];
                //                for (int kk = 0; kk < dim_array; kk++) r[kk] = 1;
                //                w = Math.PI;
                //                break;
                //            default:
                //                r = new Complex[dim_array + 1];
                //                break;
                //        }

                //        Complex[] coeff_pol_num = new Complex[dim_array + 1];

                //        coeff_pol_num = Poly(r, dim_array);

                //        Complex[] kern = new Complex[dim_array + 1];

                //        for (int kk = 0; kk < dim_array + 1; kk++) kern[kk] = Complex.Exp(-Complex.ImaginaryOne * w * kk);

                //        Complex temp_sum_I;
                //        Complex temp_sum_II;

                //        for (int kk = 0; kk < dim_array + 1; kk++)
                //        {
                //            temp_sum_I = new Complex(0.0, 0.0);
                //            temp_sum_II = new Complex(0.0, 0.0);
                //            for (int hh = 0; hh < dim_array + 1; hh++)
                //            {
                //                temp_sum_I += kern[hh] * den_filt[hh];
                //                temp_sum_II += kern[hh] * coeff_pol_num[hh];
                //            }
                //            num_filt[kk] = (coeff_pol_num[kk] * temp_sum_I / temp_sum_II).Real;
                //            save_filt_coeff[0][kk] = num_filt[kk];
                //        }
                //        return save_filt_coeff;

                //    }

                //    public void Zp2ss(int order_filt, out int temp_dim_arr_matr, out Matrix<Complex> a, out Complex[] b, out Complex[] c, out double d)
                //    {
                //        //Order the pairs of complex conjugate. The pairs with the smallest real part come first. Within pairs, the ones with the negative imaginary part comes first    
                //        Complex temp_max;
                //        Complex[] p;
                //        temp_dim_arr_matr = 0;
                //        a = Matrix<Complex>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
                //        b = new Complex[0];
                //        c = new Complex[0];

                //        //Using the selection sort algorithm to order based on the real part
                //        double order_filt_d = (double)order_filt;
                //        double temp_rem;
                //        Math.DivRem(order_filt, 2, out int res_rem);
                //        if (res_rem == 0)
                //        {
                //            temp_rem = order_filt_d;
                //        }
                //        else
                //        {
                //            temp_rem = order_filt_d - 1;
                //        }

                //        int min_real;
                //        for (int kk = 0; kk < temp_rem - 1; kk++)
                //        {
                //            min_real = kk;
                //            for (int jj = kk + 1; jj < temp_rem; jj++)
                //            {
                //                if (p[jj].Real < p[min_real].Real)
                //                {
                //                    min_real = jj;
                //                    temp_max = p[kk];
                //                    p[kk] = p[min_real];
                //                    p[min_real] = temp_max;
                //                }
                //            }
                //        }

                //        //Using the selection sort algorithm to order the values based on the imaginary part
                //        for (int kk = 0; kk < temp_rem - 1; kk += 2)
                //        {
                //            min_real = kk;
                //            if (p[kk].Imaginary > p[kk + 1].Imaginary)
                //            {
                //                min_real = kk + 1;
                //                temp_max = p[kk];
                //                p[kk] = p[min_real];
                //                p[min_real] = temp_max;
                //            }
                //        }

                //        // Initialize state - space matrices for running series
                //        d = 1;

                //        int track_index = 1;

                //        // Take care of any left over unmatched pole pairs.
                //        // H(s) = 1 / (s ^ 2 + den(2)s + den(3))
                //        Complex[] temp_poly_p;

                //        if (order_filt > 1)
                //        {
                //            double[] b1 = new double[] { 1, 0 };
                //            double[] c1 = new double[] { 0, 1 };

                //            Math.DivRem(order_filt, 2, out int rem_div);
                //            temp_rem = (double)rem_div;
                //            int order_filt_temp;
                //            int dim_matr;
                //            Matrix<Complex> temp_matrix_a;    //Temporary matrix where to save the coefficients at each interaction
                //            int coeff_numb = 3;

                //            if (temp_rem == 0)
                //            {
                //                order_filt_temp = order_filt;
                //                dim_matr = order_filt_temp / 2;
                //                temp_matrix_a = Matrix<Complex>.Build.Dense(dim_matr, coeff_numb);
                //            }
                //            else
                //            {
                //                order_filt_temp = order_filt - 1;
                //                dim_matr = order_filt_temp / 2;
                //                temp_matrix_a = Matrix<Complex>.Build.Dense(dim_matr, coeff_numb);
                //            }

                //            int track_cycles = 0;
                //            int temp_val_pos_check;
                //            int dim_poly_p = 2;
                //            temp_poly_p = new Complex[dim_poly_p];
                //            temp_dim_arr_matr = dim_poly_p + (order_filt % 2);

                //            while (track_index < order_filt_temp)
                //            {
                //                for (int rr = track_index - 1; rr < track_index + 1; rr++) temp_poly_p[rr - track_index + 1] = p[rr];
                //                Complex[] coeff_pol = new Complex[dim_poly_p + 1];
                //                coeff_pol[0] = 1;

                //                for (int ll = 0; ll < dim_poly_p; ll++)
                //                {
                //                    int yy = 0;
                //                    do
                //                    {
                //                        coeff_pol[ll + 1 - yy] = coeff_pol[ll + 1 - yy] - temp_poly_p[ll] * coeff_pol[ll - yy];
                //                        yy++;
                //                    } while (yy <= ll);
                //                }


                //                for (int qq = 0; qq < coeff_numb; qq++) temp_matrix_a[track_cycles, qq] = -(coeff_pol[qq].Real);

                //                //Update the state-space arrays/matrix
                //                track_cycles += 1;
                //                a = Matrix<Complex>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);

                //                int track_index_coeff = 0;

                //                Math.DivRem(order_filt, 2, out int div_res);
                //                if (div_res == 0)
                //                {
                //                    ///////////////////////////////////////////////////////////////////////////////////////////////
                //                    //Even number of poles
                //                    for (int kk = 0; kk < temp_dim_arr_matr; kk++)
                //                    {
                //                        for (int gg = 0; gg < temp_dim_arr_matr; gg++)
                //                        {
                //                            temp_val_pos_check = kk - gg;
                //                            switch (temp_val_pos_check)
                //                            {
                //                                case 1:
                //                                    a[kk, gg] = 1;
                //                                    break;
                //                                case 0:
                //                                    Math.DivRem(kk + 1, 2, out int div_rem_I);
                //                                    if ((div_rem_I != 0) || (kk == 0)) a[kk, gg] = temp_matrix_a[track_index_coeff, 1];
                //                                    else a[kk, gg] = 0;
                //                                    break;
                //                                case -1:
                //                                    Math.DivRem(kk + 1, 2, out int div_rem_II);
                //                                    if ((div_rem_II != 0) || (kk == 0))
                //                                    {
                //                                        a[kk, gg] = temp_matrix_a[track_index_coeff, 2];
                //                                        track_index_coeff++;
                //                                    }
                //                                    else a[kk, gg] = 0;
                //                                    break;
                //                                default:
                //                                    a[kk, gg] = 0;
                //                                    break;
                //                            }
                //                        }
                //                    }
                //                    ///////////////////////////////////////////////////////////////////////////////////////////////
                //                }
                //                else
                //                {
                //                    ///////////////////////////////////////////////////////////////////////////////////////////////
                //                    //Odd number of poles
                //                    for (int kk = 0; kk < temp_dim_arr_matr; kk++)
                //                    {
                //                        for (int gg = 0; gg < temp_dim_arr_matr; gg++)
                //                        {
                //                            temp_val_pos_check = kk - gg;

                //                            switch (temp_val_pos_check)
                //                            {
                //                                case 1:
                //                                    a[kk, gg] = 1;
                //                                    break;
                //                                case 0:
                //                                    if (kk == 0) a[kk, gg] = -1;
                //                                    else
                //                                    {
                //                                        Math.DivRem(kk + 1, 2, out int div_rem_III);
                //                                        if (div_rem_III == 0) a[kk, gg] = temp_matrix_a[track_index_coeff, 1];
                //                                        else a[kk, gg] = 0;
                //                                    }
                //                                    break;
                //                                case -1:
                //                                    Math.DivRem(kk + 1, 2, out int div_rem_IV);
                //                                    if (div_rem_IV == 0)
                //                                    {
                //                                        a[kk, gg] = temp_matrix_a[track_index_coeff, 2];
                //                                        track_index_coeff++;
                //                                    }
                //                                    else a[kk, gg] = 0;
                //                                    break;
                //                                default:
                //                                    a[kk, gg] = 0;
                //                                    break;
                //                            }
                //                        }
                //                    }
                //                }
                //                ///////////////////////////////////////////////////////////////////////////////////////////////

                //                //Initialize the vectors "b" and "c"
                //                b = new Complex[temp_dim_arr_matr];
                //                c = new Complex[temp_dim_arr_matr];

                //                for (int kk = 0; kk < temp_dim_arr_matr; kk++)
                //                {
                //                    if (kk == 0)
                //                    {
                //                        b[kk] = 1;
                //                    }
                //                    else
                //                    {
                //                        b[kk] = 0;
                //                    }
                //                }

                //                for (int kk = 0; kk < temp_dim_arr_matr; kk++)
                //                {
                //                    if (kk == temp_dim_arr_matr - 1)
                //                    {
                //                        c[kk] = 1;
                //                    }
                //                    else
                //                    {
                //                        c[kk] = 0;
                //                    }
                //                }
                //                track_index += 2;

                //                if (track_index < order_filt_temp)
                //                {
                //                    //Clean up the matrix "a" and the arrays "b" and "c", so they can be re-initialized
                //                    a.Clear();
                //                    Array.Clear(b, 0, b.Length);
                //                    Array.Clear(b, 0, c.Length);
                //                }
                //                dim_matr += 2;
                //                temp_dim_arr_matr += 2;
                //            }
                //        }
                //        else
                //        {
                //            Matrix<Complex> a = Matrix<Complex>.Build.Dense(1, 1);
                //            a[0, 0] = p[0];
                //            b = new Complex[1];
                //            b[0] = 1;
                //            c = new Complex[1];
                //            c[0] = 1;
                //        }
                //        d = 0;
                //    }


                //    //Calculate the coefficients of the characteristic polynomial (Bernard Brooks' paper (2023))
                //    public Complex[] Char_poly(Matrix<Complex> temp_matr_poly, int row_col)
                //    {
                //        Complex[] coeff_pol_ff = new Complex[row_col + 1];
                //        Matrix<Complex> temp_val = Matrix<Complex>.Build.Dense(1, 1);
                //        int num_det = 0;
                //        Matrix<Complex> temp_matr = Matrix<Complex>.Build.Dense(row_col, row_col);

                //        for (int kk = 0; kk < row_col + 1; kk++)
                //        {
                //            if (kk == 0)
                //            {
                //                coeff_pol_ff[row_col - kk] = Math.Pow(-1, row_col) * temp_matr_poly.Determinant();
                //            }
                //            else
                //            {
                //                int[][] matrix_comb;
                //                temp_val[0, 0] = 0;
                //                try
                //                {
                //                    num_det = factorial(row_col) / (factorial(row_col - kk) * factorial(kk));  //Calculate the number of combinations   
                //                }
                //                //If the DivideByZero exception is thrown, the filter is unstable. The method returns a default value of 10^10 for the coefficients.
                //                catch (System.DivideByZeroException)
                //                {
                //                    for (int ll = 0; ll < row_col + 1; ll++) coeff_pol_ff[ll] = Math.Pow(10, 10);
                //                    break;
                //                }

                //                try
                //                {
                //                    matrix_comb = new int[num_det][];
                //                }
                //                //If the numerical overflow exception is thrown, the filter is unstable. The method returns a default value of 10^10 for the coefficients.
                //                catch (System.OverflowException)
                //                {
                //                    for (int ll = 0; ll < row_col + 1; ll++) coeff_pol_ff[ll] = Math.Pow(10, 10);
                //                    break;
                //                }

                //                for (int ll = 0; ll < num_det; ll++)
                //                {
                //                    matrix_comb[ll] = new int[num_det];
                //                }

                //                // Generate the combinations 
                //                matrix_comb = Combination_method(row_col, kk, num_det);
                //                for (int mm = 0; mm < num_det; mm++)

                //                {
                //                    temp_matr = temp_matr_poly.Clone();
                //                    for (int pp = 0; pp < row_col; pp++)
                //                    {
                //                        temp_matr[matrix_comb[mm][0], pp] = 0;
                //                        temp_matr[pp, matrix_comb[mm][0]] = 0;
                //                        temp_matr[matrix_comb[mm][0], matrix_comb[mm][0]] = -1;
                //                    }

                //                    for (int nn = 1; nn < kk; nn++)
                //                    {
                //                        for (int pp = 0; pp < row_col; pp++)
                //                        {
                //                            temp_matr[matrix_comb[mm][nn], pp] = 0;
                //                            temp_matr[pp, matrix_comb[mm][nn]] = 0;
                //                            temp_matr[matrix_comb[mm][nn], matrix_comb[mm][nn]] = -1;
                //                        }
                //                    }
                //                    temp_val[0, 0] += temp_matr.Determinant();
                //                }
                //                coeff_pol_ff[row_col - kk] = Math.Pow(-1, row_col) * temp_val[0, 0];
                //            }
                //        }
                //        return coeff_pol_ff;
                //    }

                //    //Method to calculate all the possible combinations
                //    public int[][] Combination_method(int N, int K, int comb_n)
                //    {
                //        var integers = new List<int> { };
                //        for (int kk = 0; kk < N; kk++) integers.Add(kk);
                //        var matrix_comb_f_temp = new Combinations<int>(integers, K).ToArray();  //Create the combinations
                //        int[][] matrix_comb_f = new int[comb_n][];
                //        for (int hh = 0; hh < comb_n; hh++) matrix_comb_f[hh] = new int[K];
                //        //Save the combinations in a 2D array
                //        for (int hh = 0; hh < matrix_comb_f_temp.Length; hh++)
                //        {
                //            for (int kk = 0; kk < K; kk++) matrix_comb_f[hh][kk] = matrix_comb_f_temp[hh][kk];
                //        }
                //        return matrix_comb_f;
                //    }
                
                //public int main(int n, double sampling_frequency, double freq_l, double freq_u )
                //{

                //    double a = Math.Cos(Math.PI * (freq_l + freq_u) / sampling_frequency) / Math.Cos(Math.PI * (freq_l - freq_u) / sampling_frequency);
                //    double a2 = a * a;
                //    double b = Math.Tan(Math.PI * (freq_l - freq_u) / sampling_frequency);
                //    double b2 = b * b;
                //    double r;

                //    n = n / 4;
                //    double[] A = new double[n];
                //    double[] d1 = new double[n];
                //    double[] d2 = new double[n];
                //    double[] d3 = new double[n];
                //    double[] d4 = new double[n];
                //    double[] w0 = new double[n];
                //    double[] w1 = new double[n];
                //    double[] w2 = new double[n];
                //    double[] w3 = new double[n];
                //    double[] w4 = new double[n];
                //    double x;

                //    for (int i = 0; i < n; ++i)
                //    {
                //        r = Math.Sin(Math.PI * (2.0 * i + 1.0) / (4.0 * n));
                //        double s = b2 + 2.0 * b * r + 1.0;
                //        A[i] = b2 / s;
                //        d1[i] = 4.0 * a * (1.0 + b * r) / s;
                //        d2[i] = 2.0 * (b2 - 2.0 * a2 - 1.0) / s;
                //        d3[i] = 4.0 * a * (1.0 - b * r) / s;
                //        d4[i] = -(b2 - 2.0 * b * r + 1.0) / s;
                //    }

                //    //Doing the filtering

                //    while (scanf("%lf", &x) != EOF)
                //    {
                //        for (i = 0; i < n; ++i)
                //        {
                //            w0[i] = d1[i] * w1[i] + d2[i] * w2[i] + d3[i] * w3[i] + d4[i] * w4[i] + x;
                //            x = A[i] * (w0[i] - 2.0 * w2[i] + w4[i]);
                //            w4[i] = w3[i];
                //            w3[i] = w2[i];
                //            w2[i] = w1[i];
                //            w1[i] = w0[i];
                //        }
                //        printf("%lf\n", x);
                //    }

                //    return (0);
                //}
            }
        }
    }
}