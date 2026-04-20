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

using Eto.Forms;
using MathNet.Numerics.Statistics;
using Pachyderm_Acoustic.Audio;
using Pachyderm_Acoustic.Pach_Graphics;
using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Diagnostics.Eventing.Reader;
using System.Linq;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;

namespace Pachyderm_Acoustic
{
    namespace Environment
    {
        public abstract class Material
        {
            public abstract void Absorb(ref OctaveRay Ray, Hare.Geometry.Vector Normal);
            public abstract void Absorb(ref BroadRay Ray, Hare.Geometry.Vector Normal);
            public abstract void Absorb(ref OctaveRay Ray, out double cos_theta, Hare.Geometry.Vector Normal);
            public abstract void Absorb(ref BroadRay Ray, out double cos_theta, Hare.Geometry.Vector Normal);
            public abstract System.Numerics.Complex Reflection_Narrow(double frequency);
            public abstract System.Numerics.Complex Reflection_Narrow(double frequency, Hare.Geometry.Vector Dir, Hare.Geometry.Vector Normal);
            public abstract double Coefficient_A_Broad(int Octave);
            public abstract double[] Coefficient_A_Broad();
            public abstract System.Numerics.Complex[] Reflection_Spectrum(int sample_frequency, int length, Hare.Geometry.Vector Normal, Hare.Geometry.Vector Dir, int threadid);
            public abstract (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0);
            public virtual void ForceIIR(double[] a, double[] b, double fs) { }
        }

        public abstract class Scattering
        {
            public abstract double Coefficient(int octave);
            public abstract double[] Coefficient();
            public abstract void Scatter_Early(ref BroadRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, double[] Transmission = null);
            public abstract void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false);
            public abstract void Scatter_VeryLate(ref OctaveRay Ray, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false);
        }

        public class Basic_Material : Material
        {
            double[] Abs = new double[8];
            double[] Ref = new double[8];
            double[] Abs_3rd = new double[24];
            double[] Ref_3rd = new double[24];
            MathNet.Numerics.Interpolation.CubicSpline Transfer_Function;

            public Basic_Material(double[] ABS)//, double[] Phase_Delay)
            {
                //Interpolate a transfer function...
                double rt2 = Math.Sqrt(2);

                List<double> f = new List<double>();
                List<double> a = new List<double>();
                List<double> pr = new List<double>();

                f.Add(0);
                f.Add(31.25);
                for (int oct = 0; oct < 8; oct++)
                {
                    f.Add(62.5 * Math.Pow(2, oct));
                }
                f.Add(62.5 * Math.Pow(2, 8));
                f.Add(24000);

                if (Abs.Length == 8)
                {
                    Abs = ABS;

                    //F_Spectrum in...
                    a.Add(Abs[0]);
                    a.Add(Abs[0]);

                    for (int oct = 0; oct < 8; oct++)
                    {
                        a.Add(Abs[oct]);
                    }
                    a.Add(Abs[7]);
                    a.Add(Abs[7]);

                    MathNet.Numerics.Interpolation.CubicSpline Alpha = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkimaSorted(f.ToArray(), a.ToArray());

                    //F_Spectrum Out...
                    for (int i = 0; i < ABS.Length; i++) Ref[i] = 1 - ABS[i];
                    double thirdmod = Math.Pow(2, 1 / 6);

                    for (int oct = 0; oct < 24; oct++)
                    {
                        double f_center = 50 * Math.Pow(2, oct / 3.0);
                        Abs_3rd[oct] = Alpha.Interpolate(f_center);
                        Ref_3rd[oct] = 1 - Abs_3rd[oct];
                    }
                }
                else if (Abs.Length == 24)
                {
                    Abs_3rd = ABS;
                    for (int oct = 0; oct < 8; oct++)
                    {
                        int octave = oct * 3;
                        Ref_3rd[octave] = 1 - Abs_3rd[octave];
                        Ref_3rd[octave + 1] = 1 - Abs_3rd[octave + 1];
                        Ref_3rd[octave + 2] = 1 - Abs_3rd[octave + 2];
                        Abs[oct] = (Abs_3rd[octave] + Abs_3rd[octave + 1] + Abs_3rd[octave + 2]) / 3;
                        Ref[oct] = 1 - Abs[oct];
                    }
                }

                pr.Add(Math.Sqrt(1 - Abs[0]));
                pr.Add(Math.Sqrt(1 - Abs[0]));

                for (int oct = 0; oct < 8; oct++)
                {
                    pr.Add(Math.Sqrt(1 - Abs[oct]));
                }
                pr.Add(Math.Sqrt(1 - Abs[7]));
                pr.Add(Math.Sqrt(1 - Abs[7]));

                while (pr.Count < f.Count) pr.Add(1 - Abs[7]);

                Transfer_Function = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkimaSorted(f.ToArray(), pr.ToArray());
            }

            public override System.Numerics.Complex[] Reflection_Spectrum(int sample_frequency, int length, Hare.Geometry.Vector Normal, Hare.Geometry.Vector Dir, int threadid)
            {
                double[] pr = new double[8] { Math.Sqrt(Ref[0]), Math.Sqrt(Ref[1]), Math.Sqrt(Ref[2]), Math.Sqrt(Ref[3]), Math.Sqrt(Ref[4]), Math.Sqrt(Ref[5]), Math.Sqrt(Ref[6]), Math.Sqrt(Ref[7]) };

                double[] filter = Audio.Pach_SP.Magnitude_Filter(pr, sample_frequency, length, threadid);

                System.Numerics.Complex[] Ref_trns = new System.Numerics.Complex[filter.Length];
                for (int i = 0; i < filter.Length; i++) Ref_trns[i] = filter[i];

                return Ref_trns;
            }

            public override void Absorb(ref OctaveRay Ray, Hare.Geometry.Vector Normal)
            {
                Ray.Intensity *= (Ref[Ray.Octave]);
            }

            public override void Absorb(ref OctaveRay Ray, out double cos_theta, Hare.Geometry.Vector Normal)
            {
                cos_theta = Hare.Geometry.Hare_math.Dot(Normal.dx, Normal.dy, Normal.dz, Ray.dx, Ray.dy, Ray.dz);
                Ray.Intensity *= (Ref[Ray.Octave]);
            }

            public override void Absorb(ref BroadRay Ray, out double cos_theta, Hare.Geometry.Vector Normal)
            {
                cos_theta = Hare.Geometry.Hare_math.Dot(Ray.dx, Ray.dy, Ray.dz, Normal.dx, Normal.dy, Normal.dz);

                Ray.Energy[0] *= (Ref[0]);
                Ray.Energy[1] *= (Ref[1]);
                Ray.Energy[2] *= (Ref[2]);
                Ray.Energy[3] *= (Ref[3]);
                Ray.Energy[4] *= (Ref[4]);
                Ray.Energy[5] *= (Ref[5]);
                Ray.Energy[6] *= (Ref[6]);
                Ray.Energy[7] *= (Ref[7]);
            }

            public override void Absorb(ref BroadRay Ray, Hare.Geometry.Vector Normal)
            {
                Ray.Energy[0] *= (Ref[0]);
                Ray.Energy[1] *= (Ref[1]);
                Ray.Energy[2] *= (Ref[2]);
                Ray.Energy[3] *= (Ref[3]);
                Ray.Energy[4] *= (Ref[4]);
                Ray.Energy[5] *= (Ref[5]);
                Ray.Energy[6] *= (Ref[6]);
                Ray.Energy[7] *= (Ref[7]);
            }

            public override double[] Coefficient_A_Broad()
            {
                return Abs;
            }

            public override double Coefficient_A_Broad(int Octave)
            {
                return Abs[Octave];
            }

            public override System.Numerics.Complex Reflection_Narrow(double frequency)
            {
                return new System.Numerics.Complex(Transfer_Function.Interpolate(frequency), 0);
            }

            public override System.Numerics.Complex Reflection_Narrow(double frequency, Hare.Geometry.Vector Dir, Hare.Geometry.Vector Normal)
            {
                return new System.Numerics.Complex(Transfer_Function.Interpolate(frequency), 0);
            }

            int rec_order = 0;
            double[] rec_a, rec_b;
            double rec_fs = 0;
            double rec_maxfreq = -1;
            object lock_IIR = new object();
            bool _iirForced = false;

            public override void ForceIIR(double[] a, double[] b, double fs)
            {
                lock (lock_IIR)
                {
                    rec_a = a != null ? (double[])a.Clone() : null;
                    rec_b = b != null ? (double[])b.Clone() : null;
                    rec_fs = fs;
                }
            }
            public override (double[] a, double[] b) Estimate_IIR_Coefficients(
    double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            {
                int samplect = 4096;
                int K = (int)Math.Floor(samplect * max_freq / sample_frequency);
                K = Math.Max(8, Math.Min(K, samplect / 2 + 1));

                frequencies = new double[K];
                for (int i = 0; i < K; i++)
                    frequencies[i] = i * sample_frequency / samplect; // correct FFT bin mapping

                lock (lock_IIR)
                {
                    // rec_order is used here as the requested pole budget cache key
                    int requestedPoleBudget = (filter_order > 0) ? Math.Max(1, Math.Min(filter_order, 6)) : 0;

                    if (rec_a != null && rec_b != null &&
                        Math.Abs(rec_fs - sample_frequency) < 1e-9 &&
                        Math.Abs(rec_maxfreq - max_freq) < 1e-9 &&
                        rec_order == requestedPoleBudget)
                    {
                        return (rec_a, rec_b);
                    }

                    rec_fs = sample_frequency;
                    rec_maxfreq = max_freq;
                    rec_order = requestedPoleBudget;

                    // Reflection spectrum
                    Complex[] Rfull = this.Reflection_Spectrum(
                        (int)sample_frequency,
                        samplect,
                        new Hare.Geometry.Vector(0, 0, 1),
                        new Hare.Geometry.Vector(0, 0, 1),
                        0);

                    Complex[] R = new Complex[K];
                    Array.Copy(Rfull, R, K);

                    // Convert to target admittance Y = (1 - R) / (1 + R)
                    Complex[] Y = new Complex[K];
                    for (int i = 0; i < K; i++)
                    {
                        Complex denom = Complex.One + R[i];
                        if (denom.Magnitude < 1e-9) denom = new Complex(1e-9, 0);

                        Complex Yi = (Complex.One - R[i]) / denom;

                        // keep target passive
                        if (Yi.Real < 0) Yi = new Complex(0, Yi.Imaginary);

                        Y[i] = Yi;
                    }

                    // Small pole budget, better search
                    int poleBudget = (filter_order > 0) ? Math.Max(1, Math.Min(filter_order, 6)) : 4;

                    (rec_a, rec_b) = BasicMaterialPassiveFit.FitLimitedPoleBudget(
                        frequencies,
                        Y,
                        sample_frequency,
                        poleBudget,
                        fLo: 20.0,
                        fHi: Math.Min(max_freq * 0.95, 0.45 * sample_frequency),
                        poleCap: 0.95,      // keep memory short; safer at corners
                        candidateCount: 28, // search harder
                        refinePasses: 4,    // refine chosen poles
                        ridge: 1e-8,
                        iters: 1200
                    );

                    return (rec_a, rec_b);
                }
            }
            internal static class BasicMaterialPassiveFit
            {
                private static double[] Conv(double[] a, double[] b)
                {
                    double[] c = new double[a.Length + b.Length - 1];
                    for (int i = 0; i < a.Length; i++)
                        for (int j = 0; j < b.Length; j++)
                            c[i + j] += a[i] * b[j];
                    return c;
                }

                private static Complex EvalOnePole(double p, Complex z1)
                {
                    Complex den = Complex.One - p * z1;
                    if (den.Magnitude < 1e-15) den = new Complex(1e-15, 0);
                    return (1.0 - p) / den;
                }

                private static double[] BuildDen(double[] poles)
                {
                    double[] a = new double[] { 1.0 };
                    for (int i = 0; i < poles.Length; i++)
                        a = Conv(a, new double[] { 1.0, -poles[i] });
                    return a;
                }

                private static void BuildPrefixSuffix(double[] poles, out double[][] prefix, out double[][] suffix)
                {
                    int M = poles.Length;
                    prefix = new double[M + 1][];
                    suffix = new double[M + 1][];

                    prefix[0] = new double[] { 1.0 };
                    for (int i = 0; i < M; i++)
                        prefix[i + 1] = Conv(prefix[i], new double[] { 1.0, -poles[i] });

                    suffix[M] = new double[] { 1.0 };
                    for (int i = M - 1; i >= 0; i--)
                        suffix[i] = Conv(new double[] { 1.0, -poles[i] }, suffix[i + 1]);
                }

                private static double[] SolveNNLSProjected(double[,] A, double[] y, int iters, double ridge)
                {
                    int rows = A.GetLength(0);
                    int cols = A.GetLength(1);

                    double[] w = new double[cols];
                    double[] grad = new double[cols];
                    double[] Aw = new double[rows];

                    // rough Lipschitz estimate
                    double[] v = new double[cols];
                    for (int i = 0; i < cols; i++) v[i] = 1.0 / Math.Max(1, cols);

                    for (int t = 0; t < 20; t++)
                    {
                        double[] u = new double[rows];
                        for (int r = 0; r < rows; r++)
                        {
                            double s = 0;
                            for (int c = 0; c < cols; c++) s += A[r, c] * v[c];
                            u[r] = s;
                        }

                        double[] v2 = new double[cols];
                        for (int c = 0; c < cols; c++)
                        {
                            double s = 0;
                            for (int r = 0; r < rows; r++) s += A[r, c] * u[r];
                            v2[c] = s + ridge * v[c];
                        }

                        double norm = 0;
                        for (int c = 0; c < cols; c++) norm += v2[c] * v2[c];
                        norm = Math.Sqrt(Math.Max(1e-30, norm));
                        for (int c = 0; c < cols; c++) v[c] = v2[c] / norm;
                    }

                    double[] Av = new double[rows];
                    for (int r = 0; r < rows; r++)
                    {
                        double s = 0;
                        for (int c = 0; c < cols; c++) s += A[r, c] * v[c];
                        Av[r] = s;
                    }

                    double[] AtAv = new double[cols];
                    for (int c = 0; c < cols; c++)
                    {
                        double s = 0;
                        for (int r = 0; r < rows; r++) s += A[r, c] * Av[r];
                        AtAv[c] = s + ridge * v[c];
                    }

                    double num = 0, den = 0;
                    for (int c = 0; c < cols; c++)
                    {
                        num += v[c] * AtAv[c];
                        den += v[c] * v[c];
                    }

                    double L = Math.Max(1e-6, num / Math.Max(1e-30, den));
                    double step = 1.0 / L;

                    for (int it = 0; it < iters; it++)
                    {
                        for (int r = 0; r < rows; r++)
                        {
                            double s = 0;
                            for (int c = 0; c < cols; c++) s += A[r, c] * w[c];
                            Aw[r] = s;
                        }

                        for (int c = 0; c < cols; c++)
                        {
                            double s = 0;
                            for (int r = 0; r < rows; r++) s += A[r, c] * (Aw[r] - y[r]);
                            grad[c] = 2.0 * s + 2.0 * ridge * w[c];
                        }

                        for (int c = 0; c < cols; c++)
                        {
                            double wc = w[c] - step * grad[c];
                            w[c] = (wc < 0) ? 0 : wc;
                        }
                    }

                    return w;
                }

                public static (double[] a, double[] b) FitLimitedPoleBudget(
                    double[] freqs,
                    Complex[] Ytarget,
                    double fs,
                    int poleCount,
                    double fLo = 20.0,
                    double fHi = -1.0,
                    double poleCap = 0.95,
                    int candidateCount = 28,
                    int refinePasses = 4,
                    double ridge = 1e-8,
                    int iters = 1200)
                {
                    if (fHi <= 0) fHi = 0.45 * fs;
                    fHi = Math.Min(fHi, 0.49 * fs);
                    fLo = Math.Max(1.0, fLo);

                    poleCount = Math.Max(1, poleCount);
                    candidateCount = Math.Max(candidateCount, poleCount + 4);

                    double Weight(double f)
                    {
                        if (f < 20 || f > fHi) return 0.0;
                        if (f >= 125 && f <= 4000) return 1.0;
                        return 0.25;
                    }

                    double PoleFromFreq(double fk)
                    {
                        double p = Math.Exp(-2.0 * Math.PI * fk / fs);
                        p = Math.Min(p, poleCap);
                        p = Math.Max(p, 0.03);
                        return p;
                    }

                    double FreqFromPole(double p)
                    {
                        p = Math.Min(Math.Max(p, 1e-9), 0.999999);
                        return -fs * Math.Log(p) / (2.0 * Math.PI);
                    }

                    void BuildSystem(double[] poles, out double[,] A, out double[] y)
                    {
                        int K = Math.Min(freqs.Length, Ytarget.Length);
                        int rows = 2 * K;
                        int cols = 1 + poles.Length; // g0 + gk

                        A = new double[rows, cols];
                        y = new double[rows];

                        for (int i = 0; i < K; i++)
                        {
                            double f = freqs[i];
                            double wgt = Weight(f);

                            y[i] = wgt * Ytarget[i].Real;
                            y[i + K] = wgt * Ytarget[i].Imaginary;

                            A[i, 0] = wgt;
                            A[i + K, 0] = 0.0;

                            double w = 2.0 * Math.PI * f / fs;
                            Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

                            for (int k = 0; k < poles.Length; k++)
                            {
                                Complex Hk = EvalOnePole(poles[k], z1);
                                A[i, 1 + k] = wgt * Hk.Real;
                                A[i + K, 1 + k] = wgt * Hk.Imaginary;
                            }
                        }
                    }

                    double EvalError(double[] poles)
                    {
                        BuildSystem(poles, out var A, out var y);
                        double[] g = SolveNNLSProjected(A, y, iters, ridge);

                        int rows = A.GetLength(0);
                        int cols = A.GetLength(1);

                        double err2 = 0.0;
                        for (int r = 0; r < rows; r++)
                        {
                            double s = 0.0;
                            for (int c = 0; c < cols; c++) s += A[r, c] * g[c];
                            double e = s - y[r];
                            err2 += e * e;
                        }
                        return err2;
                    }

                    (double[] a, double[] b) BuildFilter(double[] poles)
                    {
                        BuildSystem(poles, out var A, out var y);
                        double[] g = SolveNNLSProjected(A, y, iters, ridge);

                        double[] a = BuildDen(poles);
                        BuildPrefixSuffix(poles, out var prefix, out var suffix);

                        double[] b = new double[a.Length];

                        // g0 * A(z)
                        for (int i = 0; i < a.Length; i++) b[i] += g[0] * a[i];

                        // sum gk (1-pk) Π_{m!=k}(1-pm z^-1)
                        for (int k = 0; k < poles.Length; k++)
                        {
                            double scale = g[1 + k] * (1.0 - poles[k]);
                            if (scale <= 0) continue;

                            double[] without = Conv(prefix[k], suffix[k + 1]);
                            for (int i = 0; i < without.Length; i++) b[i] += scale * without[i];
                        }

                        double a0 = a[0];
                        if (Math.Abs(a0) < 1e-12) a0 = 1.0;
                        if (Math.Abs(a0 - 1.0) > 1e-12)
                        {
                            for (int i = 0; i < a.Length; i++) a[i] /= a0;
                            for (int i = 0; i < b.Length; i++) b[i] /= a0;
                        }

                        return (a, b);
                    }

                    // candidate bank
                    double[] bank = new double[candidateCount];
                    for (int i = 0; i < candidateCount; i++)
                    {
                        double t = (i + 0.5) / candidateCount;
                        double fk = fLo * Math.Pow(fHi / fLo, t);
                        bank[i] = PoleFromFreq(fk);
                    }

                    // greedy selection
                    List<double> selected = new List<double>();
                    HashSet<int> used = new HashSet<int>();

                    for (int m = 0; m < poleCount; m++)
                    {
                        double bestErr = double.PositiveInfinity;
                        int bestIdx = -1;

                        for (int i = 0; i < bank.Length; i++)
                        {
                            if (used.Contains(i)) continue;

                            double[] trial = new double[selected.Count + 1];
                            for (int j = 0; j < selected.Count; j++) trial[j] = selected[j];
                            trial[selected.Count] = bank[i];

                            double err = EvalError(trial);
                            if (err < bestErr)
                            {
                                bestErr = err;
                                bestIdx = i;
                            }
                        }

                        if (bestIdx < 0) break;
                        used.Add(bestIdx);
                        selected.Add(bank[bestIdx]);
                    }

                    // local refinement
                    double[] mults = new double[] { 0.60, 0.75, 0.90, 1.00, 1.10, 1.33, 1.67 };

                    for (int pass = 0; pass < refinePasses; pass++)
                    {
                        for (int j = 0; j < selected.Count; j++)
                        {
                            double baseP = selected[j];
                            double baseF = FreqFromPole(baseP);

                            double bestP = baseP;
                            double bestErr = EvalError(selected.ToArray());

                            for (int m = 0; m < mults.Length; m++)
                            {
                                double fk = Math.Min(fHi, Math.Max(fLo, baseF * mults[m]));
                                double pk = PoleFromFreq(fk);

                                bool tooClose = false;
                                for (int q = 0; q < selected.Count; q++)
                                {
                                    if (q == j) continue;
                                    if (Math.Abs(selected[q] - pk) < 1e-3) { tooClose = true; break; }
                                }
                                if (tooClose) continue;

                                double old = selected[j];
                                selected[j] = pk;

                                double err = EvalError(selected.ToArray());
                                if (err < bestErr)
                                {
                                    bestErr = err;
                                    bestP = pk;
                                }

                                selected[j] = old;
                            }

                            selected[j] = bestP;
                        }
                    }

                    return BuildFilter(selected.ToArray());
                }
            }

            //        int rec_order = 0;
            //        double[] rec_a, rec_b;
            //        double rec_fs = 0;

            //        object lock_IIR = new object();

            //        internal static class StableAdmittanceFitter
            //        {
            //            public static double[] Convolve(double[] a, double[] b)
            //            {
            //                double[] c = new double[a.Length + b.Length - 1];
            //                for (int i = 0; i < a.Length; i++)
            //                    for (int j = 0; j < b.Length; j++)
            //                        c[i + j] += a[i] * b[j];
            //                return c;
            //            }

            //            public static Complex EvalZinvPoly(double[] c, Complex z1)
            //            {
            //                // c[0] + c[1] z^-1 + ... with z1 = e^{-j w}
            //                Complex acc = 0;
            //                Complex zp = 1;
            //                for (int k = 0; k < c.Length; k++)
            //                {
            //                    acc += c[k] * zp;
            //                    zp *= z1;
            //                }
            //                return acc;
            //            }

            //            public static double[] BuildStableDenominatorBiquads(int order, double fs, double fLo, double fHi, double Q = 1.0)
            //            {
            //                if (order < 2) order = 2;
            //                if ((order & 1) == 1) order++; // even

            //                int nsec = order / 2;

            //                fLo = Math.Max(1.0, fLo);
            //                fHi = Math.Min(0.49 * fs, Math.Max(fLo * 1.01, fHi));

            //                double[] a = new double[] { 1.0 };

            //                for (int s = 0; s < nsec; s++)
            //                {
            //                    double t = (s + 0.5) / nsec;
            //                    double fc = fLo * Math.Pow(fHi / fLo, t);

            //                    double w0 = 2.0 * Math.PI * fc / fs;

            //                    double bw = fc / Math.Max(0.25, Q);
            //                    double r = Math.Exp(-Math.PI * bw / fs);

            //                    // keep safely inside unit circle
            //                    r = Math.Min(r, 0.995);
            //                    r = Math.Max(r, 0.20);

            //                    double c0 = Math.Cos(w0);

            //                    // A_sec(z) = 1 - 2 r cos(w0) z^-1 + r^2 z^-2
            //                    double[] sec = new double[] { 1.0, -2.0 * r * c0, r * r };
            //                    a = Convolve(a, sec);
            //                }

            //                // normalize a[0]=1
            //                double a0 = a[0];
            //                if (Math.Abs(a0 - 1.0) > 1e-12 && Math.Abs(a0) > 1e-12)
            //                {
            //                    double inv = 1.0 / a0;
            //                    for (int i = 0; i < a.Length; i++) a[i] *= inv;
            //                }

            //                return a;
            //            }

            //            public static double[] SolveNumeratorLeastSquares(double[] freqs, Complex[] Y, double[] a, double fs, double ridge = 1e-10)
            //            {
            //                // Solve for real b[] minimizing ||B - Y*A|| on sampled grid.
            //                int M = a.Length; // numerator length = denominator length
            //                int K = Math.Min(freqs.Length, Y.Length);

            //                double[,] G = new double[M, M];
            //                double[] h = new double[M];

            //                for (int i = 0; i < K; i++)
            //                {
            //                    double w = 2.0 * Math.PI * freqs[i] / fs;
            //                    Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

            //                    Complex A = EvalZinvPoly(a, z1);
            //                    Complex d = Y[i] * A; // desired B(e^jw)

            //                    // basis z^{-m}
            //                    Complex bp = 1.0;

            //                    double[] xr = new double[M];
            //                    double[] xi = new double[M];

            //                    for (int m = 0; m < M; m++)
            //                    {
            //                        xr[m] = bp.Real;
            //                        xi[m] = bp.Imaginary;
            //                        bp *= z1;
            //                    }

            //                    double yr = d.Real;
            //                    double yi = d.Imaginary;

            //                    for (int m = 0; m < M; m++)
            //                    {
            //                        h[m] += xr[m] * yr + xi[m] * yi;
            //                        for (int n = 0; n < M; n++)
            //                            G[m, n] += xr[m] * xr[n] + xi[m] * xi[n];
            //                    }
            //                }

            //                for (int m = 0; m < M; m++) G[m, m] += ridge;

            //                // Gaussian elimination (M is small: <= 13)
            //                double[,] Aaug = new double[M, M + 1];
            //                for (int i = 0; i < M; i++)
            //                {
            //                    for (int j = 0; j < M; j++) Aaug[i, j] = G[i, j];
            //                    Aaug[i, M] = h[i];
            //                }

            //                for (int col = 0; col < M; col++)
            //                {
            //                    int piv = col;
            //                    double best = Math.Abs(Aaug[col, col]);
            //                    for (int r = col + 1; r < M; r++)
            //                    {
            //                        double v = Math.Abs(Aaug[r, col]);
            //                        if (v > best) { best = v; piv = r; }
            //                    }
            //                    if (best < 1e-18) break;

            //                    if (piv != col)
            //                    {
            //                        for (int j = col; j < M + 1; j++)
            //                        {
            //                            double tmp = Aaug[col, j];
            //                            Aaug[col, j] = Aaug[piv, j];
            //                            Aaug[piv, j] = tmp;
            //                        }
            //                    }

            //                    double diag = Aaug[col, col];
            //                    for (int j = col; j < M + 1; j++) Aaug[col, j] /= diag;

            //                    for (int r = 0; r < M; r++)
            //                    {
            //                        if (r == col) continue;
            //                        double f = Aaug[r, col];
            //                        if (Math.Abs(f) < 1e-18) continue;
            //                        for (int j = col; j < M + 1; j++) Aaug[r, j] -= f * Aaug[col, j];
            //                    }
            //                }

            //                double[] b = new double[M];
            //                for (int i = 0; i < M; i++) b[i] = Aaug[i, M];
            //                return b;
            //            }

            //            public static void EnforceSampledPositiveReal(double[] freqs, double fs, double[] a, ref double[] b)
            //            {
            //                double minRe = double.PositiveInfinity;

            //                for (int i = 0; i < freqs.Length; i++)
            //                {
            //                    double w = 2.0 * Math.PI * freqs[i] / fs;
            //                    Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

            //                    Complex A = EvalZinvPoly(a, z1);
            //                    Complex B = EvalZinvPoly(b, z1);
            //                    Complex H = (A.Magnitude < 1e-18) ? Complex.Zero : (B / A);

            //                    if (H.Real < minRe) minRe = H.Real;
            //                }

            //                if (minRe < -1e-9)
            //                {
            //                    // Lift by adding a small constant conductance (passive “fuse”)
            //                    b[0] += (-minRe + 1e-9);
            //                }
            //            }

            //            public static (double[] b, double[] a) FitStableAdmittanceIIR(double[] freqs, Complex[] Y, int order, double fs, double fLo, double fHi)
            //            {
            //                double[] a = BuildStableDenominatorBiquads(order, fs, fLo, fHi, Q: 1.0);
            //                double[] b = SolveNumeratorLeastSquares(freqs, Y, a, fs, ridge: 1e-10);
            //                EnforceSampledPositiveReal(freqs, fs, a, ref b);
            //                return (b, a);
            //            }
            //        }

            //        public override (double[] a, double[] b) Estimate_IIR_Coefficients(
            //double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //        {
            //            int N = 4096;
            //            double fs = sample_frequency;

            //            // bins from 0..max_freq inclusive
            //            int K = (int)Math.Floor(max_freq * N / fs);
            //            K = Math.Max(8, Math.Min(K, N / 2 + 1));

            //            frequencies = new double[K];
            //            for (int i = 0; i < K; i++)
            //                frequencies[i] = i * fs / N;

            //            lock (lock_IIR)
            //            {
            //                // NOTE: your original cache key is only fs; if max_freq changes, you may want to include it too.
            //                if (rec_fs == fs && rec_a != null && rec_b != null && (filter_order <= 0 || rec_order == filter_order))
            //                    return (rec_a, rec_b);

            //                rec_fs = fs;

            //                // 1) Reflection spectrum R(ω), take first K bins
            //                Complex[] Rfull = this.Reflection_Spectrum((int)fs, N,
            //                    new Hare.Geometry.Vector(0, 0, 1),
            //                    new Hare.Geometry.Vector(0, 0, 1),
            //                    0);

            //                Complex[] R = new Complex[K];
            //                Array.Copy(Rfull, R, K);

            //                // 2) Target admittance Y = (1 - R) / (1 + R)
            //                Complex[] Y = new Complex[K];
            //                for (int i = 0; i < K; i++)
            //                {
            //                    Complex denom = Complex.One + R[i];
            //                    if (denom.Magnitude < 1e-12) denom = new Complex(1e-12, 0);

            //                    Complex Yi = (Complex.One - R[i]) / denom;

            //                    // keep your "no active" clamp on target
            //                    if (Yi.Real < 0) Yi = new Complex(0, Yi.Imaginary);

            //                    Y[i] = Yi;
            //                }

            //                // 3) Fit: stable poles first
            //                int min_order = 2;
            //                int max_order = 12;
            //                double fLo = 20.0;
            //                double fHi = Math.Min(max_freq * 0.95, 0.45 * fs);

            //                if (filter_order > 0)
            //                {
            //                    int order = Math.Max(2, filter_order);
            //                    if ((order & 1) == 1) order++; // even

            //                    var fit = StableAdmittanceFitter.FitStableAdmittanceIIR(frequencies, Y, order, fs, fLo, fHi);

            //                    rec_order = order;
            //                    rec_a = fit.a;
            //                    rec_b = fit.b;
            //                    return (rec_a, rec_b);
            //                }
            //                else
            //                {
            //                    double bestErr = double.PositiveInfinity;
            //                    (double[] b, double[] a) best = (null, null);
            //                    int bestOrder = min_order;

            //                    for (int order = min_order; order <= max_order; order += 2)
            //                    {
            //                        var cand = StableAdmittanceFitter.FitStableAdmittanceIIR(frequencies, Y, order, fs, fLo, fHi);

            //                        // error metric: relative complex error over 20..10k (or adapt to your needs)
            //                        double errSum = 0.0;
            //                        int count = 0;
            //                        bool passive = true;

            //                        for (int i = 0; i < K; i++)
            //                        {
            //                            double f = frequencies[i];
            //                            if (f < 20 || f > 10000) continue;

            //                            double w = 2.0 * Math.PI * f / fs;
            //                            Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

            //                            Complex A = StableAdmittanceFitter.EvalZinvPoly(cand.a, z1);
            //                            Complex B = StableAdmittanceFitter.EvalZinvPoly(cand.b, z1);
            //                            Complex H = (A.Magnitude < 1e-18) ? Complex.Zero : (B / A);

            //                            if (H.Real < -1e-6) passive = false;

            //                            double denomMag = Math.Max(1e-6, Y[i].Magnitude);
            //                            errSum += (H - Y[i]).Magnitude / denomMag;
            //                            count++;
            //                        }

            //                        double avgErr = (count > 0) ? (errSum / count) : double.PositiveInfinity;
            //                        if (!passive) avgErr += 10.0;

            //                        if (avgErr < bestErr)
            //                        {
            //                            bestErr = avgErr;
            //                            best = cand;
            //                            bestOrder = order;
            //                        }
            //                    }

            //                    rec_order = bestOrder;
            //                    rec_a = best.a;
            //                    rec_b = best.b;
            //                    return (rec_a, rec_b);
            //                }
            //            }
            //        }

            //public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //{
            //    const int N = 4096;
            //    double fs = sample_frequency;

            //    int K = (int)Math.Floor(max_freq * N / fs) + 1;
            //    if (K < 8) K = Math.Min(8, N);

            //    frequencies = new double[K];
            //    for (int i = 0; i < K; i++)
            //        frequencies[i] = i * fs / N; // includes DC

            //    // Reflection spectrum at normal incidence (existing behavior)
            //    var Rfull = this.Reflection_Spectrum((int)fs, N,
            //        new Hare.Geometry.Vector(0, 0, 1),
            //        new Hare.Geometry.Vector(0, 0, 1),
            //        0);

            //    // Target admittance: Y = (1 - R)/(1 + R)
            //    System.Numerics.Complex[] Y = new System.Numerics.Complex[K];
            //    for (int i = 0; i < K; i++)
            //    {
            //        var R = Rfull[i];
            //        var denom = System.Numerics.Complex.One + R;
            //        if (denom.Magnitude < 1e-12) denom = new System.Numerics.Complex(1e-12, 0);

            //        var Yi = (System.Numerics.Complex.One - R) / denom;

            //        // Keep your “no active” pre-clamp (fine)
            //        if (Yi.Real < 0) Yi = new System.Numerics.Complex(0, Yi.Imaginary);

            //        Y[i] = Yi;
            //    }

            //    // Auto-order selection (stable poles for each candidate)
            //    if (filter_order <= 0)
            //    {
            //        int min_order = 2;
            //        int max_order = 12;

            //        double best_error = double.MaxValue;
            //        (double[] b, double[] a) best = (null, null);
            //        int best_order = min_order;

            //        for (int order = min_order; order <= max_order; order += 2)
            //        {
            //            var cand = FitStableAdmittanceIIR(frequencies, Y, order, fs, fLo: 20.0, fHi: Math.Min(max_freq * 0.95, 0.45 * fs));

            //            // error metric (complex relative error over 125–4000 Hz)
            //            double errSum = 0.0;
            //            int count = 0;

            //            for (int i = 0; i < K; i++)
            //            {
            //                double f = frequencies[i];
            //                if (f < 125 || f > 4000) continue;

            //                double w = 2.0 * Math.PI * f / fs;
            //                var z1 = System.Numerics.Complex.Exp(-System.Numerics.Complex.ImaginaryOne * w);

            //                var A = EvalZinvPoly(cand.a, z1);
            //                var B = EvalZinvPoly(cand.b, z1);
            //                var H = (A.Magnitude < 1e-18) ? System.Numerics.Complex.Zero : (B / A);

            //                double denomMag = Math.Max(1e-3, Y[i].Magnitude);
            //                double e = (H - Y[i]).Magnitude / denomMag;

            //                errSum += e;
            //                count++;
            //            }

            //            double avgErr = (count > 0) ? errSum / count : double.MaxValue;
            //            if (avgErr < best_error)
            //            {
            //                best_error = avgErr;
            //                best = cand;
            //                best_order = order;
            //            }
            //        }

            //        rec_order = best_order;
            //        rec_a = best.a;
            //        rec_b = best.b;
            //        return (rec_a, rec_b);
            //    }
            //    else
            //    {
            //        int order = Math.Max(2, filter_order);
            //        if ((order & 1) == 1) order++; // keep even

            //        var fit = FitStableAdmittanceIIR(frequencies, Y, order, fs, fLo: 20.0, fHi: Math.Min(max_freq * 0.95, 0.45 * fs));
            //        return (fit.a, fit.b);
            //    }
            //}


            ////public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            ////{
            ////    int N = 4096;
            ////    double fs = sample_frequency;

            ////    // bins from 0..max_freq inclusive
            ////    int K = (int)Math.Floor(max_freq * N / fs);
            ////    K = Math.Max(8, Math.Min(K, N / 2 + 1));

            ////    frequencies = new double[K];
            ////    for (int i = 0; i < K; i++)
            ////        frequencies[i] = i * fs / N;

            ////    lock (lock_IIR)
            ////    {
            ////        // NOTE: your original cache key is only fs; if max_freq changes, you may want to include it too.
            ////        if (rec_fs == fs && rec_a != null && rec_b != null && (filter_order <= 0 || rec_order == filter_order))
            ////            return (rec_a, rec_b);

            ////        rec_fs = fs;

            ////        // 1) Reflection spectrum R(ω), take first K bins
            ////        Complex[] Rfull = this.Reflection_Spectrum((int)fs, N,
            ////            new Hare.Geometry.Vector(0, 0, 1),
            ////            new Hare.Geometry.Vector(0, 0, 1),
            ////            0);

            ////        Complex[] R = new Complex[K];
            ////        Array.Copy(Rfull, R, K);

            ////        // 2) Target admittance Y = (1 - R) / (1 + R)
            ////        Complex[] Y = new Complex[K];
            ////        for (int i = 0; i < K; i++)
            ////        {
            ////            Complex denom = Complex.One + R[i];
            ////            if (denom.Magnitude < 1e-12) denom = new Complex(1e-12, 0);

            ////            Complex Yi = (Complex.One - R[i]) / denom;

            ////            // keep your "no active" clamp on target
            ////            if (Yi.Real < 0) Yi = new Complex(0, Yi.Imaginary);

            ////            Y[i] = Yi;
            ////        }

            ////        // 3) Fit: stable poles first
            ////        int min_order = 2;
            ////        int max_order = 12;
            ////        double fLo = 20.0;
            ////        double fHi = Math.Min(max_freq * 0.95, 0.45 * fs);

            ////        if (filter_order > 0)
            ////        {
            ////            int order = Math.Max(2, filter_order);
            ////            if ((order & 1) == 1) order++; // even

            ////            var fit = StableAdmittanceFitter.FitStableAdmittanceIIR(frequencies, Y, order, fs, fLo, fHi);

            ////            rec_order = order;
            ////            rec_a = fit.a;
            ////            rec_b = fit.b;
            ////            return (rec_a, rec_b);
            ////        }
            ////        else
            ////        {
            ////            double bestErr = double.PositiveInfinity;
            ////            (double[] b, double[] a) best = (null, null);
            ////            int bestOrder = min_order;

            ////            for (int order = min_order; order <= max_order; order += 2)
            ////            {
            ////                var cand = StableAdmittanceFitter.FitStableAdmittanceIIR(frequencies, Y, order, fs, fLo, fHi);

            ////                // error metric: relative complex error over 20..10k (or adapt to your needs)
            ////                double errSum = 0.0;
            ////                int count = 0;
            ////                bool passive = true;

            ////                for (int i = 0; i < K; i++)
            ////                {
            ////                    double f = frequencies[i];
            ////                    if (f < 20 || f > 10000) continue;

            ////                    double w = 2.0 * Math.PI * f / fs;
            ////                    Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

            ////                    Complex A = StableAdmittanceFitter.EvalZinvPoly(cand.a, z1);
            ////                    Complex B = StableAdmittanceFitter.EvalZinvPoly(cand.b, z1);
            ////                    Complex H = (A.Magnitude < 1e-18) ? Complex.Zero : (B / A);

            ////                    if (H.Real < -1e-6) passive = false;

            ////                    double denomMag = Math.Max(1e-6, Y[i].Magnitude);
            ////                    errSum += (H - Y[i]).Magnitude / denomMag;
            ////                    count++;
            ////                }

            ////                double avgErr = (count > 0) ? (errSum / count) : double.PositiveInfinity;
            ////                if (!passive) avgErr += 10.0;

            ////                if (avgErr < bestErr)
            ////                {
            ////                    bestErr = avgErr;
            ////                    best = cand;
            ////                    bestOrder = order;
            ////                }
            ////            }

            ////            rec_order = bestOrder;
            ////            rec_a = best.a;
            ////            rec_b = best.b;
            ////            return (rec_a, rec_b);
            ////        }
            ////    }
            ////}

            ////internal static class StableAdmittanceFitter
            ////{
            ////    public static double[] Convolve(double[] a, double[] b)
            ////    {
            ////        double[] c = new double[a.Length + b.Length - 1];
            ////        for (int i = 0; i < a.Length; i++)
            ////            for (int j = 0; j < b.Length; j++)
            ////                c[i + j] += a[i] * b[j];
            ////        return c;
            ////    }

            ////    public static Complex EvalZinvPoly(double[] c, Complex z1)
            ////    {
            ////        // c[0] + c[1] z^-1 + ... with z1 = e^{-j w}
            ////        Complex acc = 0;
            ////        Complex zp = 1;
            ////        for (int k = 0; k < c.Length; k++)
            ////        {
            ////            acc += c[k] * zp;
            ////            zp *= z1;
            ////        }
            ////        return acc;
            ////    }

            ////    public static double[] BuildStableDenominatorBiquads(int order, double fs, double fLo, double fHi, double Q = 1.0)
            ////    {
            ////        if (order < 2) order = 2;
            ////        if ((order & 1) == 1) order++; // even

            ////        int nsec = order / 2;

            ////        fLo = Math.Max(1.0, fLo);
            ////        fHi = Math.Min(0.49 * fs, Math.Max(fLo * 1.01, fHi));

            ////        double[] a = new double[] { 1.0 };

            ////        for (int s = 0; s < nsec; s++)
            ////        {
            ////            double t = (s + 0.5) / nsec;
            ////            double fc = fLo * Math.Pow(fHi / fLo, t);

            ////            double w0 = 2.0 * Math.PI * fc / fs;

            ////            double bw = fc / Math.Max(0.25, Q);
            ////            double r = Math.Exp(-Math.PI * bw / fs);

            ////            // keep safely inside unit circle
            ////            r = Math.Min(r, 0.995);
            ////            r = Math.Max(r, 0.20);

            ////            double c0 = Math.Cos(w0);

            ////            // A_sec(z) = 1 - 2 r cos(w0) z^-1 + r^2 z^-2
            ////            double[] sec = new double[] { 1.0, -2.0 * r * c0, r * r };
            ////            a = Convolve(a, sec);
            ////        }

            ////        // normalize a[0]=1
            ////        double a0 = a[0];
            ////        if (Math.Abs(a0 - 1.0) > 1e-12 && Math.Abs(a0) > 1e-12)
            ////        {
            ////            double inv = 1.0 / a0;
            ////            for (int i = 0; i < a.Length; i++) a[i] *= inv;
            ////        }

            ////        return a;
            ////    }

            ////    public static double[] SolveNumeratorLeastSquares(double[] freqs, Complex[] Y, double[] a, double fs, double ridge = 1e-10)
            ////    {
            ////        // Solve for real b[] minimizing ||B - Y*A|| on sampled grid.
            ////        int M = a.Length; // numerator length = denominator length
            ////        int K = Math.Min(freqs.Length, Y.Length);

            ////        double[,] G = new double[M, M];
            ////        double[] h = new double[M];

            ////        for (int i = 0; i < K; i++)
            ////        {
            ////            double w = 2.0 * Math.PI * freqs[i] / fs;
            ////            Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

            ////            Complex A = EvalZinvPoly(a, z1);
            ////            Complex d = Y[i] * A; // desired B(e^jw)

            ////            // basis z^{-m}
            ////            Complex bp = 1.0;

            ////            double[] xr = new double[M];
            ////            double[] xi = new double[M];

            ////            for (int m = 0; m < M; m++)
            ////            {
            ////                xr[m] = bp.Real;
            ////                xi[m] = bp.Imaginary;
            ////                bp *= z1;
            ////            }

            ////            double yr = d.Real;
            ////            double yi = d.Imaginary;

            ////            for (int m = 0; m < M; m++)
            ////            {
            ////                h[m] += xr[m] * yr + xi[m] * yi;
            ////                for (int n = 0; n < M; n++)
            ////                    G[m, n] += xr[m] * xr[n] + xi[m] * xi[n];
            ////            }
            ////        }

            ////        for (int m = 0; m < M; m++) G[m, m] += ridge;

            ////        // Gaussian elimination (M is small: <= 13)
            ////        double[,] Aaug = new double[M, M + 1];
            ////        for (int i = 0; i < M; i++)
            ////        {
            ////            for (int j = 0; j < M; j++) Aaug[i, j] = G[i, j];
            ////            Aaug[i, M] = h[i];
            ////        }

            ////        for (int col = 0; col < M; col++)
            ////        {
            ////            int piv = col;
            ////            double best = Math.Abs(Aaug[col, col]);
            ////            for (int r = col + 1; r < M; r++)
            ////            {
            ////                double v = Math.Abs(Aaug[r, col]);
            ////                if (v > best) { best = v; piv = r; }
            ////            }
            ////            if (best < 1e-18) break;

            ////            if (piv != col)
            ////            {
            ////                for (int j = col; j < M + 1; j++)
            ////                {
            ////                    double tmp = Aaug[col, j];
            ////                    Aaug[col, j] = Aaug[piv, j];
            ////                    Aaug[piv, j] = tmp;
            ////                }
            ////            }

            ////            double diag = Aaug[col, col];
            ////            for (int j = col; j < M + 1; j++) Aaug[col, j] /= diag;

            ////            for (int r = 0; r < M; r++)
            ////            {
            ////                if (r == col) continue;
            ////                double f = Aaug[r, col];
            ////                if (Math.Abs(f) < 1e-18) continue;
            ////                for (int j = col; j < M + 1; j++) Aaug[r, j] -= f * Aaug[col, j];
            ////            }
            ////        }

            ////        double[] b = new double[M];
            ////        for (int i = 0; i < M; i++) b[i] = Aaug[i, M];
            ////        return b;
            ////    }

            ////    public static void EnforceSampledPositiveReal(double[] freqs, double fs, double[] a, ref double[] b)
            ////    {
            ////        double minRe = double.PositiveInfinity;

            ////        for (int i = 0; i < freqs.Length; i++)
            ////        {
            ////            double w = 2.0 * Math.PI * freqs[i] / fs;
            ////            Complex z1 = Complex.Exp(-Complex.ImaginaryOne * w);

            ////            Complex A = EvalZinvPoly(a, z1);
            ////            Complex B = EvalZinvPoly(b, z1);
            ////            Complex H = (A.Magnitude < 1e-18) ? Complex.Zero : (B / A);

            ////            if (H.Real < minRe) minRe = H.Real;
            ////        }

            ////        if (minRe < -1e-9)
            ////        {
            ////            // Lift by adding a small constant conductance (passive “fuse”)
            ////            b[0] += (-minRe + 1e-9);
            ////        }
            ////    }

            ////    public static (double[] b, double[] a) FitStableAdmittanceIIR(double[] freqs, Complex[] Y, int order, double fs, double fLo, double fHi)
            ////    {
            ////        double[] a = BuildStableDenominatorBiquads(order, fs, fLo, fHi, Q: 1.0);
            ////        double[] b = SolveNumeratorLeastSquares(freqs, Y, a, fs, ridge: 1e-10);
            ////        EnforceSampledPositiveReal(freqs, fs, a, ref b);
            ////        return (b, a);
            ////    }
            ////}

            //public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //{
            //    // Number of samples for the IIR filter design
            //    int samplect = 4096;
            //    frequencies = new double[(int)(samplect * max_freq / sample_frequency)];

            //    lock (lock_IIR)
            //    {
            //        if (rec_fs == sample_frequency) return (rec_a, rec_b);

            //        rec_fs = sample_frequency;

            //        // Get reflection spectrum and calculate desired absorption coefficients
            //        Complex[] desired_Spectrum = new Complex[(int)(samplect * max_freq / sample_frequency)];
            //        Array.Copy(this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0), desired_Spectrum, desired_Spectrum.Length / 2);

            //        for (int i = 0; i < frequencies.Length; i++)
            //        {
            //            frequencies[i] = (double)((i + 1) * (sample_frequency / 2)) / samplect;
            //        }

            //        double[] desired_alpha = AbsorptionModels.Operations.Absorption_Coef(desired_Spectrum);

            //        // If filter_order is zero or negative, automatically determine the optimal order
            //        if (filter_order <= 0)
            //        {
            //            if (rec_order == 0)
            //            {
            //                int min_order = 2;
            //                int max_order = 12;
            //                double error_threshold = 0.01; // 1% average error threshold
            //                double best_error = double.MaxValue;
            //                (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);

            //                for (int order = min_order; order <= max_order; order += 2) // Try even orders
            //                {
            //                    try
            //                    {
            //                        // Design filter with current order
            //                        (double[] b, double[] a) current_filter = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
            //                            frequencies, desired_Spectrum, order, order, sample_frequency);

            //                        // Normalize the filter
            //                        double mag = 0;
            //                        Complex[] filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                            new List<double>(current_filter.b), new List<double>(current_filter.a),
            //                            frequencies, sample_frequency);

            //                        mag = filter_response.MaximumMagnitudePhase().Magnitude;
            //                        if (mag > 1)
            //                        {
            //                            for (int j = 0; j < current_filter.b.Length; j++)
            //                                current_filter.b[j] /= mag;
            //                        }
            //                        else
            //                        {
            //                            double mod = mag / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
            //                            for (int j = 0; j < current_filter.b.Length; j++)
            //                                current_filter.b[j] /= mod;
            //                        }

            //                        // Recalculate response after normalization
            //                        filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                            new List<double>(current_filter.b), new List<double>(current_filter.a),
            //                            frequencies, sample_frequency);

            //                        // Calculate error metric (average magnitude error)
            //                        double error_sum = 0;
            //                        int count = 0;

            //                        // Evaluate from 20 Hz. to 10,000 Hz. (as applicable)
            //                        for (int i = 0; i < frequencies.Length; i++)
            //                        {
            //                            if (frequencies[i] >= 20 && frequencies[i] <= 10000)
            //                            {
            //                                double mag_error = Math.Abs(filter_response[i].Magnitude - desired_Spectrum[i].Magnitude);
            //                                error_sum += mag_error / Math.Max(0.01, desired_Spectrum[i].Magnitude); // Relative error
            //                                count++;
            //                            }
            //                        }

            //                        double avg_error = count > 0 ? error_sum / count : double.MaxValue;

            //                        // Update best filter if this one is better
            //                        if (avg_error < best_error)
            //                        {
            //                            best_error = avg_error;
            //                            best_filter = current_filter;
            //                        }

            //                        // If error is below threshold, we've found a good enough match
            //                        if (avg_error < error_threshold)
            //                            break;
            //                    }
            //                    catch
            //                    {
            //                        // If filter design fails, continue to next order
            //                        continue;
            //                    }
            //                }

            //                // Use the best filter found
            //                rec_order = best_filter.a.Length - 1;
            //                rec_a = best_filter.a;
            //                rec_b = best_filter.b;
            //                return best_filter;
            //            }
            //            else
            //            {
            //                return (rec_a, rec_b); // Return previously calculated coefficients
            //            }
            //        }
            //        else
            //        {
            //            // Original code for specified filter order
            //            int numeratorOrder = filter_order;
            //            int denominatorOrder = filter_order;

            //            try
            //            {
            //                // Use the IIR filter design method from Pach_SP.IIR_Design
            //                (rec_b, rec_a) = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
            //                    frequencies, desired_Spectrum, numeratorOrder, denominatorOrder, sample_frequency);

            //                // Calculate frequency response of the fitted filter
            //                Complex[] IIR_spec = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                    new List<double>(rec_b), new List<double>(rec_a), frequencies, sample_frequency);

            //                // Normalize the filter if necessary
            //                double m = IIR_spec.MaximumMagnitudePhase().Magnitude;
            //                if (m > 1)
            //                {
            //                    for (int j = 0; j < rec_b.Length; j++) rec_b[j] /= m;
            //                }
            //                else
            //                {
            //                    double mod = m / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
            //                    for (int j = 0; j < rec_b.Length; j++) rec_b[j] /= mod;
            //                }

            //                return (rec_a, rec_b);
            //            }
            //            catch (Exception)
            //            {
            //                // If the fitting process fails, fall back to a simple filter
            //                rec_a = new double[filter_order + 1];
            //                rec_b = new double[filter_order + 1];

            //                // Create a basic filter based on octave band coefficients
            //                rec_a[0] = 1.0;
            //                for (int i = 1; i < rec_a.Length; i++) rec_a[i] = -0.7 * Math.Pow(0.8, i - 1);
            //                rec_b[0] = 1.0;

            //                return (rec_a, rec_b);
            //            }
            //        }
            //    }
            //}
            //    public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //{
            //    int samplect = 4096;
            //    frequencies = new double[(int)(samplect * max_freq / sample_frequency)];

            //    lock (lock_IIR)
            //    {
            //        if (rec_fs == sample_frequency) return (rec_a, rec_b);

            //        rec_fs = sample_frequency;

            //        // 1) Get reflection spectrum R(ω)
            //        Complex[] R = new Complex[(int)(samplect * max_freq / sample_frequency)];
            //        Array.Copy(
            //            this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0),
            //            R,
            //            R.Length
            //        );

            //        for (int i = 0; i < frequencies.Length; i++)
            //        {
            //            frequencies[i] = ((i + 1) * (sample_frequency / 2.0)) / samplect;
            //        }

            //        // 2) Convert to admittance spectrum Y = (1 - R) / (1 + R)
            //        Complex[] Y = new Complex[R.Length];
            //        for (int i = 0; i < R.Length; i++)
            //        {
            //            Complex denom = 1 + R[i];
            //            if (denom.Magnitude < 1e-9) denom = new Complex(1e-9, 0); // avoid div-by-zero
            //            Y[i] = (1 - R[i]) / denom;
            //        }

            //        // 3) Fit IIR to Y(ω)
            //        if (filter_order <= 0)
            //        {
            //            if (rec_order == 0)
            //            {
            //                int min_order = 2;
            //                int max_order = 12;
            //                double error_threshold = 0.01;
            //                double best_error = double.MaxValue;
            //                (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);

            //                for (int order = min_order; order <= max_order; order += 2)
            //                {
            //                    try
            //                    {
            //                        (double[] b, double[] a) current = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
            //                            frequencies, Y, order, order, sample_frequency);

            //                        // Normalize to a[0] = 1 to keep physical scaling consistent
            //                        double a0 = current.a[0];
            //                        if (Math.Abs(a0) > 1e-12)
            //                        {
            //                            for (int j = 0; j < current.a.Length; j++) current.a[j] /= a0;
            //                            for (int j = 0; j < current.b.Length; j++) current.b[j] /= a0;
            //                        }

            //                        Complex[] resp = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                            new List<double>(current.b), new List<double>(current.a), frequencies, sample_frequency);

            //                        double errSum = 0;
            //                        int count = 0;
            //                        for (int i = 0; i < frequencies.Length; i++)
            //                        {
            //                            if (frequencies[i] >= 20 && frequencies[i] <= 10000)
            //                            {
            //                                double magErr = Math.Abs(resp[i].Magnitude - Y[i].Magnitude);
            //                                errSum += magErr / Math.Max(1e-6, Y[i].Magnitude);
            //                                count++;
            //                            }
            //                        }
            //                        double avgErr = count > 0 ? errSum / count : double.MaxValue;

            //                        if (avgErr < best_error)
            //                        {
            //                            best_error = avgErr;
            //                            best_filter = current;
            //                        }
            //                        if (avgErr < error_threshold) break;
            //                    }
            //                    catch
            //                    {
            //                        continue;
            //                    }
            //                }

            //                rec_order = best_filter.a.Length - 1;
            //                rec_a = best_filter.a;
            //                rec_b = best_filter.b;
            //                return best_filter;
            //            }
            //            else
            //            {
            //                return (rec_a, rec_b);
            //            }
            //        }
            //        else
            //        {
            //            int numeratorOrder = filter_order;
            //            int denominatorOrder = filter_order;

            //            try
            //            {
            //                (rec_b, rec_a) = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
            //                    frequencies, Y, numeratorOrder, denominatorOrder, sample_frequency);

            //                // Normalize to a[0] = 1
            //                double a0 = rec_a[0];
            //                if (Math.Abs(a0) > 1e-12)
            //                {
            //                    for (int j = 0; j < rec_a.Length; j++) rec_a[j] /= a0;
            //                    for (int j = 0; j < rec_b.Length; j++) rec_b[j] /= a0;
            //                }

            //                return (rec_a, rec_b);
            //            }
            //            catch (Exception)
            //            {
            //                rec_a = new double[filter_order + 1];
            //                rec_b = new double[filter_order + 1];
            //                rec_a[0] = 1.0;
            //                for (int i = 1; i < rec_a.Length; i++) rec_a[i] = -0.7 * Math.Pow(0.8, i - 1);
            //                rec_b[0] = 1.0;
            //                return (rec_a, rec_b);
            //            }
            //        }
            //    }
            //}
            //public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //{
            //    int samplect = 4096;
            //    frequencies = new double[(int)(samplect * max_freq / sample_frequency)];

            //    lock (lock_IIR)
            //    {
            //        if (rec_fs == sample_frequency) return (rec_a, rec_b);

            //        rec_fs = sample_frequency;

            //        // 1) Get reflection spectrum R(ω)
            //        Complex[] R = new Complex[(int)(samplect * max_freq / sample_frequency)];
            //        Array.Copy(
            //            this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0),
            //            R,
            //            R.Length
            //        );

            //        for (int i = 0; i < frequencies.Length; i++)
            //        {
            //            frequencies[i] = ((i + 1) * (sample_frequency / 2.0)) / samplect;
            //        }

            //        // 2) Convert to admittance spectrum Y = (1 - R) / (1 + R)
            //        Complex[] Y = new Complex[R.Length];
            //        for (int i = 0; i < R.Length; i++)
            //        {
            //            Complex denom = 1 + R[i];
            //            if (denom.Magnitude < 1e-9) denom = new Complex(1e-9, 0);
            //            Y[i] = (1 - R[i]) / denom;

            //            // Passivity enforcement: real part of admittance must be >= 0
            //            // (negative real part means energy injection)
            //            if (Y[i].Real < 0)
            //                Y[i] = new Complex(0, Y[i].Imaginary);
            //        }

            //        // 3) Fit IIR to Y(ω)
            //        if (filter_order <= 0)
            //        {
            //            if (rec_order == 0)
            //            {
            //                int min_order = 2;
            //                int max_order = 12;
            //                double error_threshold = 0.01;
            //                double best_error = double.MaxValue;
            //                (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);

            //                for (int order = min_order; order <= max_order; order += 2)
            //                {
            //                    try
            //                    {
            //                        (double[] b, double[] a) current = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
            //                            frequencies, Y, order, order, sample_frequency);

            //                        // Normalize to a[0] = 1
            //                        double a0 = current.a[0];
            //                        if (Math.Abs(a0) > 1e-12)
            //                        {
            //                            for (int j = 0; j < current.a.Length; j++) current.a[j] /= a0;
            //                            for (int j = 0; j < current.b.Length; j++) current.b[j] /= a0;
            //                        }

            //                        // Verify passivity: check that the fitted filter doesn't go negative
            //                        Complex[] resp = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                            new List<double>(current.b), new List<double>(current.a), frequencies, sample_frequency);

            //                        bool passive = true;
            //                        double errSum = 0;
            //                        int count = 0;
            //                        for (int i = 0; i < frequencies.Length; i++)
            //                        {
            //                            // Check passivity: Re{Y} >= 0 at all frequencies
            //                            if (resp[i].Real < -1e-6)
            //                                passive = false;

            //                            if (frequencies[i] >= 20 && frequencies[i] <= 10000)
            //                            {
            //                                double magErr = Math.Abs(resp[i].Magnitude - Y[i].Magnitude);
            //                                errSum += magErr / Math.Max(1e-6, Y[i].Magnitude);
            //                                count++;
            //                            }
            //                        }
            //                        double avgErr = count > 0 ? errSum / count : double.MaxValue;

            //                        // Penalize non-passive filters heavily
            //                        if (!passive) avgErr += 10.0;

            //                        if (avgErr < best_error)
            //                        {
            //                            best_error = avgErr;
            //                            best_filter = current;
            //                        }
            //                        if (avgErr < error_threshold) break;
            //                    }
            //                    catch
            //                    {
            //                        continue;
            //                    }
            //                }

            //                // Final passivity check on best filter — if still non-passive,
            //                // clamp b[0] to ensure DC admittance is non-negative
            //                Complex[] bestResp = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                    new List<double>(best_filter.b), new List<double>(best_filter.a), frequencies, sample_frequency);

            //                bool bestPassive = true;
            //                for (int i = 0; i < bestResp.Length; i++)
            //                {
            //                    if (bestResp[i].Real < -1e-6) { bestPassive = false; break; }
            //                }

            //                if (!bestPassive)
            //                {
            //                    // Scale b coefficients to ensure max magnitude doesn't exceed
            //                    // the target admittance max magnitude (conservative passive bound)
            //                    double maxTarget = 0;
            //                    double maxFit = 0;
            //                    for (int i = 0; i < frequencies.Length; i++)
            //                    {
            //                        if (Y[i].Magnitude > maxTarget) maxTarget = Y[i].Magnitude;
            //                        if (bestResp[i].Magnitude > maxFit) maxFit = bestResp[i].Magnitude;
            //                    }
            //                    if (maxFit > 1e-12)
            //                    {
            //                        double scale = maxTarget / maxFit;
            //                        for (int j = 0; j < best_filter.b.Length; j++)
            //                            best_filter.b[j] *= scale;
            //                    }
            //                }

            //                rec_order = best_filter.a.Length - 1;
            //                rec_a = best_filter.a;
            //                rec_b = best_filter.b;
            //                return best_filter;
            //            }
            //            else
            //            {
            //                return (rec_a, rec_b);
            //            }
            //        }
            //        else
            //        {
            //            int numeratorOrder = filter_order;
            //            int denominatorOrder = filter_order;

            //            try
            //            {
            //                (rec_b, rec_a) = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
            //                    frequencies, Y, numeratorOrder, denominatorOrder, sample_frequency);

            //                double a0 = rec_a[0];
            //                if (Math.Abs(a0) > 1e-12)
            //                {
            //                    for (int j = 0; j < rec_a.Length; j++) rec_a[j] /= a0;
            //                    for (int j = 0; j < rec_b.Length; j++) rec_b[j] /= a0;
            //                }

            //                return (rec_a, rec_b);
            //            }
            //            catch (Exception)
            //            {
            //                rec_a = new double[filter_order + 1];
            //                rec_b = new double[filter_order + 1];
            //                rec_a[0] = 1.0;
            //                for (int i = 1; i < rec_a.Length; i++) rec_a[i] = -0.7 * Math.Pow(0.8, i - 1);
            //                rec_b[0] = 1.0;
            //                return (rec_a, rec_b);
            //            }
            //        }
            //    }
            //}
        }


        //public class Finite_Material : Material
        //{
        //    double[][][] alpha;
        //    double[] Azimuth;
        //    double[] Altitude;

        //    Smart_Material Inf_Mat;

        //    public Finite_Material(Smart_Material Mat, Rhino.Geometry.Brep Br, Rhino.Geometry.Mesh M, int face_id, Medium_Properties med)
        //    {
        //        //Strictly for the flat X,Y case - oversimplified for now.
        //        Inf_Mat = Mat;
        //        Azimuth = new double[36];
        //        Altitude = new double[Mat.Angles.Length/2];
        //        alpha = new double[Altitude.Length][][];
        //        for(int i = 0; i < Altitude.Length; i++) Altitude[i] = Mat.Angles[i].Magnitude;
        //        for(int i = 0; i < Azimuth.Length; i++) Azimuth[i] = i * 360f / Azimuth.Length;

        //            //Set up a frequency interpolated Zr for each direction individually.
        //            Rhino.Geometry.Point3d pt = M.Faces.GetFaceCenter(face_id);

        //            double[][][] ZrR = new double[Altitude.Length][][], ZrI = new double[Altitude.Length][][];
        //            double[] fr = new double[9];
        //            for (int k = 0; k < Altitude.Length; k++)
        //            {
        //                ZrR[k] = new double[Azimuth.Length][];
        //                ZrI[k] = new double[Azimuth.Length][];
        //                alpha[k] = new double[Azimuth.Length][];
        //                for (int j = 0; j < Azimuth.Length; j++)
        //                {
        //                    ZrR[k][j] = new double[9];
        //                    ZrI[k][j] = new double[9];
        //                    alpha[k][j] = new double[8];
        //                }
        //            }

        //            for (int oct = 0; oct < 9; oct++)
        //            {
        //                fr[oct] = 62.5 * Math.Pow(2, oct) / Utilities.Numerics.rt2;
        //                System.Numerics.Complex[][] Zr = AbsorptionModels.Operations.Finite_Radiation_Impedance_Rect_Longhand(pt.X, pt.Y, Br, fr[oct], Altitude, Azimuth, med.Sound_Speed(pt));

        //                for (int k = 0; k < Zr.Length; k++)
        //                {
        //                    for (int j = 0; j < Zr[k].Length; j++)
        //                    {
        //                        ZrR[k][j][oct] = Zr[k][j].Real;
        //                        ZrI[k][j][oct] = Zr[k][j].Imaginary;
        //                    }
        //                }
        //            }

        //            MathNet.Numerics.Interpolation.CubicSpline[][] Zr_r = new MathNet.Numerics.Interpolation.CubicSpline[Altitude.Length][];
        //            MathNet.Numerics.Interpolation.CubicSpline[][] Zr_i = new MathNet.Numerics.Interpolation.CubicSpline[Altitude.Length][];

        //            for (int k = 0; k < Zr_r.Length; k++)
        //            {
        //                Zr_r[k] = new MathNet.Numerics.Interpolation.CubicSpline[Azimuth.Length];
        //                Zr_i[k] = new MathNet.Numerics.Interpolation.CubicSpline[Azimuth.Length];
        //                for (int j = 0; j < Zr_r[k].Length; j++)
        //                {
        //                    //Interpolate over curve real and imaginary Zr here...
        //                    Zr_r[k][j] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(fr, ZrR[k][j]);
        //                    Zr_i[k][j] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(fr, ZrI[k][j]);
        //                }
        //            }

        //            for (int k = 0; k < Zr_r.Length; k++)
        //            {
        //                for (int j = 0; j < Zr_r[k].Length; j++)
        //                {
        //                    List<double> freq = new List<double>();
        //                    List<double> alpha_interp = new List<double>();
        //                    for (int l = 0; l < Mat.frequency.Length; l++)
        //                    {
        //                        if (Mat.frequency[l] > 10000) break;
        //                        freq.Add(Mat.frequency[l]);
        //                        alpha_interp.Add(AbsorptionModels.Operations.Finite_Unit_Absorption_Coefficient(Mat.Z[k][j], new System.Numerics.Complex(Zr_r[k][j].Interpolate(Mat.frequency[l]), Zr_i[k][j].Interpolate(Mat.frequency[l])), med.Rho(pt), med.Sound_Speed(pt)));
        //                    }
        //                    MathNet.Numerics.Interpolation.CubicSpline a = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(freq, alpha_interp);
        //                    for (int oct = 0; oct < 8; oct++)
        //                    {
        //                        alpha[k][j][oct] = 1 - a.Integrate(fr[oct], fr[oct + 1]) / (fr[oct + 1] - fr[oct]);
        //                    }
        //                }
        //            }
        //    }

        //    public override void Absorb(ref BroadRay Ray, Hare.Geometry.Vector Normal)
        //    {
        //        //Simplified for sample laid on floor...
        //        Ray.direction.Normalize();
        //        int Alt = (int)Math.Floor((Math.Acos(Hare.Geometry.Hare_math.Dot(Ray.direction, new Hare.Geometry.Vector(0, 0, -1))) * Altitude.Length) / (Math.PI / 2));
        //        if (Alt >= Altitude.Length / 2) Alt = Altitude.Length - 1;
        //        int Azi = (int)Math.Round((Math.Atan2(Ray.direction.y, Ray.direction.x) * Azimuth.Length) / Math.PI);
        //        if (Ray.direction.y < 0 && Azi < Azimuth.Length / 2) Azi = Azimuth.Length - Math.Abs(Azi);
        //        for (int oct = 0; oct < 8; oct++) Ray.Energy[oct] *= alpha[Alt][Azi][oct];
        //    }

        //    public override void Absorb(ref BroadRay Ray, out double cos_theta, Hare.Geometry.Vector Normal)
        //    {
        //        Ray.direction.Normalize();
        //        cos_theta = Math.Acos(Hare.Geometry.Hare_math.Dot(Ray.direction, new Hare.Geometry.Vector(0, 0, -1)));
        //        int Alt = (int)Math.Floor((cos_theta * Altitude.Length) / (Math.PI / 2));
        //        if (Alt >= Altitude.Length / 2) Alt = Altitude.Length - 1;
        //        int Azi = (int)Math.Round((Math.Atan2(Ray.direction.y, Ray.direction.x) * Azimuth.Length) / Math.PI);
        //        if (Ray.direction.y < 0 && Azi < Azimuth.Length / 2) Azi = Azimuth.Length - Math.Abs(Azi);

        //        for (int oct = 0; oct < 8; oct++) Ray.Energy[oct] *= alpha[Alt][Azi][oct];                
        //    }

        //    public override void Absorb(ref OctaveRay Ray, Hare.Geometry.Vector Normal)
        //    {
        //        Ray.direction.Normalize();
        //        int Alt = (int)Math.Floor((Math.Acos(Hare.Geometry.Hare_math.Dot(Ray.direction, new Hare.Geometry.Vector(0, 0, -1))) * Altitude.Length) / (Math.PI / 2));
        //        if (Alt >= Altitude.Length / 2) Alt = Altitude.Length - 1;
        //        int Azi = (int)Math.Round((Math.Atan2(Ray.direction.y, Ray.direction.x) * Azimuth.Length) / Math.PI);
        //        if (Ray.direction.y < 0 && Azi < Azimuth.Length / 2) Azi = Azimuth.Length - Math.Abs(Azi);
        //        Ray.Intensity *= alpha[Alt][Azi][Ray.Octave];
        //    }

        //    public override void Absorb(ref OctaveRay Ray, out double cos_theta, Hare.Geometry.Vector Normal)
        //    {
        //        Ray.direction.Normalize();
        //        cos_theta = Math.Acos(Hare.Geometry.Hare_math.Dot(Ray.direction, new Hare.Geometry.Vector(0, 0, -1)));
        //        int Alt = (int)Math.Floor((cos_theta * Altitude.Length) / (Math.PI / 2));
        //        if (Alt >= Altitude.Length / 2) Alt = Altitude.Length - 1;
        //        int Azi = (int)Math.Round((Math.Atan2(Ray.direction.y, Ray.direction.x) * Azimuth.Length) / Math.PI / 2);
        //        if (Ray.direction.y < 0 && Azi < Azimuth.Length / 2) Azi = Azimuth.Length - Math.Abs(Azi);
        //        if (Azi == Azimuth.Length) Azi = 0;
        //        Ray.Intensity *= alpha[Alt][Azi][Ray.Octave];
        //    }

        //    public override double[] Coefficient_A_Broad()
        //    {
        //        return Inf_Mat.Coefficient_A_Broad();
        //    }

        //    public override double Coefficient_A_Broad(int Octave)
        //    {
        //        return Inf_Mat.Coefficient_A_Broad(Octave);
        //    }

        //    public override System.Numerics.Complex Reflection_Narrow(double frequency)
        //    {
        //        return Inf_Mat.Reflection_Narrow(frequency);
        //    }

        //    public override System.Numerics.Complex Reflection_Narrow(double frequency, Hare.Geometry.Vector Dir, Hare.Geometry.Vector Normal)
        //    {
        //        return Inf_Mat.Reflection_Narrow(frequency, Dir, Normal);
        //    }

        //    public override System.Numerics.Complex[] Reflection_Spectrum(int sample_frequency, int length, Hare.Geometry.Vector Normal, Hare.Geometry.Vector Dir, int threadid)
        //    {
        //        return Inf_Mat.Reflection_Spectrum(sample_frequency, length, Normal, Dir, threadid);
        //    }
        //}

        public class Smart_Material : Material
        {
            List<AbsorptionModels.ABS_Layer> Buildup;
            int Fs;
            double rho;
            double c;
            public double[] frequency = null;
            public System.Numerics.Complex[] Angles = null;
            public System.Numerics.Complex[][] Z;
            public MathNet.Numerics.Interpolation.CubicSpline[] Transfer_FunctionR;
            public MathNet.Numerics.Interpolation.CubicSpline[] Transfer_FunctionI;
            public System.Numerics.Complex[][] Reflection_Coefficient;
            public System.Numerics.Complex[][] Trans_Loss;
            public System.Numerics.Complex[][] Trans_Coefficient;
            public double[] NI_Coef;
            public double[][] Ang_Coef_Oct;//[oct][angle]
            public double[][] Ang_tau_Oct;//[oct][angle]
            public double[] RI_Coef = new double[8];
            public double[] TI_Coef = new double[8];
            public double[] RI_Averages;
            private double angle_incr;

            public Smart_Material(bool Trans, List<AbsorptionModels.ABS_Layer> Layers, int Samplefreq, double Air_Density, double SoundSpeed, Finite_Field_Impedance Zr, double step)
            : this(Trans, Layers, Samplefreq, Air_Density, SoundSpeed, Zr, step, CancellationToken.None)
            {
            }

            public Smart_Material(bool Trans, List<AbsorptionModels.ABS_Layer> Layers, int Samplefreq, double Air_Density, double SoundSpeed, Finite_Field_Impedance Zr, double step, CancellationToken CTS = default)
            {
                Buildup = Layers;
                Fs = Samplefreq;
                rho = Air_Density;
                c = SoundSpeed;

                int min_freq = Samplefreq / 4096;
                if (Layers.Count < 1) return;

                //the current version...
                //Z = AbsorptionModels.Operations.Recursive_Transfer_Matrix(false, 10000, 343, Layers, ref frequency, ref Angles);
                //the goal...
                if (Trans)
                {
                    //Z = AbsorptionModels.Operations.Transfer_Matrix_Explicit_Tau(false, 44100, 343, Layers, ref frequency, ref Angles, ref Trans_Loss, ref Reflection_Coefficient);
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(true, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient, CTS);
                    Trans_Coefficient = new System.Numerics.Complex[Trans_Loss.Length][];
                    for (int i = 0; i < Trans_Loss.Length; i++)
                    {
                        Trans_Coefficient[i] = new System.Numerics.Complex[Trans_Loss[i].Length];
                        for (int j = 0; j < Trans_Coefficient[i].Length; j++)
                        {
                            Trans_Coefficient[i][j] = Trans_Loss[i][j] * Trans_Loss[i][j];
                        }
                    }
                }
                else
                {
                    //Z = AbsorptionModels.Operations.Transfer_Matrix_Explicit_Z(false, 44100, 343, Layers, ref frequency, ref Angles);
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(true, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient, CTS);
                    Trans_Coefficient = new System.Numerics.Complex[36][];
                    for (int i = 0; i < Trans_Coefficient.Length; i++)
                    {
                        Trans_Coefficient[i] = new System.Numerics.Complex[frequency.Length];
                        for (int j = 0; j < Trans_Coefficient[i].Length; j++) Trans_Coefficient[i][j] = 0;
                    }
                }
                //////////////////Radiation Impedance///////////////////////
                double[] a_real = new double[Angles.Length]; //prop;
                for (int i = 0; i < Angles.Length; i++) a_real[i] = Angles[i].Real;

                double[][] Angular_Absorption;

                System.Numerics.Complex[][] Zr_interp = Zr.Interpolate(frequency);

                Reflection_Coefficient = Pachyderm_Acoustic.AbsorptionModels.Operations.Reflection_Coef(Z, Zr_interp); //(Z, Air_Density, SoundSpeed); //No defined way to build a complex finite reflection coefficient.
                Angular_Absorption = Pachyderm_Acoustic.AbsorptionModels.Operations.Finite_Absorption_Coefficient(Zr_interp, Z, a_real, rho, 343);

                //if (Zf_incorp_Choice == 0)
                //{
                //Reflection_Coefficient = Pachyderm_Acoustic.AbsorptionModels.Operations.Reflection_Coef(Z, Zr_interp); //(Z, Air_Density, SoundSpeed); //No defined way to build a complex finite reflection coefficient.
                //Angular_Absorption = Pachyderm_Acoustic.AbsorptionModels.Operations.Finite_Absorption_Coefficient(Zr_interp, Z, a_real, rho, 343);
                //}
                //else if (Zf_incorp_Choice == 1)
                //{
                //    Reflection_Coefficient = Pachyderm_Acoustic.AbsorptionModels.Operations.Reflection_Coef(Z, Zr_interp, Air_Density, SoundSpeed); //No defined way to build a complex finite reflection coefficient.
                //    Angular_Absorption = Pachyderm_Acoustic.AbsorptionModels.Operations.Absorption_Coef(Reflection_Coefficient);
                //}
                //else throw new Exception("Field Impedance Incorporation choice not valid or not implemented...");

                Transfer_FunctionR = new MathNet.Numerics.Interpolation.CubicSpline[Angles.Length / 2];
                Transfer_FunctionI = new MathNet.Numerics.Interpolation.CubicSpline[Angles.Length / 2];
                for (int i = 0; i < Reflection_Coefficient.Length / 2; i++)
                {
                    List<double> real = new List<double>(), imag = new List<double>();
                    for (int j = 0; j < Reflection_Coefficient[i].Length; j++)
                    {
                        real.Add(Reflection_Coefficient[i][j].Real);
                        imag.Add(Reflection_Coefficient[i][j].Imaginary);
                    }
                    Transfer_FunctionR[Angles.Length / 2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, real);
                    Transfer_FunctionI[Angles.Length / 2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, imag);
                }

                System.Numerics.Complex[] TI_Averages;

                //if (Averaging_Choice == 0)
                //if (Zf_incorp_Choice == 0)
                //{
                RI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris(Z, Zr_interp);
                //RI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris_Finite(Angular_Absorption);
                //RI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris(Angular_Absorption);
                //RI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris_Finite(Angular_Absorption, Zr_interp, SoundSpeed * Air_Density);
                //RI_Averages = AbsorptionModels.Operations.Finite_Absorption_Coefficient(Zr_interp, Z, Angles, rho, SoundSpeed);
                TI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris(Trans_Coefficient, Zr_interp, SoundSpeed * Air_Density);
                //}
                //else
                //{
                //    RI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris(Angular_Absorption);
                //    TI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris_Finite(Trans_Coefficient);
                //}
                //else if (Averaging_Choice == 1)
                //    if (Zf_incorp_Choice == 0)
                //    {
                //        RI_Averages = AbsorptionModels.Operations.Random_Incidence_0_78(Angular_Absorption, Zr_interp, SoundSpeed * Air_Density);
                //        TI_Averages = AbsorptionModels.Operations.Random_Incidence_0_78(Trans_Coefficient, Zr_interp, SoundSpeed * Air_Density);
                //    }
                //    else
                //    {
                //        RI_Averages = AbsorptionModels.Operations.Random_Incidence_0_78(Angular_Absorption);
                //        TI_Averages = AbsorptionModels.Operations.Random_Incidence_0_78(Trans_Coefficient, Zr_interp, SoundSpeed * Air_Density);
                //    }
                //else if (Averaging_Choice == 2)
                //    if (Zf_incorp_Choice == 0)
                //    {
                //        RI_Averages = AbsorptionModels.Operations.Random_Incidence_NoWeights(Angular_Absorption, Zr_interp, SoundSpeed * Air_Density);
                //        TI_Averages = AbsorptionModels.Operations.Random_Incidence_NoWeights(Trans_Coefficient, Zr_interp, SoundSpeed * Air_Density);
                //    }
                //    else
                //    {
                //        RI_Averages = AbsorptionModels.Operations.Random_Incidence_NoWeights(Angular_Absorption);
                //        TI_Averages = AbsorptionModels.Operations.Random_Incidence_NoWeights(Trans_Coefficient);
                //    }
                //else throw new Exception("Averaging choice not valid or not implemented...");

                NI_Coef = Angular_Absorption[18];
                Ang_Coef_Oct = new double[8][];
                Ang_tau_Oct = new double[8][];

                //5 degree increments, in radians...
                angle_incr = 5 * Math.PI / 180;

                double root2 = Math.Sqrt(2);

                int f = -1;

                for (int oct = 0; oct < 8; oct++)
                {
                    double f_center = 62.5 * Math.Pow(2, oct);
                    double f_lower = (int)((Math.Floor(f_center / root2)));// - min_freq)/df);
                    double f_upper = (int)((Math.Floor(f_center * root2)));// - min_freq)/df);

                    int f_id_l = 0;//(int)Math.Floor((double)((f_lower) / 5));

                    for (int i = 0; i < frequency.Length; i++)
                    {
                        if (frequency[i] < f_lower) f_id_l = i;
                        else break;
                    }

                    int f_id_u;//(int)Math.Floor((double)((f_upper) / 5));

                    for (f_id_u = f_id_l; f_id_u < frequency.Length; f_id_u++)
                    {
                        if (frequency[f_id_u] > f_upper) break;
                    }

                    int count = 0;
                    int RI_count = 0;

                    Ang_Coef_Oct[oct] = new double[Angles.Length];
                    Ang_tau_Oct[oct] = new double[Angles.Length];
                    int[] fct = new int[Angular_Absorption.Length];

                    do
                    {
                        f++;
                        RI_count++;
                        if (f < f_id_l) { f++; continue; }
                        if (f >= frequency.Length) break;
                        RI_Coef[oct] += RI_Averages[f];
                        TI_Coef[oct] += TI_Averages[f].Real;
                        for (int a = 0; a < 19; a++)
                        {
                            if (double.IsNaN(Angular_Absorption[a][f])) continue;
                            fct[a]++;
                            count++;
                            Ang_Coef_Oct[oct][a] += Angular_Absorption[a][f];
                            Ang_tau_Oct[oct][a] += Trans_Coefficient[a][f].Real;
                        }
                        for (int a = 19; a < Angles.Length; a++)
                        {
                            if (double.IsNaN(Angular_Absorption[35 - a][f])) continue;
                            fct[a]++;
                            count++;
                            Ang_Coef_Oct[oct][a] += Angular_Absorption[35 - a][f];
                            Ang_tau_Oct[oct][a] += Trans_Coefficient[35 - a][f].Real;
                        }
                    } while (frequency[f] < f_upper);

                    for (int a = 0; a < Angles.Length; a++) Ang_Coef_Oct[oct][a] /= fct[a];
                    RI_Coef[oct] /= RI_count;
                    TI_Coef[oct] /= RI_count;
                }
            }

            public Smart_Material(bool Trans, List<AbsorptionModels.ABS_Layer> Layers, int Samplefreq, double Air_Density, double SoundSpeed)
            : this(Trans, Layers, Samplefreq, Air_Density, SoundSpeed, CancellationToken.None)
            {
            }

            public Smart_Material(bool Trans, List<AbsorptionModels.ABS_Layer> Layers, int Samplefreq, double Air_Density, double SoundSpeed, CancellationToken CTS = default)
            {
                Buildup = Layers;
                Fs = Samplefreq;
                rho = Air_Density;
                c = SoundSpeed;

                int min_freq = Samplefreq / 4096;
                int max_freq = Samplefreq / 2;

                if (Layers.Count < 1) return;

                //the current version...
                //Z = AbsorptionModels.Operations.Recursive_Transfer_Matrix(false, Samplefreq, 343, Layers, ref frequency, ref Angles);
                //the goal...
                if (Trans)
                {
                    //Z = AbsorptionModels.Operations.Transfer_Matrix_Explicit_Tau(false, 44100, 343, Layers, ref frequency, ref Angles, ref Trans_Loss, ref Reflection_Coefficient);
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(true, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient, CTS);
                    Trans_Coefficient = new System.Numerics.Complex[Trans_Loss.Length][];
                    for (int i = 0; i < Trans_Loss.Length; i++)
                    {
                        Trans_Coefficient[i] = new System.Numerics.Complex[Trans_Loss[i].Length];
                        for (int j = 0; j < Trans_Coefficient[i].Length; j++)
                        {
                            double t = Trans_Loss[i][j].Magnitude;
                            Trans_Coefficient[i][j] = t * t;
                        }
                    }
                }
                else
                {
                    //Z = AbsorptionModels.Operations.Transfer_Matrix_Explicit_Z(false, 44100, 343, Layers, ref frequency, ref Angles);
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(false, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient, CTS);
                    Trans_Coefficient = new System.Numerics.Complex[36][];
                    for (int i = 0; i < Trans_Coefficient.Length; i++) Trans_Coefficient[i] = new System.Numerics.Complex[frequency.Length];
                }

                Reflection_Coefficient = Pachyderm_Acoustic.AbsorptionModels.Operations.Reflection_Coef(Z, Air_Density, SoundSpeed);

                Transfer_FunctionR = new MathNet.Numerics.Interpolation.CubicSpline[Angles.Length / 2];
                Transfer_FunctionI = new MathNet.Numerics.Interpolation.CubicSpline[Angles.Length / 2];
                for (int i = 0; i < Reflection_Coefficient.Length / 2; i++)
                {
                    List<double> real = new List<double>(), imag = new List<double>();
                    for (int j = 0; j < Reflection_Coefficient[i].Length; j++)
                    {
                        real.Add(Reflection_Coefficient[i][j].Real);
                        imag.Add(Reflection_Coefficient[i][j].Imaginary);
                    }
                    Transfer_FunctionR[Angles.Length / 2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, real);
                    Transfer_FunctionI[Angles.Length / 2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, imag);
                }

                double[][] Angular_Absorption = Pachyderm_Acoustic.AbsorptionModels.Operations.Absorption_Coef(Reflection_Coefficient);

                NI_Coef = Angular_Absorption[18];
                double[] RI_Averages;
                System.Numerics.Complex[] TI_Averages = new System.Numerics.Complex[0];

                //if (Averaging_Choice == 0)
                //{
                RI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris(Angular_Absorption);
                TI_Averages = AbsorptionModels.Operations.Random_Incidence_Paris(Trans_Coefficient);
                //}
                //else if (Averaging_Choice == 1)
                //{
                //    RI_Averages = AbsorptionModels.Operations.Random_Incidence_0_78(Angular_Absorption);
                //    TI_Averages = AbsorptionModels.Operations.Random_Incidence_0_78(Trans_Coefficient);
                //}
                //else if (Averaging_Choice == 2)
                //{
                //    RI_Averages = AbsorptionModels.Operations.Random_Incidence_NoWeights(Angular_Absorption);
                //    TI_Averages = AbsorptionModels.Operations.Random_Incidence_NoWeights(Trans_Coefficient);
                //}
                //else throw new Exception("Averaging choice not valid or not implemented...");

                Ang_Coef_Oct = new double[8][];

                ////5 degree increments, in radians...
                //angle_incr = 5 * Math.PI / 180;

                //double root2 = Math.Sqrt(2);

                //for (int oct = 0; oct < 8; oct++)
                //{
                //    double f_center = 62.5 * Math.Pow(2, oct);
                //    int f_lower = (int)Math.Floor(f_center / root2) - min_freq;
                //    int f_upper = (int)Math.Floor(f_center * root2) - min_freq;
                //    int f_id_l = (int)Math.Floor((double)((f_lower) / 5));
                //    int f_id_u = (int)Math.Floor((double)((f_upper) / 5));
                //    int count = 0;
                //    int RI_count = 0;
                //    Ang_Coef_Oct[oct] = new double[Angles.Length];
                //    int[] fct = new int[Angular_Absorption.Length];
                //    int f = 0;

                //    do
                //    {
                //        RI_Coef[oct] += RI_Averages[f];
                //        RI_count++;
                //        for (int a = 0; a < 19; a++)
                //        {
                //            if (double.IsNaN(Angular_Absorption[a][f])) continue;
                //            fct[a]++;
                //            count++;
                //            Ang_Coef_Oct[oct][a] += Angular_Absorption[a][f];                            
                //        }
                //        for (int a = 19; a < Angles.Length; a++)
                //        {
                //            if (double.IsNaN(Angular_Absorption[35 - a][f])) continue;
                //            fct[a]++;
                //            count++;
                //            Ang_Coef_Oct[oct][a] += Angular_Absorption[35 - a][f];
                //        }
                //        f++;
                //    } while (frequency[f] < f_upper);

                //    for (int a = 0; a < Angles.Length; a++) Ang_Coef_Oct[oct][a] /= fct[a];
                //    RI_Coef[oct] /= RI_count;
                //}

                Ang_tau_Oct = new double[8][];

                //5 degree increments, in radians...
                angle_incr = 5 * Math.PI / 180;

                double root2 = Math.Sqrt(2);

                int f = -1;

                for (int oct = 0; oct < 8; oct++)
                {
                    double f_center = 62.5 * Math.Pow(2, oct);
                    double f_lower = (int)((Math.Floor(f_center / root2)));// - min_freq)/df);
                    double f_upper = (int)((Math.Floor(f_center * root2)));// - min_freq)/df);

                    int f_id_l = 0;//(int)Math.Floor((double)((f_lower) / 5));

                    for (int i = 0; i < frequency.Length; i++)
                    {
                        if (frequency[i] < f_lower) f_id_l = i;
                        else break;
                    }

                    int f_id_u;//(int)Math.Floor((double)((f_upper) / 5));

                    for (f_id_u = f_id_l; f_id_u < frequency.Length; f_id_u++)
                    {
                        if (frequency[f_id_u] > f_upper) break;
                    }

                    int count = 0;
                    int RI_count = 0;

                    Ang_Coef_Oct[oct] = new double[Angles.Length];
                    Ang_tau_Oct[oct] = new double[Angles.Length];
                    int[] fct = new int[Angular_Absorption.Length];

                    do
                    {
                        f++;
                        RI_count++;
                        if (f < f_id_l) { f++; continue; }
                        if (f >= frequency.Length) break;
                        RI_Coef[oct] += RI_Averages[f];
                        TI_Coef[oct] += TI_Averages[f].Magnitude;
                        for (int a = 0; a < 19; a++)
                        {
                            if (double.IsNaN(Angular_Absorption[a][f])) continue;
                            fct[a]++;
                            count++;
                            Ang_Coef_Oct[oct][a] += Angular_Absorption[a][f];
                            Ang_tau_Oct[oct][a] += Trans_Coefficient[a][f].Magnitude;
                        }
                        for (int a = 19; a < Angles.Length; a++)
                        {
                            if (double.IsNaN(Angular_Absorption[35 - a][f])) continue;
                            fct[a]++;
                            count++;
                            Ang_Coef_Oct[oct][a] += Angular_Absorption[35 - a][f];
                            Ang_tau_Oct[oct][a] += Trans_Coefficient[35 - a][f].Magnitude;
                        }
                    } while (frequency[f] < f_upper);

                    for (int a = 0; a < Angles.Length; a++) { Ang_Coef_Oct[oct][a] /= fct[a]; Ang_tau_Oct[oct][a] /= fct[a]; }
                    RI_Coef[oct] /= RI_count;
                    TI_Coef[oct] /= RI_count;
                }
            }

            public override double[] Coefficient_A_Broad()
            {
                return RI_Coef;
            }

            public override void Absorb(ref BroadRay Ray, Hare.Geometry.Vector Normal)
            {
                double cos_theta = Hare.Geometry.Hare_math.Dot(Ray.dx, Ray.dy, Ray.dz, Normal.dx, Normal.dy, Normal.dz);
                int index = 18 - (int)Math.Round(Math.Acos(Math.Abs(cos_theta)) / angle_incr);

                for (int oct = 0; oct < 8; oct++) Ray.Energy[oct] *= (1 - Ang_Coef_Oct[oct][index]);
            }

            public override void Absorb(ref BroadRay Ray, out double cos_theta, Hare.Geometry.Vector Normal)
            {
                cos_theta = Hare.Geometry.Hare_math.Dot(Ray.dx, Ray.dy, Ray.dz, Normal.dx, Normal.dy, Normal.dz);
                int index = 18 - (int)Math.Round(Math.Acos(Math.Abs(cos_theta)) / angle_incr);

                for (int oct = 0; oct < 8; oct++) Ray.Energy[oct] *= (1 - Ang_Coef_Oct[oct][index]);
            }

            public override void Absorb(ref OctaveRay Ray, Hare.Geometry.Vector Normal)
            {
                double cos_theta = Hare.Geometry.Hare_math.Dot(Ray.dx, Ray.dy, Ray.dz, Normal.dx, Normal.dy, Normal.dz);
                int index = 18 - (int)Math.Round(Math.Acos(Math.Abs(cos_theta)) / angle_incr);

                Ray.Intensity *= (1 - Ang_Coef_Oct[Ray.Octave][index]);
            }

            public override void Absorb(ref OctaveRay Ray, out double cos_theta, Hare.Geometry.Vector Normal)
            {
                cos_theta = Hare.Geometry.Hare_math.Dot(Ray.dx, Ray.dy, Ray.dz, Normal.dx, Normal.dy, Normal.dz);
                int index = 18 - (int)Math.Round(Math.Acos(Math.Abs(cos_theta)) / angle_incr);

                Ray.Intensity *= (1 - Ang_Coef_Oct[Ray.Octave][index]);
            }

            public override double Coefficient_A_Broad(int Octave)
            {
                return RI_Coef[Octave];
            }

            public override System.Numerics.Complex Reflection_Narrow(double frequency)
            {
                System.Numerics.Complex alpha = 0;
                for (int a = 0; a < Transfer_FunctionR.Length; a++) alpha += new System.Numerics.Complex(Transfer_FunctionR[a].Interpolate(frequency), Transfer_FunctionI[a].Interpolate(frequency));
                alpha /= Transfer_FunctionR.Length;
                return alpha;
            }

            public override System.Numerics.Complex Reflection_Narrow(double frequency, Hare.Geometry.Vector Dir, Hare.Geometry.Vector Normal)
            {
                int a = (int)(Math.Abs(Hare.Geometry.Hare_math.Dot(Dir, Normal)) * 180 / Math.PI / 18);
                return new System.Numerics.Complex(Transfer_FunctionR[a].Interpolate(frequency), Transfer_FunctionI[a].Interpolate(frequency));
            }

            public class Finite_Field_Impedance
            {
                MathNet.Numerics.Interpolation.CubicSpline[] Zr_Curves_R;
                MathNet.Numerics.Interpolation.CubicSpline[] Zr_Curves_I;

                public Finite_Field_Impedance(double Xdim, double Ydim, double freq_limit, double c_sound, double air_density, IReporting Graph)
                {
                    List<double> freq = new List<double>();
                    double f = 15.625;
                    int ct = 1;
                    while (f < freq_limit)
                    {
                        ct++;
                        f = 15.625 * Math.Pow(2, (double)ct / 3f);
                        freq.Add(f);
                    }
                    double[] anglesdeg = new double[(int)(180 / 5)];
                    anglesdeg[0] = -87.5;
                    for (int i = 1; i < anglesdeg.Length; i++) anglesdeg[i] = anglesdeg[i - 1] + 5;

                    System.Numerics.Complex[][] Zr = AbsorptionModels.Operations.Finite_Radiation_Impedance_Atalla_Rect(Xdim, Ydim, freq.ToArray(), anglesdeg, c_sound, air_density, Graph);

                    //System.Numerics.Complex[][] Zr = AbsorptionModels.Operations.Finite_Radiation_Impedance_Rect_Longhand(Xdim, Ydim, freq.ToArray(), anglesdeg, c_sound, Graph);

                    Zr_Curves_R = new MathNet.Numerics.Interpolation.CubicSpline[Zr[0].Length];
                    Zr_Curves_I = new MathNet.Numerics.Interpolation.CubicSpline[Zr[0].Length];
                    for (int a = 0; a < Zr_Curves_R.Length; a++)
                    {
                        double[] ZR = new double[freq.Count];
                        double[] ZI = new double[freq.Count];
                        for (int fr = 0; fr < freq.Count; fr++)
                        {
                            ZR[fr] = Zr[fr][a].Real;
                            ZI[fr] = Zr[fr][a].Imaginary;
                        }
                        Zr_Curves_R[a] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(freq, ZR);
                        Zr_Curves_I[a] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(freq, ZI);
                    }
                }

                public System.Numerics.Complex[][] Interpolate(double[] freq)
                {
                    System.Numerics.Complex[][] Zr = new System.Numerics.Complex[freq.Length][];

                    for (int f = 0; f < freq.Length; f++)
                    {
                        Zr[f] = new System.Numerics.Complex[Zr_Curves_R.Length];
                        for (int a = 0; a < Zr_Curves_R.Length; a++)
                        {
                            Zr[f][a] = new System.Numerics.Complex(Zr_Curves_R[a].Interpolate(freq[f]), Zr_Curves_I[a].Interpolate(freq[f]));
                        }
                    }
                    return Zr;
                }
            }

            public System.Numerics.Complex Admittance(double frequency)
            {
                System.Numerics.Complex R = new System.Numerics.Complex(Transfer_FunctionR[17].Interpolate(frequency), Transfer_FunctionI[17].Interpolate(frequency));
                return (1 - R) / (rho * c * (1 + R));
            }

            public override System.Numerics.Complex[] Reflection_Spectrum(int sample_frequency, int length, Hare.Geometry.Vector Normal, Hare.Geometry.Vector Dir, int threadid)
            {
                int a = (int)(Math.Abs(Hare.Geometry.Hare_math.Dot(Dir, Normal)) * 180 / Math.PI / 18);

                System.Numerics.Complex[] Ref_trns = new System.Numerics.Complex[length];

                for (int j = 0; j < length; j++)
                {
                    double freq = j * (sample_frequency / 2) / length;
                    Ref_trns[j] = new System.Numerics.Complex(Transfer_FunctionR[a].Interpolate(freq), Transfer_FunctionI[a].Interpolate(freq));
                }

                return Ref_trns;
            }

            // -----------------------------------------------------------------------------
            // Botts-style IIR estimation for Smart_Material
            // Replace the previous Botts block with this.
            // Assumes your class already has:
            int rec_order;
            double[] rec_a, rec_b;
            double rec_fs, rec_maxfreq;
            double[] rec_fit_freqs;
            object lock_IIR = new object();
            // and helper methods:
            //   SM_NormalizeBA(...)
            //   SM_EvalDigital(...)
            //   SM_BuildFitFrequencies(...)
            // -----------------------------------------------------------------------------

            private sealed class BottsModelSpec
            {
                public int Order;
                public int RealSections;
                public int ComplexSections;
                public int ParameterCount;
            }



            private sealed class BottsSample
            {
                public double[] Theta;
                public double LogL;
                public double[] Lower;
                public double[] Upper;

                public BottsSample(double[] theta, double logL, double[] lower, double[] upper)
                {
                    Theta = theta;
                    LogL = logL;
                    Lower = lower;
                    Upper = upper;
                }
            }

            //private static double[] SM_BuildFitFrequencies(double fs, double maxFreq)
            //{
            //    List<double> f = new List<double>();

            //    void AddLinear(double f0, double f1, int n)
            //    {
            //        if (n < 2) { f.Add(f0); return; }
            //        for (int i = 0; i < n; i++)
            //        {
            //            double t = (double)i / (n - 1);
            //            f.Add(f0 + t * (f1 - f0));
            //        }
            //    }

            //    void AddLog(double f0, double f1, int n)
            //    {
            //        if (n < 2) { f.Add(f0); return; }
            //        double l0 = Math.Log(f0);
            //        double l1 = Math.Log(f1);
            //        for (int i = 0; i < n; i++)
            //        {
            //            double t = (double)i / (n - 1);
            //            f.Add(Math.Exp(l0 + t * (l1 - l0)));
            //        }
            //    }

            //    double fHi = Math.Min(0.49 * fs, maxFreq);

            //    AddLinear(5.0, 20.0, 32);
            //    AddLinear(20.0, 60.0, 56);
            //    AddLinear(60.0, 150.0, 72);
            //    AddLog(150.0, 500.0, 40);
            //    AddLog(500.0, Math.Min(1200.0, fHi), 24);

            //    if (fHi > 1200.0)
            //        AddLog(1200.0, fHi, 12);

            //    f.Sort();

            //    List<double> outF = new List<double>();
            //    double last = -1.0;
            //    for (int i = 0; i < f.Count; i++)
            //    {
            //        if (i == 0 || Math.Abs(f[i] - last) > 1e-9)
            //        {
            //            outF.Add(f[i]);
            //            last = f[i];
            //        }
            //    }

            //    return outF.ToArray();
            //}
            private static double[] SM_BuildFitFrequencies(double fs, double maxFreq)
            {
                List<double> f = new List<double>();

                void AddLinear(double f0, double f1, int n)
                {
                    if (n < 2) { f.Add(f0); return; }
                    for (int i = 0; i < n; i++)
                    {
                        double t = (double)i / (n - 1);
                        f.Add(f0 + t * (f1 - f0));
                    }
                }

                void AddLog(double f0, double f1, int n)
                {
                    if (n < 2) { f.Add(f0); return; }
                    double l0 = Math.Log(f0);
                    double l1 = Math.Log(f1);
                    for (int i = 0; i < n; i++)
                    {
                        double t = (double)i / (n - 1);
                        f.Add(Math.Exp(l0 + t * (l1 - l0)));
                    }
                }

                double fHi = Math.Min(0.49 * fs, maxFreq);

                // Much denser low-frequency grid for narrow resonances
                AddLinear(5.0, 20.0, 48);
                AddLinear(20.0, 50.0, 72);
                AddLinear(50.0, 120.0, 96);
                AddLinear(120.0, 250.0, 96);

                AddLog(250.0, Math.Min(1000.0, fHi), 64);

                if (fHi > 1000.0)
                    AddLog(1000.0, Math.Min(4000.0, fHi), 32);

                if (fHi > 4000.0)
                    AddLog(4000.0, fHi, 16);

                f.Sort();

                List<double> outF = new List<double>();
                double last = -1.0;
                for (int i = 0; i < f.Count; i++)
                {
                    if (i == 0 || Math.Abs(f[i] - last) > 1e-9)
                    {
                        outF.Add(f[i]);
                        last = f[i];
                    }
                }

                return outF.ToArray();
            }

            private static void SM_NormalizeBA(ref double[] b, ref double[] a)
            {
                if (a == null || a.Length == 0)
                {
                    a = new double[] { 1.0 };
                    return;
                }

                double a0 = a[0];
                if (Math.Abs(a0) < 1e-14) a0 = 1.0;

                if (Math.Abs(a0 - 1.0) > 1e-12)
                {
                    for (int i = 0; i < a.Length; i++) a[i] /= a0;
                    for (int i = 0; i < b.Length; i++) b[i] /= a0;
                }
            }
            private static System.Numerics.Complex SM_EvalDigital(double[] b, double[] a, double f, double fs)
            {
                double w = 2.0 * Math.PI * f / fs;
                System.Numerics.Complex z1 = System.Numerics.Complex.Exp(-System.Numerics.Complex.ImaginaryOne * w);

                System.Numerics.Complex B = System.Numerics.Complex.Zero;
                System.Numerics.Complex A = System.Numerics.Complex.Zero;
                System.Numerics.Complex zp = System.Numerics.Complex.One;

                for (int i = 0; i < b.Length; i++)
                {
                    B += b[i] * zp;
                    zp *= z1;
                }

                zp = System.Numerics.Complex.One;
                for (int i = 0; i < a.Length; i++)
                {
                    A += a[i] * zp;
                    zp *= z1;
                }

                if (A.Magnitude < 1e-18) return System.Numerics.Complex.Zero;
                return B / A;
            }

            private static double SM_WrapPhase(double x)
            {
                while (x > Math.PI) x -= 2.0 * Math.PI;
                while (x < -Math.PI) x += 2.0 * Math.PI;
                return x;
            }

            private static double SM_LogAdd(double a, double b)
            {
                if (double.IsNegativeInfinity(a)) return b;
                if (double.IsNegativeInfinity(b)) return a;
                if (a < b) { double t = a; a = b; b = t; }
                return a + Math.Log(1.0 + Math.Exp(b - a));
            }

            private static double SM_LogSubExp(double a, double b)
            {
                // assumes a > b
                if (b >= a) return double.NegativeInfinity;
                return a + Math.Log(1.0 - Math.Exp(b - a));
            }

            private static double SM_NextGaussian(Random rng)
            {
                double u1 = Math.Max(1e-12, rng.NextDouble());
                double u2 = rng.NextDouble();
                return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            }

            private static double SM_ReflectToBounds(double x, double lo, double hi)
            {
                if (hi <= lo) return lo;
                while (x < lo || x > hi)
                {
                    if (x < lo) x = lo + (lo - x);
                    if (x > hi) x = hi - (x - hi);
                }
                return x;
            }

            private static double[] SM_SampleFromPrior(Random rng, double[] lower, double[] upper)
            {
                double[] x = new double[lower.Length];
                for (int i = 0; i < x.Length; i++)
                    x[i] = lower[i] + rng.NextDouble() * (upper[i] - lower[i]);
                return x;
            }

            private static BottsModelSpec SM_MakeSpec(int order)
            {
                int realSections = order % 2;
                int complexSections = order / 2;

                return new BottsModelSpec
                {
                    Order = order,
                    RealSections = realSections,
                    ComplexSections = complexSections,
                    ParameterCount = 1 + 2 * realSections + 4 * complexSections // 2N+1
                };
            }

            //private static void SM_GetBounds(BottsModelSpec spec, out double[] lower, out double[] upper)
            //{
            //    const double rMax = 0.9997;
            //    const double rMin = 0.7;

            //    lower = new double[spec.ParameterCount];
            //    upper = new double[spec.ParameterCount];

            //    int k = 0;

            //    // log gain
            //    lower[k] = -8.0;
            //    upper[k] = 4.0;
            //    k++;

            //    // real zero / real pole
            //    for (int i = 0; i < spec.RealSections; i++)
            //    {
            //        lower[k] = -rMax; upper[k] = rMax; k++; // q
            //        lower[k] = -rMax; upper[k] = rMax; k++; // p
            //    }

            //    // complex pair: radius and cos(angle) for zero and pole
            //    for (int i = 0; i < spec.ComplexSections; i++)
            //    {
            //        lower[k] = 0.0; upper[k] = rMax; k++; // rq
            //        lower[k] = -1.0; upper[k] = 1.0; k++; // cos phiq
            //        lower[k] = 0.0; upper[k] = rMax; k++; // rp
            //        lower[k] = -1.0; upper[k] = 1.0; k++; // cos phip
            //    }
            //}

            private static void SM_GetBounds(BottsModelSpec spec, double fs, double fMax, out double[] lower, out double[] upper)
            {
                const double rMaxReal = 0.9997;

                lower = new double[spec.ParameterCount];
                upper = new double[spec.ParameterCount];

                int k = 0;

                // log gain
                lower[k] = -8.0;
                upper[k] = 4.0;
                k++;

                // real zero / real pole
                for (int i = 0; i < spec.RealSections; i++)
                {
                    lower[k] = -rMaxReal; upper[k] = rMaxReal; k++; // q
                    lower[k] = -rMaxReal; upper[k] = rMaxReal; k++; // p
                }

                double fHi = Math.Min(0.49 * fs, Math.Max(20.0, fMax));
                double bwHi = Math.Min(2500.0, 0.35 * fs);

                // complex pair:
                // [log f_zero, log bw_zero, log f_pole, log bw_pole]
                for (int i = 0; i < spec.ComplexSections; i++)
                {
                    lower[k] = Math.Log(8.0); upper[k] = Math.Log(fHi); k++; // fq
                    lower[k] = Math.Log(1.0); upper[k] = Math.Log(bwHi); k++; // bwq
                    lower[k] = Math.Log(8.0); upper[k] = Math.Log(fHi); k++; // fp
                    lower[k] = Math.Log(0.5); upper[k] = Math.Log(bwHi); k++; // bwp
                }
            }

            private static double[] SM_Conv(double[] a, double[] b)
            {
                double[] c = new double[a.Length + b.Length - 1];
                for (int i = 0; i < a.Length; i++)
                    for (int j = 0; j < b.Length; j++)
                        c[i + j] += a[i] * b[j];
                return c;
            }

            //private static double[] SM_Conv(double[] a, double[] b)
            //{
            //    double[] c = new double[a.Length + b.Length - 1];
            //    for (int i = 0; i < a.Length; i++)
            //        for (int j = 0; j < b.Length; j++)
            //            c[i + j] += a[i] * b[j];
            //    return c;
            //}
            private static void SM_BuildReflectionCoefficients(BottsModelSpec spec, double[] theta, double fs, out double[] bR, out double[] aR)
            {
                int k = 0;
                double gain = Math.Exp(theta[k++]);

                bR = new double[] { gain };
                aR = new double[] { 1.0 };

                // First-order sections
                for (int i = 0; i < spec.RealSections; i++)
                {
                    double q = theta[k++];
                    double p = theta[k++];

                    double[] b = new double[] { 1.0, -q };
                    double[] a = new double[] { 1.0, -p };

                    bR = SM_Conv(bR, b);
                    aR = SM_Conv(aR, a);
                }

                // Second-order resonant sections parameterized by center frequency and bandwidth
                for (int i = 0; i < spec.ComplexSections; i++)
                {
                    double fq = Math.Exp(theta[k++]);
                    double bwq = Math.Exp(theta[k++]);
                    double fp = Math.Exp(theta[k++]);
                    double bwp = Math.Exp(theta[k++]);

                    fq = Math.Min(Math.Max(fq, 1.0), 0.49 * fs);
                    fp = Math.Min(Math.Max(fp, 1.0), 0.49 * fs);

                    bwq = Math.Min(Math.Max(bwq, 0.1), 0.45 * fs);
                    bwp = Math.Min(Math.Max(bwp, 0.1), 0.45 * fs);

                    double phiq = 2.0 * Math.PI * fq / fs;
                    double phip = 2.0 * Math.PI * fp / fs;

                    double rq = Math.Exp(-Math.PI * bwq / fs);
                    double rp = Math.Exp(-Math.PI * bwp / fs);

                    rq = Math.Min(Math.Max(rq, 0.0), 0.99995);
                    rp = Math.Min(Math.Max(rp, 0.0), 0.99995);

                    double[] b = new double[] { 1.0, -2.0 * rq * Math.Cos(phiq), rq * rq };
                    double[] a = new double[] { 1.0, -2.0 * rp * Math.Cos(phip), rp * rp };

                    bR = SM_Conv(bR, b);
                    aR = SM_Conv(aR, a);
                }

                SM_NormalizeBA(ref bR, ref aR);
            }
            //private static void SM_BuildReflectionCoefficients(BottsModelSpec spec, double[] theta, out double[] bR, out double[] aR)
            //{
            //    int k = 0;
            //    double gain = Math.Exp(theta[k++]);

            //    bR = new double[] { gain };
            //    aR = new double[] { 1.0 };

            //    // First-order sections
            //    for (int i = 0; i < spec.RealSections; i++)
            //    {
            //        double q = theta[k++];
            //        double p = theta[k++];

            //        double[] b = new double[] { 1.0, -q };
            //        double[] a = new double[] { 1.0, -p };

            //        bR = SM_Conv(bR, b);
            //        aR = SM_Conv(aR, a);
            //    }

            //    // Second-order conjugate-pair sections
            //    for (int i = 0; i < spec.ComplexSections; i++)
            //    {
            //        double rq = theta[k++];
            //        double cq = theta[k++];
            //        double rp = theta[k++];
            //        double cp = theta[k++];

            //        double[] b = new double[] { 1.0, -2.0 * rq * cq, rq * rq };
            //        double[] a = new double[] { 1.0, -2.0 * rp * cp, rp * rp };

            //        bR = SM_Conv(bR, b);
            //        aR = SM_Conv(aR, a);
            //    }

            //    SM_NormalizeBA(ref bR, ref aR);
            //}

            private static void SM_ConvertReflectionToAdmittance(double[] bR, double[] aR, out double[] aY, out double[] bY)
            {
                // R = B/A, so Y = (1-R)/(1+R) = (A-B)/(A+B)
                int n = Math.Max(aR.Length, bR.Length);
                aY = new double[n];
                bY = new double[n];

                for (int i = 0; i < n; i++)
                {
                    double a = (i < aR.Length) ? aR[i] : 0.0;
                    double b = (i < bR.Length) ? bR[i] : 0.0;
                    aY[i] = a + b;
                    bY[i] = a - b;
                }

                SM_NormalizeBA(ref bY, ref aY);
            }

            private static double[] SM_BuildBandpassWeights(double[] freqs, double[] alphaTarget)
            {
                int n = freqs.Length;
                double[] w = new double[n];
                if (n == 0) return w;

                double[] d1 = new double[n];
                double[] d2 = new double[n];

                for (int i = 1; i < n - 1; i++)
                {
                    double x0 = Math.Log(Math.Max(freqs[i - 1], 1e-9));
                    double x1 = Math.Log(Math.Max(freqs[i], 1e-9));
                    double x2 = Math.Log(Math.Max(freqs[i + 1], 1e-9));

                    double y0 = alphaTarget[i - 1];
                    double y1 = alphaTarget[i];
                    double y2 = alphaTarget[i + 1];

                    double dx = Math.Max(1e-9, x2 - x0);
                    d1[i] = (y2 - y0) / dx;

                    double dx1 = Math.Max(1e-9, x1 - x0);
                    double dx2 = Math.Max(1e-9, x2 - x1);
                    double s1 = (y1 - y0) / dx1;
                    double s2 = (y2 - y1) / dx2;
                    d2[i] = (s2 - s1) / Math.Max(1e-9, 0.5 * (dx1 + dx2));
                }

                double maxAlpha = Math.Max(1e-9, alphaTarget.Max());
                double maxD1 = Math.Max(1e-9, d1.Select(Math.Abs).Max());
                double maxD2 = Math.Max(1e-9, d2.Select(Math.Abs).Max());

                List<int> peaks = new List<int>();
                for (int i = 1; i < n - 1; i++)
                {
                    if (alphaTarget[i] >= alphaTarget[i - 1] &&
                        alphaTarget[i] >= alphaTarget[i + 1] &&
                        alphaTarget[i] > 0.05)
                    {
                        peaks.Add(i);
                    }
                }

                peaks = peaks
                    .OrderByDescending(i => alphaTarget[i])
                    .Take(4)
                    .ToList();

                for (int i = 0; i < n; i++)
                {
                    double f = freqs[i];

                    double lowBias;
                    if (f < 20.0) lowBias = 2.2;
                    else if (f < 50.0) lowBias = 1.8;
                    else if (f < 100.0) lowBias = 1.4;
                    else if (f < 250.0) lowBias = 1.0;
                    else lowBias = 0.5;

                    double shape =
                        0.8 * (alphaTarget[i] / maxAlpha) +
                        2.8 * (Math.Abs(d1[i]) / maxD1) +
                        2.2 * (Math.Abs(d2[i]) / maxD2);

                    double localPeak = 0.0;
                    for (int p = 0; p < peaks.Count; p++)
                    {
                        double fp = freqs[peaks[p]];
                        double d = Math.Abs(Math.Log(Math.Max(f, 1e-9)) - Math.Log(Math.Max(fp, 1e-9)));
                        localPeak += Math.Exp(-(d * d) / (2.0 * 0.10 * 0.10));
                    }

                    w[i] = lowBias + shape + 6.0 * localPeak;
                }

                return w;
            }
            //private static double[] SM_BuildBandpassWeights(double[] freqs, double[] alphaTarget)
            //{
            //    int n = freqs.Length;
            //    double[] w = new double[n];

            //    int peak = 0;
            //    for (int i = 1; i < n; i++)
            //        if (alphaTarget[i] > alphaTarget[peak]) peak = i;

            //    double f0 = freqs[peak];
            //    double logf0 = Math.Log(Math.Max(f0, 1e-9));

            //    for (int i = 0; i < n; i++)
            //    {
            //        double f = freqs[i];
            //        double local = 0.0;

            //        if (f0 > 0.0)
            //        {
            //            double d = Math.Abs(Math.Log(Math.Max(f, 1e-9)) - logf0);
            //            local = Math.Exp(-(d * d) / (2.0 * 0.45 * 0.45));
            //        }

            //        double lowBias;
            //        if (f < 20.0) lowBias = 1.6;
            //        else if (f < 50.0) lowBias = 1.3;
            //        else if (f < 100.0) lowBias = 1.1;
            //        else lowBias = 0.8;

            //        w[i] = lowBias + 2.5 * local;
            //    }
            //
            //    return w;
            //}

            private static double SM_ComputeLogLikelihoodForReflection(
                BottsModelSpec spec,
                double[] theta,
                System.Numerics.Complex[] Rtarget,
                double[] alphaTarget,
                double[] freqs,
                double[] weights,
                double fs)
            {
                const double eps = 1e-8;
                const double lambdaAlpha = 7.5;
                const double lambdaLogMag = 2.5;
                const double lambdaPhase = 0.25;

                SM_BuildReflectionCoefficients(spec, theta, fs, out double[] bR, out double[] aR);

                double E = 0.0;
                int K = 0;

                for (int i = 0; i < freqs.Length; i++)
                {
                    double wt = weights[i];
                    if (wt <= 0.0) continue;

                    System.Numerics.Complex H = SM_EvalDigital(bR, aR, freqs[i], fs);

                    // enforce passive reflection
                    if (H.Magnitude >= 0.999999 || double.IsNaN(H.Real) || double.IsNaN(H.Imaginary))
                        return double.NegativeInfinity;

                    double magFit = Math.Max(H.Magnitude, eps);
                    double magTar = Math.Max(Rtarget[i].Magnitude, eps);

                    double logMagErr = Math.Log(magFit) - Math.Log(magTar);
                    double phaseErr = SM_WrapPhase(
                        Math.Atan2(Rtarget[i].Imaginary, Rtarget[i].Real) -
                        Math.Atan2(H.Imaginary, H.Real));

                    double alphaFit = 1.0 - H.Magnitude * H.Magnitude;
                    if (alphaFit < 0.0) alphaFit = 0.0;
                    if (alphaFit > 1.0) alphaFit = 1.0;

                    double alphaErr = alphaFit - alphaTarget[i];

                    E += wt * (
                        lambdaAlpha * alphaErr * alphaErr +
                        lambdaLogMag * logMagErr * logMagErr +
                        lambdaPhase * phaseErr * phaseErr);

                    K++;
                }

                if (K == 0) return double.NegativeInfinity;
                if (E <= 1e-30) E = 1e-30;

                // Student-t likelihood, up to additive constants
                return -0.5 * K * Math.Log(E);
            }

            private static BottsSample SM_DrawPriorSample(
                Random rng,
                BottsModelSpec spec,
                System.Numerics.Complex[] Rtarget,
                double[] alphaTarget,
                double[] freqs,
                double[] weights,
                double fs)
            {
                SM_GetBounds(spec, fs, freqs.Last(), out double[] lower, out double[] upper);

                for (int tries = 0; tries < 1000; tries++)
                {
                    double[] theta = SM_SampleFromPrior(rng, lower, upper);
                    double logL = SM_ComputeLogLikelihoodForReflection(spec, theta, Rtarget, alphaTarget, freqs, weights, fs);
                    if (!double.IsNegativeInfinity(logL) && !double.IsNaN(logL))
                        return new BottsSample(theta, logL, lower, upper);
                }

                return new BottsSample(SM_SampleFromPrior(rng, lower, upper), double.NegativeInfinity, lower, upper);
            }

            private static BottsSample SM_ConstrainedRandomWalk(Random rng, BottsModelSpec spec,
                BottsSample seed, double logLmin, System.Numerics.Complex[] Rtarget, double[] alphaTarget, double[] freqs,
                double[] weights,
                double fs,
                int mhSteps)
            {
                double[] x = (double[])seed.Theta.Clone();
                double[] step = new double[x.Length];

                for (int i = 0; i < x.Length; i++)
                    step[i] = 0.08 * (seed.Upper[i] - seed.Lower[i]);

                double bestLogL = seed.LogL;
                double[] bestX = (double[])x.Clone();

                for (int t = 0; t < mhSteps; t++)
                {
                    double[] cand = new double[x.Length];
                    for (int i = 0; i < x.Length; i++)
                    {
                        cand[i] = x[i] + step[i] * SM_NextGaussian(rng);
                        cand[i] = SM_ReflectToBounds(cand[i], seed.Lower[i], seed.Upper[i]);
                    }

                    double logL = SM_ComputeLogLikelihoodForReflection(spec, cand, Rtarget, alphaTarget, freqs, weights, fs);
                    if (logL > logLmin && !double.IsNaN(logL))
                    {
                        x = cand;
                        if (logL > bestLogL)
                        {
                            bestLogL = logL;
                            bestX = (double[])cand.Clone();
                        }
                    }
                }

                if (bestLogL <= logLmin)
                {
                    for (int tries = 0; tries < 1000; tries++)
                    {
                        double[] cand = SM_SampleFromPrior(rng, seed.Lower, seed.Upper);
                        double logL = SM_ComputeLogLikelihoodForReflection(spec, cand, Rtarget, alphaTarget, freqs, weights, fs);
                        if (logL > logLmin && !double.IsNaN(logL))
                            return new BottsSample(cand, logL, seed.Lower, seed.Upper);
                    }
                }

                return new BottsSample(bestX, bestLogL, seed.Lower, seed.Upper);
            }

            private static double[] SM_BestThetaFromSamples(List<BottsSample> dead, List<BottsSample> live)
            {
                BottsSample best = null;

                for (int i = 0; i < dead.Count; i++)
                    if (best == null || dead[i].LogL > best.LogL) best = dead[i];

                for (int i = 0; i < live.Count; i++)
                    if (best == null || live[i].LogL > best.LogL) best = live[i];

                return (best != null) ? (double[])best.Theta.Clone() : null;
            }

            private static (double logZ, double[] thetaBest, int iterations) SM_RunNestedSampling(
                Random rng,
                BottsModelSpec spec,
                System.Numerics.Complex[] Rtarget,
                double[] alphaTarget,
                double[] freqs,
                double[] weights,
                double fs,
                int nLive = 96,
                int mhSteps = 100,
                int maxIter = 12000,
                double stopDeltaLogZ = 1e-7)
            {
                List<BottsSample> live = new List<BottsSample>();
                for (int i = 0; i < nLive; i++)
                    live.Add(SM_DrawPriorSample(rng, spec, Rtarget, alphaTarget, freqs, weights, fs));

                List<BottsSample> dead = new List<BottsSample>();
                List<double> deadLogWt = new List<double>();

                double logZ = double.NegativeInfinity;
                double logXPrev = 0.0;
                int iter;

                for (iter = 0; iter < maxIter; iter++)
                {
                    live.Sort((a, b) => a.LogL.CompareTo(b.LogL));

                    BottsSample worst = live[0];
                    double logXNew = -(double)(iter + 1) / nLive;
                    double logWt = SM_LogSubExp(logXPrev, logXNew) + worst.LogL;

                    dead.Add(worst);
                    deadLogWt.Add(logWt);
                    logZ = SM_LogAdd(logZ, logWt);

                    double maxLiveLogL = live[live.Count - 1].LogL;
                    double logRemaining = logXNew + maxLiveLogL;
                    if ((logRemaining - logZ) < stopDeltaLogZ)
                    {
                        logXPrev = logXNew;
                        break;
                    }

                    int seedIndex = 1 + rng.Next(Math.Max(1, live.Count - 1));
                    BottsSample seed = live[seedIndex];
                    BottsSample repl = SM_ConstrainedRandomWalk(
                        rng, spec, seed, worst.LogL, Rtarget, alphaTarget, freqs, weights, fs, mhSteps);

                    live[0] = repl;
                    logXPrev = logXNew;
                }

                double[] thetaBest = SM_BestThetaFromSamples(dead, live);
                return (logZ, thetaBest, iter + 1);
            }

            public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 6)
            {
                double fs = sample_frequency;

                lock (lock_IIR)
                {
                    int requestedOrder = (filter_order > 0) ? filter_order : 0;

                    if (rec_a != null && rec_b != null &&
                        Math.Abs(rec_fs - fs) < 1e-9 &&
                        Math.Abs(rec_maxfreq - max_freq) < 1e-9 &&
                        rec_order == requestedOrder &&
                        rec_fit_freqs != null)
                    {
                        frequencies = (double[])rec_fit_freqs.Clone();
                        return (rec_a, rec_b);
                    }

                    frequencies = SM_BuildFitFrequencies(fs, max_freq);

                    // Keep this aligned with your current Admittance(double frequency) convention.
                    int idx = 17;
                    if (Transfer_FunctionR == null || Transfer_FunctionR.Length == 0) idx = 0;
                    idx = Math.Max(0, Math.Min(idx, Transfer_FunctionR.Length - 1));

                    System.Numerics.Complex[] Rtarget = new System.Numerics.Complex[frequencies.Length];
                    double[] alphaTarget = new double[frequencies.Length];

                    for (int i = 0; i < frequencies.Length; i++)
                    {
                        double f = frequencies[i];

                        System.Numerics.Complex R = new System.Numerics.Complex(
                            Transfer_FunctionR[idx].Interpolate(f),
                            Transfer_FunctionI[idx].Interpolate(f));

                        double mag = R.Magnitude;
                        if (mag > 0.999) R *= 0.999 / mag;

                        Rtarget[i] = R;

                        double r2 = R.Real * R.Real + R.Imaginary * R.Imaginary;
                        alphaTarget[i] = Math.Max(0.0, Math.Min(1.0, 1.0 - r2));
                    }

                    double[] weights = SM_BuildBandpassWeights(frequencies, alphaTarget);

                    List<int> orders = new List<int>();
                    if (filter_order > 0)
                    {
                        orders.Add(Math.Max(1, filter_order));
                    }
                    else
                    {
                        for (int n = 2; n <= 10; n++) orders.Add(n);
                    }

                    double bestLogZ = double.NegativeInfinity;
                    double[] bestTheta = null;
                    BottsModelSpec bestSpec = null;

                    object bestLock = new object();

                    double[] freqs = frequencies.Clone() as double[];

                    Parallel.ForEach(
                        orders,
                        new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, System.Environment.ProcessorCount - 1) },
                        order =>
                        {
                            Random rngLocal = new Random(Guid.NewGuid().GetHashCode());

                            BottsModelSpec spec = SM_MakeSpec(order);

                            var ns = SM_RunNestedSampling(
                                rngLocal,
                                spec,
                                Rtarget,
                                alphaTarget,
                                freqs,
                                weights,
                                fs,
                                nLive: 96,
                                mhSteps: 100,
                                maxIter: 12000,
                                stopDeltaLogZ: 1e-7);

                            if (ns.thetaBest == null) return;

                            lock (bestLock)
                            {
                                if (ns.logZ > bestLogZ)
                                {
                                    bestLogZ = ns.logZ;
                                    bestTheta = ns.thetaBest;
                                    bestSpec = spec;
                                }
                            }
                        });

                    if (bestTheta == null || bestSpec == null)
                    {
                        rec_fs = fs;
                        rec_maxfreq = max_freq;
                        rec_order = requestedOrder;
                        rec_a = new double[] { 1.0 };
                        rec_b = new double[] { 0.0 };
                        rec_fit_freqs = (double[])frequencies.Clone();
                        return (rec_a, rec_b);
                    }

                    SM_BuildReflectionCoefficients(bestSpec, bestTheta, fs, out double[] bR, out double[] aR);
                    SM_ConvertReflectionToAdmittance(bR, aR, out double[] aY, out double[] bY);

                    rec_fs = fs;
                    rec_maxfreq = max_freq;
                    rec_order = (requestedOrder > 0) ? requestedOrder : bestSpec.Order;
                    rec_a = aY;
                    rec_b = bY;
                    rec_fit_freqs = (double[])frequencies.Clone();

                    return (rec_a, rec_b);
                }
            }
            ////public override (double[] a, double[] b) Estimate_IIR_Coefficients(
            ////    double sample_frequency,
            ////    double max_freq,
            ////    out double[] frequencies,
            ////    int filter_order = 0)
            ////{
            ////    double fs = sample_frequency;

            ////    lock (lock_IIR)
            ////    {
            ////        int requestedOrder = (filter_order > 0) ? filter_order : 0;

            ////        if (rec_a != null && rec_b != null &&
            ////            Math.Abs(rec_fs - fs) < 1e-9 &&
            ////            Math.Abs(rec_maxfreq - max_freq) < 1e-9 &&
            ////            rec_order == requestedOrder &&
            ////            rec_fit_freqs != null)
            ////        {
            ////            frequencies = (double[])rec_fit_freqs.Clone();
            ////            return (rec_a, rec_b);
            ////        }

            ////        // Keep aligned with your current Admittance(double frequency) convention.
            ////        int idx = 17;
            ////        if (Transfer_FunctionR == null || Transfer_FunctionR.Length == 0) idx = 0;
            ////        idx = Math.Max(0, Math.Min(idx, Transfer_FunctionR.Length - 1));

            ////        System.Numerics.Complex[] BuildTarget(double[] freqGrid)
            ////        {
            ////            var Y = new System.Numerics.Complex[freqGrid.Length];

            ////            for (int i = 0; i < freqGrid.Length; i++)
            ////            {
            ////                double f = freqGrid[i];

            ////                System.Numerics.Complex R = new System.Numerics.Complex(
            ////                    Transfer_FunctionR[idx].Interpolate(f),
            ////                    Transfer_FunctionI[idx].Interpolate(f));

            ////                System.Numerics.Complex denom = System.Numerics.Complex.One + R;
            ////                if (denom.Magnitude < 1e-12) denom = new System.Numerics.Complex(1e-12, 0.0);

            ////                System.Numerics.Complex Yn = (System.Numerics.Complex.One - R) / denom;

            ////                // Small passive clamp
            ////                if (Yn.Real < 0.0) Yn = new System.Numerics.Complex(0.0, Yn.Imaginary);

            ////                Y[i] = Yn;
            ////            }

            ////            return Y;
            ////        }

            ////        double[] BuildValueWeights(double[] freqGrid, System.Numerics.Complex[] Ytarget)
            ////        {
            ////            double[] wts = new double[freqGrid.Length];
            ////            double fHiUse = Math.Min(0.49 * fs, max_freq);

            ////            for (int i = 0; i < freqGrid.Length; i++)
            ////            {
            ////                double f = freqGrid[i];
            ////                double w;

            ////                if (f < 5.0 || f > fHiUse)
            ////                {
            ////                    wts[i] = 0.0;
            ////                    continue;
            ////                }

            ////                if (f < 12.0) w = 5.0;
            ////                else if (f < 25.0) w = 4.2;
            ////                else if (f < 50.0) w = 3.4;
            ////                else if (f < 100.0) w = 2.4;
            ////                else if (f < 200.0) w = 1.7;
            ////                else if (f < 1000.0) w = 1.0;
            ////                else if (f < 2500.0) w = 0.35;
            ////                else w = 0.12;

            ////                wts[i] = w / Math.Max(0.02, Ytarget[i].Magnitude);
            ////            }

            ////            return wts;
            ////        }

            ////        double[] BuildSlopeWeights(double[] freqGrid, System.Numerics.Complex[] Ytarget, System.Numerics.Complex[] dYdLog)
            ////        {
            ////            double[] wts = new double[freqGrid.Length];
            ////            double fHiUse = Math.Min(0.49 * fs, max_freq);

            ////            const double lambdaSlope = 0.18;

            ////            for (int i = 0; i < freqGrid.Length; i++)
            ////            {
            ////                double f = freqGrid[i];
            ////                double w;

            ////                if (f < 5.0 || f > fHiUse)
            ////                {
            ////                    wts[i] = 0.0;
            ////                    continue;
            ////                }

            ////                if (f < 12.0) w = 5.0;
            ////                else if (f < 25.0) w = 4.0;
            ////                else if (f < 50.0) w = 3.0;
            ////                else if (f < 100.0) w = 2.0;
            ////                else if (f < 200.0) w = 1.2;
            ////                else if (f < 1000.0) w = 0.55;
            ////                else if (f < 2500.0) w = 0.20;
            ////                else w = 0.08;

            ////                // Keep slope rows from dominating where target is tiny or very flat.
            ////                double scale = Math.Max(0.05, Ytarget[i].Magnitude + 0.5 * dYdLog[i].Magnitude);
            ////                wts[i] = lambdaSlope * w / scale;
            ////            }

            ////            return wts;
            ////        }

            ////        (double[] bestA, double[] bestB, int bestBudget, double bestErr) RunSolve(
            ////            double[] freqGrid,
            ////            System.Numerics.Complex[] Ytarget,
            ////            System.Numerics.Complex[] dYdLog,
            ////            double[] valueWeights,
            ////            double[] slopeWeights)
            ////        {
            ////            // Positive-real dictionary
            ////            double[] hpHzBase = new double[] { 3, 6, 10, 16, 25, 40, 63, 100, 160, 250 };
            ////            double[] bpHzBase = new double[] { 12, 18, 25, 35, 50, 70, 100, 140, 200, 280, 400, 560, 800 };
            ////            double[] qBase = new double[] { 0.7, 1.4, 3.0, 6.0 };

            ////            double[] hpShiftCandidates = new double[] { 0.85, 1.0, 1.2 };
            ////            double[] bpShiftCandidates = new double[] { 0.85, 1.0, 1.15 };

            ////            int[] budgets = (filter_order > 0)
            ////                ? new int[] { Math.Max(4, Math.Min(filter_order, 10)) }
            ////                : new int[] { 5, 7, 9 };

            ////            double fHiUse = Math.Min(0.45 * fs, max_freq);

            ////            double bestErrLocal = double.MaxValue;
            ////            double[] bestALocal = null;
            ////            double[] bestBLocal = null;
            ////            int bestBudgetLocal = 0;

            ////            foreach (int budget in budgets)
            ////            {
            ////                foreach (double hpShift in hpShiftCandidates)
            ////                {
            ////                    foreach (double bpShift in bpShiftCandidates)
            ////                    {
            ////                        List<double> hpHz = new List<double>();
            ////                        for (int i = 0; i < hpHzBase.Length; i++)
            ////                        {
            ////                            double f = hpHzBase[i] * hpShift;
            ////                            if (f >= 2.0 && f <= Math.Min(fHiUse, 400.0))
            ////                                hpHz.Add(f);
            ////                        }

            ////                        List<double> bpHz = new List<double>();
            ////                        List<double> bpQ = new List<double>();

            ////                        for (int i = 0; i < bpHzBase.Length; i++)
            ////                        {
            ////                            double fc = bpHzBase[i] * bpShift;
            ////                            if (fc >= 8.0 && fc <= Math.Min(fHiUse, 1200.0))
            ////                            {
            ////                                for (int q = 0; q < qBase.Length; q++)
            ////                                {
            ////                                    bpHz.Add(fc);
            ////                                    bpQ.Add(qBase[q]);
            ////                                }
            ////                            }
            ////                        }

            ////                        int nHP = hpHz.Count;
            ////                        int nBP = bpHz.Count;
            ////                        int cols = 1 + nHP + nBP; // constant + hp + bp
            ////                        int K = freqGrid.Length;
            ////                        int rows = 4 * K;         // value real/imag + slope real/imag

            ////                        double[,] A = new double[rows, cols];
            ////                        double[] y = new double[rows];
            ////                        double[] colNorm2 = new double[cols];

            ////                        for (int i = 0; i < K; i++)
            ////                        {
            ////                            double f = freqGrid[i];
            ////                            double w = 2.0 * Math.PI * f;
            ////                            System.Numerics.Complex s = new System.Numerics.Complex(0.0, w);

            ////                            int rv = i;
            ////                            int iv = i + K;
            ////                            int rd = i + 2 * K;
            ////                            int id = i + 3 * K;

            ////                            y[rv] = valueWeights[i] * Ytarget[i].Real;
            ////                            y[iv] = valueWeights[i] * Ytarget[i].Imaginary;
            ////                            y[rd] = slopeWeights[i] * dYdLog[i].Real;
            ////                            y[id] = slopeWeights[i] * dYdLog[i].Imaginary;

            ////                            int col = 0;

            ////                            // Constant conductance atom
            ////                            A[rv, col] = valueWeights[i];
            ////                            A[iv, col] = 0.0;
            ////                            A[rd, col] = 0.0;
            ////                            A[id, col] = 0.0;
            ////                            colNorm2[col] += valueWeights[i] * valueWeights[i];
            ////                            col++;

            ////                            // High-pass atoms: H = s/(s+sigma), dH/dlnw = s*sigma/(s+sigma)^2
            ////                            for (int h = 0; h < nHP; h++)
            ////                            {
            ////                                double sigma = SM_WarpHzToAnalogRad(hpHz[h], fs);
            ////                                System.Numerics.Complex H = s / (s + sigma);
            ////                                System.Numerics.Complex dH = s * sigma / ((s + sigma) * (s + sigma));

            ////                                double vr = valueWeights[i] * H.Real;
            ////                                double vi = valueWeights[i] * H.Imaginary;
            ////                                double dr = slopeWeights[i] * dH.Real;
            ////                                double di = slopeWeights[i] * dH.Imaginary;

            ////                                A[rv, col] = vr;
            ////                                A[iv, col] = vi;
            ////                                A[rd, col] = dr;
            ////                                A[id, col] = di;

            ////                                colNorm2[col] += vr * vr + vi * vi + dr * dr + di * di;
            ////                                col++;
            ////                            }

            ////                            // Band-pass atoms:
            ////                            // H = (bw s)/(s^2 + bw s + w0^2)
            ////                            // dH/dlnw = s * bw * (w0^2 - s^2) / D^2
            ////                            for (int bpi = 0; bpi < nBP; bpi++)
            ////                            {
            ////                                double w0 = SM_WarpHzToAnalogRad(bpHz[bpi], fs);
            ////                                double bw = w0 / bpQ[bpi];

            ////                                System.Numerics.Complex D = s * s + bw * s + w0 * w0;
            ////                                if (D.Magnitude < 1e-18) D = new System.Numerics.Complex(1e-18, 0.0);

            ////                                System.Numerics.Complex H = (bw * s) / D;
            ////                                System.Numerics.Complex dH = s * bw * (w0 * w0 - s * s) / (D * D);

            ////                                double vr = valueWeights[i] * H.Real;
            ////                                double vi = valueWeights[i] * H.Imaginary;
            ////                                double dr = slopeWeights[i] * dH.Real;
            ////                                double di = slopeWeights[i] * dH.Imaginary;

            ////                                A[rv, col] = vr;
            ////                                A[iv, col] = vi;
            ////                                A[rd, col] = dr;
            ////                                A[id, col] = di;

            ////                                colNorm2[col] += vr * vr + vi * vi + dr * dr + di * di;
            ////                                col++;
            ////                            }
            ////                        }

            ////                        double[] coeff = SM_SolveNNLSProjected(A, y, iters: 1800, ridge: 1e-6);

            ////                        // Keep only the strongest atoms
            ////                        int dynamicCount = cols - 1;
            ////                        int keepDynamic = Math.Min(Math.Max(3, budget), dynamicCount);

            ////                        List<(double score, int col)> ranked = new List<(double score, int col)>();
            ////                        for (int col = 1; col < cols; col++)
            ////                        {
            ////                            double score = coeff[col] * Math.Sqrt(Math.Max(0.0, colNorm2[col]));
            ////                            ranked.Add((score, col));
            ////                        }

            ////                        ranked.Sort((x, y2) => y2.score.CompareTo(x.score));

            ////                        HashSet<int> keepCols = new HashSet<int>();
            ////                        for (int i = 0; i < Math.Min(keepDynamic, ranked.Count); i++)
            ////                            if (ranked[i].score > 1e-12) keepCols.Add(ranked[i].col);

            ////                        double[] bTot = new double[] { 0.0 };
            ////                        double[] aTot = new double[] { 1.0 };

            ////                        double g0 = Math.Max(0.0, coeff[0]);
            ////                        if (g0 > 0.0)
            ////                        {
            ////                            double[] bc = new double[] { g0 };
            ////                            double[] ac = new double[] { 1.0 };
            ////                            SM_AddParallelSection(ref bTot, ref aTot, bc, ac);
            ////                        }

            ////                        int colHP0 = 1;
            ////                        for (int h = 0; h < nHP; h++)
            ////                        {
            ////                            int col = colHP0 + h;
            ////                            if (!keepCols.Contains(col)) continue;

            ////                            double gain = Math.Max(0.0, coeff[col]);
            ////                            if (gain <= 0.0) continue;

            ////                            double sigma = SM_WarpHzToAnalogRad(hpHz[h], fs);
            ////                            SM_AddDigitalHighPassAtom(ref bTot, ref aTot, gain, sigma, fs);
            ////                        }

            ////                        int colBP0 = 1 + nHP;
            ////                        for (int bpi = 0; bpi < nBP; bpi++)
            ////                        {
            ////                            int col = colBP0 + bpi;
            ////                            if (!keepCols.Contains(col)) continue;

            ////                            double gain = Math.Max(0.0, coeff[col]);
            ////                            if (gain <= 0.0) continue;

            ////                            double w0 = SM_WarpHzToAnalogRad(bpHz[bpi], fs);
            ////                            double Q = bpQ[bpi];
            ////                            SM_AddDigitalBandPassAtom(ref bTot, ref aTot, gain, w0, Q, fs);
            ////                        }

            ////                        double gScale = SM_BestPositiveScale(bTot, aTot, freqGrid, Ytarget, valueWeights, fs);

            ////                        // Apply scalar to numerator only
            ////                        if (Math.Abs(gScale - 1.0) > 1e-12)
            ////                        {
            ////                            for (int ii = 0; ii < bTot.Length; ii++) bTot[ii] *= gScale;
            ////                        }

            ////                        double err = SM_DigitalComplexError(bTot, aTot, freqGrid, Ytarget, valueWeights, fs);

            ////                        if (err < bestErrLocal)
            ////                        {
            ////                            bestErrLocal = err;
            ////                            bestALocal = aTot;
            ////                            bestBLocal = bTot;
            ////                            bestBudgetLocal = budget;
            ////                        }
            ////                    }
            ////                }
            ////            }

            ////            return (bestALocal, bestBLocal, bestBudgetLocal, bestErrLocal);
            ////        }

            ////        // Pass 1
            ////        double[] freqPass1 = SM_BuildFitFrequencies(fs, max_freq);
            ////        System.Numerics.Complex[] targetPass1 = BuildTarget(freqPass1);
            ////        System.Numerics.Complex[] dTargetPass1 = SM_BuildDerivativeLogTarget(freqPass1, targetPass1);
            ////        double[] valueWeightsPass1 = BuildValueWeights(freqPass1, targetPass1);
            ////        double[] slopeWeightsPass1 = BuildSlopeWeights(freqPass1, targetPass1, dTargetPass1);

            ////        var pass1 = RunSolve(freqPass1, targetPass1, dTargetPass1, valueWeightsPass1, slopeWeightsPass1);

            ////        // Pass 2: adaptive refinement around worst-fit frequencies
            ////        double[] freqPass2;
            ////        if (pass1.bestA != null && pass1.bestB != null)
            ////        {
            ////            freqPass2 = SM_RefineFitFrequencies(
            ////                freqPass1,
            ////                pass1.bestB,
            ////                pass1.bestA,
            ////                targetPass1,
            ////                fs,
            ////                extraPerPeak: 6,
            ////                peakCount: 10);
            ////        }
            ////        else
            ////        {
            ////            freqPass2 = freqPass1;
            ////        }

            ////        System.Numerics.Complex[] targetPass2 = BuildTarget(freqPass2);
            ////        System.Numerics.Complex[] dTargetPass2 = SM_BuildDerivativeLogTarget(freqPass2, targetPass2);
            ////        double[] valueWeightsPass2 = BuildValueWeights(freqPass2, targetPass2);
            ////        double[] slopeWeightsPass2 = BuildSlopeWeights(freqPass2, targetPass2, dTargetPass2);

            ////        var pass2 = RunSolve(freqPass2, targetPass2, dTargetPass2, valueWeightsPass2, slopeWeightsPass2);

            ////        double[] bestA = null;
            ////        double[] bestB = null;
            ////        int bestBudget = 0;

            ////        if (pass2.bestA != null && pass2.bestB != null &&
            ////            (pass1.bestA == null || pass2.bestErr <= pass1.bestErr))
            ////        {
            ////            bestA = pass2.bestA;
            ////            bestB = pass2.bestB;
            ////            bestBudget = pass2.bestBudget;
            ////            frequencies = freqPass2;
            ////        }
            ////        else if (pass1.bestA != null && pass1.bestB != null)
            ////        {
            ////            bestA = pass1.bestA;
            ////            bestB = pass1.bestB;
            ////            bestBudget = pass1.bestBudget;
            ////            frequencies = freqPass1;
            ////        }
            ////        else
            ////        {
            ////            frequencies = freqPass2;
            ////        }

            ////        // Guaranteed non-rigid fallback: constant positive conductance
            ////        if (bestA == null || bestB == null)
            ////        {
            ////            double[] freqFallback = frequencies;
            ////            System.Numerics.Complex[] targetFallback = BuildTarget(freqFallback);
            ////            double[] weightsFallback = BuildValueWeights(freqFallback, targetFallback);

            ////            double num = 0.0;
            ////            double den = 0.0;
            ////            for (int i = 0; i < freqFallback.Length; i++)
            ////            {
            ////                if (weightsFallback[i] <= 0.0) continue;
            ////                num += weightsFallback[i] * Math.Max(0.0, targetFallback[i].Real);
            ////                den += weightsFallback[i];
            ////            }

            ////            double g = (den > 0.0) ? (num / den) : 0.05;
            ////            if (g < 1e-6) g = 0.05;

            ////            bestA = new double[] { 1.0 };
            ////            bestB = new double[] { g };
            ////            bestBudget = 1;
            ////        }

            ////        rec_fs = fs;
            ////        rec_maxfreq = max_freq;
            ////        rec_order = (requestedOrder > 0) ? requestedOrder : bestBudget;
            ////        rec_a = bestA;
            ////        rec_b = bestB;
            ////        rec_fit_freqs = (double[])frequencies.Clone();

            ////        return (rec_a, rec_b);
            ////    }
            ////}

            //private static double[] SM_Convolve(double[] a, double[] b)
            //{
            //    double[] c = new double[a.Length + b.Length - 1];
            //    for (int i = 0; i < a.Length; i++)
            //        for (int j = 0; j < b.Length; j++)
            //            c[i + j] += a[i] * b[j];
            //    return c;
            //}

            //private static System.Numerics.Complex SM_EvalZinvPoly(double[] c, System.Numerics.Complex z1)
            //{
            //    // c[0] + c[1] z^-1 + c[2] z^-2 + ...
            //    System.Numerics.Complex acc = 0;
            //    System.Numerics.Complex zp = 1;
            //    for (int k = 0; k < c.Length; k++)
            //    {
            //        acc += c[k] * zp;
            //        zp *= z1;
            //    }
            //    return acc;
            //}

            //private static double[] SM_BuildStableDenominatorBiquads(int order, double fs, double fLo, double fHi, double poleCap, double Q)
            //{
            //    if (order < 2) order = 2;
            //    if ((order & 1) == 1) order++; // even
            //    int nsec = order / 2;

            //    fLo = Math.Max(1.0, fLo);
            //    fHi = Math.Min(0.49 * fs, Math.Max(fLo * 1.01, fHi));

            //    double[] a = new double[] { 1.0 };

            //    for (int s = 0; s < nsec; s++)
            //    {
            //        double t = (s + 0.5) / nsec;
            //        double fc = fLo * Math.Pow(fHi / fLo, t);
            //        double w0 = 2.0 * Math.PI * fc / fs;

            //        double bw = fc / Math.Max(0.25, Q);
            //        double r = Math.Exp(-Math.PI * bw / fs);
            //        r = Math.Min(r, poleCap);
            //        r = Math.Max(r, 0.10);

            //        double c0 = Math.Cos(w0);

            //        // 1 - 2 r cos(w0) z^-1 + r^2 z^-2
            //        double[] sec = new double[] { 1.0, -2.0 * r * c0, r * r };
            //        a = SM_Convolve(a, sec);
            //    }

            //    double a0 = a[0];
            //    if (Math.Abs(a0) > 1e-14 && Math.Abs(a0 - 1.0) > 1e-12)
            //    {
            //        double inv = 1.0 / a0;
            //        for (int i = 0; i < a.Length; i++) a[i] *= inv;
            //    }

            //    return a;
            //}

            //private static double[] SM_SolveNumeratorLSComplexWeighted(
            //    double[] freqs,
            //    System.Numerics.Complex[] Ytarget,
            //    double[] a,
            //    double fs,
            //    double[] weights,
            //    double ridge = 1e-10)
            //{
            //    int M = a.Length;
            //    int K = Math.Min(freqs.Length, Ytarget.Length);

            //    double[,] G = new double[M, M];
            //    double[] h = new double[M];

            //    for (int i = 0; i < K; i++)
            //    {
            //        double wgt = weights[i];
            //        if (wgt <= 0) continue;

            //        double w = 2.0 * Math.PI * freqs[i] / fs;
            //        System.Numerics.Complex z1 = System.Numerics.Complex.Exp(-System.Numerics.Complex.ImaginaryOne * w);

            //        System.Numerics.Complex A = SM_EvalZinvPoly(a, z1);
            //        System.Numerics.Complex d = Ytarget[i] * A;

            //        System.Numerics.Complex bp = System.Numerics.Complex.One;
            //        double[] xr = new double[M];
            //        double[] xi = new double[M];

            //        for (int m = 0; m < M; m++)
            //        {
            //            xr[m] = bp.Real;
            //            xi[m] = bp.Imaginary;
            //            bp *= z1;
            //        }

            //        double yr = d.Real;
            //        double yi = d.Imaginary;
            //        double ww = wgt * wgt;

            //        for (int m = 0; m < M; m++)
            //        {
            //            h[m] += ww * (xr[m] * yr + xi[m] * yi);
            //            for (int n = 0; n < M; n++)
            //                G[m, n] += ww * (xr[m] * xr[n] + xi[m] * xi[n]);
            //        }
            //    }

            //    for (int m = 0; m < M; m++) G[m, m] += ridge;

            //    // Gaussian elimination
            //    double[,] Aaug = new double[M, M + 1];
            //    for (int i = 0; i < M; i++)
            //    {
            //        for (int j = 0; j < M; j++) Aaug[i, j] = G[i, j];
            //        Aaug[i, M] = h[i];
            //    }

            //    for (int col = 0; col < M; col++)
            //    {
            //        int piv = col;
            //        double best = Math.Abs(Aaug[col, col]);
            //        for (int r = col + 1; r < M; r++)
            //        {
            //            double v = Math.Abs(Aaug[r, col]);
            //            if (v > best) { best = v; piv = r; }
            //        }
            //        if (best < 1e-18) break;

            //        if (piv != col)
            //        {
            //            for (int j = col; j < M + 1; j++)
            //            {
            //                double tmp = Aaug[col, j];
            //                Aaug[col, j] = Aaug[piv, j];
            //                Aaug[piv, j] = tmp;
            //            }
            //        }

            //        double diag = Aaug[col, col];
            //        for (int j = col; j < M + 1; j++) Aaug[col, j] /= diag;

            //        for (int r = 0; r < M; r++)
            //        {
            //            if (r == col) continue;
            //            double f = Aaug[r, col];
            //            if (Math.Abs(f) < 1e-18) continue;
            //            for (int j = col; j < M + 1; j++) Aaug[r, j] -= f * Aaug[col, j];
            //        }
            //    }

            //    double[] b = new double[M];
            //    for (int i = 0; i < M; i++) b[i] = Aaug[i, M];
            //    return b;
            //}

            //private static bool SM_DenomIsStable(double[] a, double radius = 0.999999)
            //{
            //    int n = a.Length - 1;
            //    if (n <= 0) return true;

            //    double a0 = a[0];
            //    if (Math.Abs(a0) < 1e-14) return false;

            //    var C = MathNet.Numerics.LinearAlgebra.Double.Matrix.Build.Dense(n, n, 0.0);
            //    for (int j = 0; j < n; j++)
            //        C[0, j] = -a[j + 1] / a0;

            //    for (int i = 1; i < n; i++)
            //        C[i, i - 1] = 1.0;

            //    var evd = C.Evd();
            //    foreach (var p in evd.EigenValues)
            //        if (p.Magnitude >= radius) return false;

            //    return true;
            //}

            //private static double SM_MinRealY(double[] freqs, double[] a, double[] b, double fs, double[] weights)
            //{
            //    double minRe = double.PositiveInfinity;

            //    for (int i = 0; i < freqs.Length; i++)
            //    {
            //        if (weights[i] <= 0) continue;

            //        double w = 2.0 * Math.PI * freqs[i] / fs;
            //        System.Numerics.Complex z1 = System.Numerics.Complex.Exp(-System.Numerics.Complex.ImaginaryOne * w);

            //        System.Numerics.Complex A = SM_EvalZinvPoly(a, z1);
            //        System.Numerics.Complex B = SM_EvalZinvPoly(b, z1);
            //        System.Numerics.Complex Yfit = (A.Magnitude < 1e-18) ? System.Numerics.Complex.Zero : (B / A);

            //        if (Yfit.Real < minRe) minRe = Yfit.Real;
            //    }

            //    return minRe;
            //}

            //private static double SM_ComplexErrorY(
            //    double[] freqs,
            //    System.Numerics.Complex[] Ytarget,
            //    double[] a,
            //    double[] b,
            //    double fs,
            //    double[] weights)
            //{
            //    double errSum = 0.0;
            //    double wSum = 0.0;

            //    for (int i = 0; i < freqs.Length; i++)
            //    {
            //        double wgt = weights[i];
            //        if (wgt <= 0) continue;

            //        double w = 2.0 * Math.PI * freqs[i] / fs;
            //        System.Numerics.Complex z1 = System.Numerics.Complex.Exp(-System.Numerics.Complex.ImaginaryOne * w);

            //        System.Numerics.Complex A = SM_EvalZinvPoly(a, z1);
            //        System.Numerics.Complex B = SM_EvalZinvPoly(b, z1);
            //        System.Numerics.Complex Yfit = (A.Magnitude < 1e-18) ? System.Numerics.Complex.Zero : (B / A);

            //        double denom = Math.Max(0.05, Ytarget[i].Magnitude);
            //        double e = (Yfit - Ytarget[i]).Magnitude / denom;

            //        errSum += wgt * e;
            //        wSum += wgt;
            //    }

            //    return (wSum > 0) ? errSum / wSum : double.MaxValue;
            //}

            //// -----------------------------------------------------------------------------
            //// Replacement estimator for Smart_Material
            //// -----------------------------------------------------------------------------

            //public override (double[] a, double[] b) Estimate_IIR_Coefficients(
            //    double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //{
            //    const int N = 4096;
            //    double fs = sample_frequency;

            //    int K = (int)Math.Floor(max_freq * N / fs) + 1;
            //    K = Math.Max(16, Math.Min(K, N / 2 + 1));

            //    frequencies = new double[K];
            //    for (int i = 0; i < K; i++) frequencies[i] = i * fs / N;

            //    lock (lock_IIR)
            //    {
            //        int requestedOrder = (filter_order > 0) ? filter_order : 0;

            //        if (rec_a != null && rec_b != null &&
            //            Math.Abs(rec_fs - fs) < 1e-9 &&
            //            Math.Abs(rec_maxfreq - max_freq) < 1e-9 &&
            //            rec_order == requestedOrder)
            //        {
            //            return (rec_a, rec_b);
            //        }

            //        // Use the same index convention as your current Admittance(double frequency) method.
            //        int idx = 17;
            //        if (Transfer_FunctionR == null || Transfer_FunctionR.Length == 0) idx = 0;
            //        idx = Math.Max(0, Math.Min(idx, Transfer_FunctionR.Length - 1));

            //        System.Numerics.Complex[] Ytarget = new System.Numerics.Complex[K];

            //        for (int i = 0; i < K; i++)
            //        {
            //            double f = frequencies[i];

            //            System.Numerics.Complex R = new System.Numerics.Complex(
            //                Transfer_FunctionR[idx].Interpolate(f),
            //                Transfer_FunctionI[idx].Interpolate(f));

            //            System.Numerics.Complex denom = System.Numerics.Complex.One + R;
            //            if (denom.Magnitude < 1e-12) denom = new System.Numerics.Complex(1e-12, 0);

            //            // Normalized admittance
            //            System.Numerics.Complex Y = (System.Numerics.Complex.One - R) / denom;

            //            // Small passive fuse
            //            if (Y.Real < 0) Y = new System.Numerics.Complex(0, Y.Imaginary);

            //            Ytarget[i] = Y;
            //        }

            //        double[] weights = new double[K];
            //        for (int i = 0; i < K; i++)
            //        {
            //            double f = frequencies[i];

            //            double bandWeight = 0.0;
            //            if (f >= 20 && f <= Math.Min(0.45 * fs, max_freq))
            //                bandWeight = (f >= 125 && f <= 4000) ? 1.0 : 0.25;

            //            // Relative-error style weighting so small admittances are not ignored
            //            weights[i] = bandWeight / Math.Max(0.05, Ytarget[i].Magnitude);
            //        }

            //        int[] orders = (filter_order > 0)
            //            ? new int[] { ((filter_order & 1) == 1) ? filter_order + 1 : filter_order }
            //            : new int[] { 2, 4, 6 };

            //        double[] poleCaps = new double[] { 0.85, 0.90, 0.93, 0.95 };
            //        double[] Qs = new double[] { 0.6, 0.8, 1.0 };

            //        double bestErr = double.MaxValue;
            //        double[] bestA = null;
            //        double[] bestB = null;
            //        int bestOrder = 0;

            //        foreach (int order0 in orders)
            //        {
            //            int order = Math.Max(2, order0);
            //            if ((order & 1) == 1) order++;

            //            foreach (double poleCap in poleCaps)
            //            {
            //                foreach (double Q in Qs)
            //                {
            //                    double[] aY = SM_BuildStableDenominatorBiquads(
            //                        order,
            //                        fs,
            //                        20.0,
            //                        Math.Min(max_freq * 0.95, 0.45 * fs),
            //                        poleCap,
            //                        Q);

            //                    double[] bY = SM_SolveNumeratorLSComplexWeighted(
            //                        frequencies,
            //                        Ytarget,
            //                        aY,
            //                        fs,
            //                        weights,
            //                        1e-10);

            //                    if (!SM_DenomIsStable(aY)) continue;
            //                    if (SM_MinRealY(frequencies, aY, bY, fs, weights) < -1e-4) continue;

            //                    double err = SM_ComplexErrorY(frequencies, Ytarget, aY, bY, fs, weights);

            //                    if (err < bestErr)
            //                    {
            //                        bestErr = err;
            //                        bestA = aY;
            //                        bestB = bY;
            //                        bestOrder = order;
            //                    }
            //                }
            //            }
            //        }

            //        if (bestA == null || bestB == null)
            //        {
            //            // rigid fallback: normalized admittance = 0
            //            rec_fs = fs;
            //            rec_maxfreq = max_freq;
            //            rec_order = requestedOrder;
            //            rec_a = new double[] { 1.0 };
            //            rec_b = new double[] { 0.0 };
            //            return (rec_a, rec_b);
            //        }

            //        rec_fs = fs;
            //        rec_maxfreq = max_freq;
            //        rec_order = requestedOrder;
            //        rec_a = bestA;
            //        rec_b = bestB;

            //        return (rec_a, rec_b);
            //    }
            //}

            //public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //{
            //    const int N = 4096; // design FFT length / frequency grid basis
            //    double fs = sample_frequency;

            //    // Number of frequency points up to max_freq using bin spacing fs/N
            //    int K = (int)Math.Floor(max_freq * N / fs) + 1;
            //    if (K < 4) K = Math.Min(4, N);

            //    frequencies = new double[K];
            //    for (int i = 0; i < K; i++)
            //        frequencies[i] = i * fs / N; // 0 .. max_freq

            //    // Compute reflection spectrum at normal incidence (as before)
            //    Complex[] Rfull = this.Reflection_Spectrum((int)fs, N,
            //        new Hare.Geometry.Vector(0, 0, 1),
            //        new Hare.Geometry.Vector(0, 0, 1),
            //        0);

            //    // Build admittance target Y = (1 - R) / (1 + R)
            //    Complex[] Y = new Complex[K];
            //    for (int i = 0; i < K; i++)
            //    {
            //        Complex R = Rfull[i];
            //        Complex denom = Complex.One + R;
            //        if (denom.Magnitude < 1e-12) denom = new Complex(1e-12, 0); // avoid R ~ -1 singularity
            //        Y[i] = (Complex.One - R) / denom;

            //        // Optional: crude passivity pre-clamp (helps prevent active fits)
            //        if (Y[i].Real < 0) Y[i] = new Complex(0, Y[i].Imaginary);
            //    }

            //    // If filter_order is zero/negative, pick an order automatically (cached)
            //    if (filter_order <= 0)
            //    {
            //        if (rec_order > 0 && rec_a != null && rec_b != null)
            //            return (rec_a, rec_b);

            //        int min_order = 2;
            //        int max_order = 12;
            //        double best_error = double.MaxValue;
            //        (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);
            //        int best_order = min_order;

            //        for (int order = min_order; order <= max_order; order += 2)
            //        {
            //            try
            //            {
            //                // Fit IIR to admittance spectrum Y(ω)
            //                (double[] b, double[] a) current =
            //                    Pach_SP.IIR_Design.FitIIRToFrequencyResponse(frequencies, Y, order, order, fs);

            //                // Evaluate response (your helper)
            //                Complex[] H =
            //                    Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                        new List<double>(current.b),
            //                        new List<double>(current.a),
            //                        frequencies,
            //                        fs);

            //                // Error metric: compare COMPLEX (magnitude+phase), weighted to 125–4k like you did
            //                double errSum = 0.0;
            //                int count = 0;

            //                for (int i = 0; i < K; i++)
            //                {
            //                    double f = frequencies[i];
            //                    if (f >= 125 && f <= 4000)
            //                    {
            //                        // Relative complex error
            //                        double denomMag = Math.Max(1e-3, Y[i].Magnitude);
            //                        double e = (H[i] - Y[i]).Magnitude / denomMag;
            //                        errSum += e;
            //                        count++;
            //                    }
            //                }

            //                double avgErr = (count > 0) ? errSum / count : double.MaxValue;

            //                // Keep best
            //                if (avgErr < best_error)
            //                {
            //                    best_error = avgErr;
            //                    best_filter = current;
            //                    best_order = order;
            //                }
            //            }
            //            catch
            //            {
            //                continue;
            //            }
            //        }

            //        // Cache and return best
            //        rec_order = best_order;
            //        rec_a = best_filter.a;
            //        rec_b = best_filter.b;
            //        return (rec_a, rec_b);
            //    }
            //    else
            //    {
            //        // Fixed specified order
            //        int order = Math.Max(2, filter_order);
            //        (double[] b, double[] a) =
            //            Pach_SP.IIR_Design.FitIIRToFrequencyResponse(frequencies, Y, order, order, fs);

            //        return (a, b);
            //    }
            //}
            //    public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            //    {
            //        int samplect = 4096; // Number of samples for the IIR filter design
            //        frequencies = new double[(int)(samplect * max_freq / sample_frequency)];

            //        // Get reflection spectrum and calculate desired absorption coefficients
            //        Complex[] desired_Spectrum = new Complex[(int)(samplect * max_freq / sample_frequency)];
            //        Array.Copy(this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0), desired_Spectrum, desired_Spectrum.Length);

            //        for (int i = 0; i < desired_Spectrum.Length; i++)
            //        {
            //            frequencies[i] = (double)((i + 1) * (sample_frequency / 2)) / 4096;
            //        }

            //        double[] desired_alpha = AbsorptionModels.Operations.Absorption_Coef(desired_Spectrum);

            //        // If filter_order is zero or negative, automatically determine the optimal order
            //        if (filter_order <= 0)
            //        {
            //            if (rec_order == 0)
            //            {
            //                // Use a range of filter orders to find the best fit

            //                int min_order = 2;
            //                int max_order = 12;
            //                double error_threshold = 0.01; // 1% average error threshold
            //                double best_error = double.MaxValue;
            //                (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);

            //                for (int order = min_order; order <= max_order; order += 2) // Try even orders
            //                {
            //                    try
            //                    {
            //                        // Design filter with current order
            //                        (double[] b, double[] a) current_filter = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(frequencies, desired_Spectrum, order, order);

            //                        // Normalize the filter
            //                        Complex[] filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                            new List<double>(current_filter.b), new List<double>(current_filter.a),
            //                            frequencies, sample_frequency);

            //                        double mag = filter_response.MaximumMagnitudePhase().Magnitude;
            //                        if (mag > 1)
            //                        {
            //                            for (int j = 0; j < current_filter.b.Length; j++)
            //                                current_filter.b[j] /= mag;
            //                        }
            //                        else
            //                        {
            //                            double mod = mag / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
            //                            for (int j = 0; j < current_filter.b.Length; j++)
            //                                current_filter.b[j] /= mod;
            //                        }

            //                        // Recalculate response after normalization
            //                        filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
            //                            new List<double>(current_filter.b), new List<double>(current_filter.a),
            //                            frequencies, sample_frequency);

            //                        // Calculate error metric focusing on important frequency range
            //                        double error_sum = 0;
            //                        int count = 0;

            //                        for (int i = 0; i < frequencies.Length; i++)
            //                        {
            //                            if (frequencies[i] >= 125 && frequencies[i] <= 4000)
            //                            {
            //                                double mag_error = Math.Abs(filter_response[i].Magnitude - desired_Spectrum[i].Magnitude);
            //                                error_sum += mag_error / Math.Max(0.01, desired_Spectrum[i].Magnitude); // Relative error
            //                                count++;
            //                            }
            //                        }

            //                        double avg_error = count > 0 ? error_sum / count : double.MaxValue;

            //                        // Update best filter if this one is better
            //                        if (avg_error < best_error)
            //                        {
            //                            best_error = avg_error;
            //                            best_filter = current_filter;
            //                        }

            //                        // If error is below threshold, we've found a good enough match
            //                        if (avg_error < error_threshold)
            //                        {
            //                            // We found a filter with acceptable error, no need to try higher orders
            //                            break;
            //                        }
            //                    }
            //                    catch
            //                    {
            //                        // If filter design fails, continue to next order
            //                        continue;
            //                    }
            //                }

            //                // Use the best filter found
            //                return (best_filter.a, best_filter.b);
            //            }
            //            else
            //            {
            //                // Use previously calculated coefficients
            //                return (rec_a, rec_b);
            //            }
            //        }
            //        else
            //        {
            //            // Original code for specified filter order
            //            int numeratorOrder = filter_order;
            //            int denominatorOrder = filter_order;

            //            // Use Modified Yule-Walker method to fit the IIR filter
            //            (double[] b, double[] a) = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(frequencies, desired_Spectrum, numeratorOrder, denominatorOrder);

            //            // Calculate frequency response of the fitted filter
            //            Complex[] IIR_spec = Audio.Pach_SP.IIR_Design.AB_FreqResponse(new List<double>(b), new List<double>(a), frequencies, 44100);

            //            // Normalize the filter if necessary
            //            double m = IIR_spec.MaximumMagnitudePhase().Magnitude;
            //            if (m > 1)
            //            {
            //                for (int j = 0; j < b.Length; j++) b[j] /= m;
            //            }
            //            else
            //            {
            //                double mod = m / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
            //                for (int j = 0; j < b.Length; j++) b[j] /= mod;
            //            }

            //            return (a, b);
            //        }
            //    }
        }

        public class Lambert_Scattering : Scattering
        {
            double[,] Scattering_Coefficient;
            public Lambert_Scattering(double[] Scattering)
            {
                Scattering_Coefficient = new double[8, 3];
                for (int oct = 0; oct < 8; oct++)
                {
                    Scattering_Coefficient[oct, 1] = Scattering[oct];

                    if (Scattering[oct] < 0.25)
                    {
                        Scattering_Coefficient[oct, 0] = 0.1;
                        Scattering_Coefficient[oct, 2] = 0.4;
                    }
                    else
                    {
                        double Mod = Math.Abs(1 - Scattering[oct]) / 5; //((Scattering[oct] < (1 - Scattering[oct])) ? (Scattering[oct] * SplitRatio / 2) : ((1 - Scattering[oct]) * SplitRatio / 2));
                        Scattering_Coefficient[oct, 0] = Scattering_Coefficient[oct, 1] - Mod;
                        Scattering_Coefficient[oct, 2] = Scattering_Coefficient[oct, 1] + Mod;
                    }
                }
            }
            public override double[] Coefficient()
            {
                double[] Scat = new double[8];
                for (int oct = 0; oct < 8; oct++) Scat[oct] = Scattering_Coefficient[oct, 1];
                return Scat;
            }

            public override double Coefficient(int octave)
            {
                return Scattering_Coefficient[octave, 1];
            }

            public override void Scatter_Early(ref BroadRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, double[] Transmission = null)
            {
                if (Cos_Theta > 0)
                {
                    Normal *= -1;
                    Cos_Theta *= -1;
                }

                for (int oct = 0; oct < 8; oct++)
                {
                    if (Ray.Energy[oct] == 0) continue;
                    // 3. Apply Scattering.
                    //// a. Create new source for scattered energy (E * Scattering).
                    //// b. Modify E (E * 1 - Scattering).
                    OctaveRay R = Ray.SplitRay(oct, Scattering_Coefficient[oct, 1]);

                    Hare.Geometry.Vector diffx;
                    Hare.Geometry.Vector diffy;
                    Hare.Geometry.Vector diffz;
                    double proj;
                    //Check that the ray and the normal are both on the same side...
                    diffz = Normal;
                    diffx = new Hare.Geometry.Vector(0, 0, 1);
                    proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz, diffx));

                    if (0.99 < proj && 1.01 > proj) diffx = new Hare.Geometry.Vector(1, 0, 0);
                    diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
                    diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
                    diffx.Normalize();
                    diffy.Normalize();
                    diffz.Normalize();

                    double u1;
                    double u2;
                    double x;
                    double y;
                    double z;
                    Hare.Geometry.Vector vect;
                    u1 = 2.0 * Math.PI * rand.NextDouble();
                    // random azimuth
                    double Scat_Mod = rand.NextDouble();
                    u2 = Math.Acos(Scat_Mod);
                    // random zenith (elevation)
                    x = Math.Cos(u1) * Math.Sin(u2);
                    y = Math.Sin(u1) * Math.Sin(u2);
                    z = Math.Cos(u2);

                    vect = (diffx * x) + (diffy * y) + (diffz * z);
                    vect.Normalize();

                    //Return the new direction
                    R.dx = vect.dx;
                    R.dy = vect.dy;
                    R.dz = vect.dz;
                    R.SetScattered();

                    if (Transmission != null && Transmission[oct] > 0)
                    {
                        OctaveRay tr = Ray.SplitRay(oct, Transmission[oct]);
                        Rays.Enqueue(tr);
                        tr.SetSpecular();
                        OctaveRay td = R.SplitRay(Transmission[oct]);
                        td.Reverse();
                        td.SetScattered();
                        Rays.Enqueue(td);
                    }

                    Rays.Enqueue(R);
                }
                Ray.dx -= Normal.dx * Cos_Theta * 2;
                Ray.dy -= Normal.dy * Cos_Theta * 2;
                Ray.dz -= Normal.dz * Cos_Theta * 2;
            }

            public override void Scatter_VeryLate(ref OctaveRay Ray, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false)
            {
                if (rand.NextDouble() < Scattering_Coefficient[Ray.Octave, 1])
                {
                    Hare.Geometry.Vector diffx;
                    Hare.Geometry.Vector diffy;
                    Hare.Geometry.Vector diffz;
                    double proj;
                    //Check that the ray and the normal are both on the same side...
                    if (Cos_Theta > 0) Normal *= -1;
                    diffz = Normal;
                    diffx = new Hare.Geometry.Vector(0, 0, 1);
                    proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz, diffx));

                    if (0.99 < proj && 1.01 > proj) diffx = new Hare.Geometry.Vector(1, 0, 0);
                    diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
                    diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
                    diffx.Normalize();
                    diffy.Normalize();
                    diffz.Normalize();

                    double u1;
                    double u2;
                    double x;
                    double y;
                    double z;
                    Hare.Geometry.Vector vect;
                    u1 = 2.0 * Math.PI * rand.NextDouble();
                    // random azimuth
                    double Scat_Mod = rand.NextDouble();
                    u2 = Math.Acos(Scat_Mod);
                    // random zenith (elevation)
                    x = Math.Cos(u1) * Math.Sin(u2);
                    y = Math.Sin(u1) * Math.Sin(u2);
                    z = Scat_Mod; //Math.Cos(u2);

                    vect = (diffx * x) + (diffy * y) + (diffz * z);
                    vect.Normalize();

                    //Return the new direction
                    Ray.dx = vect.dx;
                    Ray.dy = vect.dy;
                    Ray.dz = vect.dz;
                    Ray.SetScattered();
                }
                else
                {
                    //Specular Reflection
                    Ray.dx -= Normal.dx * Cos_Theta * 2;
                    Ray.dy -= Normal.dy * Cos_Theta * 2;
                    Ray.dz -= Normal.dz * Cos_Theta * 2;
                    Ray.SetSpecular();
                }

                if (Transmit) Ray.Reverse();
            }

            public override void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false)
            {
                double scat_sel = rand.NextDouble();
                if (scat_sel > Scattering_Coefficient[Ray.Octave, 2])
                {
                    // Specular Reflection
                    Ray.dx -= Normal.dx * Cos_Theta * 2;
                    Ray.dy -= Normal.dy * Cos_Theta * 2;
                    Ray.dz -= Normal.dz * Cos_Theta * 2;
                    if (Transmit) Ray.Reverse();
                    Ray.SetSpecular();
                    return;
                }
                else if (scat_sel > Scattering_Coefficient[Ray.Octave, 0])
                {
                    //Only for a certain portion of high benefit cases--
                    //// a. Create new source for scattered energy (E * Scattering).
                    //// b. Modify E (E * 1 - Scattering).
                    //Create a new ray...
                    OctaveRay tr = Ray.SplitRay(1 - Scattering_Coefficient[Ray.Octave, 1]);
                    tr.SetSpecular();
                    // this is the specular reflection. Save it for later.
                    tr.dx -= Normal.dx * Cos_Theta * 2;
                    tr.dy -= Normal.dy * Cos_Theta * 2;
                    tr.dz -= Normal.dz * Cos_Theta * 2;
                    if (Transmit) tr.Reverse();
                    tr.SetSpecular();
                    Rays.Enqueue(tr);
                }

                //If we are here, the Original ray needs a scattered direction:
                Hare.Geometry.Vector diffx;
                Hare.Geometry.Vector diffy;
                Hare.Geometry.Vector diffz;
                double proj;
                //Check that the ray and the normal are both on the same side...
                if (Cos_Theta > 0) Normal *= -1;
                diffz = Normal;
                diffx = new Hare.Geometry.Vector(0, 0, 1);
                proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz, diffx));

                if (0.99 < proj && 1.01 > proj) diffx = new Hare.Geometry.Vector(1, 0, 0);
                diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
                diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
                diffx.Normalize();
                diffy.Normalize();
                diffz.Normalize();

                double u1;
                double u2;
                double x;
                double y;
                double z;
                Hare.Geometry.Vector vect;
                u1 = 2.0 * Math.PI * rand.NextDouble();
                // random azimuth
                double Scat_Mod = rand.NextDouble();
                u2 = Math.Acos(Scat_Mod);
                // random zenith (elevation)
                x = Math.Cos(u1) * Math.Sin(u2);
                y = Math.Sin(u1) * Math.Sin(u2);
                z = Math.Cos(u2);

                vect = (diffx * x) + (diffy * y) + (diffz * z);
                vect.Normalize();

                //Return the new direction
                Ray.dx = vect.dx;
                Ray.dy = vect.dy;
                Ray.dz = vect.dz;
                Ray.SetScattered();
                if (Transmit) Ray.Reverse();
            }
        }

        //public class Lambert_Scattering : Scattering
        //{
        //    double[,] Scattering_Coefficient;
        //    public Lambert_Scattering(double[] Scattering)
        //    {
        //        Scattering_Coefficient = new double[8, 3];
        //        for (int oct = 0; oct < 8; oct++)
        //        {
        //            double Mod = Math.Abs(Scattering[oct] - 0.5); //((Scattering[oct] < (1 - Scattering[oct])) ? (Scattering[oct] * SplitRatio / 2) : ((1 - Scattering[oct]) * SplitRatio / 2));
        //            Scattering_Coefficient[oct, 1] = Scattering[oct];
        //            Scattering_Coefficient[oct, 0] = Scattering_Coefficient[oct, 1] - Mod;
        //            Scattering_Coefficient[oct, 2] = Scattering_Coefficient[oct, 1] + Mod;
        //        }
        //    }
        //    public override double[] Coefficient()
        //    {
        //        double[] Scat = new double[8];
        //        for (int oct = 0; oct < 8; oct++) Scat[oct] = Scattering_Coefficient[oct, 1];
        //        return Scat;
        //    }

        //    public override double Coefficient(int octave)
        //    {
        //        return Scattering_Coefficient[octave, 1];
        //    }

        //    public override void Scatter_Early(ref BroadRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, double[] Transmission = null)
        //    {
        //        if (Cos_Theta > 0)
        //        {
        //            Normal *= -1;
        //            Cos_Theta *= -1;
        //        }

        //        for(int oct = 0; oct < 8; oct++)
        //        {
        //            if (Ray.Energy[oct] == 0) continue;
        //            // 3. Apply Scattering.
        //            //// a. Create new source for scattered energy (E * Scattering).
        //            //// b. Modify E (E * 1 - Scattering).
        //            OctaveRay R = Ray.SplitRay(oct, Scattering_Coefficient[oct, 1]);

        //            Hare.Geometry.Vector diffx;
        //            Hare.Geometry.Vector diffy;
        //            Hare.Geometry.Vector diffz;
        //            double proj;
        //            //Check that the ray and the normal are both on the same side...
        //            diffz = Normal;
        //            diffx = new Hare.Geometry.Vector(0, 0, 1);
        //            proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz, diffx));

        //            if (0.99 < proj && 1.01 > proj) diffx = new Hare.Geometry.Vector(1, 0, 0);
        //            diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
        //            diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
        //            diffx.Normalize();
        //            diffy.Normalize();
        //            diffz.Normalize();

        //            double u1;
        //            double u2;
        //            double x;
        //            double y;
        //            double z;
        //            Hare.Geometry.Vector vect;
        //            u1 = 2.0 * Math.PI * rand.NextDouble();
        //            // random azimuth
        //            double Scat_Mod = rand.NextDouble();
        //            u2 = Math.Acos(Scat_Mod);
        //            // random zenith (elevation)
        //            x = Math.Cos(u1) * Math.Sin(u2);
        //            y = Math.Sin(u1) * Math.Sin(u2);
        //            z = Math.Cos(u2);

        //            vect = (diffx * x) + (diffy * y) + (diffz * z);
        //            vect.Normalize();

        //            //Return the new direction
        //            R.dx = vect.dx;
        //            R.dy = vect.dy;
        //            R.dz = vect.dz;
        //            R.SetScattered();

        //            if (Transmission != null && Transmission[oct] > 0)
        //            {
        //                OctaveRay tr = Ray.SplitRay(oct, Transmission[oct]);
        //                Rays.Enqueue(tr);
        //                tr.SetSpecular();
        //                OctaveRay td = R.SplitRay(Transmission[oct]);
        //                td.Reverse();
        //                td.SetScattered();
        //                Rays.Enqueue(td);
        //            }

        //            Rays.Enqueue(R);
        //        }
        //        Ray.dx -= Normal.dx * Cos_Theta * 2;
        //        Ray.dy -= Normal.dy * Cos_Theta * 2;
        //        Ray.dz -= Normal.dz * Cos_Theta * 2;
        //    }

        //    public override void Scatter_VeryLate(ref OctaveRay Ray, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false)
        //    {
        //        if (rand.NextDouble() < Scattering_Coefficient[Ray.Octave, 1])
        //        {
        //            Hare.Geometry.Vector diffx;
        //            Hare.Geometry.Vector diffy;
        //            Hare.Geometry.Vector diffz;
        //            double proj;
        //            //Check that the ray and the normal are both on the same side...
        //            if (Cos_Theta > 0) Normal *= -1;
        //            diffz = Normal;
        //            diffx = new Hare.Geometry.Vector(0, 0, 1);
        //            proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz, diffx));

        //            if (0.99 < proj && 1.01 > proj) diffx = new Hare.Geometry.Vector(1, 0, 0);
        //            diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
        //            diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
        //            diffx.Normalize();
        //            diffy.Normalize();
        //            diffz.Normalize();

        //            double u1;
        //            double u2;
        //            double x;
        //            double y;
        //            double z;
        //            Hare.Geometry.Vector vect;
        //            u1 = 2.0 * Math.PI * rand.NextDouble();
        //            // random azimuth
        //            double Scat_Mod = rand.NextDouble();
        //            u2 = Math.Acos(Scat_Mod);
        //            // random zenith (elevation)
        //            x = Math.Cos(u1) * Math.Sin(u2);
        //            y = Math.Sin(u1) * Math.Sin(u2);
        //            z = Scat_Mod; //Math.Cos(u2);

        //            vect = (diffx * x) + (diffy * y) + (diffz * z);
        //            vect.Normalize();

        //            //Return the new direction
        //            Ray.dx = vect.dx;
        //            Ray.dy = vect.dy;
        //            Ray.dz = vect.dz;
        //            Ray.SetScattered();
        //        }
        //        else
        //        {
        //            //Specular Reflection
        //            Ray.dx -= Normal.dx * Cos_Theta * 2;
        //            Ray.dy -= Normal.dy * Cos_Theta * 2;
        //            Ray.dz -= Normal.dz * Cos_Theta * 2;
        //            Ray.SetSpecular();
        //        }

        //        if (Transmit) Ray.Reverse();
        //    }

        //    public override void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false)
        //    {
        //        double scat_sel = rand.NextDouble();
        //        if (scat_sel > Scattering_Coefficient[Ray.Octave, 2])
        //        {
        //            // Specular Reflection
        //            Ray.dx -= Normal.dx * Cos_Theta * 2;
        //            Ray.dy -= Normal.dy * Cos_Theta * 2;
        //            Ray.dz -= Normal.dz * Cos_Theta * 2;
        //            if (Transmit) Ray.Reverse();
        //            Ray.SetSpecular();
        //            return;
        //        }
        //        else if (scat_sel > Scattering_Coefficient[Ray.Octave, 0])
        //        {
        //            //Only for a certain portion of high benefit cases--
        //            //// a. Create new source for scattered energy (E * Scattering).
        //            //// b. Modify E (E * 1 - Scattering).
        //            //Create a new ray...
        //            OctaveRay tr = Ray.SplitRay(1 - Scattering_Coefficient[Ray.Octave,1]);
        //            tr.SetSpecular();
        //            // this is the specular reflection. Save it for later.
        //            tr.dx -= Normal.dx * Cos_Theta * 2;
        //            tr.dy -= Normal.dy * Cos_Theta * 2;
        //            tr.dz -= Normal.dz * Cos_Theta * 2;
        //            if (Transmit) tr.Reverse();
        //            tr.SetSpecular();
        //            Rays.Enqueue(tr);
        //        }

        //        //If we are here, the Original ray needs a scattered direction:
        //        Hare.Geometry.Vector diffx;
        //        Hare.Geometry.Vector diffy;
        //        Hare.Geometry.Vector diffz;
        //        double proj;
        //        //Check that the ray and the normal are both on the same side...
        //        if (Cos_Theta > 0) Normal *= -1;
        //        diffz = Normal;
        //        diffx = new Hare.Geometry.Vector(0, 0, 1);
        //        proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz, diffx));

        //        if (0.99 < proj && 1.01 > proj) diffx = new Hare.Geometry.Vector(1, 0, 0);
        //        diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
        //        diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
        //        diffx.Normalize();
        //        diffy.Normalize();
        //        diffz.Normalize();

        //        double u1;
        //        double u2;
        //        double x;
        //        double y;
        //        double z;
        //        Hare.Geometry.Vector vect;
        //        u1 = 2.0 * Math.PI * rand.NextDouble();
        //        // random azimuth
        //        double Scat_Mod = rand.NextDouble();
        //        u2 = Math.Acos(Scat_Mod);
        //        // random zenith (elevation)
        //        x = Math.Cos(u1) * Math.Sin(u2);
        //        y = Math.Sin(u1) * Math.Sin(u2);
        //        z = Math.Cos(u2);

        //        vect = (diffx * x) + (diffy * y) + (diffz * z);
        //        vect.Normalize();

        //        //Return the new direction
        //        Ray.dx = vect.dx;
        //        Ray.dy = vect.dy;
        //        Ray.dz = vect.dz;
        //        Ray.SetScattered();
        //        if (Transmit) Ray.Reverse();
        //    }
        //}
    }
}