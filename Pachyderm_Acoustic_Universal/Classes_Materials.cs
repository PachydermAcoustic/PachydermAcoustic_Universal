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
            public virtual void ForceIIR(double[] a, double[] b, double fs, double max_freq = 0) { }
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

            public override void ForceIIR(double[] a, double[] b, double fs, double max_freq = 0)
            {
                lock (lock_IIR)
                {
                    rec_a = a != null ? (double[])a.Clone() : null;
                    rec_b = b != null ? (double[])b.Clone() : null;
                    rec_fs = fs;
                    if (max_freq > 0) rec_maxfreq = max_freq;
                    rec_order = 0;
                    _iirForced = true;
                }
            }
            public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            {
                int samplect = 4096;
                int K = (int)Math.Floor(samplect * max_freq / sample_frequency);
                K = Math.Max(8, Math.Min(K, samplect / 2));

                frequencies = new double[K];
                for (int i = 0; i < K; i++)
                    frequencies[i] = i * sample_frequency / samplect; // correct FFT bin mapping

                if (_iirForced && rec_a != null && rec_b != null && Math.Abs(rec_fs - sample_frequency) < 1e-9 && Math.Abs(rec_maxfreq - max_freq) < 1)
                {
                    return (rec_a, rec_b);
                }

                lock (lock_IIR)
                {
                    // rec_order is used here as the requested pole budget cache key
                    //int requestedPoleBudget = (filter_order > 0) ? Math.Max(1, Math.Min(filter_order, 6) : 0;
                    int requestedPoleBudget = filter_order;

                    if (rec_a != null && rec_b != null && Math.Abs(rec_fs - sample_frequency) < 1e-9 && Math.Abs(rec_maxfreq - max_freq) < 1e-9 && rec_order == requestedPoleBudget)
                    {
                        return (rec_a, rec_b);
                    }

                    rec_fs = sample_frequency;
                    rec_maxfreq = max_freq;
                    rec_order = requestedPoleBudget;

                    // Reflection spectrum
                    Complex[] Rfull = this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0);

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

                    // Small pole budget, better search.
                    // The NNLS solver operates on the (tiny) normal-equation system
                    // instead of the full frequency-sampled design matrix, which is
                    // roughly rows/cols times cheaper per iteration (rows ~ thousands of
                    // frequency bins, cols ~ a handful of pole/gain parameters). Search
                    // parameters are kept at their previously validated, conservative
                    // levels (higher pole budgets/candidate density were found to make
                    // the greedy pole search unstable/degenerate for some materials);
                    // the performance headroom from the cheaper solver is realized as
                    // faster evaluation at the same fit quality/reliability, rather than
                    // being spent on a riskier search.
                    int poleBudget = (filter_order > 0) ? Math.Max(1, Math.Min(filter_order, 6)) : 4;

                    (rec_a, rec_b) = BasicMaterialPassiveFit.FitProfileAware(frequencies, Y, sample_frequency, poleBudget, fLo: 20.0, fHi: Math.Min(max_freq * 0.95, 0.45 * sample_frequency), poleCap: 0.95
                    );

                    // Failsafe: never surface a degenerate (NaN/Infinity) fit. If the
                    // requested pole budget produced an unusable result, fall back to a
                    // minimal single-pole fit, which is numerically the most stable
                    // configuration and always well-conditioned.
                    if (!IsFinite(rec_a) || !IsFinite(rec_b))
                    {
                        (rec_a, rec_b) = BasicMaterialPassiveFit.FitLimitedPoleBudget(frequencies, Y, sample_frequency, poleCount: 1, fLo: 20.0, fHi: Math.Min(max_freq * 0.95, 0.45 * sample_frequency), poleCap: 0.95, candidateCount: 12, refinePasses: 2, ridge: 1e-6, iters: 1200);
                    }

                    return (rec_a, rec_b);
                }
            }

            private static bool IsFinite(double[] values)
            {
                if (values == null || values.Length == 0) return false;
                for (int i = 0; i < values.Length; i++)
                    if (double.IsNaN(values[i]) || double.IsInfinity(values[i])) return false;
                return true;
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

                // Builds the (small) normal-equation matrices AtA (cols x cols) and Atb (cols)
                // from the (large) design matrix A (rows x cols) and target y (rows).
                // Solving the NNLS problem via these normal equations is mathematically
                // equivalent to operating on A/y directly (same gradient: 2*Aᵀ(Aw-y) == 2*(AtA*w - Atb)),
                // but reduces the cost of each solver iteration from O(rows*cols) to O(cols^2),
                // which matters enormously here since rows (~thousands of frequency bins) is
                // far larger than cols (a handful of pole/gain parameters).
                private static void BuildNormalEquations(double[,] A, double[] y, out double[,] AtA, out double[] Atb)
                {
                    int rows = A.GetLength(0);
                    int cols = A.GetLength(1);

                    AtA = new double[cols, cols];
                    Atb = new double[cols];

                    for (int r = 0; r < rows; r++)
                    {
                        for (int c1 = 0; c1 < cols; c1++)
                        {
                            double a1 = A[r, c1];
                            if (a1 == 0.0) continue;
                            Atb[c1] += a1 * y[r];
                            for (int c2 = c1; c2 < cols; c2++)
                            {
                                AtA[c1, c2] += a1 * A[r, c2];
                            }
                        }
                    }

                    // mirror symmetric entries
                    for (int c1 = 0; c1 < cols; c1++)
                        for (int c2 = c1 + 1; c2 < cols; c2++)
                            AtA[c2, c1] = AtA[c1, c2];
                }

                private static double[] SolveNNLSFromNormalEquations(double[,] AtA, double[] Atb, int iters, double ridge)
                {
                    int cols = AtA.GetLength(0);

                    double[] w = new double[cols];
                    double[] grad = new double[cols];

                    // rough Lipschitz estimate via power iteration on the (already tiny) AtA matrix
                    double[] v = new double[cols];
                    for (int i = 0; i < cols; i++) v[i] = 1.0 / Math.Max(1, cols);

                    for (int t = 0; t < 20; t++)
                    {
                        double[] v2 = new double[cols];
                        for (int c1 = 0; c1 < cols; c1++)
                        {
                            double s = ridge * v[c1];
                            for (int c2 = 0; c2 < cols; c2++) s += AtA[c1, c2] * v[c2];
                            v2[c1] = s;
                        }

                        double norm = 0;
                        for (int c = 0; c < cols; c++) norm += v2[c] * v2[c];
                        norm = Math.Sqrt(Math.Max(1e-30, norm));
                        for (int c = 0; c < cols; c++) v[c] = v2[c] / norm;
                    }

                    double[] AtAv = new double[cols];
                    for (int c1 = 0; c1 < cols; c1++)
                    {
                        double s = ridge * v[c1];
                        for (int c2 = 0; c2 < cols; c2++) s += AtA[c1, c2] * v[c2];
                        AtAv[c1] = s;
                    }

                    double num = 0, den = 0;
                    for (int c = 0; c < cols; c++)
                    {
                        num += v[c] * AtAv[c];
                        den += v[c] * v[c];
                    }

                    double L = Math.Max(1e-6, num / Math.Max(1e-30, den));
                    double step = 1.0 / L;

                    const double tol = 1e-12;

                    for (int it = 0; it < iters; it++)
                    {
                        double gradNorm2 = 0;
                        for (int c1 = 0; c1 < cols; c1++)
                        {
                            double s = -Atb[c1] + ridge * w[c1];
                            for (int c2 = 0; c2 < cols; c2++) s += AtA[c1, c2] * w[c2];
                            grad[c1] = 2.0 * s;
                            gradNorm2 += grad[c1] * grad[c1];
                        }

                        for (int c = 0; c < cols; c++)
                        {
                            double wc = w[c] - step * grad[c];
                            w[c] = (wc < 0) ? 0 : wc;
                        }

                        // Early exit once the gradient is negligible; a handful of unknowns
                        // typically converges in well under 1200 iterations, so this avoids
                        // burning the full iteration budget on every candidate evaluation.
                        if (gradNorm2 < tol) break;
                    }

                    return w;
                }

                private static double[] SolveNNLSProjected(double[,] A, double[] y, int iters, double ridge)
                {
                    BuildNormalEquations(A, y, out double[,] AtA, out double[] Atb);
                    return SolveNNLSFromNormalEquations(AtA, Atb, iters, ridge);
                }

                private enum Profile
                {
                    Flat,
                    LowPass,
                    PorousRising,
                    Resonant
                }

                private struct ResonantSection
                {
                    public double[] Denominator;
                    public double[] Numerator;
                }

                public static (double[] a, double[] b) FitProfileAware(
                    double[] freqs,
                    Complex[] Ytarget,
                    double fs,
                    int poleBudget,
                    double fLo,
                    double fHi,
                    double poleCap)
                {
                    Profile profile = ClassifyProfile(freqs, Ytarget, fLo, fHi, out List<double> peakFrequencies, out double meanAdmittance);

                    if (profile == Profile.Flat)
                    {
                        // Keep a first-order, exactly cancelling numerator/denominator pair
                        // instead of a zero-order gain so the FDTD boundary always has state.
                        const double pole = 0.25;
                        return (new double[] { 1.0, -pole }, new double[] { meanAdmittance, -pole * meanAdmittance });
                    }

                    if (profile == Profile.Resonant)
                    {
                        int maxSections = Math.Max(1, Math.Min(3, poleBudget > 0 ? poleBudget / 2 : 2));
                        if (peakFrequencies.Count > maxSections)
                            peakFrequencies = peakFrequencies.Take(maxSections).ToList();

                        (double[] a, double[] b) = FitResonantSections(freqs, Ytarget, fs, peakFrequencies, fLo, fHi, poleCap);
                        if (IsFinite(a) && IsFinite(b)) return (a, b);
                    }

                    if (profile == Profile.PorousRising)
                    {
                        (double[] a, double[] b) = FitPorousRising(freqs, Ytarget, fs, poleBudget, fLo, fHi, poleCap);
                        if (IsFinite(a) && IsFinite(b)) return (a, b);
                    }

                    return FitLimitedPoleBudget(
                        freqs,
                        Ytarget,
                        fs,
                        poleBudget,
                        fLo,
                        fHi,
                        poleCap,
                        candidateCount: 28,
                        refinePasses: 4,
                        ridge: 1e-8,
                        iters: 1200);
                }

                private static Profile ClassifyProfile(
                    double[] freqs,
                    Complex[] Ytarget,
                    double fLo,
                    double fHi,
                    out List<double> peakFrequencies,
                    out double meanAdmittance)
                {
                    peakFrequencies = new List<double>();
                    List<int> valid = new List<int>();
                    double min = double.PositiveInfinity;
                    double max = double.NegativeInfinity;
                    double sum = 0.0;
                    double lowMin = double.PositiveInfinity;
                    double lowMax = double.NegativeInfinity;
                    double lowSum = 0.0;
                    int lowCount = 0;
                    double highSum = 0.0;
                    int highCount = 0;

                    int count = Math.Min(freqs.Length, Ytarget.Length);
                    for (int i = 0; i < count; i++)
                    {
                        if (freqs[i] < fLo || freqs[i] > fHi) continue;
                        double value = Math.Max(0.0, Ytarget[i].Real);
                        if (double.IsNaN(value) || double.IsInfinity(value)) continue;
                        valid.Add(i);
                        min = Math.Min(min, value);
                        max = Math.Max(max, value);
                        sum += value;

                        if (freqs[i] <= Math.Min(250.0, fHi))
                        {
                            lowMin = Math.Min(lowMin, value);
                            lowMax = Math.Max(lowMax, value);
                            lowSum += value;
                            lowCount++;
                        }

                        if (freqs[i] >= 1000.0 && freqs[i] <= Math.Min(4000.0, fHi))
                        {
                            highSum += value;
                            highCount++;
                        }
                    }

                    if (valid.Count == 0)
                    {
                        meanAdmittance = 0.0;
                        return Profile.Flat;
                    }

                    meanAdmittance = sum / valid.Count;
                    double range = max - min;

                    if (lowCount >= 3 && highCount >= 3)
                    {
                        double lowMean = lowSum / lowCount;
                        double highMean = highSum / highCount;
                        double rise = highMean - lowMean;

                        // For porous materials: low-band variation must be small relative to the rise
                        // Allow up to 50% of the rise as "flat enough" variation in the low band
                        double lowVariation = lowMax - lowMin;
                        bool lowIsRelativelyFlat = lowVariation <= Math.Max(0.02, 0.5 * rise);

                        // Rise must be significant (at least 5% of low-band level or 0.01 absolute)
                        bool risesAtHighFrequency = rise >= Math.Max(0.01, 0.05 * Math.Max(lowMean, 0.05));

                        if (lowIsRelativelyFlat && risesAtHighFrequency) return Profile.PorousRising;
                    }

                    if (range <= Math.Max(0.01, 0.12 * Math.Max(meanAdmittance, 0.05)))
                        return Profile.Flat;

                    double lowPeakLimit = Math.Min(fHi, 1200.0);
                    double minProminence = Math.Max(0.01, 0.15 * range);
                    for (int i = 1; i < count - 1; i++)
                    {
                        if (freqs[i] < fLo || freqs[i] > lowPeakLimit) continue;

                        double current = Math.Max(0.0, Ytarget[i].Real);
                        double left = Math.Max(0.0, Ytarget[i - 1].Real);
                        double right = Math.Max(0.0, Ytarget[i + 1].Real);
                        if (current < left || current < right || current - Math.Max(left, right) < minProminence) continue;

                        bool distinct = true;
                        for (int p = 0; p < peakFrequencies.Count; p++)
                        {
                            double ratio = freqs[i] / peakFrequencies[p];
                            if (ratio > 0.72 && ratio < 1.38)
                            {
                                distinct = false;
                                break;
                            }
                        }

                        if (distinct) peakFrequencies.Add(freqs[i]);
                    }

                    return peakFrequencies.Count > 0 ? Profile.Resonant : Profile.LowPass;
                }

                private static (double[] a, double[] b) FitResonantSections(
                    double[] freqs,
                    Complex[] Ytarget,
                    double fs,
                    List<double> peakFrequencies,
                    double fLo,
                    double fHi,
                    double poleCap)
                {
                    int sectionCount = peakFrequencies.Count;
                    if (sectionCount == 0) return (null, null);

                    ResonantSection[] sections = new ResonantSection[sectionCount];
                    for (int i = 0; i < sectionCount; i++)
                    {
                        double center = Math.Min(fHi, Math.Max(fLo, peakFrequencies[i]));
                        double bandwidth = Math.Max(80.0, 0.35 * center);
                        double radius = Math.Exp(-Math.PI * bandwidth / fs);
                        radius = Math.Min(Math.Min(poleCap, 0.95), Math.Max(0.70, radius));
                        double omega = 2.0 * Math.PI * center / fs;

                        sections[i] = new ResonantSection
                        {
                            Denominator = new double[] { 1.0, -2.0 * radius * Math.Cos(omega), radius * radius },
                            Numerator = new double[] { 1.0 - radius * radius }
                        };
                    }

                    int kCount = Math.Min(freqs.Length, Ytarget.Length);
                    int rows = 2 * kCount;
                    int cols = sectionCount + 1;
                    double[,] matrix = new double[rows, cols];
                    double[] target = new double[rows];

                    for (int i = 0; i < kCount; i++)
                    {
                        double weight = freqs[i] >= 125.0 && freqs[i] <= 4000.0 ? 1.0 : 0.25;
                        if (freqs[i] < fLo || freqs[i] > fHi) weight = 0.0;
                        target[i] = weight * Ytarget[i].Real;
                        target[i + kCount] = weight * Ytarget[i].Imaginary;
                        matrix[i, 0] = weight;

                        double omega = 2.0 * Math.PI * freqs[i] / fs;
                        Complex z1 = Complex.Exp(-Complex.ImaginaryOne * omega);
                        for (int section = 0; section < sectionCount; section++)
                        {
                            Complex response = EvalTransfer(sections[section].Numerator, sections[section].Denominator, z1);
                            matrix[i, section + 1] = weight * response.Real;
                            matrix[i + kCount, section + 1] = weight * response.Imaginary;
                        }
                    }

                    double[] gains = SolveNNLSProjected(matrix, target, iters: 1200, ridge: 1e-7);
                    double[][] prefix = new double[sectionCount + 1][];
                    double[][] suffix = new double[sectionCount + 1][];
                    prefix[0] = new double[] { 1.0 };
                    for (int i = 0; i < sectionCount; i++) prefix[i + 1] = Conv(prefix[i], sections[i].Denominator);
                    suffix[sectionCount] = new double[] { 1.0 };
                    for (int i = sectionCount - 1; i >= 0; i--) suffix[i] = Conv(sections[i].Denominator, suffix[i + 1]);

                    double[] a = prefix[sectionCount];
                    double[] b = new double[a.Length];
                    for (int i = 0; i < a.Length; i++) b[i] = gains[0] * a[i];

                    for (int section = 0; section < sectionCount; section++)
                    {
                        double[] withoutSection = Conv(prefix[section], suffix[section + 1]);
                        double[] contribution = Conv(sections[section].Numerator, withoutSection);
                        for (int i = 0; i < contribution.Length; i++) b[i] += gains[section + 1] * contribution[i];
                    }

                    return (a, b);
                }

                private static (double[] a, double[] b) FitPorousRising(
                    double[] freqs,
                    Complex[] Ytarget,
                    double fs,
                    int poleBudget,
                    double fLo,
                    double fHi,
                    double poleCap)
                {
                    int kCount = Math.Min(freqs.Length, Ytarget.Length);
                    double lowSum = 0.0;
                    int lowCount = 0;
                    double highSum = 0.0;
                    int highCount = 0;
                    double lowUpper = Math.Min(250.0, fHi);
                    double highLower = Math.Max(1000.0, fLo);
                    for (int i = 0; i < kCount; i++)
                    {
                        double value = Math.Max(0.0, Ytarget[i].Real);
                        if (freqs[i] >= fLo && freqs[i] <= lowUpper)
                        {
                            lowSum += value;
                            lowCount++;
                        }
                        if (freqs[i] >= highLower && freqs[i] <= Math.Min(4000.0, fHi))
                        {
                            highSum += value;
                            highCount++;
                        }
                    }

                    if (lowCount == 0 || highCount == 0) return (null, null);

                    double lowLevel = lowSum / lowCount;
                    double highLevel = highSum / highCount;
                    double rise = highLevel - lowLevel;
                    // For admittance, use a small absolute threshold since typical porous
                    // admittance values are O(0.1-1.0), much smaller than absorption (0-1)
                    if (rise < Math.Max(0.01, 0.05 * Math.Max(lowLevel, 0.05))) return (null, null);

                    double halfLevel = lowLevel + 0.5 * rise;
                    double corner = Math.Sqrt(lowUpper * highLower);
                    for (int i = 0; i < kCount; i++)
                    {
                        // Search for the transition frequency in the band between low and high
                        if (freqs[i] <= lowUpper || freqs[i] >= highLower) continue;
                        if (Ytarget[i].Real >= halfLevel)
                        {
                            corner = freqs[i];
                            break;
                        }
                    }

                    double pole = Math.Exp(-2.0 * Math.PI * corner / fs);
                    pole = Math.Min(Math.Min(poleCap, 0.95), Math.Max(0.20, pole));
                    double highReference = Math.Min(Math.Max(highLower, 1500.0), Math.Min(4000.0, fHi));
                    Complex shelfAtHigh = EvalTransfer(
                        new double[] { 1.0 - pole, pole - 1.0 },
                        new double[] { 1.0, -pole },
                        Complex.Exp(-Complex.ImaginaryOne * 2.0 * Math.PI * highReference / fs));
                    double shelfReal = Math.Max(1e-6, shelfAtHigh.Real);
                    double shelfGain = rise / shelfReal;

                    double[] a = new double[] { 1.0, -pole };
                    double[] b = new double[]
                    {
                        lowLevel - shelfGain * (pole - 1.0),
                        -lowLevel * pole + shelfGain * (pole - 1.0)
                    };

                    return (a, b);
                }

                private static Complex EvalTransfer(double[] numerator, double[] denominator, Complex z1)
                {
                    Complex num = Complex.Zero;
                    Complex den = Complex.Zero;
                    Complex zPower = Complex.One;
                    int length = Math.Max(numerator.Length, denominator.Length);
                    for (int i = 0; i < length; i++)
                    {
                        if (i < numerator.Length) num += numerator[i] * zPower;
                        if (i < denominator.Length) den += denominator[i] * zPower;
                        zPower *= z1;
                    }

                    if (den.Magnitude < 1e-12) den = new Complex(1e-12, 0.0);
                    return num / den;
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
        }

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

            // Botts-style IIR estimation for Smart_Material0
            int rec_order;
            double[] rec_a, rec_b;
            double rec_fs, rec_maxfreq;
            object lock_IIR = new object();
            bool _iirForced = false;

            public override void ForceIIR(double[] a, double[] b, double fs, double max_freq = 0)
            {
                lock (lock_IIR)
                {
                    rec_a = a != null ? (double[])a.Clone() : null;
                    rec_b = b != null ? (double[])b.Clone() : null;
                    rec_fs = fs;
                    if (max_freq > 0) rec_maxfreq = max_freq;
                    rec_order = 0;
                    _iirForced = true;
                }
            }

            private sealed class BottsModelSpec
            {
                public int Order;
                public int RealSections;
                public int ComplexSections;
                public int ParameterCount;
            }

            private sealed class BottsTargetResonance
            {
                public double CenterLog;
                public double WidthLog;
                public double WindowLoLog;
                public double WindowHiLog;
                public double Prominence;
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
                double[] prominence = new double[n];

                for (int i = 1; i < n - 1; i++)
                {
                    if (alphaTarget[i] < alphaTarget[i - 1] || alphaTarget[i] < alphaTarget[i + 1]) continue;

                    // Determine how far this peak rises above its nearby response.
                    // Use approximately a third-octave neighborhood on each side.
                    double xi = Math.Log(Math.Max(freqs[i], 1e-9));
                    double leftMin = alphaTarget[i];
                    double rightMin = alphaTarget[i];

                    for (int j = i - 1; j >= 0; j--)
                    {
                        if (xi - Math.Log(Math.Max(freqs[j], 1e-9)) > 0.24) break;
                        leftMin = Math.Min(leftMin, alphaTarget[j]);
                    }

                    for (int j = i + 1; j < n; j++)
                    {
                        if (Math.Log(Math.Max(freqs[j], 1e-9)) - xi > 0.24) break;
                        rightMin = Math.Min(rightMin, alphaTarget[j]);
                    }

                    double baseline = 0.5 * (leftMin + rightMin);
                    prominence[i] = Math.Max(0.0, alphaTarget[i] - baseline);

                    // Don't promote tiny numerical wiggles to resonances.
                    if (prominence[i] > 0.015 || -d2[i] / maxD2 > 0.10) peaks.Add(i);
                }

                double maxProminence = Math.Max(1e-9, prominence.Max());

                // Favor the strongest and sharpest resonances rather than simply
                // choosing the four locations having the highest absolute alpha.
                peaks = peaks.OrderByDescending(i => 0.65 * prominence[i] / maxProminence + 0.35 * Math.Max(0.0, -d2[i]) / maxD2).Take(6).ToList();

                for (int i = 0; i < n; i++)
                {
                    double f = freqs[i];

                    double lowBias;
                    if (f < 20.0) lowBias = 2.2;
                    else if (f < 50.0) lowBias = 1.8;
                    else if (f < 100.0) lowBias = 1.4;
                    else if (f < 250.0) lowBias = 1.0;
                    else lowBias = 0.5;

                    // Slope matters around the sides of a resonance, but curvature
                    // matters more when distinguishing a narrow peak from a broad hump.
                    double shape = 0.5 * (alphaTarget[i] / maxAlpha) + 2.0 * (Math.Abs(d1[i]) / maxD1) + 4.0 * (Math.Abs(d2[i]) / maxD2);

                    double localPeak = 0.0;

                    for (int p = 0; p < peaks.Count; p++)
                    {
                        int pi = peaks[p];
                        double fp = freqs[pi];

                        double d = Math.Abs(Math.Log(Math.Max(f, 1e-9)) - Math.Log(Math.Max(fp, 1e-9)));

                        // Narrower than the old 0.10-log-frequency weighting.
                        // This makes the optimizer care about reproducing the peak
                        // itself rather than merely producing energy in its vicinity.
                        double peakShape = Math.Exp(-(d * d) / (2.0 * 0.065 * 0.065));
                        double strength = 0.5 + 1.5 * prominence[pi] / maxProminence + 0.75 * Math.Max(0.0, -d2[pi]) / maxD2;
                        localPeak += strength * peakShape;
                    }

                    w[i] = lowBias + shape + 8.0 * localPeak;
                }

                return w;
            }

            private static double SM_ComputeLogLikelihoodForReflection(BottsModelSpec spec, double[] theta, System.Numerics.Complex[] Rtarget, double[] alphaTarget, double[] freqs, double[] weights, double fs, IList<BottsTargetResonance> resonanceTargets = null)
            {
                const double eps = 1e-8;
                const double lambdaAlpha = 7.5;
                const double lambdaLogMag = 2.5;
                const double lambdaPhase = 0.25;

                // Shape terms are deliberately modest. They should distinguish a
                // narrow resonance from a broad hump, not overwhelm magnitude/phase.
                const double lambdaSlope = 0.10;
                const double lambdaCurvature = 0.15;

                SM_BuildReflectionCoefficients(spec, theta, fs, out double[] bR, out double[] aR);

                int n = freqs.Length;
                double[] alphaFit = new double[n];

                double E = 0.0;
                int K = 0;

                for (int i = 0; i < n; i++)
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

                    double phaseErr = SM_WrapPhase(Math.Atan2(Rtarget[i].Imaginary, Rtarget[i].Real) - Math.Atan2(H.Imaginary, H.Real));

                    alphaFit[i] = 1.0 - H.Magnitude * H.Magnitude;
                    if (alphaFit[i] < 0.0) alphaFit[i] = 0.0;
                    if (alphaFit[i] > 1.0) alphaFit[i] = 1.0;

                    double alphaErr = alphaFit[i] - alphaTarget[i];

                    E += wt * (lambdaAlpha * alphaErr * alphaErr + lambdaLogMag * logMagErr * logMagErr + lambdaPhase * phaseErr * phaseErr);

                    K++;
                }

                if (K == 0) return double.NegativeInfinity;

                // Compare the local shape of the target and fitted alpha curves in
                // log-frequency coordinates. This is especially important for narrow
                // resonant absorbers, where pointwise error alone tends to reward a
                // broader hump.
                if (n > 2)
                {
                    double[] targetSlope = new double[n];
                    double[] fitSlope = new double[n];
                    double[] targetCurvature = new double[n];
                    double[] fitCurvature = new double[n];

                    double maxTargetSlope = 0.0;
                    double maxTargetCurvature = 0.0;

                    for (int i = 1; i < n - 1; i++)
                    {
                        double x0 = Math.Log(Math.Max(freqs[i - 1], 1e-9));
                        double x1 = Math.Log(Math.Max(freqs[i], 1e-9));
                        double x2 = Math.Log(Math.Max(freqs[i + 1], 1e-9));

                        double dx = Math.Max(1e-9, x2 - x0);
                        double dx1 = Math.Max(1e-9, x1 - x0);
                        double dx2 = Math.Max(1e-9, x2 - x1);

                        targetSlope[i] = (alphaTarget[i + 1] - alphaTarget[i - 1]) / dx;
                        fitSlope[i] = (alphaFit[i + 1] - alphaFit[i - 1]) / dx;

                        double targetS1 = (alphaTarget[i] - alphaTarget[i - 1]) / dx1;
                        double targetS2 = (alphaTarget[i + 1] - alphaTarget[i]) / dx2;
                        double fitS1 = (alphaFit[i] - alphaFit[i - 1]) / dx1;
                        double fitS2 = (alphaFit[i + 1] - alphaFit[i]) / dx2;

                        double dxc = Math.Max(1e-9, 0.5 * (dx1 + dx2));

                        targetCurvature[i] = (targetS2 - targetS1) / dxc;
                        fitCurvature[i] = (fitS2 - fitS1) / dxc;

                        maxTargetSlope = Math.Max(maxTargetSlope, Math.Abs(targetSlope[i]));
                        maxTargetCurvature = Math.Max(maxTargetCurvature, Math.Abs(targetCurvature[i]));
                    }

                    // Prevent essentially flat targets from magnifying tiny numerical
                    // slope/curvature differences.
                    maxTargetSlope = Math.Max(0.05, maxTargetSlope);
                    maxTargetCurvature = Math.Max(0.05, maxTargetCurvature);

                    for (int i = 1; i < n - 1; i++)
                    {
                        double targetSlopeStrength = Math.Abs(targetSlope[i]) / maxTargetSlope;
                        double targetCurvatureStrength = Math.Abs(targetCurvature[i]) / maxTargetCurvature;

                        // Only invoke the shape penalty where the target actually has
                        // meaningful structure.
                        double structure = Math.Min(1.0, targetSlopeStrength + targetCurvatureStrength);
                        if (structure < 0.05) continue;

                        double slopeErr = (fitSlope[i] - targetSlope[i]) / maxTargetSlope;
                        double curvatureErr = (fitCurvature[i] - targetCurvature[i]) / maxTargetCurvature;

                        // The ordinary point weights are intentionally capped here:
                        // SM_BuildBandpassWeights already boosts resonances strongly.
                        double shapeWeight = Math.Min(3.0, weights[i]) * structure;

                        E += shapeWeight * (lambdaSlope * slopeErr * slopeErr + lambdaCurvature * curvatureErr * curvatureErr);
                    }
                }

                if (resonanceTargets != null && resonanceTargets.Count > 0)
                {
                    const double lambdaResonanceCenter = 0.25;
                    const double lambdaResonanceWidth = 0.50;

                    double CrossingLog(int i0, int i1, double level)
                    {
                        double x0 = Math.Log(Math.Max(freqs[i0], 1e-9));
                        double x1 = Math.Log(Math.Max(freqs[i1], 1e-9));
                        double y0 = alphaFit[i0];
                        double y1 = alphaFit[i1];

                        if (Math.Abs(y1 - y0) < 1e-12)
                            return 0.5 * (x0 + x1);

                        double t = (level - y0) / (y1 - y0);
                        t = Math.Max(0.0, Math.Min(1.0, t));

                        return x0 + t * (x1 - x0);
                    }

                    foreach (BottsTargetResonance target in resonanceTargets)
                    {
                        int first = -1;
                        int last = -1;
                        int peak = -1;
                        int localCount = 0;
                        double peakAlpha = double.NegativeInfinity;

                        for (int i = 0; i < freqs.Length; i++)
                        {
                            double x = Math.Log(Math.Max(freqs[i], 1e-9));

                            if (x < target.WindowLoLog || x > target.WindowHiLog)
                                continue;

                            if (first < 0) first = i;
                            last = i;
                            localCount++;

                            if (alphaFit[i] > peakAlpha)
                            {
                                peakAlpha = alphaFit[i];
                                peak = i;
                            }
                        }

                        if (first < 0 || last <= first || peak < 0)
                            continue;

                        double leftMin = peakAlpha;
                        double rightMin = peakAlpha;

                        for (int i = first; i <= peak; i++)
                            leftMin = Math.Min(leftMin, alphaFit[i]);

                        for (int i = peak; i <= last; i++)
                            rightMin = Math.Min(rightMin, alphaFit[i]);

                        double baseline = 0.5 * (leftMin + rightMin);
                        double fitProminence = Math.Max(0.0, peakAlpha - baseline);

                        // If the fitted response has effectively lost the resonance,
                        // give it a substantial but finite penalty.
                        if (fitProminence < Math.Max(0.01, 0.10 * target.Prominence))
                        {
                            E += Math.Max(8.0, 0.35 * localCount) *
                                (4.0 * lambdaResonanceCenter + 4.0 * lambdaResonanceWidth);

                            continue;
                        }

                        double halfLevel = baseline + 0.5 * fitProminence;

                        int left = peak;
                        while (left > first && alphaFit[left] > halfLevel)
                            left--;

                        int right = peak;
                        while (right < last && alphaFit[right] > halfLevel)
                            right++;

                        double leftLog;
                        if (left == first && alphaFit[left] > halfLevel)
                            leftLog = target.WindowLoLog;
                        else
                            leftLog = CrossingLog(left, Math.Min(left + 1, peak), halfLevel);

                        double rightLog;
                        if (right == last && alphaFit[right] > halfLevel)
                            rightLog = target.WindowHiLog;
                        else
                            rightLog = CrossingLog(Math.Max(peak, right - 1), right, halfLevel);

                        double fitCenterLog = Math.Log(Math.Max(freqs[peak], 1e-9));
                        double fitWidthLog = Math.Max(1e-4, rightLog - leftLog);

                        // Express center error in units of approximately one target
                        // half-bandwidth, so shifts matter more for narrow resonances.
                        double centerScale = Math.Max(0.025, 0.5 * target.WidthLog);

                        double centerErr = (fitCenterLog - target.CenterLog) / centerScale;

                        // A ratio is more meaningful for bandwidth than an absolute
                        // difference. 2x too broad and 2x too narrow are symmetric.
                        double widthErr =
                            Math.Log(fitWidthLog / Math.Max(1e-4, target.WidthLog));

                        centerErr = Math.Max(-4.0, Math.Min(4.0, centerErr));
                        widthErr = Math.Max(-4.0, Math.Min(4.0, widthErr));

                        // Scale like a local group of frequency samples so this descriptor
                        // has enough influence relative to the existing pointwise errors.
                        double resonanceWeight = Math.Max(8.0, 0.35 * localCount);

                        E += resonanceWeight * (lambdaResonanceCenter * centerErr * centerErr + lambdaResonanceWidth * widthErr * widthErr);
                    }
                }

                if (E <= 1e-30) E = 1e-30;

                return -0.5 * K * Math.Log(E);
            }

            private static BottsSample SM_DrawPriorSample(Random rng, BottsModelSpec spec, System.Numerics.Complex[] Rtarget, double[] alphaTarget, double[] freqs, double[] weights, double fs)
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
                    BottsSample repl = SM_ConstrainedRandomWalk(rng, spec, seed, worst.LogL, Rtarget, alphaTarget, freqs, weights, fs, mhSteps);

                    live[0] = repl;
                    logXPrev = logXNew;
                }

                double[] thetaBest = SM_BestThetaFromSamples(dead, live);
                return (logZ, thetaBest, iter + 1);
            }

            public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            {
                double fs = sample_frequency;
                frequencies = SM_BuildFitFrequencies(fs, max_freq);

                double[] fitFrequencies = (double[])frequencies.Clone();

                if (_iirForced && rec_a != null && rec_b != null && Math.Abs(rec_fs - sample_frequency) < 1e-9 && Math.Abs(rec_maxfreq - max_freq) < 1)
                {
                    return (rec_a, rec_b);
                }

                lock (lock_IIR)
                {
                    int requestedOrder = filter_order;

                    rec_fs = sample_frequency;
                    rec_maxfreq = max_freq;
                    rec_order = requestedOrder;

                    // Keep this aligned with your current Admittance(double frequency) convention.
                    int idx = 17;
                    if (Transfer_FunctionR == null || Transfer_FunctionR.Length == 0) idx = 0;
                    idx = Math.Max(0, Math.Min(idx, Transfer_FunctionR.Length - 1));

                    System.Numerics.Complex[] Rtarget = new System.Numerics.Complex[frequencies.Length];
                    double[] alphaTarget = new double[frequencies.Length];

                    for (int i = 0; i < frequencies.Length; i++)
                    {
                        double f = frequencies[i];

                        System.Numerics.Complex R = new System.Numerics.Complex(Transfer_FunctionR[idx].Interpolate(f), Transfer_FunctionI[idx].Interpolate(f));

                        double mag = R.Magnitude;
                        if (mag > 0.999) R *= 0.999 / mag;

                        Rtarget[i] = R;

                        double r2 = R.Real * R.Real + R.Imaginary * R.Imaginary;
                        alphaTarget[i] = Math.Max(0.0, Math.Min(1.0, 1.0 - r2));
                    }

                    double[] weights = SM_BuildBandpassWeights(frequencies, alphaTarget);

                    // Build resonance descriptors once. These are used only for final
                    // fit refinement, not for the nested-sampling evidence calculation.
                    List<BottsTargetResonance> resonanceTargets = new List<BottsTargetResonance>();

                    double TargetCrossingLog(int i0, int i1, double level)
                    {
                        double x0 = Math.Log(Math.Max(fitFrequencies[i0], 1e-9));
                        double x1 = Math.Log(Math.Max(fitFrequencies[i1], 1e-9));
                        double y0 = alphaTarget[i0];
                        double y1 = alphaTarget[i1];

                        if (Math.Abs(y1 - y0) < 1e-12) return 0.5 * (x0 + x1);

                        double t = (level - y0) / (y1 - y0);
                        t = Math.Max(0.0, Math.Min(1.0, t));

                        return x0 + t * (x1 - x0);
                    }

                    for (int peakIndex = 1; peakIndex < alphaTarget.Length - 1; peakIndex++)
                    {
                        if (alphaTarget[peakIndex] < alphaTarget[peakIndex - 1] || alphaTarget[peakIndex] < alphaTarget[peakIndex + 1]) continue;

                        double peakF = fitFrequencies[peakIndex];
                        double windowLo = peakF / Math.Sqrt(2.0);
                        double windowHi = peakF * Math.Sqrt(2.0);

                        int first = peakIndex;
                        int last = peakIndex;

                        while (first > 0 && fitFrequencies[first] > windowLo) first--;
                        while (last < fitFrequencies.Length - 1 && fitFrequencies[last] < windowHi) last++;

                        double leftMin = alphaTarget[peakIndex];
                        double rightMin = alphaTarget[peakIndex];

                        for (int i = first; i <= peakIndex; i++) leftMin = Math.Min(leftMin, alphaTarget[i]);

                        for (int i = peakIndex; i <= last; i++) rightMin = Math.Min(rightMin, alphaTarget[i]);

                        double baseline = 0.5 * (leftMin + rightMin);
                        double peakProminence = alphaTarget[peakIndex] - baseline;

                        if (peakProminence < 0.03) continue;

                        double halfLevel = baseline + 0.5 * peakProminence;

                        int left = peakIndex;
                        while (left > first && alphaTarget[left] > halfLevel) left--;

                        int right = peakIndex;
                        while (right < last && alphaTarget[right] > halfLevel) right++;

                        double leftLog = left == first && alphaTarget[left] > halfLevel ? Math.Log(Math.Max(windowLo, 1e-9)) : TargetCrossingLog(left, Math.Min(left + 1, peakIndex), halfLevel);
                        double rightLog = right == last && alphaTarget[right] > halfLevel ? Math.Log(Math.Max(windowHi, 1e-9)) : TargetCrossingLog(Math.Max(peakIndex, right - 1), right, halfLevel);

                        double widthLog = rightLog - leftLog;
                        if (widthLog <= 1e-4) continue;

                        resonanceTargets.Add(new BottsTargetResonance{CenterLog = Math.Log(Math.Max(peakF, 1e-9)), WidthLog = widthLog, WindowLoLog = Math.Log(Math.Max(windowLo, 1e-9)), WindowHiLog = Math.Log(Math.Max(windowHi, 1e-9)), Prominence = peakProminence});
                    }

                    resonanceTargets = resonanceTargets.OrderByDescending(r => r.Prominence).Take(6).ToList();


                    double ScoreTheta(BottsModelSpec spec, double[] theta)
                    {
                        return SM_ComputeLogLikelihoodForReflection(spec, theta, Rtarget, alphaTarget, fitFrequencies, weights, fs, resonanceTargets);
                    }


                    // Refine any valid Botts starting point. This contains both the guided
                    // resonant-section adjustment and the joint Nelder-Mead optimization.
                    double[] RefineTheta(BottsModelSpec spec, double[] startTheta)
                    {
                        if (startTheta == null) return null;

                        double[] refined = (double[])startTheta.Clone();
                        double refinedLogL = ScoreTheta(spec, refined);

                        SM_GetBounds(spec, fs, fitFrequencies.Last(), out double[] lower, out double[] upper);

                        // First point the available complex sections toward the strongest
                        // target resonances. Every candidate must improve the complete score.
                        if (spec.ComplexSections > 0 && resonanceTargets.Count > 0)
                        {
                            bool[] usedSection = new bool[spec.ComplexSections];

                            double[] centerScale = { 0.94, 0.97, 1.00, 1.03, 1.06 };
                            double[] zeroBWScale = { 0.10, 0.20, 0.35, 0.50, 0.75 };
                            double[] poleBWScale = { 0.50, 0.75, 1.00, 1.50, 2.00, 3.00 };
                            double[] gainDb = { -1.5, 0.0, 1.5 };

                            int complexStart = 1 + 2 * spec.RealSections;

                            for (int r = 0; r < Math.Min(resonanceTargets.Count, spec.ComplexSections); r++)
                            {
                                BottsTargetResonance resonance = resonanceTargets[r];

                                double peakF = Math.Exp(resonance.CenterLog);

                                // Convert the measured log-frequency half-prominence width
                                // to an approximate bandwidth in Hz.
                                double halfWidthLog = 0.5 * resonance.WidthLog;
                                double width = peakF * (Math.Exp(halfWidthLog) - Math.Exp(-halfWidthLog));
                                width = Math.Max(1.0, width);

                                double localBestLogL = refinedLogL;
                                double[] localBestTheta = null;
                                int localBestSection = -1;

                                for (int section = 0; section < spec.ComplexSections; section++)
                                {
                                    if (usedSection[section]) continue;

                                    int k = complexStart + 4 * section;

                                    for (int cf = 0; cf < centerScale.Length; cf++)
                                    {
                                        double f = peakF * centerScale[cf];
                                        f = Math.Max(Math.Exp(lower[k]), Math.Min(Math.Exp(upper[k]), f));

                                        for (int z = 0; z < zeroBWScale.Length; z++)
                                        {
                                            double bwZero = Math.Max(1.0, width * zeroBWScale[z]);

                                            for (int p = 0; p < poleBWScale.Length; p++)
                                            {
                                                double bwPole = Math.Max(0.5, width * poleBWScale[p]);

                                                for (int g = 0; g < gainDb.Length; g++)
                                                {
                                                    double[] candidate = (double[])refined.Clone();

                                                    // Keep zero/pole center together during this coarse
                                                    // guided placement. Nelder-Mead may separate them later.
                                                    candidate[k] = Math.Log(f);
                                                    candidate[k + 1] = Math.Log(bwZero);
                                                    candidate[k + 2] = Math.Log(f);
                                                    candidate[k + 3] = Math.Log(bwPole);

                                                    candidate[k + 1] = Math.Max(lower[k + 1], Math.Min(upper[k + 1], candidate[k + 1]));
                                                    candidate[k + 3] = Math.Max(lower[k + 3], Math.Min(upper[k + 3], candidate[k + 3]));

                                                    candidate[0] = refined[0] + gainDb[g] * Math.Log(10.0) / 20.0;
                                                    candidate[0] = Math.Max(lower[0], Math.Min(upper[0], candidate[0]));

                                                    double logL = ScoreTheta(spec, candidate);

                                                    if (logL > localBestLogL)
                                                    {
                                                        localBestLogL = logL;
                                                        localBestTheta = candidate;
                                                        localBestSection = section;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                if (localBestTheta != null)
                                {
                                    refined = localBestTheta;
                                    refinedLogL = localBestLogL;

                                    if (localBestSection >= 0) usedSection[localBestSection] = true;
                                }
                            }
                        }

                        // Joint local optimization: allow every pole, zero, bandwidth and
                        // gain parameter to rebalance simultaneously.
                        double[] perturbation = new double[refined.Length];
                        int pk = 0;

                        perturbation[pk++] = 0.5 * Math.Log(10.0) / 20.0;

                        for (int i = 0; i < spec.RealSections; i++)
                        {
                            perturbation[pk++] = 0.025;
                            perturbation[pk++] = 0.025;
                        }

                        for (int i = 0; i < spec.ComplexSections; i++)
                        {
                            perturbation[pk++] = Math.Log(1.04);
                            perturbation[pk++] = Math.Log(1.20);
                            perturbation[pk++] = Math.Log(1.04);
                            perturbation[pk++] = Math.Log(1.20);
                        }

                        var initialVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(refined);

                        var perturbationVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(perturbation);

                        var objective = MathNet.Numerics.Optimization.ObjectiveFunction.Value(v =>
                        {
                            double[] candidate = v.ToArray();
                            for (int i = 0; i < candidate.Length; i++) candidate[i] = SM_ReflectToBounds(candidate[i], lower[i], upper[i]);

                            double logL = ScoreTheta(spec, candidate);
                            if (double.IsNaN(logL) || double.IsInfinity(logL)) return 1e100;
                            return -logL;
                        });

                        try
                        {
                            var nm = MathNet.Numerics.Optimization.NelderMeadSimplex.Minimum(objective, initialVector, perturbationVector, 1e-6, 1000);

                            double[] nmTheta = nm.MinimizingPoint.ToArray();
                            for (int i = 0; i < nmTheta.Length; i++) nmTheta[i] = SM_ReflectToBounds(nmTheta[i], lower[i], upper[i]);

                            double nmLogL = ScoreTheta(spec, nmTheta);

                            // Never let local refinement make the starting solution worse.
                            if (!double.IsNaN(nmLogL) && !double.IsInfinity(nmLogL) && nmLogL > refinedLogL)  refined = nmTheta;
                        }
                        catch
                        {
                            // Retain the valid pre-Nelder-Mead solution.
                        }

                        return refined;
                    }

                    // Embed a lower-order solution into the next model of the same parity.
                    // The added pole/zero section is exactly neutral because its numerator
                    // and denominator are identical.
                    double[] ExpandTheta(BottsModelSpec oldSpec, double[] oldTheta, BottsModelSpec newSpec)
                    {
                        if (oldTheta == null || oldSpec.RealSections != newSpec.RealSections || newSpec.ComplexSections != oldSpec.ComplexSections + 1) return null;

                        double[] theta = new double[newSpec.ParameterCount];
                        Array.Copy(oldTheta, theta, oldTheta.Length);

                        SM_GetBounds(newSpec, fs, fitFrequencies.Last(), out double[] lower, out double[] upper);

                        int k = oldTheta.Length;
                        double center = resonanceTargets.Count > 0 ? Math.Exp(resonanceTargets[0].CenterLog) : Math.Min(125.0, fitFrequencies.Last());
                        center = Math.Max(Math.Exp(lower[k]), Math.Min(Math.Exp(upper[k]), center));

                        double bandwidth = Math.Max(1.0, 0.25 * center);
                        bandwidth = Math.Max(Math.Exp(lower[k + 1]), Math.Min(Math.Exp(upper[k + 1]), bandwidth));

                        theta[k] = Math.Log(center);
                        theta[k + 1] = Math.Log(bandwidth);
                        theta[k + 2] = Math.Log(center);
                        theta[k + 3] = Math.Log(bandwidth);

                        return theta;
                    }

                    // For an explicit order, calculate the preceding same-parity models as
                    // stepping stones. Auto evaluates the complete 2..10 family.
                    List<int> orders = new List<int>();
                    double[] bestTheta = null;
                    BottsModelSpec bestSpec = null;

                    if (requestedOrder > 0)
                    {
                        if (requestedOrder <= 2)
                        {
                            orders.Add(requestedOrder);
                        }
                        else
                        {
                            int firstOrder = requestedOrder % 2 == 0 ? 2 : 3;
                            for (int n = firstOrder; n <= requestedOrder; n += 2) orders.Add(n);
                        }
                    }
                    else
                    {
                        for (int n = 2; n <= 10; n++) orders.Add(n);
                    }

                    // First run the ordinary Botts nested sampler for every required order.
                    // Keep these searches independent and deterministic.
                    var nestedResults = new Dictionary<int, (double logZ, double[] theta, BottsModelSpec spec)>();

                    object nestedLock = new object();
                    double[] freqs = (double[])fitFrequencies.Clone();

                    Parallel.ForEach(
                        orders,
                        new ParallelOptions
                        {
                            MaxDegreeOfParallelism = Math.Max(1, System.Environment.ProcessorCount - 1)
                        },
                        order =>
                        {
                            Random rngLocal = new Random(unchecked(7919 + order * 104729));
                            BottsModelSpec spec = SM_MakeSpec(order);
                            var ns = SM_RunNestedSampling(rngLocal, spec, Rtarget, alphaTarget, freqs, weights, fs, nLive: 96, mhSteps: 100, maxIter: 12000, stopDeltaLogZ: 1e-7);
                            lock (nestedLock) nestedResults[order] = (ns.logZ, ns.thetaBest, spec);
                        });


                    // Now walk each parity chain. A higher-order model gets both:
                    // 1. its independently discovered Botts solution, and
                    // 2. the already-refined lower-order solution plus a neutral section.
                    //
                    // Whichever refines to the better complete fit becomes the solution
                    // inherited by the next model.
                    var finalResults = new Dictionary<int, (double logZ, double logL, double[] theta, BottsModelSpec spec)>();

                    void ProcessChain(IEnumerable<int> chain)
                    {
                        double[] previousTheta = null;
                        BottsModelSpec previousSpec = null;

                        foreach (int order in chain)
                        {
                            if (!nestedResults.TryGetValue(order, out var nested)) continue;

                            double[] bestOrderTheta = null;
                            double bestOrderLogL = double.NegativeInfinity;

                            // Independent nested-sampling candidate.
                            if (nested.theta != null)
                            {
                                double[] candidate = RefineTheta(nested.spec, nested.theta);

                                if (candidate != null)
                                {
                                    double logL = ScoreTheta(nested.spec, candidate);

                                    if (logL > bestOrderLogL)
                                    {
                                        bestOrderLogL = logL;
                                        bestOrderTheta = candidate;
                                    }
                                }
                            }

                            // Inherited candidate from the preceding order of the same parity.
                            if (previousTheta != null && previousSpec != null)
                            {
                                double[] inherited = ExpandTheta(previousSpec, previousTheta, nested.spec);

                                if (inherited != null)
                                {
                                    // Before optimization this produces exactly the old response;
                                    // the added section cancels identically.
                                    double[] candidate = RefineTheta(nested.spec, inherited);

                                    if (candidate != null)
                                    {
                                        double logL = ScoreTheta(nested.spec, candidate);

                                        if (logL > bestOrderLogL)
                                        {
                                            bestOrderLogL = logL;
                                            bestOrderTheta = candidate;
                                        }
                                    }
                                }
                            }

                            if (bestOrderTheta == null)
                                continue;

                            finalResults[order] =
                                (nested.logZ, bestOrderLogL, bestOrderTheta, nested.spec);

                            previousTheta = bestOrderTheta;
                            previousSpec = nested.spec;
                        }
                    }

                    if (requestedOrder > 0)
                    {
                        ProcessChain(orders);
                    }
                    else
                    {
                        ProcessChain(new int[] { 2, 4, 6, 8, 10 });
                        ProcessChain(new int[] { 3, 5, 7, 9 });
                    }

                    if (finalResults.Count == 0)
                    {
                        rec_a = new double[] { 1.0 };
                        rec_b = new double[] { 0.0 };
                        return (rec_a, rec_b);
                    }

                    // Explicit order means exactly that order.
                    if (requestedOrder > 0)
                    {
                        if (!finalResults.TryGetValue(requestedOrder, out var selectedExplicit))
                        {
                            rec_a = new double[] { 1.0 };
                            rec_b = new double[] { 0.0 };
                            return (rec_a, rec_b);
                        }

                        bestTheta = selectedExplicit.theta;
                        bestSpec = selectedExplicit.spec;
                    }
                    else
                    {
                        double bestFitLogL = finalResults.Values.Max(r => r.logL);

                        int K = Math.Max(1, weights.Count(w => w > 0.0));
                        double fitTolerance = 0.5 * K * Math.Log(1.01);

                        var acceptable = finalResults.Where(r => r.Value.logL >= bestFitLogL - fitTolerance).ToList();
                        var selectedAuto = acceptable.OrderByDescending(r => r.Value.logZ).ThenBy(r => r.Key).First();

                        bestTheta = selectedAuto.Value.theta;
                        bestSpec = selectedAuto.Value.spec;
                    }
                    SM_BuildReflectionCoefficients(bestSpec, bestTheta, fs, out double[] bR, out double[] aR);
                    SM_ConvertReflectionToAdmittance(bR, aR, out double[] aY, out double[] bY);

                    rec_a = aY;
                    rec_b = bY;

                    return (rec_a, rec_b);
                }
            }
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
                    z = Math.Cos(u2);

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
    }
}