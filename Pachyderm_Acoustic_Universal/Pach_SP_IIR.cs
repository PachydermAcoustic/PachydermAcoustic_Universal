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

using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

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
                    double[] H = new double[length/2];
                    double w_l = Utilities.Numerics.PiX2 * freq_l;
                    double w_u = Utilities.Numerics.PiX2 * freq_u;
                    double max = 0;

                    for (int i = 0; i < H.Length; i++)
                    {
                        double w_f = Utilities.Numerics.PiX2 * i * sampling_freq / length;
                        H[i] = 1 / (Math.Sqrt(1 + Math.Pow(w_f / w_l, 2 * order)) * Math.Sqrt(1 + Math.Pow(w_u / w_f, 2 * order)));
                        max = max > H[i]? max : H[i];
                    }
                    for (int i = 0; i < H.Length; i++) H[i] /= max;

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

                /// <summary>
                /// Fits IIR filter coefficients to match a given frequency response
                public static (double[] b, double[] a) FitIIRToFrequencyResponse(double[] frequencies, Complex[] desiredResponse, int numeratorOrder, int denominatorOrder, double fs)
                {
                    if (frequencies == null || desiredResponse == null)
                        throw new ArgumentNullException();
                    if (frequencies.Length != desiredResponse.Length)
                        throw new ArgumentException("frequencies and desiredResponse must have same length");
                    if (frequencies.Length < 4)
                        throw new ArgumentException("Need at least a few frequency points");

                    // Weights: emphasize low-ish frequencies but avoid DC singularities and runaway dominance
                    double[] weights = new double[frequencies.Length];
                    for (int i = 0; i < frequencies.Length; i++)
                    {
                        double f = Math.Max(frequencies[i], 20.0); // avoid DC blow-up
                        double w = 1.0 / Math.Sqrt(f);
                        weights[i] = Math.Min(w, 0.25);            // cap
                    }

                    int rows = frequencies.Length * 2;
                    int cols = numeratorOrder + denominatorOrder + 1;

                    var A = Matrix<double>.Build.Dense(rows, cols);
                    var rhs = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(rows);

                    for (int i = 0; i < frequencies.Length; i++)
                    {
                        double omega = 2.0 * Math.PI * frequencies[i] / fs;
                        double weight = weights[i];

                        double Dr = desiredResponse[i].Real;
                        double Di = desiredResponse[i].Imaginary;

                        // Real equation row
                        for (int j = 0; j <= numeratorOrder; j++)
                            A[i, j] = weight * Math.Cos(-j * omega);

                        for (int j = 0; j <= denominatorOrder - 1; j++)
                        {
                            double c = Math.Cos(-(j + 1) * omega);
                            double s = Math.Sin(-(j + 1) * omega);
                            A[i, numeratorOrder + 1 + j] = -weight * Dr * c + weight * Di * s;
                        }

                        rhs[i] = weight * Dr;

                        // Imag equation row
                        int ii = i + frequencies.Length;
                        for (int j = 0; j <= numeratorOrder; j++)
                            A[ii, j] = weight * Math.Sin(-j * omega);

                        for (int j = 0; j <= denominatorOrder - 1; j++)
                        {
                            double c = Math.Cos(-(j + 1) * omega);
                            double s = Math.Sin(-(j + 1) * omega);
                            A[ii, numeratorOrder + 1 + j] = -weight * Dr * s - weight * Di * c;
                        }

                        rhs[ii] = weight * Di;
                    }

                    var x = A.QR().Solve(rhs);

                    double[] numerator = new double[numeratorOrder + 1];
                    double[] denominator = new double[denominatorOrder + 1];

                    for (int i = 0; i <= numeratorOrder; i++)
                        numerator[i] = x[i];

                    denominator[0] = 1.0;
                    int k = numeratorOrder + 1;
                    for (int i = 1; i <= denominatorOrder; i++)
                        denominator[i] = x[k++];

                    EnsureStability(ref denominator);

                    return (numerator, denominator);
                }


                /// <summary>
                /// Ensures filter stability by reflecting poles outside the unit circle.
                /// Uses a self-contained QR eigenvalue solver to avoid MathNet root-finding failures.
                /// </summary>
                public static void EnsureStability(ref double[] a)
                {
                    if (a == null || a.Length < 2) return;

                    // Check for NaN or Infinity in coefficients
                    for (int i = 0; i < a.Length; i++)
                    {
                        if (double.IsNaN(a[i]) || double.IsInfinity(a[i]))
                        {
                            for (int j = 0; j < a.Length; j++)
                                a[j] = (j == 0) ? 1.0 : 0.0; 
                            return;
                        }
                    }

                    int n = a.Length - 1; // polynomial order
                    if (n < 1) return;

                    // Normalize so a[0] = 1
                    if (Math.Abs(a[0]) < 1e-15)
                    {
                        for (int i = 0; i < a.Length; i++)
                            a[i] = (i == 0) ? 1.0 : 0.0;
                        return;
                    }
                    double inv_a0 = 1.0 / a[0];
                    for (int i = 0; i < a.Length; i++) a[i] *= inv_a0;

                    // Build companion matrix for z^n + a[1]*z^(n-1) + ... + a[n] = 0
                    double[,] companion = new double[n, n];
                    for (int i = 0; i < n - 1; i++)
                        companion[i, i + 1] = 1.0;
                    for (int i = 0; i < n; i++)
                        companion[n - 1, i] = -a[n - i];

                    Complex[] poles = EigenvaluesHessenbergQR(companion, n);
                    if (poles == null)
                    {
                        for (int i = 0; i < a.Length; i++)
                            a[i] = (i == 0) ? 1.0 : 0.0;
                        return;
                    }

                    // Check for invalid poles
                    for (int i = 0; i < poles.Length; i++)
                    {
                        if (double.IsNaN(poles[i].Real) || double.IsNaN(poles[i].Imaginary) ||
                            double.IsInfinity(poles[i].Real) || double.IsInfinity(poles[i].Imaginary))
                        {
                            for (int j = 0; j < a.Length; j++)
                                a[j] = (j == 0) ? 1.0 : 0.0; return;
                        }
                    }

                    // Reflect any poles outside the unit circle
                    bool modified = false;
                    for (int i = 0; i < poles.Length; i++)
                    {
                        if (poles[i].Magnitude > 1.0)
                        {
                            modified = true;
                            poles[i] = 1.0 / Complex.Conjugate(poles[i]);
                        }
                    }

                    if (!modified) return;

                    // Reconstruct a(z) from stabilized poles, pairing conjugates
                    bool[] used = new bool[poles.Length];
                    double[] result = new double[] { 1.0 };

                    for (int i = 0; i < poles.Length; i++)
                    {
                        if (used[i]) continue;
                        used[i] = true;

                        if (Math.Abs(poles[i].Imaginary) > 1e-10)
                        {
                            // Find conjugate partner
                            int conj = -1;
                            for (int j = i + 1; j < poles.Length; j++)
                            {
                                if (!used[j] &&
                                    Math.Abs(poles[j].Real - poles[i].Real) < 1e-10 &&
                                    Math.Abs(poles[j].Imaginary + poles[i].Imaginary) < 1e-10)
                                {
                                    conj = j;
                                    break;
                                }
                            }

                            if (conj >= 0)
                            {
                                used[conj] = true;
                                double mag2 = poles[i].Real * poles[i].Real + poles[i].Imaginary * poles[i].Imaginary;
                                result = ConvolvePolynomials(result, new double[] { 1.0, -2.0 * poles[i].Real, mag2 });
                            }
                            else
                            {
                                // Orphan complex pole — treat as real
                                result = ConvolvePolynomials(result, new double[] { 1.0, -poles[i].Real });
                            }
                        }
                        else
                        {
                            result = ConvolvePolynomials(result, new double[] { 1.0, -poles[i].Real });
                        }
                    }

                    // Validate result
                    for (int i = 0; i < result.Length; i++)
                    {
                        if (double.IsNaN(result[i]) || double.IsInfinity(result[i]))
                        {
                            for (int j = 0; j < a.Length; j++)
                                a[j] = (j == 0) ? 1.0 : 0.0; return;
                        }
                    }

                    // Copy back
                    for (int i = 0; i < a.Length; i++)
                        a[i] = (i < result.Length) ? result[i] : 0.0;
                }

                private static double[] ConvolvePolynomials(double[] p, double[] q)
                {
                    double[] r = new double[p.Length + q.Length - 1];
                    for (int i = 0; i < p.Length; i++)
                        for (int j = 0; j < q.Length; j++)
                            r[i + j] += p[i] * q[j];
                    return r;
                }

                /// <summary>
                /// QR iteration on a real matrix to find eigenvalues (poles).
                /// Reduces to upper Hessenberg form first, then iterates with Wilkinson shift.
                /// Self-contained — no MathNet dependency.
                /// </summary>
                private static Complex[] EigenvaluesHessenbergQR(double[,] M, int n)
                {
                    double[,] H = (double[,])M.Clone();

                    // Reduce to upper Hessenberg via Householder
                    for (int k = 0; k < n - 2; k++)
                    {
                        double norm = 0;
                        for (int i = k + 1; i < n; i++) norm += H[i, k] * H[i, k];
                        norm = Math.Sqrt(norm);
                        if (norm < 1e-15) continue;

                        double sign = (H[k + 1, k] >= 0) ? 1.0 : -1.0;
                        double[] v = new double[n];
                        v[k + 1] = H[k + 1, k] + sign * norm;
                        for (int i = k + 2; i < n; i++) v[i] = H[i, k];

                        double vv = 0;
                        for (int i = k + 1; i < n; i++) vv += v[i] * v[i];
                        if (vv < 1e-30) continue;
                        double scale = 2.0 / vv;

                        for (int j = k; j < n; j++)
                        {
                            double dot = 0;
                            for (int i = k + 1; i < n; i++) dot += v[i] * H[i, j];
                            dot *= scale;
                            for (int i = k + 1; i < n; i++) H[i, j] -= dot * v[i];
                        }
                        for (int i = 0; i < n; i++)
                        {
                            double dot = 0;
                            for (int j = k + 1; j < n; j++) dot += H[i, j] * v[j];
                            dot *= scale;
                            for (int j = k + 1; j < n; j++) H[i, j] -= dot * v[j];
                        }
                    }

                    // QR iteration with Wilkinson shift
                    var eigenvalues = new Complex[n];
                    int nn = n;
                    int maxIter = 200 * n;
                    int iter = 0;

                    while (nn > 0 && iter < maxIter)
                    {
                        iter++;

                        if (nn == 1)
                        {
                            eigenvalues[nn - 1] = H[nn - 1, nn - 1];
                            nn--;
                            continue;
                        }

                        double tol = 1e-14 * (Math.Abs(H[nn - 2, nn - 2]) + Math.Abs(H[nn - 1, nn - 1]));
                        if (tol < 1e-30) tol = 1e-30;

                        if (Math.Abs(H[nn - 1, nn - 2]) < tol)
                        {
                            eigenvalues[nn - 1] = H[nn - 1, nn - 1];
                            nn--;
                            continue;
                        }

                        // 2×2 block deflation
                        if (nn == 2 || (nn >= 3 && Math.Abs(H[nn - 2, nn - 3]) <
                            1e-14 * (Math.Abs(H[nn - 3, nn - 3]) + Math.Abs(H[nn - 2, nn - 2]))))
                        {
                            double a11 = H[nn - 2, nn - 2], a12 = H[nn - 2, nn - 1];
                            double a21 = H[nn - 1, nn - 2], a22 = H[nn - 1, nn - 1];
                            double tr = a11 + a22;
                            double det = a11 * a22 - a12 * a21;
                            double disc = tr * tr - 4 * det;
                            if (disc >= 0)
                            {
                                double sqd = Math.Sqrt(disc);
                                eigenvalues[nn - 1] = (tr + sqd) / 2;
                                eigenvalues[nn - 2] = (tr - sqd) / 2;
                            }
                            else
                            {
                                double sqd = Math.Sqrt(-disc);
                                eigenvalues[nn - 1] = new Complex(tr / 2, sqd / 2);
                                eigenvalues[nn - 2] = new Complex(tr / 2, -sqd / 2);
                            }
                            nn -= 2;
                            continue;
                        }

                        // Wilkinson shift
                        double p = H[nn - 2, nn - 2], q = H[nn - 2, nn - 1];
                        double r = H[nn - 1, nn - 2], ss = H[nn - 1, nn - 1];
                        double shift = ss;
                        double dd = (p - ss) / 2.0;
                        double sqr = dd * dd + q * r;
                        if (sqr >= 0)
                        {
                            double sqrtSqr = Math.Sqrt(sqr);
                            shift = ss - r * q / (dd + (dd >= 0 ? sqrtSqr : -sqrtSqr));
                        }

                        // Givens rotations
                        for (int i = 0; i < nn - 1; i++)
                        {
                            double x = (i == 0) ? H[0, 0] - shift : H[i, i - 1];
                            double z = (i == 0) ? H[1, 0] : H[i + 1, i - 1];
                            double rr = Math.Sqrt(x * x + z * z);
                            if (rr < 1e-30) continue;
                            double cs = x / rr, sn = z / rr;

                            for (int j = Math.Max(0, i - 1); j < nn; j++)
                            {
                                double t1 = H[i, j], t2 = H[i + 1, j];
                                H[i, j] = cs * t1 + sn * t2;
                                H[i + 1, j] = -sn * t1 + cs * t2;
                            }
                            for (int j = 0; j <= Math.Min(i + 2, nn - 1); j++)
                            {
                                double t1 = H[j, i], t2 = H[j, i + 1];
                                H[j, i] = cs * t1 + sn * t2;
                                H[j, i + 1] = -sn * t1 + cs * t2;
                            }
                        }
                    }

                    return (iter >= maxIter) ? null : eigenvalues;
                }

                public static double CompareFunctions(Complex[] f1, Complex[] f2, int[] sampled)
                {
                    if (f1.Length != f2.Length) throw new Exception("f1 and f2 must be equal in length");

                    double Corr = 0;

                    for (int i = 0; i < f1.Length; i++)
                    {
                        Corr += Math.Pow(f1[i].Real - f2[i].Real, 2) + Math.Pow(f1[i].Imaginary - f2[i].Imaginary, 2);
                    }

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
            }
        }
    }
}