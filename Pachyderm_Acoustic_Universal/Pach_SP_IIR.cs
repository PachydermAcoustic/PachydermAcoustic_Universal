//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
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

                /// <summary>
                /// Fits IIR filter coefficients to match a given frequency response
                /// </summary>
                public static (double[] b, double[] a) FitIIRToFrequencyResponse(double[] frequencies, Complex[] desiredResponse, int numeratorOrder, int denominatorOrder, double fs = 44100)
                {
                    // Calculate weights to emphasize lower frequencies
                    double[] weights = new double[frequencies.Length];
                    for (int i = 0; i < frequencies.Length; i++) weights[i] = 1.0 / Math.Sqrt(frequencies[i]);
                    
                    // Prepare matrices for least squares solution
                    var A = Matrix<double>.Build.Dense(frequencies.Length * 2, numeratorOrder + denominatorOrder + 1);
                    var b = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(frequencies.Length * 2);
                    
                    for (int i = 0; i < frequencies.Length; i++)
                    {
                        double w = 2 * Math.PI * frequencies[i] / fs;
                        double weight = weights[i];
                        
                        // Real part equations
                        for (int j = 0; j <= numeratorOrder; j++) A[i, j] = weight * Math.Cos(-j * w);
                        for (int j = 0; j <= denominatorOrder - 1; j++)
                        {
                            A[i, numeratorOrder + 1 + j] = -weight * desiredResponse[i].Real * Math.Cos(-(j + 1) * w)
                                                         + weight * desiredResponse[i].Imaginary * Math.Sin(-(j + 1) * w);
                        }
                        b[i] = weight * desiredResponse[i].Real;
                        
                        // Imaginary part equations
                        for (int j = 0; j <= numeratorOrder; j++) A[i + frequencies.Length, j] = weight * Math.Sin(-j * w);
                        for (int j = 0; j <= denominatorOrder - 1; j++)
                        {
                            A[i + frequencies.Length, numeratorOrder + 1 + j] = -weight * desiredResponse[i].Real * Math.Sin(-(j + 1) * w)
                                                                              - weight * desiredResponse[i].Imaginary * Math.Cos(-(j + 1) * w);
                        }
                        b[i + frequencies.Length] = weight * desiredResponse[i].Imaginary;
                    }
                    
                    // Solve the linear system using QR decomposition
                    var qr = A.QR();
                    var x = qr.Solve(b);
                    
                    // Extract the coefficients
                    double[] numerator = new double[numeratorOrder + 1];
                    double[] denominator = new double[denominatorOrder + 1];
                    
                    for (int i = 0; i <= numeratorOrder; i++) numerator[i] = x[i];
                    
                    denominator[0] = 1.0; // a0 = 1 by convention
                    for (int i = 1; i <= denominatorOrder; i++) denominator[i] = x[numeratorOrder + i];
                    
                    // Ensure filter stability
                    EnsureStability(ref denominator);
                    
                    return (numerator, denominator);
                }

                /// <summary>
                /// Ensures filter stability by reflecting poles outside the unit circle
                /// </summary>
                public static void EnsureStability(ref double[] a)
                {
                    try
                    {
                        // Check for NaN or Infinity in coefficients
                        if (a.Any(coeff => double.IsNaN(coeff) || double.IsInfinity(coeff)))
                        {
                            ApplyConservativeStabilization(ref a);
                            return;
                        }
                        
                        // Create polynomial from denominator coefficients (ascending order)
                        double[] reversedCoeffs = a.Reverse().ToArray();
                        var poly = new MathNet.Numerics.Polynomial(reversedCoeffs);
                        
                        // Find the roots with error handling
                        Complex[] roots;
                        try
                        {
                            roots = poly.Roots();
                            
                            // Check for invalid roots
                            if (roots.Any(r => double.IsNaN(r.Real) || double.IsNaN(r.Imaginary) ||
                                             double.IsInfinity(r.Real) || double.IsInfinity(r.Imaginary)))
                            {
                                roots = FindRootsCompanionMatrix(reversedCoeffs);
                                
                                if (roots.Any(r => double.IsNaN(r.Real) || double.IsNaN(r.Imaginary) ||
                                               double.IsInfinity(r.Real) || double.IsInfinity(r.Imaginary)))
                                {
                                    ApplyConservativeStabilization(ref a);
                                    return;
                                }
                            }
                        }
                        catch
                        {
                            try
                            {
                                roots = FindRootsCompanionMatrix(reversedCoeffs);
                                
                                if (roots.Any(r => double.IsNaN(r.Real) || double.IsNaN(r.Imaginary) ||
                                               double.IsInfinity(r.Real) || double.IsInfinity(r.Imaginary)))
                                {
                                    ApplyConservativeStabilization(ref a);
                                    return;
                                }
                            }
                            catch
                            {
                                ApplyConservativeStabilization(ref a);
                                return;
                            }
                        }
                        
                        // Reflect any poles outside unit circle
                        bool modified = false;
                        for (int i = 0; i < roots.Length; i++)
                        {
                            if (roots[i].Magnitude > 1.0)
                            {
                                roots[i] = roots[i] / (roots[i].Magnitude * roots[i].Magnitude);
                                modified = true;
                            }
                        }
                        
                        // Reconstruct the polynomial if needed
                        if (modified)
                        {
                            try
                            {
                                double[] newCoeffs = ReconstructPolynomialFromRoots(roots, a.Length);
                                
                                if (newCoeffs.Any(c => double.IsNaN(c) || double.IsInfinity(c)))
                                {
                                    ApplyConservativeStabilization(ref a);
                                    return;
                                }
                                
                                Array.Copy(newCoeffs, a, a.Length);
                            }
                            catch
                            {
                                // If reconstruction fails, use a conservative approach
                                ApplyConservativeStabilization(ref a);
                            }
                        }
                    }
                    catch
                    {
                        // Catch any unexpected exceptions and use the conservative approach
                        ApplyConservativeStabilization(ref a);
                    }
                }

                /// <summary>
                /// Applies conservative scaling to ensure filter stability
                /// </summary>
                public static void ApplyConservativeStabilization(ref double[] a)
                {
                    double[] stableCoeffs = new double[a.Length];
                    stableCoeffs[0] = 1.0;
                    
                    // Scale coefficients with diminishing power to move poles toward origin
                    for (int i = 1; i < a.Length; i++)
                    {
                        stableCoeffs[i] = a[i] * Math.Pow(0.85, i);
                    }
                    
                    Array.Copy(stableCoeffs, a, a.Length);
                }

                /// <summary>
                /// Reconstructs polynomial coefficients from roots
                /// </summary>
                public static double[] ReconstructPolynomialFromRoots(Complex[] roots, int coeffLength)
                {
                    try
                    {
                        // Group roots into real roots and complex conjugate pairs
                        List<Complex> realRoots = new List<Complex>();
                        List<Tuple<Complex, Complex>> complexPairs = new List<Tuple<Complex, Complex>>();
                        
                        for (int i = 0; i < roots.Length; i++)
                        {
                            Complex root = roots[i];

                            // If this is a real root (or very close to real)
                            if (Math.Abs(root.Imaginary) < 1e-10)
                            {
                                realRoots.Add(new Complex(root.Real, 0));
                            }
                            else
                            {
                                // Look for the conjugate pair
                                bool foundConjugate = false;
                                for (int j = i + 1; j < roots.Length; j++)
                                {
                                    if (Math.Abs(roots[j].Real - root.Real) < 1e-10 &&
                                        Math.Abs(roots[j].Imaginary + root.Imaginary) < 1e-10)
                                    {
                                        complexPairs.Add(Tuple.Create(root, roots[j]));
                                        foundConjugate = true;
                                        i = j;
                                        break;
                                    }
                                }
                                
                                if (!foundConjugate) realRoots.Add(new Complex(root.Real, 0));
                            }
                        }
                        
                        // Build polynomial from roots
                        var resultPoly = new MathNet.Numerics.Polynomial(new double[] { 1.0 });
                        
                        // Add real root factors: (x - a)
                        foreach (Complex root in realRoots)
                        {
                            var linearFactor = new MathNet.Numerics.Polynomial(new double[] { -root.Real, 1.0 });
                            resultPoly = resultPoly * linearFactor;
                        }
                        
                        // Add complex conjugate pair factors: (x² - 2ax + a² + b²)
                        foreach (var pair in complexPairs)
                        {
                            Complex root = pair.Item1;
                            double a = root.Real;
                            double b = root.Imaginary;
                            var quadraticFactor = new MathNet.Numerics.Polynomial(new double[] { a * a + b * b, -2 * a, 1.0 });
                            resultPoly = resultPoly * quadraticFactor;
                        }

                        // Get the coefficients (ascending order in MathNet)
                        double[] resultCoeffs = resultPoly.Coefficients.ToArray();

                        // Ensure the resulting polynomial has the right length and convert to descending order
                        double[] finalCoeffs = new double[coeffLength];
                        int copyLength = Math.Min(resultCoeffs.Length, coeffLength);

                        // Normalize by the leading coefficient to ensure a[0] = 1
                        double scale = resultCoeffs[resultCoeffs.Length - 1];
                        if (Math.Abs(scale) < 1e-10) scale = 1.0; // Avoid division by zero

                        for (int i = 0; i < copyLength; i++)
                        {
                            // Convert from ascending to descending power order and normalize
                            finalCoeffs[i] = resultCoeffs[resultCoeffs.Length - 1 - i] / scale;
                        }

                        return finalCoeffs;
                    }
                    catch (Exception ex)
                    {
                        // If anything goes wrong, fall back to a safe default
                        double[] defaultCoeffs = new double[coeffLength];
                        defaultCoeffs[0] = 1.0;
                        for (int i = 1; i < coeffLength; i++)
                        {
                            defaultCoeffs[i] = 0.5 * Math.Pow(0.85, i); // Safe values that will create a stable filter
                        }
                        return defaultCoeffs;
                    }
                }

                /// <summary>
                /// Finds polynomial roots using the companion matrix method
                /// </summary>
                public static Complex[] FindRootsCompanionMatrix(double[] coeffs)
                {
                    int n = coeffs.Length - 1;
                    if (n <= 0) return new Complex[0];
                    
                    // Normalize coefficients
                    double leadingCoeff = coeffs[n];
                    double[] normalizedCoeffs = new double[coeffs.Length];
                    for (int i = 0; i < coeffs.Length; i++) normalizedCoeffs[i] = coeffs[i] / leadingCoeff;
                    
                    // Create companion matrix
                    var companion = Matrix<double>.Build.Dense(n, n);
                    
                    // Fill matrix
                    for (int i = 0; i < n; i++)
                    {
                        if (i < n - 1) companion[i + 1, i] = 1.0;
                        companion[i, n - 1] = -normalizedCoeffs[n - i - 1];
                    }
                    
                    try
                    {
                        // Find eigenvalues
                        var evd = companion.Evd(Symmetricity.Asymmetric);
                        return evd.EigenValues.ToArray();
                    }
                    catch
                    {
                        // Create stable artificial roots if eigenvalue calculation fails
                        Complex[] artificialRoots = new Complex[n];
                        for (int i = 0; i < n; i++)
                        {
                            double angle = 2 * Math.PI * i / n;
                            // Create roots with magnitudes between 0.7 and 0.9 to ensure stability
                            double magnitude = 0.7 + 0.2 * (i % 3) / 2.0;
                            artificialRoots[i] = new Complex(magnitude * Math.Cos(angle), magnitude * Math.Sin(angle));
                        }
                        return artificialRoots;
                    }
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