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

using Eto.Forms;
using MathNet.Numerics.Statistics;
using Pachyderm_Acoustic.Audio;
using Pachyderm_Acoustic.Pach_Graphics;
using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Diagnostics.Eventing.Reader;
using System.Numerics;

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
                        double f_center = 50 * Math.Pow(2, oct/3.0);
                        Abs_3rd[oct] = Alpha.Interpolate(f_center);
                        Ref_3rd[oct] = 1 - Abs_3rd[oct];
                    }
                }
                else if(Abs.Length == 24)
                {
                    Abs_3rd = ABS;
                    for(int oct = 0; oct < 8; oct++)
                    {
                        int octave = oct * 3;
                        Ref_3rd[octave] = 1 - Abs_3rd[octave];
                        Ref_3rd[octave+1] = 1 - Abs_3rd[octave+1];
                        Ref_3rd[octave+2] = 1 - Abs_3rd[octave+2];
                        Abs[oct] = (Abs_3rd[octave] + Abs_3rd[octave + 1]+ Abs_3rd[octave + 2])/3;
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
                double[] pr = new double[8] { Math.Sqrt(Ref[0]), Math.Sqrt(Ref[1]), Math.Sqrt(Ref[2]), Math.Sqrt(Ref[3]), Math.Sqrt(Ref[4]), Math.Sqrt(Ref[5]), Math.Sqrt(Ref[6]), Math.Sqrt(Ref[7])};

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

            public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            {
                // Number of samples for the IIR filter design
                int samplect = 4096;
                frequencies = new double[(int)(samplect * max_freq / sample_frequency)];

                // Get reflection spectrum and calculate desired absorption coefficients
                Complex[] desired_Spectrum = new Complex[(int)(samplect * max_freq / sample_frequency)];
                Array.Copy(this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0), desired_Spectrum, desired_Spectrum.Length);

                for (int i = 0; i < frequencies.Length; i++)
                {
                    frequencies[i] = (double)((i + 1) * (sample_frequency / 2)) / samplect;
                }

                double[] desired_alpha = AbsorptionModels.Operations.Absorption_Coef(desired_Spectrum);

                // If filter_order is zero or negative, automatically determine the optimal order
                if (filter_order <= 0)
                {
                    if (rec_order == 0)
                    {
                        int min_order = 2;
                        int max_order = 12;
                        double error_threshold = 0.01; // 1% average error threshold
                        double best_error = double.MaxValue;
                        (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);

                        for (int order = min_order; order <= max_order; order += 2) // Try even orders
                        {
                            try
                            {
                                // Design filter with current order
                                (double[] b, double[] a) current_filter = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
                                    frequencies, desired_Spectrum, order, order, sample_frequency);

                                // Normalize the filter
                                double mag = 0;
                                Complex[] filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
                                    new List<double>(current_filter.b), new List<double>(current_filter.a),
                                    frequencies, sample_frequency);

                                mag = filter_response.MaximumMagnitudePhase().Magnitude;
                                if (mag > 1)
                                {
                                    for (int j = 0; j < current_filter.b.Length; j++)
                                        current_filter.b[j] /= mag;
                                }
                                else
                                {
                                    double mod = mag / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
                                    for (int j = 0; j < current_filter.b.Length; j++)
                                        current_filter.b[j] /= mod;
                                }

                                // Recalculate response after normalization
                                filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
                                    new List<double>(current_filter.b), new List<double>(current_filter.a),
                                    frequencies, sample_frequency);

                                // Calculate error metric (average magnitude error)
                                double error_sum = 0;
                                int count = 0;

                                // Focus on the 125 Hz to 4 kHz range (most important for acoustics)
                                for (int i = 0; i < frequencies.Length; i++)
                                {
                                    if (frequencies[i] >= 125 && frequencies[i] <= 4000)
                                    {
                                        double mag_error = Math.Abs(filter_response[i].Magnitude - desired_Spectrum[i].Magnitude);
                                        error_sum += mag_error / Math.Max(0.01, desired_Spectrum[i].Magnitude); // Relative error
                                        count++;
                                    }
                                }

                                double avg_error = count > 0 ? error_sum / count : double.MaxValue;

                                // Update best filter if this one is better
                                if (avg_error < best_error)
                                {
                                    best_error = avg_error;
                                    best_filter = current_filter;
                                }

                                // If error is below threshold, we've found a good enough match
                                if (avg_error < error_threshold)
                                    break;
                            }
                            catch
                            {
                                // If filter design fails, continue to next order
                                continue;
                            }
                        }

                        // Use the best filter found
                        rec_order = best_filter.a.Length - 1;
                        rec_a = best_filter.a;
                        rec_b = best_filter.b;
                        return best_filter;
                    }
                    else 
                    {
                        return (rec_a, rec_b); // Return previously calculated coefficients
                    }
                }
                else
                {
                    // Original code for specified filter order
                    int numeratorOrder = filter_order;
                    int denominatorOrder = filter_order;

                    try
                    {
                        // Use the IIR filter design method from Pach_SP.IIR_Design
                        (double[] b, double[] a) = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(
                            frequencies, desired_Spectrum, numeratorOrder, denominatorOrder, sample_frequency);

                        // Calculate frequency response of the fitted filter
                        Complex[] IIR_spec = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
                            new List<double>(b), new List<double>(a), frequencies, sample_frequency);

                        // Normalize the filter if necessary
                        double m = IIR_spec.MaximumMagnitudePhase().Magnitude;
                        if (m > 1)
                        {
                            for (int j = 0; j < b.Length; j++) b[j] /= m;
                        }
                        else
                        {
                            double mod = m / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
                            for (int j = 0; j < b.Length; j++) b[j] /= mod;
                        }

                        return (a, b);
                    }
                    catch (Exception)
                    {
                        // If the fitting process fails, fall back to a simple filter
                        double[] a = new double[filter_order + 1];
                        double[] b = new double[filter_order + 1];

                        // Create a basic filter based on octave band coefficients
                        a[0] = 1.0;
                        for (int i = 1; i < a.Length; i++) a[i] = -0.7 * Math.Pow(0.8, i - 1);
                        b[0] = 1.0;


                        return (a, b);
                    }
                }
            }
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
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(true, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient);
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
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(true, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient);
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

                System.Numerics.Complex [][] Zr_interp = Zr.Interpolate(frequency);

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
                for(int i = 0; i < Reflection_Coefficient.Length / 2; i++)
                {
                    List<double> real = new List<double>(), imag = new List<double>();
                    for(int j = 0; j < Reflection_Coefficient[i].Length; j++)
                    {
                        real.Add(Reflection_Coefficient[i][j].Real);
                        imag.Add(Reflection_Coefficient[i][j].Imaginary);
                    }
                    Transfer_FunctionR[Angles.Length/2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, real);
                    Transfer_FunctionI[Angles.Length/2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, imag);
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
                    
                    for(int i = 0; i < frequency.Length; i++)
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
                            Ang_tau_Oct[oct][a] += Trans_Coefficient[35-a][f].Real;
                        }
                    } while (frequency[f] < f_upper);

                    for (int a = 0; a < Angles.Length; a++) Ang_Coef_Oct[oct][a] /= fct[a];
                    RI_Coef[oct] /= RI_count;
                    TI_Coef[oct] /= RI_count;
                }
            }

            public Smart_Material(bool Trans, List<AbsorptionModels.ABS_Layer> Layers, int Samplefreq, double Air_Density, double SoundSpeed)
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
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(true, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient);
                    Trans_Coefficient = new System.Numerics.Complex[Trans_Loss.Length][];
                    for (int i = 0; i < Trans_Loss.Length; i++)
                    {
                        Trans_Coefficient[i] = new System.Numerics.Complex[Trans_Loss[i].Length];
                        for (int j = 0; j < Trans_Coefficient[i].Length; j++)
                        {
                            double t = Trans_Loss[i][j].Magnitude;
                            Trans_Coefficient[i][j] = t*t;
                        }
                    }
                }
                else
                {
                    //Z = AbsorptionModels.Operations.Transfer_Matrix_Explicit_Z(false, 44100, 343, Layers, ref frequency, ref Angles);
                    Z = AbsorptionModels.Operations.Transfer_Matrix_Divisible(false, false, 44100, 343, Layers, ref frequency, ref Angles, out Trans_Loss, out Reflection_Coefficient);
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
                    Transfer_FunctionR[Angles.Length/2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, real);
                    Transfer_FunctionI[Angles.Length/2 - i - 1] = MathNet.Numerics.Interpolation.CubicSpline.InterpolateAkima(frequency, imag);
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

                for(int oct = 0; oct < 8; oct++) Ray.Energy[oct] *= (1 - Ang_Coef_Oct[oct][index]);
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
                for(int a = 0; a < Transfer_FunctionR.Length; a++) alpha += new System.Numerics.Complex(Transfer_FunctionR[a].Interpolate(frequency), Transfer_FunctionI[a].Interpolate(frequency)); 
                alpha /= Transfer_FunctionR.Length;
                return alpha;
            }

            public override System.Numerics.Complex Reflection_Narrow(double frequency, Hare.Geometry.Vector Dir, Hare.Geometry.Vector Normal)
            {
                int a = (int)(Math.Abs(Hare.Geometry.Hare_math.Dot(Dir, Normal))*180/Math.PI / 18);
                return new System.Numerics.Complex(Transfer_FunctionR[a].Interpolate(frequency), Transfer_FunctionI[a].Interpolate(frequency));
            }

            public  class Finite_Field_Impedance
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
                    for (int i = 1; i < anglesdeg.Length; i++) anglesdeg[i] = anglesdeg[i-1] + 5;

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

            public System.Numerics.Complex Admittance (double frequency)
            {
                System.Numerics.Complex R = new System.Numerics.Complex (Transfer_FunctionR[17].Interpolate(frequency), Transfer_FunctionI[17].Interpolate(frequency));
                return (1 - R) / (rho * c * (1 + R));
            }

            public override System.Numerics.Complex[] Reflection_Spectrum(int sample_frequency, int length, Hare.Geometry.Vector Normal, Hare.Geometry.Vector Dir, int threadid)
            {
                int a = (int)(Math.Abs(Hare.Geometry.Hare_math.Dot(Dir, Normal))*180/Math.PI / 18);

                System.Numerics.Complex[] Ref_trns = new System.Numerics.Complex[length];

                for (int j = 0; j < length; j++)
                {
                    double freq = j * (sample_frequency / 2) / length;
                    Ref_trns[j] = new System.Numerics.Complex(Transfer_FunctionR[a].Interpolate(freq), Transfer_FunctionI[a].Interpolate(freq));
                }

                return Ref_trns;
            }

            int rec_order = 0;
            double[] rec_a, rec_b;

            public override (double[] a, double[] b) Estimate_IIR_Coefficients(double sample_frequency, double max_freq, out double[] frequencies, int filter_order = 0)
            {
                int samplect = 4096; // Number of samples for the IIR filter design
                frequencies = new double[(int)(samplect * max_freq / sample_frequency)];

                // Get reflection spectrum and calculate desired absorption coefficients
                Complex[] desired_Spectrum = new Complex[(int)(samplect * max_freq / sample_frequency)];
                Array.Copy(this.Reflection_Spectrum((int)sample_frequency, samplect, new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1), 0), desired_Spectrum, desired_Spectrum.Length);

                for (int i = 0; i < desired_Spectrum.Length; i++)
                {
                    frequencies[i] = (double)((i + 1) * (sample_frequency / 2)) / 4096;
                }

                double[] desired_alpha = AbsorptionModels.Operations.Absorption_Coef(desired_Spectrum);

                // If filter_order is zero or negative, automatically determine the optimal order
                if (filter_order <= 0)
                {
                    if (rec_order == 0)
                    {
                        // Use a range of filter orders to find the best fit

                        int min_order = 2;
                        int max_order = 12;
                        double error_threshold = 0.01; // 1% average error threshold
                        double best_error = double.MaxValue;
                        (double[] b, double[] a) best_filter = (new double[min_order + 1], new double[min_order + 1]);

                        for (int order = min_order; order <= max_order; order += 2) // Try even orders
                        {
                            try
                            {
                                // Design filter with current order
                                (double[] b, double[] a) current_filter = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(frequencies, desired_Spectrum, order, order);

                                // Normalize the filter
                                Complex[] filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
                                    new List<double>(current_filter.b), new List<double>(current_filter.a),
                                    frequencies, sample_frequency);

                                double mag = filter_response.MaximumMagnitudePhase().Magnitude;
                                if (mag > 1)
                                {
                                    for (int j = 0; j < current_filter.b.Length; j++)
                                        current_filter.b[j] /= mag;
                                }
                                else
                                {
                                    double mod = mag / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
                                    for (int j = 0; j < current_filter.b.Length; j++)
                                        current_filter.b[j] /= mod;
                                }

                                // Recalculate response after normalization
                                filter_response = Audio.Pach_SP.IIR_Design.AB_FreqResponse(
                                    new List<double>(current_filter.b), new List<double>(current_filter.a),
                                    frequencies, sample_frequency);

                                // Calculate error metric focusing on important frequency range
                                double error_sum = 0;
                                int count = 0;

                                for (int i = 0; i < frequencies.Length; i++)
                                {
                                    if (frequencies[i] >= 125 && frequencies[i] <= 4000)
                                    {
                                        double mag_error = Math.Abs(filter_response[i].Magnitude - desired_Spectrum[i].Magnitude);
                                        error_sum += mag_error / Math.Max(0.01, desired_Spectrum[i].Magnitude); // Relative error
                                        count++;
                                    }
                                }

                                double avg_error = count > 0 ? error_sum / count : double.MaxValue;

                                // Update best filter if this one is better
                                if (avg_error < best_error)
                                {
                                    best_error = avg_error;
                                    best_filter = current_filter;
                                }

                                // If error is below threshold, we've found a good enough match
                                if (avg_error < error_threshold)
                                {
                                    // We found a filter with acceptable error, no need to try higher orders
                                    break;
                                }
                            }
                            catch
                            {
                                // If filter design fails, continue to next order
                                continue;
                            }
                        }

                        // Use the best filter found
                        return (best_filter.a, best_filter.b);
                    }
                    else
                    {
                        // Use previously calculated coefficients
                        return (rec_a, rec_b);
                    }
                }
                else
                {
                    // Original code for specified filter order
                    int numeratorOrder = filter_order;
                    int denominatorOrder = filter_order;

                    // Use Modified Yule-Walker method to fit the IIR filter
                    (double[] b, double[] a) = Pach_SP.IIR_Design.FitIIRToFrequencyResponse(frequencies, desired_Spectrum, numeratorOrder, denominatorOrder);

                    // Calculate frequency response of the fitted filter
                    Complex[] IIR_spec = Audio.Pach_SP.IIR_Design.AB_FreqResponse(new List<double>(b), new List<double>(a), frequencies, 44100);

                    // Normalize the filter if necessary
                    double m = IIR_spec.MaximumMagnitudePhase().Magnitude;
                    if (m > 1)
                    {
                        for (int j = 0; j < b.Length; j++) b[j] /= m;
                    }
                    else
                    {
                        double mod = m / desired_Spectrum.MaximumMagnitudePhase().Magnitude;
                        for (int j = 0; j < b.Length; j++) b[j] /= mod;
                    }

                    return (a, b);
                }
            }
        }

        public class Lambert_Scatteringex : Scattering
        {
            double[,] Scattering_Coefficient;
            public Lambert_Scatteringex(double[] Scattering, double SplitRatio)
            {
                Scattering_Coefficient = new double[8, 3];
                for (int oct = 0; oct < 8; oct++)
                {
                    double Mod = Math.Abs(Scattering[oct] - 0.5); //((Scattering[oct] < (1 - Scattering[oct])) ? (Scattering[oct] * SplitRatio / 2) : ((1 - Scattering[oct]) * SplitRatio / 2));
                    Scattering_Coefficient[oct, 1] = Scattering[oct];
                    Scattering_Coefficient[oct, 0] = Scattering_Coefficient[oct, 1] - Mod;
                    Scattering_Coefficient[oct, 2] = Scattering_Coefficient[oct, 1] + Mod;
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

                foreach (int oct in Ray.Freq_Bands)
                {
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

                    if (0.99 < proj && 1.01 > proj)// diffx = new Hare.Geometry.Vector(1, 0, 0);
                        diffy = Hare.Geometry.Hare_math.Cross(diffz.dx, diffz.dy, diffz.dz, 1, 0, 0);
                    else
                        diffy = Hare.Geometry.Hare_math.Cross(diffz.dx, diffz.dy, diffz.dz, 0, 0, 1);

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

                    if (Transmission != null && Transmission[oct] > 0)
                    {
                        OctaveRay tr = Ray.SplitRay(oct, Transmission[oct]);
                        Rays.Enqueue(tr);

                        OctaveRay td = R.SplitRay(Transmission[oct]);
                        td.Reverse();
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
                    //Hare.Geometry.Vector diffx;
                    //Hare.Geometry.Vector diffy;
                    //Hare.Geometry.Vector diffz;
                    double proj;
                    //Check that the ray and the normal are both on the same side...
                    if (Cos_Theta > 0) Normal *= -1;
                    //diffz = Normal;
                    double diffz_x = Normal.dx;
                    double diffz_y = Normal.dy;
                    double diffz_z = Normal.dz;
                    //diffx = new Hare.Geometry.Vector(0, 0, 1);
                    double diffx_x = 0, diffx_y = 0, diffx_z = 0;
                    proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz_x, diffz_y, diffz_z, 0, 0, 1));

                    if (0.99 < proj && 1.01 > proj)
                    { //diffx = new Hare.Geometry.Vector(1, 0, 0);
                        diffx_x = 1;
                    }
                    else
                    {
                        diffx_z = 1;
                    }
                    //diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
                    double diffy_x = diffz_y * diffx_z - diffz_z * diffx_y;
                    double diffy_y = -(diffz_x * diffx_z - diffz_z * diffx_x);
                    double diffy_z = diffz_x * diffx_y - diffz_y * diffx_x;

                    //diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
                    diffx_x = diffy_y * diffz_z - diffy_z * diffz_y;
                    diffx_y = -(diffy_x * diffz_z - diffy_z * diffz_x);
                    diffx_z = diffy_x * diffz_y - diffy_y * diffz_x;

                    //diffx.Normalize();
                    //diffy.Normalize();
                    //diffz.Normalize();
                    Hare.Geometry.Hare_math.Normalize(ref diffx_x, ref diffx_y, ref diffx_z);
                    Hare.Geometry.Hare_math.Normalize(ref diffy_x, ref diffy_y, ref diffy_z);
                    Hare.Geometry.Hare_math.Normalize(ref diffz_x, ref diffz_y, ref diffz_z);

                    double u1;
                    double u2;
                    double x;
                    double y;
                    double z;
                    //Hare.Geometry.Vector vect;
                    double vectx, vecty, vectz;
                    u1 = 2.0 * Math.PI * rand.NextDouble();
                    // random azimuth
                    double Scat_Mod = rand.NextDouble();
                    u2 = Math.Acos(Scat_Mod);
                    // random zenith (elevation)
                    x = Math.Cos(u1) * Math.Sin(u2);
                    y = Math.Sin(u1) * Math.Sin(u2);
                    z = Scat_Mod;

                    //vect = (diffx * x) + (diffy * y) + (diffz * z);
                    vectx = diffx_x * x + diffx_y * x - diffx_z * x;
                    vecty = diffy_x * y + diffy_y * y - diffy_z * y;
                    vectz = diffz_x * z + diffz_y * z - diffz_z * z;
                    Hare.Geometry.Hare_math.Normalize(ref vectx, ref vecty, ref vectz);

                    //Return the new direction
                    Ray.dx = vectx;
                    Ray.dy = vecty;
                    Ray.dz = vectz;
                }
                else
                {
                    //Specular Reflection
                    Ray.dx -= Normal.dx * Cos_Theta * 2;
                    Ray.dy -= Normal.dy * Cos_Theta * 2;
                    Ray.dz -= Normal.dz * Cos_Theta * 2;
                }

                if (Transmit) { Ray.Reverse(); }
            }

            public override void Scatter_Late(ref OctaveRay Ray, ref Queue<OctaveRay> Rays, ref Random rand, Hare.Geometry.Vector Normal, double Cos_Theta, bool Transmit = false)
            {
                //Approximation of a logarithmic function to split proportionally with logarithmic decay, but faster.
                double percent =  Ray.Intensity * Ray.Decimation_threshold;
                percent *= 1E9;
                //if (percent < 0.0001) percent = percent * 1250;
                //else if (percent < 0.001) percent = percent * 125 + 0.1;
                //else if (percent < 0.01) percent = percent * 12.5 + 0.2;
                //else if (percent < 0.1) percent = percent * 1.25 + 0.3;
                //else if (percent < 1) percent = percent * .125 + 0.4;

                if (percent < 0.00001) percent = percent * 25000;
                else if (percent < 0.0001) percent = percent * 2500 + 0.1;
                else if (percent < 0.001) percent = percent * 250 + 0.2;
                else if (percent < 0.01) percent = percent * 25 + 0.3;
                else if (percent < .1) percent = percent * 2.5 + 0.4;

                double scat_sel = rand.NextDouble();
                if (scat_sel > percent + 0.5)
                {
                    // Specular Reflection
                    Ray.dx -= Normal.dx * Cos_Theta * 2;
                    Ray.dy -= Normal.dy * Cos_Theta * 2;
                    Ray.dz -= Normal.dz * Cos_Theta * 2;
                    if (Transmit) Ray.Reverse();
                    return;
                }
                else if (scat_sel >= 0.5 - percent)
                {
                    //Only for a certain portion of high benefit cases--
                    //// a. Create new source for scattered energy (E * Scattering).
                    //// b. Modify E (E * 1 - Scattering).
                    //Create a new ray...
                    OctaveRay tr = Ray.SplitRay(1 - Scattering_Coefficient[Ray.Octave, 1]);
                    // this is the specular reflection. Save it for later.
                    tr.dx -= Normal.dx * Cos_Theta * 2;
                    tr.dy -= Normal.dy * Cos_Theta * 2;
                    tr.dz -= Normal.dz * Cos_Theta * 2;
                    if (Transmit) tr.Reverse();

                    Rays.Enqueue(tr);
                }

                //If we are here, the Original ray needs a scattered direction:
                //Hare.Geometry.Vector diffx;
                //Hare.Geometry.Vector diffy;
                //Hare.Geometry.Vector diffz;
                double proj;
                //Check that the ray and the normal are both on the same side...
                if (Cos_Theta > 0) Normal *= -1;
                //diffz = Normal;
                double diffz_x = Normal.dx;
                double diffz_y = Normal.dy;
                double diffz_z = Normal.dz;
                //diffx = new Hare.Geometry.Vector(0, 0, 1);
                double diffx_x= 0, diffx_y = 0, diffx_z = 0;
                proj = Math.Abs(Hare.Geometry.Hare_math.Dot(diffz_x, diffz_y, diffz_z, 0, 0, 1));

                if (0.99 < proj && 1.01 > proj) 
                { //diffx = new Hare.Geometry.Vector(1, 0, 0);
                    diffx_x = 1;
                }
                else
                {
                    diffx_z = 1;
                }
                //diffy = Hare.Geometry.Hare_math.Cross(diffz, diffx);
                double diffy_x = diffz_y * diffx_z - diffz_z * diffx_y;
                double diffy_y = -(diffz_x * diffx_z - diffz_z * diffx_x);
                double diffy_z = diffz_x * diffx_y - diffz_y * diffx_x;

                //diffx = Hare.Geometry.Hare_math.Cross(diffy, diffz);
                diffx_x = diffy_y * diffz_z - diffy_z * diffz_y;
                diffx_y = -(diffy_x * diffz_z - diffy_z * diffz_x);
                diffx_z = diffy_x * diffz_y - diffy_y * diffz_x;

                //diffx.Normalize();
                //diffy.Normalize();
                //diffz.Normalize();
                Hare.Geometry.Hare_math.Normalize(ref diffx_x, ref diffx_y, ref diffx_z);
                Hare.Geometry.Hare_math.Normalize(ref diffy_x, ref diffy_y, ref diffy_z);
                Hare.Geometry.Hare_math.Normalize(ref diffz_x, ref diffz_y, ref diffz_z);

                double u1;
                double u2;
                double x;
                double y;
                double z;
                //Hare.Geometry.Vector vect;
                double vectx, vecty, vectz;
                u1 = 2.0 * Math.PI * rand.NextDouble();
                // random azimuth
                double Scat_Mod = rand.NextDouble();
                u2 = Math.Acos(Scat_Mod);
                // random zenith (elevation)
                x = Math.Cos(u1) * Math.Sin(u2);
                y = Math.Sin(u1) * Math.Sin(u2);
                z = Scat_Mod;//Math.Cos(u2);

                //vect = (diffx * x) + (diffy * y) + (diffz * z);
                vectx = diffx_x * x + diffx_y * x - diffx_z * x;
                vecty = diffy_x * y + diffy_y * y - diffy_z * y;
                vectz = diffz_x * z + diffz_y * z - diffz_z * z;
                Hare.Geometry.Hare_math.Normalize(ref vectx, ref vecty, ref vectz);

                //Return the new direction
                Ray.dx = vectx;
                Ray.dy = vecty;
                Ray.dz = vectz;
                if (Transmit) Ray.Reverse();
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
                    double Mod = Math.Abs(Scattering[oct] - 0.5); //((Scattering[oct] < (1 - Scattering[oct])) ? (Scattering[oct] * SplitRatio / 2) : ((1 - Scattering[oct]) * SplitRatio / 2));
                    Scattering_Coefficient[oct, 1] = Scattering[oct];
                    Scattering_Coefficient[oct, 0] = Scattering_Coefficient[oct, 1] - Mod;
                    Scattering_Coefficient[oct, 2] = Scattering_Coefficient[oct, 1] + Mod;
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
                    
                for(int oct = 0; oct < 8; oct++)
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
                    OctaveRay tr = Ray.SplitRay(1 - Scattering_Coefficient[Ray.Octave,1]);
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