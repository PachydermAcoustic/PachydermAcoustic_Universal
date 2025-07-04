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

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using System;
using System.ComponentModel;
using System.Numerics;

namespace Pachyderm_Acoustic
{
    namespace AbsorptionModels
    {
        public static class Explicit_TMM
        {
            public static SparseMatrix PorousMatrix(bool Rigid, double d, Complex kt, double freq, double porosity, double tortuosity, double YoungsModulus, double PoissonRatio, double flow_resistivity, double FrameDensity, double AmbientMeanPressure, double Viscous_Characteristic_Length = 0 , double Viscous_Permeability = 0, double Thermal_Characteristic_Length = 0, double Thermal_Permeability = 0)
            {
                //Viscous_Characteristic_Length = 0.56E-4;
                //Thermal_Characteristic_Length = 1.1E-4;
                //Get the Viscous Parameters...
                if (Viscous_Characteristic_Length == 0) Viscous_Characteristic_Length = System.Math.Sqrt(8 * AbsorptionModels.Biot_Porous_Absorbers.Shear_Viscosity * tortuosity / (flow_resistivity * porosity));
                if (Viscous_Permeability == 0) Viscous_Permeability = Biot_Porous_Absorbers.Viscous_Permeability(flow_resistivity);
                ////Get the Thermal Parameters...
                if (Thermal_Characteristic_Length == 0) Thermal_Characteristic_Length = Biot_Porous_Absorbers.Thermal_Characteristic_Length(Viscous_Characteristic_Length);
                if (Thermal_Permeability == 0) Thermal_Permeability = Viscous_Permeability;//;porosity * Thermal_Characteristic_Length * Thermal_Characteristic_Length / 8;
                

                double w = Utilities.Numerics.PiX2 * freq;
                Complex FrameShear = AbsorptionModels.Biot_Porous_Absorbers.Shear_Modulus(YoungsModulus, PoissonRatio) * new Complex(1,.1);
                Complex Gw = Biot_Porous_Absorbers.G_w(tortuosity, porosity, Viscous_Permeability, Viscous_Characteristic_Length, w);
                //Complex Gw = Biot_Porous_Absorbers.Gp_w(1, tortuosity, Viscous_Permeability, porosity, Viscous_Characteristic_Length, w);
                Complex Gpw = Biot_Porous_Absorbers.G_w_prime(porosity, Thermal_Permeability, Thermal_Characteristic_Length, w);
                Complex kb = 2 * FrameShear * (PoissonRatio + 1) / (3 * (1 - 2 * PoissonRatio));
                Complex Kf = Biot_Porous_Absorbers.BulkMod_Fluid(w, AmbientMeanPressure, porosity,Thermal_Permeability, Gpw);//AmbientMeanPressure / (1 - (gamma - 1) / (gamma * alpha));
                double rho12 = Biot_Porous_Absorbers.rho12(porosity, tortuosity);
                double rhoa = Biot_Porous_Absorbers.rhoA(rho12);
                
                Complex P, Q, R;
                if (!Rigid)
                {
                    //Universal (Limp) Frame Case:
                    double Ks = AbsorptionModels.Biot_Porous_Absorbers.BulkMod_Solid(YoungsModulus, PoissonRatio);
                    P = ((1 - porosity) * (1 - kb / Ks) * Ks + porosity * Ks * kb / Kf) / (1 - porosity - kb / Ks + porosity * Ks / Kf);
                    Q = (1 - porosity - kb / Ks) * porosity * Ks / (1 - porosity - kb / Ks + porosity * Ks / Kf);
                    R = porosity * porosity * Ks / (1 - porosity - kb / Ks + porosity * Ks / Kf);
                }
                else
                {
                    //Rigid Frame Case:
                    P = 4 * FrameShear / 3 + kb + ((1 - porosity) * (1 - porosity)) * Kf / porosity;
                    R = porosity * Kf;
                    Q = Kf * (1 - porosity);
                }

                Complex rho12eff = Biot_Porous_Absorbers.rho12eff(rhoa, porosity, flow_resistivity, Gw, w);
                Complex rho22eff = Biot_Porous_Absorbers.rho22eff(rhoa, porosity, flow_resistivity, Gw, w);
                Complex rho11eff = Biot_Porous_Absorbers.rho11eff(FrameDensity, rhoa, porosity, flow_resistivity, Gw, w);

                Complex test = w * Complex.Sqrt(rho22eff / R);

                Complex D = (P * rho22eff + R * rho11eff - 2 * Q * rho12eff);
                Complex rtDELTA = Complex.Sqrt(D * D - 4 * (P * R - Q * Q) * (rho11eff * rho22eff - rho12eff * rho12eff));
                Complex delta = w * w / (2 * (P * R - Q * Q));
                Complex drho = P * rho22eff + R * rho11eff - Q * rho12eff;
                Complex delta21 = delta * (drho - rtDELTA);
                Complex delta22 = delta * (drho + rtDELTA);
                Complex delta23 = (w * w / FrameShear) * ((rho11eff * rho22eff - rho12eff * rho12eff) / rho22eff);

                //Complex kt = k * sin_theta;
                Complex k13 = Complex.Sqrt(delta21 - kt * kt);
                Complex k23 = Complex.Sqrt(delta22 - kt * kt);
                Complex k33 = Complex.Sqrt(delta23 - kt * kt);
                Complex Mu1 = P * delta21 - w * w * rho11eff / (w * w * rho12eff - Q * delta21);
                Complex Mu2 = P * delta22 - w * w * rho11eff / (w * w * rho12eff - Q * delta22);
                Complex Mu3 = FrameShear * delta23 - w * w * rho11eff / (w * w * rho22eff);

                SparseMatrix GH = GammaH_P(kt, w, d, FrameShear, P, Q, R, k13, k23, k33, Mu1, Mu2, Mu3);
                SparseMatrix G0T = GammaH_P(kt, w, 0, FrameShear, P, Q, R, k13, k23, k33, Mu1, Mu2, Mu3).Inverse() as SparseMatrix;
                //SparseMatrix G0T = Gamma0T_P(kt, w, FrameShear, P, Q, R, k13, k23, k33, Mu1, Mu2, Mu3);

                return GH * G0T;
            }

            public static SparseMatrix GammaH_P(Complex kt, double w, double h, Complex ShearModulus, Complex P, Complex Q, Complex R, Complex k13, Complex k23, Complex k33, Complex Mu1, Complex Mu2, Complex Mu3)
            {
                SparseMatrix M = new SparseMatrix(6, 6);
                h *= -1;

                //double w = Utilities.Numerics.PiX2 * freq;
                //Complex LameL = FrameElasticity * PoissonRatio / ((1 + PoissonRatio) * (1 - 2 * PoissonRatio));
                //Complex LameMu = FrameElasticity / (2 * (1 + PoissonRatio));
                //Complex delta21 = w * w * FrameDensity;
                //Complex delta22 = w * w * FrameDensity;
                //Complex delta23 = delta21 / LameMu;
                // delta21 /= (LameL + 2 * LameMu);
                //delta22 /= (LameL + 2 * LameMu);
                ////Complex k1 = k * Complex.Cos(theta);
                ////Complex k3 = k * Complex.Sin(theta);
                //Complex k13 = Complex.Sqrt(delta21 - kt * kt);
                //Complex k23 = Complex.Sqrt(delta22 - kt * kt);
                //Complex k33 = Complex.Sqrt(delta23 - kt * kt);
                ////Complex Damping1 = LameMu * (k13 * k13 - k0 * k0);
                ////Complex Damping2 = 2 * LameMu * k0;

                //double rho12 = Biot_Porous_Absorbers.rho12(porosity, tortuosity);
                //double rhoa = Biot_Porous_Absorbers.rhoA(rho12);
                //double Viscous_Permeability = Biot_Porous_Absorbers.Viscous_Permeability(flow_resistivity);
                //double v = Biot_Porous_Absorbers.v();
                //Complex Gw = Biot_Porous_Absorbers.G_w(tortuosity, porosity, Viscous_Permeability, characteristic_length, freq, v);
                //Complex rho12eff = Biot_Porous_Absorbers.rho12eff(rhoa, porosity, flow_resistivity, Gw, freq);
                //Complex rho22eff = Biot_Porous_Absorbers.rho22eff(rhoa, porosity, flow_resistivity, Gw, freq);
                //Complex rho11eff = Biot_Porous_Absorbers.rho11eff(FrameDensity, rhoa, porosity, flow_resistivity, Gw, freq);

                //Complex Kf = AmbientMeanPressure / (1 - (gamma - 1) / (gamma * alpha));
                //Complex P = 4 * FrameShear / 3 + kb + (phi_1 * phi_1) * Kf / porosity;
                //Complex R = porosity * Kf;
                //Complex Q = Kf * (1 - porosity);
                //Complex Mu1 = Q * delta21 - w * w * rho11eff / (w * w *rho22eff - R *delta21);
                //Complex Mu2 = Q * delta22 - w * w * rho11eff / (w * w * rho22eff - R * delta22);
                //Complex Mu3 = FrameShear * delta23 - w * w * rho11eff / (w * w * rho22eff);

                Complex D1 = (P + Q * Mu1) * (kt * kt + k13 * k13) - 2 * ShearModulus * kt * kt;
                Complex D2 = (P + Q * Mu2) * (kt * kt + k23 * k23) - 2 * ShearModulus * kt * kt;
                Complex E1 = (R * Mu1 + Q) * (kt * kt + k13 * k13);
                Complex E2 = (R * Mu2 + Q) * (kt * kt + k23 * k23);

                //Complex Mu1 =                         

                //Taken from Lauriks, et. al 1990.
                M[0, 0] = w * kt * Complex.Cos(k13 * h);
                M[1, 0] = -Complex.ImaginaryOne * w * k13 * Complex.Sin(k13 * h);
                M[2, 0] = -Complex.ImaginaryOne * w * k13 * Mu1 * Complex.Sin(k13 * h);
                M[3, 0] = -D1 * Complex.Cos(k13 * h);
                M[4, 0] = 2 * Complex.ImaginaryOne * ShearModulus * kt * k13 * Complex.Sin(k13 * h);
                M[5, 0] = -E1 * Complex.Cos(k13 * h);
                M[0, 1] = -Complex.ImaginaryOne * w * kt * Complex.Sin(k13 * h);
                M[1, 1] = w * k13 * Complex.Cos(k13 * h);
                M[2, 1] = w * k13 * Mu1 * Complex.Cos(k13 * h);
                M[3, 1] = Complex.ImaginaryOne * D1 * Complex.Sin(k13 * h);
                M[4, 1] = -2 * ShearModulus * kt * k13 * Complex.Cos(k13 * h);
                M[5, 1] = Complex.ImaginaryOne * E1 * Complex.Sin(k13 * h);

                M[0, 2] = w * kt * Complex.Cos(k23 * h);
                M[1, 2] = -Complex.ImaginaryOne * w * k23 * Complex.Sin(k23 * h);
                M[2, 2] = -Complex.ImaginaryOne * w * k23 * Mu2 * Complex.Sin(k23 * h);
                M[3, 2] = -D2 * Complex.Cos(k23 * h);
                M[4, 2] = 2 * Complex.ImaginaryOne * ShearModulus * kt * k23 * Complex.Sin(k23 * h);
                M[5, 2] = -E2 * Complex.Cos(k23 * h);
                M[0, 3] = -Complex.ImaginaryOne * w * kt * Complex.Sin(k23 * h);
                M[1, 3] = w * k23 * Complex.Cos(k23 * h);
                M[2, 3] = w * k23 * Mu2 * Complex.Cos(k23 * h);
                M[3, 3] = Complex.ImaginaryOne * D2 * Complex.Sin(k23 * h);
                M[4, 3] = -2 * ShearModulus * kt * k23 * Complex.Cos(k23 * h);
                M[5, 3] = Complex.ImaginaryOne * E2 * Complex.Sin(k23 * h);

                M[0, 4] = Complex.ImaginaryOne * w * k33 * Complex.Sin(k33 * h);
                M[1, 4] = w * kt * Complex.Cos(k33 * h);
                M[2, 4] = w * kt * Mu3 * Complex.Cos(k33 * h);
                M[3, 4] = 2 * Complex.ImaginaryOne * ShearModulus * k33 * kt * Complex.Sin(k33 * h);
                M[4, 4] = ShearModulus * (k33 * k33 - kt * kt) * Complex.Cos(k33 * h);
                M[5, 4] = 0;
                M[0, 5] = -w * k33 * Complex.Cos(k33 * h);
                M[1, 5] = -Complex.ImaginaryOne * w * kt * Complex.Sin(k33 * h);
                M[2, 5] = -Complex.ImaginaryOne * w * kt * Mu3 * Complex.Sin(k33 * h);
                M[3, 5] = -2 * ShearModulus * k33 * kt * Complex.Cos(k33 * h);
                M[4, 5] = -Complex.ImaginaryOne * ShearModulus * (k33 * k33 - kt * kt) * Complex.Sin(k33 * h);
                M[5, 5] = 0;

                return M;
            }

            public static SparseMatrix Gamma0T_P(Complex kt, double w, Complex ShearModulus, Complex P, Complex Q, Complex R, Complex k13, Complex k23, Complex k33, Complex Mu1, Complex Mu2, Complex Mu3)
            {
                SparseMatrix M = new SparseMatrix(6, 6);

                Complex A1 = -(P - 2 * ShearModulus + Mu1 * Q) * kt * kt - 2 * ShearModulus * k13 * k13;
                Complex A2 = -(P - 2 * ShearModulus + Mu2 * Q) * kt * kt - 2 * ShearModulus * k23 * k23;
                Complex B1 = -(Q + R * Mu1) * kt * kt;
                Complex B2 = -(Q + R * Mu2) * kt * kt;
                Complex X = B1 * A2 - B2 * A1;
                Complex Y = X - 2 * ShearModulus * w * w * kt * kt * (B1 - B2);
                Complex Mu21 = Mu2 - Mu1;
                Complex Mu31 = Mu3 - Mu1;

                //  gamma = kt?

                M[0, 0] = -B2 * 2 * ShearModulus * w * kt / (Complex.ImaginaryOne * Y);
                //M[0,1] = 
                //M[0,2] = 
                M[0, 3] = -B2 * (1 / X + 2 * ShearModulus * w * w * kt * kt * (B1 - B2)) / (X * Y);
                //M[0,4] = 
                M[0, 5] = 1 + B2 * ((A1 / X) - (1 - A1 * (B1 - B2) / X) * (2 * ShearModulus * w * w * kt * kt) / Y) / B1;
                //M[1,0] = 
                M[1, 1] = -((1 + Mu1 / Mu21) - 2 * w * w * kt * kt * (1 - Mu31 / Mu21) / (kt * kt)) / (Complex.ImaginaryOne * k13);
                M[1, 2] = 1 / (k13 * Mu21);
                //M[1,3] = 
                M[1, 4] = -w * kt * (1 - (Mu31 / Mu21)) / (k13 * ShearModulus * kt * kt);
                //M[1,5] = 
                M[2, 0] = (B1 * 2 * ShearModulus * w * kt) / (Complex.ImaginaryOne * Y);
                //M[2,1] = 
                //M[2,2] = 
                M[2, 3] = B1 * (1 / X + (B1 - B2) * (2 * ShearModulus * w * w * kt * kt) / (X * Y));
                //M[2,4] = 
                M[2, 5] = -(A1 / X - 2 * ShearModulus * w * w * kt * kt * (1 - A1 * (B1 - B2) / X) / Y);
                //M[3,0] = 
                M[3, 1] = (Mu1 / Mu21 + 2 * (w * w * kt * kt * Mu31) / (kt * kt * Mu21)) / (Complex.ImaginaryOne * k23); ///Shown in paper as simply alpha... alpha2 = k23 in Atalla's work... 
                M[3, 2] = -1 / (Complex.ImaginaryOne * k23 * Mu21);
                //M[3,3] = 
                M[3, 4] = -(w * w * kt * kt * Mu31) / (k23 * ShearModulus * kt * kt * Mu21);
                //M[3,5] =                                                 M[0,0] = 
                //M[0,4] =
                M[4, 1] = -(2 * w * w * kt * kt) / (Complex.ImaginaryOne * kt * kt);
                //M[4,2] = 
                //M[4,3] = 
                M[4, 4] = 1 / (ShearModulus * kt * kt);
                //M[4,5] =
                M[5, 0] = X / (Complex.ImaginaryOne * k33 * Y);
                //M[5,1] = 
                //M[5,2] = 
                M[5, 3] = w * kt * (B1 - B2) / (k33 * Y);
                //M[5,4] = 
                M[5, 5] = w * kt * (1 - A1 * (B1 - B2) / X) / (k33 * Y * B1);

                return M;
            }

            public static SparseMatrix FluidMatrix(double h, Complex kz, Complex K, Complex Zc)
            {
                SparseMatrix M = new SparseMatrix(2, 2);
                Complex kzh = kz * h;
                Complex sinkzh = Complex.Sin(kzh);
                Complex coskzh = Complex.Cos(kzh);
                //Complex pr2 = w * K * Zc / w;
                Complex pr2 = K * Zc;

                M[0, 0] = 1d * coskzh;
                M[0, 1] = (Complex.ImaginaryOne * pr2 / kz) * sinkzh;
                M[1, 0] = (Complex.ImaginaryOne * kz / pr2) * sinkzh;
                M[1, 1] = 1d * coskzh;

                return M;
            }

            public static SparseMatrix AirMatrix(double h, Complex kz, Complex K, Complex Zc)
            {
                SparseMatrix M = new SparseMatrix(2, 2);
                Complex kzh = kz * h;
                Complex sinkzh = Complex.Sin(kzh);
                Complex coskzh = Complex.Cos(kzh);

                M[0, 0] = 1d * coskzh;
                M[0, 1] = (Complex.ImaginaryOne * Zc / kz) * sinkzh;
                M[1, 0] = (Complex.ImaginaryOne * kz / Zc) * sinkzh;
                M[1, 1] = 1d * coskzh;

                return M;
            }

            public static SparseMatrix Solid_Matrix(Complex kt, double h, double freq, double density, double Youngs_Modulus, double Poisson_Ratio)
            {
                // Reference: Allard & Atalla "Propagation of Sound in Porous Media", 2nd ed. with corrections and augmentations
                h *= -1; // Maintain original sign convention

                // Calculate Lamé parameters from Young's modulus and Poisson's ratio
                double w = Utilities.Numerics.PiX2 * freq;

                (Complex Mu, Complex Lambda) = Solids.ComplexLameParameters(Youngs_Modulus, Poisson_Ratio, density, freq);

                // Calculate wave numbers for longitudinal and transverse waves
                Complex d21 = w * w * density / (Lambda + 2 * Mu);
                Complex d23 = w * w * density / Mu;

                Complex k13, k33;

                if ((d21 - kt * kt).Real < 0 && Math.Abs((d21 - kt * kt).Real) < 1e-10)
                    k13 = Complex.Sqrt(new Complex(0, (d21 - kt * kt).Imaginary));
                else
                    k13 = Complex.Sqrt(d21 - kt * kt);

                if ((d23 - kt * kt).Real < 0 && Math.Abs((d23 - kt * kt).Real) < 1e-10)
                    k33 = Complex.Sqrt(new Complex(0, (d23 - kt * kt).Imaginary));
                else
                    k33 = Complex.Sqrt(d23 - kt * kt);

                // Avoid small imaginary parts that should be zero
                if (Math.Abs(k13.Imaginary) < 1e-10)
                    k13 = new Complex(k13.Real, 0);
                if (Math.Abs(k33.Imaginary) < 1e-10)
                    k33 = new Complex(k33.Real, 0);

                // Ensure non-zero wavenumbers for numerical stability
                if (k13.Magnitude < 1e-10)
                    k13 = new Complex(1e-10, 1e-10);
                if (k33.Magnitude < 1e-10)
                    k33 = new Complex(1e-10, 1e-10);
                // Avoid small imaginary parts that should be zero
                if (Math.Abs(k13.Imaginary) < 1e-10) k13 = new Complex(k13.Real, 0);
                if (Math.Abs(k33.Imaginary) < 1e-10) k33 = new Complex(k33.Real, 0);

                // Trigonometric functions
                Complex cos_k13_h = Complex.Cos(k13 * h);
                Complex sin_k13_h = Complex.Sin(k13 * h);
                Complex cos_k33_h = Complex.Cos(k33 * h);
                Complex sin_k33_h = Complex.Sin(k33 * h);

                // Create matrix D(h) at depth h according to Allard & Atalla
                SparseMatrix Dh = new SparseMatrix(4, 4);
                //double k0 = w / 343.0;
                //Allart and Atalla equation... likely incorrect.
                //Complex D1 = Lambda * (k0 * k0 + k13 * k13) + 2 * Mu * k13 * k13;
                //Complex D2 = 2 * Mu * k0;
                Complex D1 = Lambda * (kt * kt + k13 * k13) + 2 * Mu * k13 * k13;
                Complex D2 = 2 * Mu * kt;

                //// Fill according to equations in Allard & Atalla
                Dh[0, 0] = w * kt * cos_k13_h;                                 // u_x for longitudinal wave, first potential
                Dh[0, 1] = -w * Complex.ImaginaryOne * kt * sin_k13_h;          // u_x for longitudinal wave, second potential
                Dh[0, 2] = w * Complex.ImaginaryOne * k33 * sin_k33_h;        // u_x for transverse wave, first potential
                Dh[0, 3] = -w * k33 * cos_k33_h;                              // u_x for transverse wave, second potential

                Dh[1, 0] = -w * Complex.ImaginaryOne * k13 * sin_k13_h;                               // u_z for longitudinal wave, first potential
                Dh[1, 1] = w * k13 * cos_k13_h;       // u_z for longitudinal wave, second potential
                Dh[1, 2] = w * kt * cos_k33_h;                                 // u_z for transverse wave, first potential
                Dh[1, 3] = -w * Complex.ImaginaryOne * kt * sin_k33_h;          // u_z for transverse wave, second potential

                Dh[2, 0] = -D1 * cos_k13_h; // σ_xx for longitudinal wave, first potential
                Dh[2, 1] = Complex.ImaginaryOne * D1 * sin_k13_h; // σ_xx for longitudinal wave, second potential
                Dh[2, 2] = Complex.ImaginaryOne * D2 * k33 * sin_k33_h; // σ_xx for transverse wave, first potential
                Dh[2, 3] = -D2 * k33 * cos_k33_h; // σ_xx for transverse wave, second potential

                Dh[3, 0] = Complex.ImaginaryOne * D2 * k13 * sin_k13_h; // σ_zz for longitudinal wave, first potential
                Dh[3, 1] = -D2 * k13 * cos_k13_h; // σ_zz for longitudinal wave, second potential
                Dh[3, 2] = D1 * cos_k33_h; // σ_zz for transverse wave, first potential
                Dh[3, 3] = -Complex.ImaginaryOne * D1 * sin_k33_h; // σ_zz for transverse wave, second potential

                //This is what is in the Allard Atalla book for an explicit version of the inverted D(0) matrix, but it leads to negative absorption values.
                //SparseMatrix D0_Inverse = new SparseMatrix(4, 4);
                //D0_Inverse[0, 0] = 2 * kt / (w * d23);
                //D0_Inverse[0, 2] = -1 / (Mu * d23);
                //D0_Inverse[1, 1] = (k33 * k33 - kt * kt) / (w * k13 * d23);
                //D0_Inverse[1, 3] = -kt / (Mu * k13 * d23);
                //D0_Inverse[2, 1] = kt / (w * d23);
                //D0_Inverse[2, 3] = 1 / (Mu * d23);
                //D0_Inverse[3, 0] = (k33 * k33 - kt * kt) / (w * k33 * d23);
                //D0_Inverse[3, 2] = -kt / (Mu * k33 * d23);

                //return D0_Inverse * Dh;

                //Instead, building the D(0) matrix and inverting seems to work better.
                SparseMatrix D0 = new SparseMatrix(4, 4);
                // Build D0 matrix at z=0
                D0[0, 0] = w * kt;                      // w * kt * cos(0)
                D0[0, 1] = 0;                          // -w * i * kt * sin(0) = 0
                D0[0, 2] = 0;                          // w * i * k33 * sin(0) = 0
                D0[0, 3] = -w * k33;                   // -w * k33 * cos(0)

                D0[1, 0] = 0;                          // -w * i * k13 * sin(0) = 0
                D0[1, 1] = w * k13;                    // w * k13 * cos(0)
                D0[1, 2] = w * kt;                     // w * kt * cos(0)
                D0[1, 3] = 0;                          // -w * i * kt * sin(0) = 0

                D0[2, 0] = -D1;                        // -D1 * cos(0)
                D0[2, 1] = 0;                          // i * D1 * sin(0) = 0
                D0[2, 2] = 0;                          // i * D2 * k33 * sin(0) = 0
                D0[2, 3] = -D2 * k33;                  // -D2 * k33 * cos(0)

                D0[3, 0] = 0;                          // i * D2 * k13 * sin(0) = 0
                D0[3, 1] = -D2 * k13;                  // -D2 * k13 * cos(0)
                D0[3, 2] = D1;                         // D1 * cos(0)
                D0[3, 3] = 0;                          // -i * D1 * sin(0) = 0

                // Use the analytical inverse from Allard & Atalla if available
                // Otherwise use numerical inverse with regularization
                try
                {
                    return Dh * (D0.Inverse() as SparseMatrix);
                }
                catch
                {
                    // Add small regularization for numerical stability
                    var identity = SparseMatrix.CreateIdentity(4);
                    var regularizedD0 = D0 + Complex.ImaginaryOne * 1e-12 * identity;
                    return Dh * (regularizedD0.Inverse() as SparseMatrix);
                }
            }

            public static SparseMatrix InterfacePP(double porosity_of_1, double porosity_of_2)
            //    : base(6, 6)
            {
                SparseMatrix M = new SparseMatrix(6, 6);
                double phi2_1 = porosity_of_2 / porosity_of_1;
                double phi1_2 = porosity_of_1 / porosity_of_2;
                M[0, 0] = 1;
                M[1, 1] = 1;
                M[2, 1] = 1 - phi2_1;
                M[2, 2] = phi2_1;
                M[3, 3] = 1;
                M[3, 5] = 1 - phi1_2;
                M[4, 4] = 1;
                M[5, 5] = phi1_2;
                return M;
            }

            public static SparseMatrix InterfaceSF_Solid()
            //: base(3, 4)
            {
                SparseMatrix M = new SparseMatrix(3, 4);
                M[0, 1] = 1;
                M[1, 2] = 1;
                M[2, 3] = 1;
                return M;
            }

            public static SparseMatrix InterfaceSF_Fluid()
            //: base(3, 2)
            {
                SparseMatrix M = new SparseMatrix(3, 2);
                M[0, 1] = -1;
                M[1, 0] = 1;
                return M;
            }

            public static SparseMatrix Interfacepf_Porous(double porosity)
            //: base(4, 6)
            {
                SparseMatrix M = new SparseMatrix(4, 6);
                M[0, 1] = 1 - porosity;
                M[0, 2] = porosity;
                M[1, 3] = 1;
                M[2, 4] = 1;
                M[3, 5] = 1;
                return M;
            }

            public static SparseMatrix Interfacepf_Fluid(double porosity)
            //: base(4, 2)
            {
                SparseMatrix M = new SparseMatrix(4, 2);
                M[0, 1] = -1;
                M[1, 0] = 1 - porosity;
                M[3, 0] = porosity;
                return M;
            }

            public static SparseMatrix Interfacesp_Solid()
            //: base(5, 4)
            {
                SparseMatrix M = new SparseMatrix(5, 4);
                M[0, 0] = 1;
                M[1, 1] = 1;
                M[2, 1] = 1;
                M[3, 2] = 1;
                M[4, 3] = 1;
                return M;
            }

            public static SparseMatrix Interfacesp_Porous()
            {
                SparseMatrix M = new SparseMatrix(5, 6);
                M[0, 0] = 1;
                M[1, 1] = 1;
                M[2, 2] = 1;
                M[3, 3] = 1;
                M[3, 5] = 1;
                M[4, 4] = 1;
                return M;
            }

            public static SparseMatrix Interfacepi_Porous()
            //: base(4, 6)
            {
                SparseMatrix M = new SparseMatrix(4, 6);
                M[0, 1] = 1;
                M[1, 2] = 1;
                M[2, 3] = 1;
                M[2, 5] = 1;
                M[3, 4] = 1;
                return M;
            }

            public static SparseMatrix Interfacepi_Thinplate()
            {
                SparseMatrix M = new SparseMatrix(4, 2);
                M[0, 1] = -1;
                M[1, 1] = -1;
                M[2, 0] = 1;
                return M;
            }

            public static SparseMatrix RigidTerminationP()
            {
                SparseMatrix M = new SparseMatrix(3, 6);
                M[0, 0] = 1;
                M[1, 1] = 1;
                M[2, 2] = 1;
                return M;
            }

            public static SparseMatrix RigidTerminationF()
            {
                SparseMatrix M = new SparseMatrix(1, 2);
                M[0, 1] = 1;
                return M;
            }

            public static SparseMatrix RigidTerminationS()
            {
                SparseMatrix M = new SparseMatrix(2, 4);
                M[0, 0] = 1;
                M[1, 1] = 1;
                return M;
            }
        }
    }
}