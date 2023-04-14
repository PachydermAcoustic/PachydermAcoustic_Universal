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
//using System.Windows;
using Hare.Geometry;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.Providers.LinearAlgebra;
using MathNet.Spatial.Euclidean;
//using MathNet.Numerics.Providers.LinearAlgebra;

namespace Pachyderm_Acoustic
{
    public abstract class Edge
    {
        public static int[] Rigid = new int[] { 1, 1, 1, 1 };
        public static int[] Soft = new int[] { -1, 1, 1, -1 };
        protected List<EdgeSource> Sources = new List<EdgeSource>();
        protected int Plane1;
        protected int Plane2;

        public Edge (int p1, int p2)
        {
            Plane1 = p1;
            Plane2 = p2;
        }
        public class EdgeSource
        {
            const double sincos45 = 0.70710678118654752440084436210485;
            double[] v;
            double[] v_4pi;
            double[] v_2pi;
            public Hare.Geometry.Point[] Z_limits;
            public Hare.Geometry.Point Z_mid;
            public Hare.Geometry.Vector Z_Norm;
            double Z_Range;
            double Z_Range_2;
            double Z_dot;//Z_Range squared.
            Hare.Geometry.Vector[] Normal = new Hare.Geometry.Vector[2];
            Hare.Geometry.Vector[] Tangent = new Hare.Geometry.Vector[2];
            int[] attr;
            public int Poly_1;
            public int Poly_2;
            public bool isthinplate;

            public EdgeSource(int[] attr_in, Hare.Geometry.Point PtZ0, Hare.Geometry.Point _PtZ, Hare.Geometry.Vector[] _Tangents, int poly1, int poly2 = -1)
            {
                isthinplate = poly2 < 0; 
                attr = attr_in;
                Tangent = _Tangents;
                Z_Norm = _PtZ - PtZ0;
                Poly_1 = poly1;
                Poly_2 = poly2;

                Z_Range = Z_Norm.Length();
                Z_dot = Z_Range * Z_Range;//Hare_math.Dot(Z_Norm, Z_Norm);
                Z_Norm /= Z_Range;
                Z_Range_2 = Z_Range / 2;
                Z_mid = (PtZ0 + _PtZ) / 2;
                Hare.Geometry.Vector Bisector = (Tangent[0] + Tangent[1]) / 2;
                Bisector.Normalize();
                double BisectAngle = Math.Acos(Hare_math.Dot(Tangent[0], Bisector));
                if (BisectAngle == 0)
                {
                    BisectAngle = 1E-12;
                    v = new double[2] { Math.PI / (Utilities.Numerics.PiX2 - 2 * BisectAngle), Math.PI / (Utilities.Numerics.PiX2 - 2 * BisectAngle) };
                }
                else//v = new double[2] { Math.PI / (Utilities.Numerics.PiX2 - 2 * BisectAngle), Math.PI / (2 * BisectAngle) };
                {
                    v = new double[2] { Math.PI / (2 * BisectAngle), Math.PI / (Utilities.Numerics.PiX2 - 2 * BisectAngle) };
                }
                v_4pi = new double[2] { v[0] / (4 * Math.PI), v[1] / (4 * Math.PI) };
                v_2pi = new double[2] { v[0] / (2 * Math.PI), v[1] / (2 * Math.PI) };
                Normal[0] = Hare_math.Cross(_Tangents[0], Z_Norm);
                Normal[1] = Hare_math.Cross(_Tangents[1], Z_Norm * -1);

                if (Hare_math.Dot(Normal[0], Bisector) >= 0) Normal[0] *= -1;
                if (Hare_math.Dot(Normal[0], Bisector) >= 0) Normal[0] *= -1;
                Z_limits = new Hare.Geometry.Point[2] { PtZ0, _PtZ };

            }

            public EdgeSource(int[] attr_in, Hare.Geometry.Point Z_mid_in, double Delta_Z, Hare.Geometry.Vector[] _Tangents, int poly1, int poly2 = -1)
            {
                attr = attr_in;
                Tangent = _Tangents;
                Poly_1 = poly1;
                Poly_2 = poly2;
                Z_Norm = Hare.Geometry.Hare_math.Cross(_Tangents[0], _Tangents[1]);
                Z_Range = Delta_Z;
                Z_dot = Z_Range * Z_Range;//Hare_math.Dot(Z_Norm, Z_Norm);
                Z_Range_2 = Z_Range / 2;
                Z_Norm.Normalize();
                Z_mid = Z_mid_in;
                Hare.Geometry.Vector Bisector = (Tangent[0] + Tangent[1]) / 2;
                Bisector.Normalize();
                double BisectAngle = Math.Acos(Hare_math.Dot(Tangent[0], Bisector));
                if (BisectAngle == 0) BisectAngle = 1E-12;
                v = new double[2] { Math.PI / (2 * BisectAngle), Math.PI / (Utilities.Numerics.PiX2 - 2 * BisectAngle) };
                v_4pi = new double[2] { v[0] / (4 * Math.PI), v[1] / (4 * Math.PI) };
                v_2pi = new double[2] { v[0] / (2 * Math.PI), v[1] / (2 * Math.PI) };
                Normal[0] = Hare_math.Cross(_Tangents[0], Z_Norm);
                Normal[1] = Hare_math.Cross(_Tangents[1], Z_Norm * -1);

                if (Hare_math.Dot(Normal[0], Bisector) > 0)
                {
                    Normal[0] *= -1;
                    Normal[1] *= -1;
                }
            }

            public bool Cyl_Coord(Hare.Geometry.Point S, Hare.Geometry.Point R, ref double rs, ref double thetas, ref double zs, ref double rr, ref double thetar, ref double zr, out int Obtuse_Side)//, out double[] tm, out double[] tl)
            {
                //diffx = Tangent;
                //diffy = Normal;
                //diffz = Z_Norm;
                Hare.Geometry.Vector S_D = S - Z_mid;
                Hare.Geometry.Vector R_D = R - Z_mid;

                zs = Hare_math.Dot(S_D, Z_Norm);//S_Coord.z;
                zr = Hare_math.Dot(R_D, Z_Norm);//R_Coord.z;

                Hare.Geometry.Point S_p = S - (Z_mid + zs * Z_Norm);
                Hare.Geometry.Point R_p = R - (Z_mid + zr * Z_Norm);

                rs = Math.Sqrt(S_p.x * S_p.x + S_p.y * S_p.y + S_p.z * S_p.z);//Math.Sqrt(S_Coord.x * S_Coord.x + S_Coord.y * S_Coord.y + S_Coord.z * S_Coord.z);
                rr = Math.Sqrt(R_p.x * R_p.x + R_p.y * R_p.y + R_p.z * R_p.z);//Math.Sqrt(R_Coord.x * R_Coord.x + R_Coord.y * R_Coord.y + R_Coord.z * R_Coord.z);

                Hare.Geometry.Vector S_Norm = new Vector(S_p.x, S_p.y, S_p.z);
                S_Norm.Normalize();
                Hare.Geometry.Vector R_Norm = new Vector(R_p.x, R_p.y, R_p.z);
                R_Norm.Normalize();

                uint SDIR = 0;

                if (isthinplate)
                {
                    if (Hare_math.Dot(S_Norm, Tangent[0]) < 0) SDIR = 1;
                }
                else
                {
                    double S0 = Math.Abs(Hare_math.Dot(S_Norm, Tangent[0]));
                    double S1 = Math.Abs(Hare_math.Dot(S_Norm, Tangent[1]));
                    if (S1 > S0) SDIR = 1;
                }
                Obtuse_Side = (Hare_math.Dot(S_Norm, Normal[SDIR]) > 0 ? 1 : 0);

                Vector AltSNorm = Hare_math.Cross(S_Norm, Z_Norm);
                Vector AltRNorm = Hare_math.Cross(R_Norm, Z_Norm);
                thetas = Math.Acos(Hare_math.Dot(AltSNorm, Normal[SDIR]));//Math.Atan2(S_Coord.y, S_Coord.x);

                double rdt = Hare_math.Dot(Normal[SDIR] * (Obtuse_Side == 1 ? 1 : -1), AltRNorm);

                if (rdt > 0)
                {
                    thetar = Utilities.Numerics.PiX2 - Math.Acos(Math.Abs(Hare_math.Dot(AltRNorm, Normal[SDIR])));//Math.Atan2(R_Coord.y, R_Coord.x);
                }
                else
                {
                    thetar = Math.Acos(Math.Abs(Hare_math.Dot(AltRNorm, Normal[SDIR])));
                }
                return true;
            }

            public double[] Phi(double thetaS, double thetaR)
            {
                return new double[] {
                        Math.PI + thetaS + thetaR,
                        Math.PI + thetaS - thetaR,
                        Math.PI - thetaS + thetaR,
                        Math.PI - thetaS - thetaR };
            }

            public double[] Beta(double thetaS, double thetaR, double Aux_N, int obtuse)
            {
                double[] phi = Phi(thetaS, thetaR);
                double[] B = new double[4];
                double v_A = v[obtuse] * Aux_N;
                double v_P = v[obtuse] * phi[0];
                B[0] = attr[0] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                v_P = v[obtuse] * phi[1];
                B[1] = attr[1] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                v_P = v[obtuse] * phi[2];
                B[2] = attr[2] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                v_P = v[obtuse] * phi[3];
                B[3] = attr[3] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                return B;
            }

            public double Beta_i(double thetaS, double thetaR, double Aux_N, int obtuse)
            {
                double[] phi = Phi(thetaS, thetaR);
                double v_A = v[obtuse] * Aux_N;
                double v_P = v[obtuse] * phi[0];
                double B = attr[0] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                v_P = v[obtuse] * phi[1];
                B += attr[1] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                v_P = v[obtuse] * phi[2];
                B += attr[2] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                v_P = v[obtuse] * phi[3];
                B += attr[3] * Math.Sin(v_P) / (Math.Cosh(v_A) - Math.Cos(v_P));
                return B;
            }

            public double Z_apex(double Z_S, double Z_R, double r_S, double r_R)
            {
                return (Z_R * r_S + Z_S * r_R) / (r_R + r_S);
            }

            public bool Is_Apex(Hare.Geometry.Point Src, Hare.Geometry.Point Rec)
            {
                //Hare.Geometry.Vector T_Comp = Hare.Geometry.Hare_math.Cross(this.Z_mid - Src, this.Z_mid - Rec);
                //return (Z_Norm.x * T_Comp.x + Z_Norm.y * T_Comp.y + Z_Norm.z * T_Comp.z) > .95;
                //Do this as line closest point operation between two points. t > 1 or t < 0 will result in better certainty that the apex lies within this segment.
                return false;
            }

            public double Calc_Pressure(double B, double m_in, double l_out, bool obtuse)
            {
                return (-v_4pi[Convert.ToInt32(obtuse)] * B) / (m_in * l_out);
            }

            public double Aux_n(double d_S, double z_S, double z_R, double r_S, double r_R, double thetaS, double thetaR)
            {
                double z_rel = (z_R - z_S);
                double Z = (d_S * d_S - (r_S * r_S + r_R * r_R + z_rel * z_rel)) / (2 * r_S * r_R);
                if (Z < 1)
                    return 0;
                return Math.Log(Z + Math.Sqrt(Z * Z - 1));
            }

            public double Aux_n(double alpha_in, double gamma_out, double thetaS, double thetaR)
            {
                double Z = (1 + Math.Sin(alpha_in) * Math.Sin(gamma_out)) / (Math.Cos(alpha_in) * Math.Cos(gamma_out));
                return Math.Log(Z + Math.Sqrt(Z * Z - 1));
            }

            public double R_o(double zr, double zs, double rr, double rs)
            {
                double r = rs + rr;
                double z = zr - zs;
                return Math.Sqrt(r * r + z * z);
            }

            public double Rho(double rr, double rs)
            {
                return rr / rs;
            }

            public double sinpsi(double rs, double rr, double R0)
            {
                return (rs + rr) / R0;
            }

            public double cospsi(double zs, double zr, double R0)
            {
                return (zr + zs) / R0;
            }

            public double[] B_0(double R0, double rho, double[] phi, double rs, double rr, int obtuse)
            {
                double[] B = new double[4];
                double rt1_2 = (1 + rho);
                rt1_2 *= rt1_2;
                double sinpsi = (rs + rr) / R0;
                double mod = 4 * R0 * R0 * rho * rho * rho / (v[obtuse] * v[obtuse] * rt1_2 * rt1_2 * (rt1_2 * sinpsi * sinpsi - 2 * rho));

                for (int i = 0; i < 4; i++)
                {
                    B[i] = attr[i] * Math.Sin(v[obtuse] * phi[i]) * mod;
                }

                return B;
            }

            public double B_0(double R0, double rho, double[] phi, double rr, double rs, double sinpsi, int obtuse)
            {
                double B = 0;
                double rt1_2 = (1 + rho);
                rt1_2 *= rt1_2;
                double mod = 4 * R0 * R0 * rho * rho * rho / (v[obtuse] * v[obtuse] * rt1_2 * rt1_2 * (rt1_2 * sinpsi * sinpsi - 2 * rho));

                for (int i = 0; i < 4; i++)
                {
                    B += attr[i] * Math.Sin(v[obtuse] * phi[i]);
                }

                return B * mod;
            }

            public double B_1(double R0, double rho, double[] phi, int obtuse)
            {
                double B = 0;
                double v_2 = v[obtuse] * 0.5;
                double rt1_2 = (1 + rho);
                rt1_2 *= rt1_2;
                for (int i = 0; i < 4; i++)
                {
                    double svp_2 = Math.Sin(phi[i] * v_2);
                    B += attr[i] * svp_2 * svp_2;
                }
                return B * 4 * R0 * R0 * rho * rho / (v[obtuse] * v[obtuse] * rt1_2 * rt1_2);
            }

            public double B_2(double R0, double rho, double sinpsi, double cospsi)
            {
                double rt1 = (1 + rho);
                return 2 * R0 * (1 - rho) * rho * cospsi / (rt1 * (rt1 * rt1 * sinpsi * sinpsi - 2 * rho));
            }

            public double B_3(double R0, double rho, double sinpsi)
            {
                double rt1_2 = (1 + rho);
                rt1_2 *= rt1_2;
                return 2 * R0 * R0 * rho * rho / (rt1_2 * (rt1_2 * sinpsi * sinpsi - 2 * rho));
            }

            public double B_4(double[] phi, int obtuse)
            {
                double B4 = 0;
                for (int i = 0; i < 4; i++)
                {
                    B4 += attr[i] * Math.Sin(v[obtuse] * phi[i]);
                }
                return B4 / (2 * v[obtuse] * v[obtuse]);
            }

            public double B_4(double B0, double B3)
            {
                return B0 / B3;
            }

            public double B_5(double B0, double B2)
            {
                return B0 / B2;
            }

            public double B_5(double R0, double rho, double[] phi, double cospsi, int obtuse)
            {
                double B5 = 0;
                for (int i = 0; i < 4; i++)
                {
                    B5 += attr[i] * Math.Sin(v[obtuse] * phi[i]);
                }
                double rt1 = rho + 1;
                return B5 * 2 * R0 * rho * rho / (v[obtuse] * v[obtuse] * rt1 * rt1 * rt1 * (1 - rho) * cospsi);
            }

            public double B_6(double B3, double B2)
            {
                return B3 / B2;
            }

            public double B_6(double R0, double rho, double cospsi)
            {
                return R0 * rho / ((1 - rho * rho) * cospsi);
            }

            public double F(double B2, double B3)
            {
                double q = 4 * B3 - B2 * B2;
                if (q < 0)
                {
                    double rtQ = Math.Sqrt(-q);
                    double ZR2 = 2 * Z_Range;
                    double rtQB2 = rtQ * B2;
                    return Math.Log((ZR2 + B2 - rtQB2 + rtQ) / (ZR2 + B2 + rtQB2 - rtQ), Math.E) / rtQ;
                }
                else if (q > 0)
                {
                    double rtQ = 1 / Math.Sqrt(q);
                    double ZR2 = 2 * Z_Range;
                    return 2 * rtQ * (Math.Atan((ZR2 + B2) * rtQ) - Math.Atan(B2 * rtQ));
                }

                return (4 * Z_Range) / (B2 * (2 * Z_Range + B2));
            }

            public double Apex_Solve(double Za, double zs, double zr, double rs, double rr, double thetas, double thetar, int obtuse, out double m, out double l)
            {
                double rho = Rho(rr, rs);
                int side = Convert.ToInt32(obtuse);
                m = Math.Sqrt(rs * rs + zs * zs);
                l = Math.Sqrt(rr * rr + zr * zr);
                double R0 = m + l;//R_o(zr, zs, rr, rs);
                double[] phi = Phi(thetas, thetar);
                double sinpsi = this.sinpsi(rs, rr, R0);
                double B1 = B_1(R0, rho, phi, obtuse);
                double sqrt_B1_inv = 1 / Math.Sqrt(B1);
                if (rho == 1 || zs == zr)
                {
                    //symmetrical case
                    if ((zr - zr) / R0 == sincos45)
                    {
                        double B4 = B_4(phi, obtuse);
                        return -v_2pi[side] * (B4 * sqrt_B1_inv) * Math.Atan(Z_Range * sqrt_B1_inv);
                    }
                    else
                    {
                        double B3 = B_3(R0, rho, sinpsi);
                        double sqrt_B3_inv = 1 / Math.Sqrt(B3);
                        double B0 = B_0(R0, rho, phi, rr, rs, sinpsi, obtuse);
                        return -v_2pi[side] * (B0 / (B3 - B1)) * (sqrt_B1_inv * Math.Atan(Z_Range * sqrt_B1_inv) - sqrt_B3_inv * Math.Atan(Z_Range * sqrt_B3_inv));
                    }
                }
                //asymmetrical case
                double rt1 = 1 + rho;
                double cospsi = this.cospsi(zs, zr, R0);
                if (sinpsi * sinpsi == 2 * rho / (rt1 * rt1))
                {
                    double B5 = B_5(R0, rho, phi, cospsi, obtuse);
                    double B6 = B_6(R0, rho, cospsi);
                    double B6sqr = B6 * B6;
                    double ZRB6 = Z_Range + B6;
                    return v_2pi[side] * (B5 * B5 / (B1 + B6sqr)) * (.5 * Math.Log(Math.Abs(B6sqr * (Z_dot + B1) / (B1 * (ZRB6 * ZRB6))), Math.E) - B6 * sqrt_B1_inv * Math.Atan(Z_Range * sqrt_B1_inv));
                }
                else
                {
                    double B0 = B_0(R0, rho, phi, rs, rr, sinpsi, obtuse);
                    double B2 = B_2(R0, rho, sinpsi, cospsi);
                    double B3 = B_3(R0, rho, sinpsi);
                    double B2_2 = B2 * B2;
                    double B1_B3 = B1 - B3;
                    double F = this.F(B2, B3);
                    return v_2pi[side] * (B0 * B2 / (B1 * B2_2 + B1_B3 * B1_B3)) * (0.5 * Math.Log(Math.Abs((B3 * (Z_dot + B1)) / B1 * (Z_dot + B2 * Z_Range + B3)), Math.E) + ((2 * B1_B3 - B2_2) * sqrt_B1_inv / B2) * Math.Atan(Z_Range * sqrt_B1_inv) + ((-2 * B1_B3 - B2_2) / (2 * B2)) * F);
                }
            }

            //public double Apex_Solve(Hare.Geometry.Point src, Hare.Geometry.Point rec, ref double m, ref double l)
            //{
            //    double rr = 0, rs = 0, zr = 0, zs = 0, thetar = 0, thetas = 0;
            //    int obtuse;
            //    if (!Cyl_Coord(src, rec, ref rs, ref thetas, ref zs, ref rr, ref thetar, ref zr, out obtuse, out )) return 0;
            //    m = Math.Sqrt(rs * rs + zs * zs);
            //    l = Math.Sqrt(rr * rr + zr * zr);
            //    double rho = Rho(rr, rs);
            //    double R0 = m + l;//R_o(zr, zs, rr, rs);
            //    double[] phi = Phi(thetas, thetar);
            //    double sinpsi = this.sinpsi(rs, rr, R0);
            //    double B1 = B_1(R0, rho, phi, obtuse);
            //    double sqrt_B1_inv = 1 / Math.Sqrt(B1);
            //    if (rho == 1 || zs == zr)
            //    {
            //        //symmetrical case
            //        if ((zr - zr) / R0 == sincos45)
            //        {
            //            double B4 = B_4(phi, obtuse);
            //            return -v_2pi[obtuse] * (B4 * sqrt_B1_inv) * Math.Atan(Z_Range * sqrt_B1_inv);
            //        }
            //        else
            //        {
            //            double B3 = B_3(R0, rho, sinpsi);
            //            double sqrt_B3_inv = 1 / Math.Sqrt(B3);
            //            double B0 = B_0(R0, rho, phi, rr, rs, sinpsi, obtuse);
            //            return -v_2pi[obtuse] * (B0 / (B3 - B1)) * (sqrt_B1_inv * Math.Atan(Z_Range * sqrt_B1_inv) - sqrt_B3_inv * Math.Atan(Z_Range * sqrt_B3_inv));
            //        }
            //    }
            //    //asymmetrical case
            //    double rt1 = 1 + rho;
            //    double cospsi = this.cospsi(zs, zr, R0);
            //    if (sinpsi * sinpsi == 2 * rho / (rt1 * rt1))
            //    {
            //        double B5 = B_5(R0, rho, phi, cospsi, obtuse);
            //        double B6 = B_6(R0, rho, cospsi);
            //        double B6sqr = B6 * B6;
            //        double ZRB6 = Z_Range + B6;
            //        return v_2pi[obtuse] * (B5 * B5 / (B1 + B6sqr)) * (.5 * Math.Log(Math.Abs(B6sqr * (Z_dot + B1) / (B1 * (ZRB6 * ZRB6))), Math.E) - B6 * sqrt_B1_inv * Math.Atan(Z_Range * sqrt_B1_inv));
            //    }
            //    else
            //    {
            //        double B0 = B_0(R0, rho, phi, rs, rr, sinpsi, obtuse);
            //        double B2 = B_2(R0, rho, sinpsi, cospsi);
            //        double B3 = B_3(R0, rho, sinpsi);
            //        double B2_2 = B2 * B2;
            //        double B1_B3 = B1 - B3;
            //        double F = this.F(B2, B3);
            //        return v_2pi[obtuse] * (B0 * B2 / (B1 * B2_2 + B1_B3 * B1_B3)) * (0.5 * Math.Log(Math.Abs((B3 * (Z_dot + B1)) / B1 * (Z_dot + B2 * Z_Range + B3)), Math.E) + ((2 * B1_B3 - B2_2) * sqrt_B1_inv / B2) * Math.Atan(Z_Range * sqrt_B1_inv) + ((-2 * B1_B3 - B2_2) / (2 * B2)) * F);
            //    }
            //}

            public double Flex_Solve(Hare.Geometry.Point src, Hare.Geometry.Point rec, ref double m, ref double l, ref double[] dM, ref double[] dL)
            {
                double zr = 0, zs = 0, thetar = 0, thetas = 0;
                int obtuse;
                double rr = 0, rs = 0;
                if (!Cyl_Coord(src, rec, ref rs, ref thetas, ref zs, ref rr, ref thetar, ref zr, out obtuse)) return 0;
                double Za = Z_apex(zs, zr, rs, rr);

                ////get the range of times that this sample occupies./////////////////////////////////
                //double z0 = Math.Abs(zs);
                //double zl = (z0 - Z_Range_2;
                //double zu = z0 + Z_Range_2;
                //tm = new double[2] { Math.Sqrt(rs * rs + zl * zl), Math.Sqrt(rs * rs + zu * zu) };
                //z0 = Math.Abs(zr);
                //zl = z0 - Z_Range_2;
                //zu = z0 + Z_Range_2;
                //tl = new double[2] { Math.Sqrt(rr * rr + zl * zl), Math.Sqrt(rr * rr + zu * zu) };
                //////////////////////////////////////////////////////////////////////////////////////

                //////////////////////////////////
                //double zm = Za - zs, zl = Za - zr;
                double zm = Math.Abs(zs);
                double zl = Math.Abs(zr);
                double zml = zm - Z_Range_2;
                double zmu = zm + Z_Range_2;
                double zll = zl - Z_Range_2;
                double zlu = zl + Z_Range_2;
                double rs_2 = rs * rs, rr_2 = rr * rr;
                dM = new double[] { Math.Sqrt(zml * zml + rs_2), Math.Sqrt(zmu * zmu + rs_2) };
                dL = new double[] { Math.Sqrt(zll * zll + rr_2), Math.Sqrt(zlu * zlu + rr_2) };
                //////////////////////////////////
                if (Math.Abs(Za) < Z_Range_2)
                {
                    double B = Apex_Solve(Za, zs, zr, rs, rr, thetas, thetar, obtuse, out m, out l);
                    dM[0] = m;
                    dL[0] = l;
                    return B * Z_Range;
                }
                return Gen_Solve(Za, zs, zr, rs, rr, thetas, thetar, obtuse, ref m, ref l) * Z_Range;
            }

            //public double Gen_Solve(Hare.Geometry.Point src, Hare.Geometry.Point rec, ref double m, ref double l)
            //{
            //    double rr = 0, rs = 0, zr = 0, zs = 0, thetar = 0, thetas = 0;
            //    int obtuse;
            //    if (!Cyl_Coord(src, rec, ref rs, ref thetas, ref zs, ref rr, ref thetar, ref zr, out obtuse)) return 0;
            //    double Za = Z_apex(zs, zr, rs, rr);
            //    double zl = Za - zs;
            //    double zm = Za - zr;
            //    m = Math.Sqrt((rs * rs) + (zl * zl));
            //    l = Math.Sqrt((rr * rr) + (zm * zm));
            //    double N = this.Aux_n(m + l, zs, zr, rs, rr, thetas, thetar);
            //    double  Beta = Beta_i(thetas, thetar, N, obtuse);
            //    return -v_4pi[obtuse] * Beta / (m * l);
            //}

            public double Gen_Solve(double Za, double zs, double zr, double rs, double rr, double thetas, double thetar, int obtuse, ref double m, ref double l)
            {
                //double zl = Za - zs;
                //double zm = Za - zr;
                m = Math.Sqrt((rs * rs) + (zs * zs));
                l = Math.Sqrt((rr * rr) + (zr * zr));
                double N = this.Aux_n(m + l, zs, zr, rs, rr, thetas, thetar);
                double Beta = Beta_i(thetas, thetar, N, obtuse);
                double B = -v_4pi[obtuse] * Beta; /// / (m * l)
                return B;
            }

        }

        /// <summary>
        /// Returns the list of Apexes for this source/receiver combination.
        /// </summary>
        /// <param name="src">source point</param>
        /// <param name="rec">receiver point</param>
        /// <returns></returns>
        public List<int> Find_Apex(Hare.Geometry.Point src, Hare.Geometry.Point rec)
        {
            List<int> Apex = new List<int>();
            for (int i = 0; i < Sources.Count; i++)
            {
                Vector T_Comp = Hare.Geometry.Hare_math.Cross(Sources[i].Z_mid - src, Sources[i].Z_mid - rec);
                if ((Sources[i].Z_Norm.x * T_Comp.x + Sources[i].Z_Norm.y * T_Comp.y + Sources[i].Z_Norm.z * T_Comp.z) > .95) Apex.Add(i);
            }
            return Apex;
        }

        public int find_poly_id(Environment.Polygon_Scene Room, Point c, int Plane_ID)
        {
            double d = double.PositiveInfinity;
            int p1 = 0;
                foreach (int poly in Room.ObjectMembers[Plane_ID])
                {
                    double d_temp = Room.Hare_Data.Polys[poly].SqDistanceToEdges(c);
                    if (d_temp < d)
                    {
                        p1 = poly;
                        d = d_temp;
                    } 
                }
            return p1;
        }

        public List<EdgeSource> EdgeSources
        {
            get
            {
                return Sources;
            }
        }

        public int[] SurfaceIDs
        {
            get { return new int[2] { Plane1, Plane2 }; }
        }
    }


public class Edge_Straight : Edge
    {
        Hare.Geometry.Point PointA;
        Hare.Geometry.Point PointB;

        public Edge_Straight(Environment.Polygon_Scene Room, IEnumerable<Hare.Geometry.Point> SPT, IEnumerable<Hare.Geometry.Point> RPT, Environment.Medium_Properties Att_Props, bool isSoft, Hare.Geometry.Point[] PointS, MathNet.Numerics.Interpolation.CubicSpline[] TangentA, MathNet.Numerics.Interpolation.CubicSpline[] TangentB, int plane_1, int plane_2 = -1)
        : base(plane_1, plane_2)
        {
            PointA = PointS[0];
            PointB = PointS[PointS.Length-1];

            //Get the open wedge angle. This needs the source location, and occlusion info...
            //Assuming tangent angles are always correctly oriented...
            Vector Z_Dir = PointA - PointB;
            double length = Z_Dir.Length();
            Z_Dir.Normalize();
            List<Hare.Geometry.Point> Dpt = new List<Hare.Geometry.Point>();

            //Find the secondary source spacing DeltaZ
            Dpt.AddRange(RPT);
            Dpt.AddRange(SPT);

            //Two options.
            //1. Use a constant step
            double MinAngle = double.PositiveInfinity;
            for (int j = 0; j < Dpt.Count; j++)
            {
                double angle = Math.Abs(Hare_math.Dot(PointA - Dpt[j], Z_Dir));
                if (angle < MinAngle) MinAngle = angle;
                angle = Math.Abs(Hare_math.Dot(PointB - Dpt[j], Z_Dir));
                if (angle < MinAngle) MinAngle = angle;
            }
            double fs = 176400; //Hz.
            double DeltaZ = Att_Props.Sound_Speed(this.PointA) / (fs * MinAngle);//TODO: Adjust depending on distance from source to receiver... (nearest, farthest?)

            double El_Ct = Math.Ceiling(length / DeltaZ);

            Random r = new Random();

            Hare.Geometry.Point Pt1 = new Hare.Geometry.Point(PointA.x, PointA.y, PointA.z);

            for (int i = 1; i < El_Ct; i++)
            {
                Hare.Geometry.Point Pt2 = PointA - i * Z_Dir * DeltaZ;
                Point c = (Pt1 + Pt2) / 2; 
                Vector[] HTangents = new Vector[2] { new Vector(TangentA[0].Interpolate(i), TangentA[1].Interpolate(i), TangentA[2].Interpolate(i)), new Vector(TangentB[0].Interpolate(i), TangentB[1].Interpolate(i), TangentB[2].Interpolate(i)) };
                Sources.Add(new EdgeSource(Edge.Rigid, Pt1, Pt2, HTangents, find_poly_id(Room, c, Plane1), Plane2 < 0? -1 : find_poly_id(Room, c, Plane2)));
                Pt1 = Pt2;
            }

            //2. Modify DeltaZ per change in velocity.
            //Point pt1 = PointA;
            //double DeltaZ = 0;
            ////for (int i = 1; i < El_Ct; i++)
            //for(double i = 0; i < length; i += DeltaZ)
            //{
            //    double MinAngle = double.PositiveInfinity;
            //    foreach (Point pt in Dpt)
            //    {
            //        double angle = Math.Abs(Hare_math.Dot(PointA - pt, Z_Dir));
            //        if (angle < MinAngle) MinAngle = angle;
            //    }
            //    DeltaZ = Att_Props.Sound_Speed(pt1) / (fs * MinAngle);
            //    Hare.Geometry.Point Pt2 = PointA - Z_Dir * i;
            //    Sources.Add(new EdgeSource(Edge.Rigid, Pt1, Pt2, new Vector[2] { new Vector(TangentA[0].Interpolate(i), TangentA[1].Interpolate(i), TangentA[2].Interpolate(i)), new Vector(TangentB[0].Interpolate(i), TangentB[1].Interpolate(i), TangentB[2].Interpolate(i)) }));
            //    Pt1 = Pt2;
            //}
        }
    }

    public class Edge_Curved : Edge
    {
        public Edge_Curved(Environment.Polygon_Scene Room, Hare.Geometry.Point From, Point To, bool issoft, double EdgeLength, Environment.Medium_Properties Att_Props, MathNet.Numerics.Interpolation.CubicSpline[] Edge_Crv, MathNet.Numerics.Interpolation.CubicSpline[] TangentA, int Plane_1, int Plane_2 = -1)//MathNet.Numerics.Interpolation.CubicSpline[] TangentB,
        : base(Plane_1, Plane_2)
        {
            //Find the secondary source spacing DeltaZ
            double fs = 176400; //Hz.
            double DeltaZ = 0;
            
            Point pt1 = new Point(Edge_Crv[0].Interpolate(0), Edge_Crv[1].Interpolate(0), Edge_Crv[2].Interpolate(0));
            Vector Z_Dir = pt1 - new Point(Edge_Crv[0].Interpolate(0.001), Edge_Crv[1].Interpolate(0.001), Edge_Crv[2].Interpolate(0.001));
            Z_Dir.Normalize();
            List<Hare.Geometry.Point> Dpt = new List<Hare.Geometry.Point>();

            //Find the secondary source spacing DeltaZ
            double i = 0;
            do
            {
                double MinAngle = double.NegativeInfinity;
                Vector D1 = (pt1 - To);
                Vector D2 = (pt1 - From);
                D1.Normalize(); D2.Normalize();
                double angle = Math.Abs(Hare_math.Dot((D2 + D1) / 2, Z_Dir));
                if (angle > MinAngle) MinAngle = angle;
                DeltaZ = Math.Min(Att_Props.Sound_Speed(pt1) / (2 * fs * MinAngle), 0.1);

                i += DeltaZ;
                if (i > EdgeLength) break;
                Hare.Geometry.Point pt2 = new Point(Edge_Crv[0].Interpolate(i), Edge_Crv[1].Interpolate(i), Edge_Crv[2].Interpolate(i));
                Point c = (pt1 + pt2) / 2;
                Sources.Add(new EdgeSource(Edge.Rigid, pt1, pt2, new Vector[2] { new Vector(TangentA[0].Interpolate(i), TangentA[1].Interpolate(i), TangentA[2].Interpolate(i)), new Vector(TangentA[3].Interpolate(i), TangentA[4].Interpolate(i), TangentA[5].Interpolate(i)) }, find_poly_id(Room, c, Plane1), Plane2 < 0 ? -1 : find_poly_id(Room, c, Plane_2) ));
                pt1 = pt2;
            } while (true);
        }

        public Edge_Curved(Environment.Polygon_Scene Room, IEnumerable<Hare.Geometry.Point> SPT, IEnumerable<Hare.Geometry.Point> RPT, bool issoft, double EdgeLength, Environment.Medium_Properties Att_Props, MathNet.Numerics.Interpolation.CubicSpline[] Edge_Crv, MathNet.Numerics.Interpolation.CubicSpline[] TangentA, int Plane_1, int Plane_2 = -1)// MathNet.Numerics.Interpolation.CubicSpline[] TangentB,
        : base(Plane_1, Plane_2)
        {
            //Find the secondary source spacing DeltaZ
            double fs = 176400; //Hz.
            double DeltaZ = 0;
            Point pt1 = new Point(Edge_Crv[0].Interpolate(0), Edge_Crv[1].Interpolate(0), Edge_Crv[2].Interpolate(0));
            Vector Z_Dir = pt1 - new Point(Edge_Crv[0].Interpolate(0.001), Edge_Crv[1].Interpolate(0.001), Edge_Crv[2].Interpolate(0.001));
            Z_Dir.Normalize();
            List<Hare.Geometry.Point> Dpt = new List<Hare.Geometry.Point>();

            //Find the secondary source spacing DeltaZ
            double i = 0;
            do
            {
                double MinAngle = double.NegativeInfinity;

                foreach (Point spt in SPT)
                {
                    foreach (Point pt in RPT)
                    {
                        Vector D1 = (pt1 - pt);
                        Vector D2 = (pt1 - spt);
                        D1.Normalize(); D2.Normalize();
                        double angle = Math.Abs(Hare_math.Dot((D2 + D1) / 2, Z_Dir));
                        //double angle = 1 - Math.Abs(Hare_math.Dot(Hare_math.Cross(D2,D1), Z_Dir));
                        if (angle > MinAngle) MinAngle = angle;
                    } 
                }
                DeltaZ = Math.Min(Att_Props.Sound_Speed(pt1) / (2 * fs * MinAngle), 0.1);

                i += DeltaZ;
                if (i > EdgeLength) break;
                Hare.Geometry.Point pt2 = new Point(Edge_Crv[0].Interpolate(i), Edge_Crv[1].Interpolate(i), Edge_Crv[2].Interpolate(i));
                //find the particular polygons that apply to the center of this node.
                Point c = (pt1 + pt2) / 2;
                int p1 = find_poly_id(Room, c, Plane_1);
                int p2 = Plane2 >= 0 ? find_poly_id(Room, c, Plane_2) : -1;

                Sources.Add(new EdgeSource(Edge.Rigid, pt1, pt2, new Vector[2] { new Vector(TangentA[0].Interpolate(i), TangentA[1].Interpolate(i), TangentA[2].Interpolate(i)), new Vector(TangentA[3].Interpolate(i), TangentA[4].Interpolate(i), TangentA[5].Interpolate(i)) }, p1, p2));
                Z_Dir = pt2 - pt1;
                Z_Dir.Normalize();
                pt1 = pt2;
            } while (true);
        }
    }
}