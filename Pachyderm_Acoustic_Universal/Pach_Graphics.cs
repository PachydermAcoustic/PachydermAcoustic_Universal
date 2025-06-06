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

using Pachyderm_Acoustic.Utilities;
using Eto.Drawing;

namespace Pachyderm_Acoustic
{
    namespace Pach_Graphics
    {
        public interface IReporting
        {
            void Populate(double[] arrayx, double[] arrayy, string xdom, string ydom, int WaitDuration, int Prog_Percent);
            void Populate(double[] arrayx, double[] arrayy, double[] arrayx2, double[] arrayy2, string xdom, string ydom, int WaitDuration, int Prog_Percent);
            void Populate(double[] arrayx, double[] arrayy, double[] arrayx2, double[] arrayy2, double[] arrayx3, double[] arrayy3, string xdom, string ydom, int WaitDuration, int Prog_Percent);
            void Populate(double[] arrayx, double[] arrayy, double[] arrayx2, double[] arrayy2, double[] arrayx3, double[] arrayy3, double[] arrayx4, double[] arrayy4, string xdom, string ydom, int WaitDuration, int Prog_Percent);
            void Fill(double[] arrayx, double[] arrayy, double[] arrayx2, double[] arrayy2, double[] arrayx3, double[] arrayy3, double[] arrayx4, double[] arrayy4, string xdom, string ydom, int WaitDuration, int Prog_Percent);
        }

        public interface IProgressFeedback
        {
            void Report(int progress_percent);
            void change_title(string msg);
        }

        /// <summary>
        /// Colorscale base class.
        /// </summary>
        public abstract class Colorscale
        {
            public Bitmap PIC;
            public abstract Color GetValue(double Value, double VMin, double VMax);
        }

        /// <summary>
        /// Hue Saturation Value based colorscale.
        /// </summary>
        public class HSV_colorscale:Colorscale
        {            
            double H_OFFSET;
            double H_BREADTH;
            double S_OFFSET;
            double S_BREADTH;
            double V_OFFSET;
            double V_BREADTH;
            bool Discretized;
            int Steps;

            public HSV_colorscale(int Ht, int Wd, double H_OFFSETin, double H_BREADTHin, double S_OFFSETin, double S_BREADTHin, double V_OFFSETin, double V_BREADTHin, bool Discretize, int Stepsin)
            {
                H_OFFSET = H_OFFSETin;
                H_BREADTH = H_BREADTHin;
                S_OFFSET = S_OFFSETin;
                S_BREADTH = S_BREADTHin;
                V_OFFSET = V_OFFSETin;
                V_BREADTH = V_BREADTHin;
                Discretized = Discretize;
                Steps = Stepsin;
                PIC = new Bitmap(new Size(Wd, Ht), PixelFormat.Format32bppRgb);
                Graphics g = new Graphics(PIC);
                    
                if (Discretized)
                {
                    for (int y = 0; y < Ht; y++)
                    {
                        double loc = System.Math.Floor(Steps * (double)y / (double)Ht);
                        Color c = PachTools.HsvColor((((loc / Steps) * System.Math.PI) * H_BREADTH) + H_OFFSET, (loc / Steps * S_BREADTH) + S_OFFSET, (loc / Steps * V_BREADTH) + V_OFFSET);
                        g.DrawLine(new Pen(c, 1), 0, y, Ht, y);
                    }
                }
                else
                {
                    for (int y = 0; y < Ht; y++)
                    {
                        double loc = y / (double)Ht;
                        Color c = PachTools.HsvColor(((loc * System.Math.PI) * H_BREADTH) + H_OFFSET, (loc * S_BREADTH) + S_OFFSET, (loc * V_BREADTH) + V_OFFSET);
                        g.DrawLine(new Pen(c, 1), 0, y, Ht, y);
                    }
                }

                g.DrawImage(PIC, Wd, Ht);

                g.Dispose();
            }

            /// <summary>
            /// Assigns a color based on the value's location in the color scale.
            /// </summary>
            /// <param name="Value">the value</param>
            /// <param name="VMin">scale min</param>
            /// <param name="VMax">scale max</param>
            /// <returns></returns>
            public override Color GetValue(double Value, double VMin, double VMax)
            {
                Color c;

                if (Discretized)
                {
                    if (Value < VMin)
                    {
                        c = PachTools.HsvColor(H_BREADTH * System.Math.PI + H_OFFSET, S_BREADTH + S_OFFSET, V_OFFSET + V_BREADTH);
                    }
                    else if (Value > VMax)
                    {
                        c = PachTools.HsvColor(H_OFFSET, S_OFFSET, V_OFFSET);
                    }
                    else
                    {
                        double loc = System.Math.Floor(Steps * (1 - (Value - VMin) / (VMax - VMin)));
                        c = PachTools.HsvColor((loc / Steps * System.Math.PI) * H_BREADTH + H_OFFSET, loc * S_BREADTH + S_OFFSET, loc * V_BREADTH + V_OFFSET);
                    }
                }
                else
                {
                    if (Value < VMin)
                    {
                        c = PachTools.HsvColor(H_BREADTH * System.Math.PI + H_OFFSET, S_BREADTH + S_OFFSET, V_OFFSET + V_BREADTH);
                    }
                    else if (Value > VMax)
                    {
                        c = PachTools.HsvColor(H_OFFSET, S_OFFSET, V_OFFSET);
                    }
                    else
                    {
                        double loc = (1 - (Value - VMin) / (VMax - VMin));
                        c = PachTools.HsvColor((loc * System.Math.PI) * H_BREADTH + H_OFFSET, loc * S_BREADTH + S_OFFSET, loc * V_BREADTH + V_OFFSET);
                    }
                }
                return c;
            }
        }

        //For best results, only use this without a need for false color interpolation...
        public class Explicit_colorscale : Colorscale
        {
            Color[] C_Bounds;
            double[] R_f;
            double[] R_i;
            double[] G_f;
            double[] G_i;
            double[] B_f;
            double[] B_i;
            double[] CBR;
            double[] CBG;
            double[] CBB;
            double[] V_Bounds;
            double Scale_Breadth;

            public Explicit_colorscale(Color[] Colors, double[] Values)
            {
                C_Bounds = Colors;
                V_Bounds = Values;

                R_f = new double[C_Bounds.Length];
                R_i = new double[C_Bounds.Length];
                G_f = new double[C_Bounds.Length];
                G_i = new double[C_Bounds.Length];
                B_f = new double[C_Bounds.Length];
                B_i = new double[C_Bounds.Length];
                CBR = new double[C_Bounds.Length];
                CBG = new double[C_Bounds.Length];
                CBB = new double[C_Bounds.Length];


                for (int i = 0; i < Colors.Length; i++)
                {
                    R_f[i] = C_Bounds[i + 1].R;
                    R_i[i] = C_Bounds[i].R;
                    G_f[i] = C_Bounds[i + 1].G;
                    G_i[i] = C_Bounds[i].G;
                    B_f[i] = C_Bounds[i + 1].B;
                    B_i[i] = C_Bounds[i].B;
                    CBR[i] = R_f[i] - R_i[i];
                    CBG[i] = G_f[i] - G_i[i];
                    CBB[i] = B_f[i] - B_i[i];

                    Scale_Breadth = V_Bounds[1] - V_Bounds[0];

                    for (int y = 0; y < PIC.Height; y++)
                    {
                        int Current_Band = (int)System.Math.Floor((double)C_Bounds.Length * y / PIC.Height);
                        int R_r = (int)(((y - R_i[Current_Band]) * CBR[Current_Band] / (PIC.Height / C_Bounds.Length)) + R_i[Current_Band]);
                        int G_r = (int)(((y - B_i[Current_Band]) * CBG[Current_Band] / (PIC.Height / C_Bounds.Length)) + R_i[Current_Band]);
                        int B_r = (int)(((y - B_i[Current_Band]) * CBB[Current_Band] / (PIC.Height / C_Bounds.Length)) + R_i[Current_Band]);
                        Color.FromArgb(R_r, G_r, B_r);
                    }
                }
            }

            public override Color GetValue(double Value, double VMin, double VMax)
            {
                if (Value < V_Bounds[0])
                {
                    return C_Bounds[0];
                }
                else if (Value > V_Bounds[1])
                {
                    return C_Bounds[C_Bounds.Length - 1];
                }
                else
                {
                    int Current_Band = (int)System.Math.Floor((Value - V_Bounds[0]) / (C_Bounds.Length * Value / Scale_Breadth));
                    int R_r = (int)(((R_f[Current_Band] - R_i[Current_Band]) * ((Value - V_Bounds[0]) - (Current_Band * CBR[0])) / CBR[0]) + R_i[Current_Band]);
                    int G_r = (int)(((G_f[Current_Band] - G_i[Current_Band]) * ((Value - V_Bounds[0]) - (Current_Band * CBG[0])) / CBG[0]) + G_i[Current_Band]);
                    int B_r = (int)(((B_f[Current_Band] - B_i[Current_Band]) * ((Value - V_Bounds[0]) - (Current_Band * CBB[0])) / CBB[0]) + B_i[Current_Band]);
                    return Color.FromArgb(R_r, G_r, B_r);
                }
            }
        }
    }
}