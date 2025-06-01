using Hare.Geometry;
using Pachyderm_Acoustic.Environment;
using Pachyderm_Acoustic.Utilities;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Net.Mail;

namespace Pachyderm_Acoustic
{
    public class Hemisphere_Plot
    {
        Hare.Geometry.Point Ctr;
        Topology hemisphere;

        public Hemisphere_Plot(Hare.Geometry.Point Center)
        {
            hemisphere = Utilities.Geometry.GeoHemiSphere(5, 1);
            Ctr = Center;
        }

        public Hare.Geometry.Topology Output(IEnumerable<double> magnitude, double Min, double Max, double Diameter)
        {
            if (magnitude.Count() != hemisphere.Vertex_Count) throw new Exception("Invalid data input to spherical plot...");
            Hare.Geometry.Point[] points = new Hare.Geometry.Point[hemisphere.Vertex_Count];
            for(int i = 0; i < magnitude.Count(); i++)
            {
                double mag = (magnitude.ElementAt(i));
                if (double.IsInfinity(mag)) mag = 0;
                mag = Math.Max(mag, Min);
                mag = Math.Min(mag, Max);
                mag -= Min;
                mag /= (Max - Min) * Diameter;
                points[i] = mag * hemisphere[i] + Ctr;
            }
            Hare.Geometry.Topology T = Utilities.Geometry.GeoHemiSphere(5,1);

            for (int i = 0; i < T.Vertex_Count; i++) T.Set_Vertex(i, points[i]);
            T.Finish_Topology();
            return T;
        }

        public Point[] Vertices 
        {
            get
            {
                return hemisphere.Polygon_Vertices(0);
            }
        }
    }

    public class Sphere_Plot
    {
        Hare.Geometry.Point Ctr;
        Voxel_Grid Sphere;
        double[] alt, azi;

        public Sphere_Plot(Hare.Geometry.Point Center)
        {
            Sphere = Utilities.Geometry.GeoSphere(3);
            Ctr = Center;
            alt = new double[Sphere.Model[0].Vertex_Count];
            azi = new double[Sphere.Model[0].Vertex_Count];
            for(int i = 0; i < Sphere.Model[0].Vertex_Count; i++)
            {
                alt[i] = Math.Asin(Sphere.Model[0][i].z);
                azi[i] = Math.Atan2(Sphere.Model[0][i].y, Sphere.Model[0][i].x);
            }
        }

        public Hare.Geometry.Topology Output(IEnumerable<double> magnitude, double Min = double.PositiveInfinity, double Max = 0, double Diameter = .3)
        {
            int vert_ct = Sphere.Model[0].Vertex_Count;
            if (magnitude.Count() != vert_ct) throw new Exception("Invalid data input to spherical plot...");
            Hare.Geometry.Point[] points = new Hare.Geometry.Point[vert_ct];

            if (Max == 0) Max = magnitude.Max();
            if (Min == double.PositiveInfinity) Min = Max - 30;

            for (int i = 0; i < magnitude.Count(); i++)
            {
                double mag = (magnitude.ElementAt(i));
                if (double.IsInfinity(mag)) mag = 0;
                mag = Math.Max(mag, Min);
                mag = Math.Min(mag, Max);
                mag -= Min;
                mag /= (Max - Min);
                mag = Math.Max(0, mag);
                mag *= Diameter;
                points[i] = mag * Sphere.Model[0][i] + Ctr;
            }
            Hare.Geometry.Topology T = Utilities.Geometry.GeoSphere(3).Model[0];

            for (int i = 0; i < T.Vertex_Count; i++) T.Set_Vertex(i, points[i]);
            T.Finish_Topology();
            return T;
        }
        
        public IEnumerable<double> SPL_From_IR(int receiver_id, int octave, int sample_start, int sample_end, Direct_Sound[] Ds, ImageSourceData[] IS = null, Receiver_Bank[] R = null)
        {
            double[] values = new double[Sphere.Model[0].Vertex_Count];

            if (receiver_id < 0) return null;

            for (int s = 0; s < Ds.Length; s++)
            {
                if (Ds != null)
                {
                    Ctr = Ds[s].Rec_Origin.ElementAt(receiver_id);
                    int dsstart = (int)(Ds[s].Time_Pt[receiver_id] * 44100);
                    Vector[] detc = Ds[s].Dir_Energy(octave, receiver_id);
                    int dsend = dsstart + detc.Length;
                    if (!(dsend < sample_start || dsstart > sample_end))
                    {
                        int start = Math.Max(sample_start - dsstart, 0);
                        int end = Math.Min(detc.Length, sample_end - sample_start);

                        for (int i = 0; i < values.Length; i++)
                        {
                            Vector dir = new Vector(Sphere.Model[0][i]);

                            for (int t = start; t < end; t++)
                            {
                                double length = detc[t].Length();
                                Vector d = detc[t] / length;
                                double dot = Hare_math.Dot(dir, d);
                                if (dot > 0)
                                {
                                    double I = Math.Pow(dot, 16) * length;
                                    values[i - start] += double.IsNaN(I) ? 0 : I;
                                }
                            }
                        }
                    }
                }
                if (IS != null && IS[s] != null)
                {
                    for (int p = 0; p < IS[s].Paths[receiver_id].Count; p++)
                    {
                        int dsstart = (int)(IS[s].Paths[receiver_id][p].TravelTime * 44100);
                        Vector[] d = IS[s].Paths[receiver_id][p].Dir_Energy(octave);
                        int dsend = (int)(dsstart + d.Length);
                        if (!(dsend < sample_start || dsstart > sample_end))
                        {
                            int end = Math.Min(d.Length, sample_end - sample_start);
                            int start = Math.Max(sample_start - dsstart, 0);

                            for (int i = 0; i < values.Length; i++)
                            {
                                Vector dir = new Vector(Sphere.Model[0][i]);

                                for (int t = start; t < end; t++)
                                {
                                    double length = d[t].Length();
                                    Vector dsamp = d[t] / -length;
                                    double dot = Hare_math.Dot(dsamp, dir);
                                    if (dot > 0)
                                    {
                                        double I = Math.Pow(dot, 16) * length;
                                        values[i] += double.IsNaN(I) || I < 0 ? 0 : I;
                                    }
                                }
                            }
                        }
                    }
                }
                if (R != null && R[s] != null)
                {
                    for (int i = 0; i < values.Length; i++)
                    {
                        double[][] V = new double[6][];
                        int samples = sample_end - sample_start;
                        for (int d = 0; d < 6; d++) V[d] = new double[samples];

                        for (int t = 0; t < samples; t++)
                        {
                            if (t > R[s].CO_Time * 44100) continue;

                            int sample = t + sample_start;

                            Hare.Geometry.Vector pos = R[s].Rec_List[receiver_id].Directions_Pos(octave, sample);
                            Hare.Geometry.Vector neg = R[s].Rec_List[receiver_id].Directions_Neg(octave, sample);

                            V[0][t] = pos.dx;
                            V[1][t] = neg.dx;
                            V[2][t] = pos.dy;
                            V[3][t] = neg.dy;
                            V[4][t] = pos.dz;
                            V[5][t] = neg.dz;
                        }

                        V = PachTools.Rotate_Vector_Rose(V, -azi[i], -alt[i], false);

                        for (int t = 0; t < samples; t++)
                        {
                            //Hare.Geometry.Vector dir = new Hare.Geometry.Vector(Sphere.Model[0][i]);

                            double V2 = Math.Abs(V[t][2]), V3 = Math.Abs(V[t][3]), V4 = Math.Abs(V[t][4]), V5 = Math.Abs(V[t][5]);
                            Vector comp = new Vector(V[t][0], V2 > V3 ? V2 : -V3, V4 > V5 ? V4 : -V5);
                            //Vector comp = new Vector(dir.dx > 0 ? Math.Abs(pos.dx) : -Math.Abs(neg.dx), dir.dy > 0 ? Math.Abs(pos.dy) : -Math.Abs(neg.dy), dir.dz > 0 ? Math.Abs(pos.dz) : -Math.Abs(neg.dz));

                            double length = comp.Length();
                            V[t][0] /= length;

                            double I = Math.Pow(V[t][0], 16) * length;

                            values[i] += double.IsNaN(I) ? 0 : I;
                        }
                    }
                }
            }
            return AcousticalMath.SPL_Intensity_Signal(values);
        }
    }
}