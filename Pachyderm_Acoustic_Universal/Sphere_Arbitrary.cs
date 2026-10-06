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
        
        // Window indices use the 44.1 kHz display timebase; the end is exclusive.
        public IEnumerable<double> SPL_From_IR(int receiver_id, int octave, int sample_start, int sample_end, Direct_Sound[] Ds, ImageSourceData[] IS = null, Receiver_Bank[] R = null)
        {
            double[] values = new double[Sphere.Model[0].Vertex_Count];
            sample_start = Math.Max(0, sample_start);
            if (receiver_id < 0 || octave < 0 || octave > 8 || sample_end <= sample_start)
                return AcousticalMath.SPL_Intensity_Signal(values);

            int sourceCount = Math.Max(Ds == null ? 0 : Ds.Length,
                Math.Max(IS == null ? 0 : IS.Length, R == null ? 0 : R.Length));
            for (int s = 0; s < sourceCount; s++)
            {
                Direct_Sound direct = Ds != null && s < Ds.Length ? Ds[s] : null;
                if (direct != null && receiver_id < direct.Rec_Origin.Count())
                {
                    Ctr = direct.Rec_Origin.ElementAt(receiver_id);
                    double rate = direct.SampleRate > 0 ? direct.SampleRate : 44100;
                    int arrival = (int)(direct.Time_Pt[receiver_id] * rate);
                    Vector[] energy = direct.Dir_Energy(octave, receiver_id);
                    int start = Math.Max(0, (int)(sample_start * rate / 44100) - arrival);
                    int end = Math.Min(energy.Length, (int)(sample_end * rate / 44100) - arrival);
                    AccumulateDirectionalEnergy(values, energy, start, end, false);
                }

                ImageSourceData images = IS != null && s < IS.Length ? IS[s] : null;
                if (images != null && images.Paths != null && receiver_id < images.Paths.Length
                    && images.Paths[receiver_id] != null)
                {
                    foreach (Deterministic_Reflection path in images.Paths[receiver_id])
                    {
                        if (path == null) continue;
                        int arrival = (int)(path.TravelTime * 44100);
                        Vector[] energy = path.Dir_Energy(octave);
                        int start = Math.Max(0, sample_start - arrival);
                        int end = Math.Min(energy.Length, sample_end - arrival);
                        AccumulateDirectionalEnergy(values, energy, start, end, true);
                    }
                }

                Receiver_Bank bank = R != null && s < R.Length ? R[s] : null;
                if (bank == null || bank.SampleRate <= 0 || bank.Rec_List == null
                    || receiver_id >= bank.Rec_List.Length || bank.Rec_List[receiver_id] == null)
                    continue;

                var receiver = bank.Rec_List[receiver_id];
                int startSample = Math.Max(0, (int)(sample_start * (double)bank.SampleRate / 44100));
                int endSample = Math.Min(bank.Duration(), (int)(sample_end * (double)bank.SampleRate / 44100));
                // Use the stored histogram length as well as the bank's declared duration.
                endSample = Math.Min(endSample, receiver.Recs.Energy[octave < 8 ? octave : 0].Length);
                int samples = endSample - startSample;
                if (samples <= 0) continue;

                double[][] directions = new double[6][];
                for (int d = 0; d < 6; d++) directions[d] = new double[samples];
                for (int t = 0; t < samples; t++)
                {
                    int sample = startSample + t;
                    Vector pos = receiver.Directions_Pos(octave, sample);
                    Vector neg = receiver.Directions_Neg(octave, sample);
                    directions[0][t] = pos.dx;
                    directions[1][t] = neg.dx;
                    directions[2][t] = pos.dy;
                    directions[3][t] = neg.dy;
                    directions[4][t] = pos.dz;
                    directions[5][t] = neg.dz;
                }

                for (int i = 0; i < values.Length; i++)
                {
                    // Rotation returns one six-component row per sample.
                    double[][] rotated = PachTools.Rotate_Vector_Rose(directions, -azi[i], -alt[i], false);
                    for (int t = 0; t < samples; t++)
                    {
                        double yPos = Math.Abs(rotated[t][2]), yNeg = Math.Abs(rotated[t][3]);
                        double zPos = Math.Abs(rotated[t][4]), zNeg = Math.Abs(rotated[t][5]);
                        Vector comp = new Vector(rotated[t][0],
                            yPos > yNeg ? yPos : -yNeg, zPos > zNeg ? zPos : -zNeg);
                        double length = comp.Length();
                        if (length <= 0 || double.IsNaN(length) || double.IsInfinity(length)) continue;
                        double intensity = Math.Pow(rotated[t][0] / length, 16) * length;
                        if (!double.IsNaN(intensity) && !double.IsInfinity(intensity))
                            values[i] += intensity;
                    }
                }
            }
            return AcousticalMath.SPL_Intensity_Signal(values);
        }

        private void AccumulateDirectionalEnergy(double[] values, Vector[] energy, int start, int end, bool reverse)
        {
            for (int t = start; t < end; t++)
            {
                double length = energy[t].Length();
                if (length <= 0 || double.IsNaN(length) || double.IsInfinity(length)) continue;
                Vector direction = energy[t] / (reverse ? -length : length);
                for (int i = 0; i < values.Length; i++)
                {
                    double dot = Hare_math.Dot(new Vector(Sphere.Model[0][i]), direction);
                    if (dot <= 0) continue;
                    double intensity = Math.Pow(dot, 16) * length;
                    if (!double.IsNaN(intensity) && !double.IsInfinity(intensity))
                        values[i] += intensity;
                }
            }
        }

    }
}