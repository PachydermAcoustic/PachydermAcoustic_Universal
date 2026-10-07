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

        public Sphere_Plot(Hare.Geometry.Point Center)
        {
            Sphere = Utilities.Geometry.GeoSphere(3);
            Ctr = Center;
        }

        public Hare.Geometry.Topology Output(IEnumerable<double> magnitude, double Min = double.PositiveInfinity, double Max = 0, double Diameter = .3)
        {
            if (magnitude == null) throw new ArgumentNullException(nameof(magnitude));
            double[] levels = magnitude.ToArray();
            int count = Sphere.Model[0].Vertex_Count;
            if (levels.Length != count) throw new ArgumentException("Invalid data input to spherical plot...", nameof(magnitude));
            if (Max == 0)
                Max = levels.Where(level => !double.IsNaN(level) && !double.IsInfinity(level)).DefaultIfEmpty(0).Max();
            if (double.IsPositiveInfinity(Min)) Min = Max - 30;

            Hare.Geometry.Topology output = Utilities.Geometry.GeoSphere(3).Model[0];
            for (int i = 0; i < count; i++)
            {
                double radius = 0;
                if (!double.IsNaN(levels[i]) && !double.IsInfinity(levels[i]) && Max > Min)
                    radius = Math.Max(0, Math.Min(1, (levels[i] - Min) / (Max - Min))) * Diameter;
                output.Set_Vertex(i, radius * Sphere.Model[0][i] + Ctr);
            }
            output.Finish_Topology();
            return output;
        }


        // Window indices use the 44.1 kHz display timebase; the end is exclusive.
        public IEnumerable<double> SPL_From_IR(int receiver_id, int octave, int sample_start, int sample_end, Direct_Sound[] Ds, ImageSourceData[] IS = null, Receiver_Bank[] R = null)
        {
            double[] values = new double[Sphere.Model[0].Vertex_Count];
            sample_start = Math.Max(0, sample_start);
            if (receiver_id < 0 || octave < 0 || octave > 8 || sample_end <= sample_start)
                return DirectionalLevels(values);

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

                // The histogram stores signed component sums, not individual rays.
                // Estimate octant lobes directly in world coordinates. Independent
                // sign weights retain opposing arrivals without cancelling them.
                for (int t = startSample; t < endSample; t++)
                {
                    Vector pos = receiver.Directions_Pos(octave, t);
                    Vector neg = receiver.Directions_Neg(octave, t);
                    double xPos = Math.Abs(pos.dx), xNeg = Math.Abs(neg.dx);
                    double yPos = Math.Abs(pos.dy), yNeg = Math.Abs(neg.dy);
                    double zPos = Math.Abs(pos.dz), zNeg = Math.Abs(neg.dz);
                    double x = xPos + xNeg, y = yPos + yNeg, z = zPos + zNeg;
                    double length = new Vector(x, y, z).Length();
                    if (length <= 0 || double.IsNaN(length) || double.IsInfinity(length)) continue;
                    double intensity = receiver.Energy(t, octave);
                    if (double.IsNaN(intensity) || double.IsInfinity(intensity)) continue;
                    if (intensity <= 0) intensity = length;

                    for (int octant = 0; octant < 8; octant++)
                    {
                        bool nx = (octant & 1) != 0, ny = (octant & 2) != 0, nz = (octant & 4) != 0;
                        double wx = x > 0 ? (nx ? xNeg : xPos) / x : (nx ? 0 : 1);
                        double wy = y > 0 ? (ny ? yNeg : yPos) / y : (ny ? 0 : 1);
                        double wz = z > 0 ? (nz ? zNeg : zPos) / z : (nz ? 0 : 1);
                        double weight = wx * wy * wz;
                        if (weight <= 0) continue;
                        Vector direction = new Vector(nx ? -x : x, ny ? -y : y, nz ? -z : z) / length;
                        AccumulateDirection(values, direction, intensity * weight);
                    }
                }

            }
            return DirectionalLevels(values);
        }

        private static double[] DirectionalLevels(double[] values)
        {
            // Silence must stay below the relative display floor, even when
            // the largest valid level is below 30 dB or below 0 dB.
            return values.Select(value => value > 0 ? 10 * Math.Log10(value / 1E-12) : double.NegativeInfinity).ToArray();
        }


        private void AccumulateDirectionalEnergy(double[] values, Vector[] energy, int start, int end, bool reverse)
        {
            for (int t = start; t < end; t++)
            {
                double length = energy[t].Length();
                if (length <= 0 || double.IsNaN(length) || double.IsInfinity(length)) continue;
                Vector direction = energy[t] / (reverse ? -length : length);
                AccumulateDirection(values, direction, length);
            }
        }

        private void AccumulateDirection(double[] values, Vector direction, double energy)
        {
            for (int i = 0; i < values.Length; i++)
            {
                double dot = Hare_math.Dot(new Vector(Sphere.Model[0][i]), direction);
                if (dot <= 0) continue;
                double intensity = Math.Pow(Math.Min(1, dot), 16) * energy;
                if (!double.IsNaN(intensity) && !double.IsInfinity(intensity))
                    values[i] += intensity;
            }
        }


    }
}