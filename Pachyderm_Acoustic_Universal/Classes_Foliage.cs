// 'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)
// Copyright (c) 2008-2026, Open Research in Acoustical Science and Education, Inc.
// Distributed under the GNU General Public License, version 3 or later.

using System;
using System.Collections.Generic;
using Hare.Geometry;
using Pachyderm_Acoustic.Environment;

namespace Pachyderm_Acoustic.Environment
{
    /// <summary>Excess loss in a closed, outward-oriented stand envelope, in metres.</summary>
    public sealed class Foliage_Volume
    {
        internal readonly Topology Shape;
        internal readonly AABB Bounds;
        internal readonly double[] Coefficients;
        internal readonly double Tolerance;
        internal readonly bool Estimated;
        internal readonly double LeafWidth, LeafAreaDensity;
        internal readonly double[] LeafWidths, CanopyShares;

        /// <summary>
        /// Mean one-sided leaf area per canopy volume (m2/m3), from trees per hectare,
        /// typical circular crown diameter (m), crown LAI (m2 leaf / m2 projected crown), and canopy depth (m).
        /// Crowns share these assumptions; leaf areas add, including where crowns overlap.
        /// This geometric conversion does not provide species-specific LAI or acoustic calibration.
        /// </summary>
        public static double EstimateLeafAreaDensity(double treesPerHectare, double crownDiameter, double crownLeafAreaIndex, double canopyDepth)
        {
            if (!(treesPerHectare >= 0) || double.IsInfinity(treesPerHectare)) throw new ArgumentOutOfRangeException(nameof(treesPerHectare));
            if (!(crownDiameter >= 0) || double.IsInfinity(crownDiameter) || (treesPerHectare > 0 && crownDiameter == 0)) throw new ArgumentOutOfRangeException(nameof(crownDiameter));
            if (!(crownLeafAreaIndex >= 0) || double.IsInfinity(crownLeafAreaIndex)) throw new ArgumentOutOfRangeException(nameof(crownLeafAreaIndex));
            if (!(canopyDepth > 0) || double.IsInfinity(canopyDepth)) throw new ArgumentOutOfRangeException(nameof(canopyDepth));
            double density = treesPerHectare * (Math.PI * .25 * crownDiameter * crownDiameter) * crownLeafAreaIndex / (10000 * canopyDepth);
            if (double.IsNaN(density) || double.IsInfinity(density)) throw new ArgumentOutOfRangeException(nameof(treesPerHectare), "Estimated density exceeds the numeric range.");
            return density;
        }

        /// <summary>
        /// Legacy constant-coefficient conversion; retained for saved v3 estimates and existing callers.
        /// Cumulative volumes are prepared with the scene sound speed by Foliage_Field.
        /// Approximate Aylor (1972), JASA 51, 411-414, Fig. 3 (doi:10.1121/1.1912852).
        /// Width is in metres, leaf area density in m2/m3, reference depth in metres, sound speed in m/s.
        /// Returns equivalent octave dB/m at that depth, NOT the paper's nonlinear law at every depth.
        /// Corn/reed direct-beam scattering data are extrapolated to trees; ground, trunks and branches are excluded.
        /// </summary>
        public static double[] Estimate_Aylor(double leafWidth, double leafAreaDensity, double referenceDepth, double soundSpeed = 343)
        {
            if (!(leafWidth > 0) || double.IsInfinity(leafWidth)) throw new ArgumentOutOfRangeException(nameof(leafWidth));
            if (!(leafAreaDensity >= 0) || double.IsInfinity(leafAreaDensity)) throw new ArgumentOutOfRangeException(nameof(leafAreaDensity));
            if (!(referenceDepth > 0) || double.IsInfinity(referenceDepth)) throw new ArgumentOutOfRangeException(nameof(referenceDepth));
            if (!(soundSpeed > 0) || double.IsInfinity(soundSpeed)) throw new ArgumentOutOfRangeException(nameof(soundSpeed));
            // Coarse readings of the drawn trend, not published tabulated constants or a fitted species model.
            // h(ka) = excess dB / sqrt(F*L), ka = 2*pi*f*width/c. Linear interpolation;
            // below ka=.3 use zero, above 6.5 hold the drawn plateau. These extensions are assumptions.
            double[] ka = { 0, .3, .6, 1, 1.5, 2, 2.7, 3.5, 4.5, 5.5, 6.5 };
            double[] trend = { 0, 0, .5, .95, 1.45, 1.8, 2.15, 2.45, 2.7, 2.8, 2.8 };
            double[] result = Utilities.AcousticalMath.OctaveCenters();
            double scale = Math.Sqrt(leafAreaDensity) / Math.Sqrt(referenceDepth);
            for (int band = 0; band < result.Length; band++)
            {
                double x = 2 * Math.PI * result[band] * leafWidth / soundSpeed;
                int i = 1;
                while (i < ka.Length - 1 && x > ka[i]) i++;
                double h = x >= ka[ka.Length - 1] ? trend[trend.Length - 1] : trend[i - 1] + (trend[i] - trend[i - 1]) * (x - ka[i - 1]) / (ka[i] - ka[i - 1]);
                result[band] = h * scale;
                if (double.IsNaN(result[band]) || double.IsInfinity(result[band])) throw new ArgumentOutOfRangeException(nameof(leafAreaDensity), "Estimated loss exceeds the numeric range.");
            }
            return result;
        }

        /// <summary>Canopy-area weighted leaf-width responses at unit square-root distance; an uncalibrated mixed-stand approximation.</summary>
        public static double[] Estimate_Aylor(double[] leafWidths, double[] canopyShares, double leafAreaDensity, double soundSpeed)
        {
            ValidateMixture(leafWidths, canopyShares, leafAreaDensity);
            double[] result = new double[8];
            for (int species = 0; species < leafWidths.Length; species++)
            {
                double[] response = Estimate_Aylor(leafWidths[species], leafAreaDensity, 1, soundSpeed);
                for (int band = 0; band < result.Length; band++) result[band] += canopyShares[species] * response[band];
            }
            return result;
        }

        /// <summary>Interpolate dB loss coefficients on log frequency, clamping the outer thirds to the nearest octave.</summary>
        public static double[] InterpolateThirdOctaves(double[] octave)
        {
            if (octave == null || octave.Length != 8) throw new ArgumentException("Provide eight octave coefficients.", nameof(octave));
            foreach (double value in octave) if (!(value >= 0) || double.IsInfinity(value)) throw new ArgumentOutOfRangeException(nameof(octave));
            double[] x = Array.ConvertAll(Utilities.AcousticalMath.OctaveCenters(), frequency => Math.Log(frequency, 2));
            var interpolation = MathNet.Numerics.Interpolation.LinearSpline.InterpolateSorted(x, octave);
            double[] result = Utilities.AcousticalMath.ThirdOctaveCenters();
            for (int band = 0; band < result.Length; band++) result[band] = Math.Max(0, interpolation.Interpolate(Math.Max(x[0], Math.Min(x[7], Math.Log(result[band], 2)))));
            return result;
        }

        /// <summary>Aylor leaf-width estimate; the owning field receives the scene sound speed before tracing.</summary>
        public Foliage_Volume(Topology shape, double leafWidth, double leafAreaDensity, double tolerance) : this(shape, new[] { leafWidth }, new[] { 1.0 }, leafAreaDensity, tolerance) { }

        public Foliage_Volume(Topology shape, double[] leafWidths, double[] canopyShares, double leafAreaDensity, double tolerance) : this(shape, new double[8], tolerance)
        {
            ValidateMixture(leafWidths, canopyShares, leafAreaDensity);
            LeafWidths = (double[])leafWidths.Clone(); CanopyShares = (double[])canopyShares.Clone();
            LeafWidth = LeafWidths[0]; LeafAreaDensity = leafAreaDensity; Estimated = true;
        }

        static void ValidateMixture(double[] leafWidths, double[] canopyShares, double leafAreaDensity)
        {
            if (leafWidths == null || canopyShares == null || leafWidths.Length == 0 || leafWidths.Length != canopyShares.Length) throw new ArgumentException("Provide matching widths and canopy shares.");
            double sum = 0;
            for (int i = 0; i < leafWidths.Length; i++)
            {
                if (!(leafWidths[i] > 0) || double.IsInfinity(leafWidths[i])) throw new ArgumentOutOfRangeException(nameof(leafWidths));
                if (!(canopyShares[i] >= 0) || double.IsInfinity(canopyShares[i])) throw new ArgumentOutOfRangeException(nameof(canopyShares));
                sum += canopyShares[i];
            }
            if (Math.Abs(sum - 1) > 1e-9) throw new ArgumentException("Canopy shares must total one.", nameof(canopyShares));
            if (!(leafAreaDensity >= 0) || double.IsInfinity(leafAreaDensity)) throw new ArgumentOutOfRangeException(nameof(leafAreaDensity));
        }

        public Foliage_Volume(Topology shape, double[] attenuation, double tolerance)
        {
            if (shape == null || shape.Polys.Count == 0) throw new ArgumentException("A closed foliage mesh is required.", nameof(shape));
            if (attenuation == null || attenuation.Length != 8) throw new ArgumentException("Provide eight octave-band losses in dB/m.", nameof(attenuation));
            foreach (double value in attenuation) if (value < 0 || double.IsNaN(value) || double.IsInfinity(value)) throw new ArgumentOutOfRangeException(nameof(attenuation));
            if (!(tolerance > 0) || double.IsInfinity(tolerance)) throw new ArgumentOutOfRangeException(nameof(tolerance));
            Shape = shape; Coefficients = (double[])attenuation.Clone(); Tolerance = tolerance;
            Bounds = new AABB(shape.Min, shape.Max);
        }
    }

    /// <summary>Scene data constructed at setup. Overlapping stand assignments contribute additive excess loss.</summary>
    public sealed class Foliage_Field
    {
        readonly Foliage_Volume[] volumes;
        readonly double[][] coefficients, thirdOctaveCoefficients;
        bool prepared = true;
        public Foliage_Field(IEnumerable<Foliage_Volume> volumes)
        {
            this.volumes = new List<Foliage_Volume>(volumes).ToArray();
            coefficients = new double[this.volumes.Length][]; thirdOctaveCoefficients = new double[this.volumes.Length][];
            for (int i = 0; i < this.volumes.Length; i++)
            {
                coefficients[i] = this.volumes[i].Coefficients;
                thirdOctaveCoefficients[i] = this.volumes[i].Estimated ? null : Foliage_Volume.InterpolateThirdOctaves(coefficients[i]);
                prepared &= !this.volumes[i].Estimated;
            }
        }

        /// <summary>Prepare estimated spectra from the scene medium, once before tracing. Volumes may be shared; spectra belong to this field.</summary>
        public void SetSoundSpeed(double soundSpeed)
        {
            if (!(soundSpeed > 0) || double.IsInfinity(soundSpeed)) throw new ArgumentOutOfRangeException(nameof(soundSpeed));
            for (int i = 0; i < volumes.Length; i++)
                if (volumes[i].Estimated)
                {
                    coefficients[i] = Foliage_Volume.Estimate_Aylor(volumes[i].LeafWidths, volumes[i].CanopyShares, volumes[i].LeafAreaDensity, soundSpeed);
                    thirdOctaveCoefficients[i] = Foliage_Volume.InterpolateThirdOctaves(coefficients[i]);
                }
            prepared = true;
        }

        internal sealed class Path
        {
            internal struct Interval
            {
                internal double Start, End, Before, RootBefore, FullWeight, CachedDistance, CachedWeight;
                internal Foliage_Volume Volume;
                internal double[] Coefficients, ThirdOctaveCoefficients;
            }
            internal readonly List<Interval> Intervals = new List<Interval>();
            // Completed histories are immutable, so clones/split rays can share them without copying.
            internal Dictionary<Foliage_Volume, double> History;
            internal double Gain(int band, bool thirdOctave, double distance, bool pressure = false)
            {
                double loss = 0;
                for (int i = 0; i < Intervals.Count; i++)
                {
                    Interval interval = Intervals[i];
                    double length = Math.Max(0, Math.Min(distance, interval.End) - interval.Start), weight = length;
                    if (interval.Volume.Estimated)
                    {
                        if (distance >= interval.End) weight = interval.FullWeight;
                        else if (interval.CachedDistance == distance) weight = interval.CachedWeight;
                        else
                        {
                            // Stable sqrt(D+L)-sqrt(D); cache the receiver's partial distance across bands.
                            weight = length == 0 ? 0 : length / (Math.Sqrt(interval.Before + length) + interval.RootBefore);
                            interval.CachedDistance = distance; interval.CachedWeight = weight; Intervals[i] = interval;
                        }
                    }
                    loss += weight * (thirdOctave ? interval.ThirdOctaveCoefficients[band] : interval.Coefficients[band]);
                }
                return Math.Pow(10, (pressure ? -.05 : -.1) * loss);
            }

            internal Path Complete()
            {
                Dictionary<Foliage_Volume, double> next = null;
                foreach (Interval interval in Intervals)
                    if (interval.Volume.Estimated)
                    {
                        if (next == null) next = History == null ? new Dictionary<Foliage_Volume, double>() : new Dictionary<Foliage_Volume, double>(History);
                        next[interval.Volume] = interval.Before + (interval.End - interval.Start);
                    }
                History = next ?? History;
                Intervals.Clear();
                return History == null ? null : this;
            }
        }

        internal Path Trace(Point start, Point end, Path previous = null)
        {
            if (!prepared) throw new InvalidOperationException("Set foliage sound speed from the scene before tracing.");
            Path path = new Path { History = previous?.History };
            Vector direction = end - start;
            double length = direction.Length();
            if (!(length > 0)) return path;
            direction /= length;
            Ray ray = new Ray(start, direction, 0, 0);
            for (int volumeIndex = 0; volumeIndex < volumes.Length; volumeIndex++)
            {
                Foliage_Volume volume = volumes[volumeIndex];
                double entry = 0;
                Point point = default(Point);
                if (!volume.Bounds.Intersect(ray, ref entry, ref point) || entry > length) continue;
                // One intersection pass per candidate volume; receipts reuse the resulting intervals.
                var hits = new List<(double Distance, double Sign)>();
                foreach (Polygon polygon in volume.Shape.Polys)
                    if (polygon.Intersect(ray, polygon.Points, out point, out double u, out double v, out double distance, out int id))
                        hits.Add((distance, Hare_math.Dot(polygon.Normal, direction)));
                hits.Sort((a, b) => a.Distance.CompareTo(b.Distance));
                double begin = 0, traversed = 0, previousDistance = 0;
                if (volume.Estimated && path.History != null) path.History.TryGetValue(volume, out previousDistance);
                bool inside = false;
                for (int i = 0; i < hits.Count;)
                {
                    double distance = hits[i].Distance;
                    bool positive = false, negative = false;
                    do
                    {
                        positive |= hits[i].Sign > 0; negative |= hits[i].Sign < 0;
                        i++;
                    } while (i < hits.Count && hits[i].Distance - distance <= volume.Tolerance);
                    // Shared triangle edges count once; an entry and exit at a tangent cancel.
                    if (positive == negative) continue;
                    if (inside)
                    {
                        double a = Math.Max(0, begin), b = Math.Min(length, distance);
                        if (b > a)
                        {
                            double before = volume.Estimated ? previousDistance + traversed : 0, root = volume.Estimated ? Math.Sqrt(before) : 0, lengthInside = b - a;
                            path.Intervals.Add(new Path.Interval { Start = a, End = b, Before = before, RootBefore = root, FullWeight = volume.Estimated ? lengthInside / (Math.Sqrt(before + lengthInside) + root) : lengthInside, CachedDistance = -1, Volume = volume, Coefficients = coefficients[volumeIndex], ThirdOctaveCoefficients = thirdOctaveCoefficients[volumeIndex] });
                            traversed += lengthInside;
                        }
                    }
                    inside = !inside;
                    begin = distance;
                }
            }
            return path;
        }

        /// <summary>Attenuate one complete deterministic path, retaining each stand's accumulated distance across its legs.</summary>
        public void Attenuate(double[] spectrum, Point[] vertices, bool pressure = false)
        {
            if (spectrum.Length != 8 && spectrum.Length != 24) throw new ArgumentException("Expected eight or twenty-four bands.", nameof(spectrum));
            Path history = null;
            for (int leg = 1; leg < vertices.Length; leg++)
            {
                Path path = Trace(vertices[leg - 1], vertices[leg], history);
                for (int band = 0; band < spectrum.Length; band++) spectrum[band] *= path.Gain(band, spectrum.Length == 24, double.PositiveInfinity, pressure);
                history = path.Complete();
            }
        }

        public void Attenuate(double[] spectrum, Point start, Point end, bool pressure = false)
        {
            if (spectrum.Length != 8 && spectrum.Length != 24) throw new ArgumentException("Expected eight or twenty-four bands.", nameof(spectrum));
            Path path = Trace(start, end);
            for (int band = 0; band < spectrum.Length; band++) spectrum[band] *= path.Gain(band, spectrum.Length == 24, double.PositiveInfinity, pressure);
        }
    }
}

namespace Pachyderm_Acoustic.Environment
{
    [Serializable]
    internal sealed class Foliage_Receiver : Spherical_Receiver
    {
        internal Foliage_Receiver(Point p, Point source, double[] attenuation, double speed, double rho, int rate, double cutoff) : base(p, source, attenuation, speed, rho, rate, cutoff) { }
        protected override double MediumLoss(BroadRay ray, int band, double distance) { return ray.Foliage?.Gain(band, ray.ThirdOctave, distance) ?? 1; }
        protected override double MediumLoss(OctaveRay ray, int band, double distance) { return ray.Foliage?.Gain(band, ray.ThirdOctave, distance) ?? 1; }
    }

    [Serializable]
    internal sealed class Foliage_Expanding_Receiver : Expanding_Receiver
    {
        internal Foliage_Expanding_Receiver(Point p, Point source, int rays, double[] attenuation, double speed, double rho, int rate, double cutoff) : base(p, source, rays, attenuation, speed, rho, rate, cutoff) { }
        protected override double MediumLoss(BroadRay ray, int band, double distance) { return ray.Foliage?.Gain(band, ray.ThirdOctave, distance) ?? 1; }
        protected override double MediumLoss(OctaveRay ray, int band, double distance) { return ray.Foliage?.Gain(band, ray.ThirdOctave, distance) ?? 1; }
    }

    [Serializable]
    internal sealed class Foliage_Map_Receiver : PachMapReceiver.Map_Receiver
    {
        internal Foliage_Map_Receiver(double diameter, int index, Point p, Source source, double speed, double rho, double[] attenuation, int rate, double cutoff, uint boxes, bool onePoint, bool directional) : base(diameter, index, p, source, speed, rho, attenuation, rate, cutoff, boxes, onePoint, directional) { }
        protected override double MediumLoss(BroadRay ray, int band, double distance) { return ray.Foliage?.Gain(band, ray.ThirdOctave, distance) ?? 1; }
        protected override double MediumLoss(OctaveRay ray, int band, double distance) { return ray.Foliage?.Gain(band, ray.ThirdOctave, distance) ?? 1; }
    }
}
