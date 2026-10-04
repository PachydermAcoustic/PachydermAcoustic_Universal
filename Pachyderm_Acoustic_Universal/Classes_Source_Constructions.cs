using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System;
using System.Collections.Generic;
using Hare.Geometry;

namespace Pachyderm_Acoustic
{
    namespace Source_Constructions
    {
        public interface ICabinet_Diffraction
        {
            System.Numerics.Complex Pressure(Point driver, Hare.Geometry.Vector direction, int octave, double driverDiameter, double distance);
            string[] Driver_Balloon(Point driver, double driverDiameter);
        }

        /// <summary>
        /// Experimental exterior Neumann radiation model using interior fundamental sources.
        /// A boundary-integral projection imposes zero normal velocity on rigid faces
        /// and uniform velocity on a finite circular piston (Neumann radiation problem).
        /// Default continuation uses this field through 500 Hz, a complex blend at 1 kHz,
        /// and the existing DED field from 2 kHz upward. Disable continuation for research only.
        /// Edges and vertices belong to this whole-cabinet solution, not additive point sources.
        /// </summary>
        public class Cabinet_Diffraction_Numerical : ICabinet_Diffraction
        {
            private readonly double Width, Height, Depth, SoundSpeed;
            private readonly int SamplesPerWavelength;
            private readonly bool HighFrequencyContinuation;
            private readonly Cabinet_Diffraction DED;
            private readonly Band[] Bands = new Band[8];

            private sealed class Band
            {
                internal Point[] Sources;
                internal Matrix<Complex> Matrix;
                internal LU<Complex> Factor;
                internal Field LastField;
            }

            private sealed class Field
            {
                internal Point Driver;
                internal double Diameter, K;
                internal Point[] Sources;
                internal Complex[] Strength;

                internal Complex At(Point receiver)
                {
                    Complex pressure = Complex.Zero;
                    for (int j = 0; j < Sources.Length; j++)
                    {
                        double r = (receiver - Sources[j]).Length();
                        pressure += Strength[j] * Complex.FromPolarCoordinates(1.0 / r, -K * r);
                    }
                    if (!Finite(pressure)) throw new InvalidOperationException("Numerical cabinet pressure is not finite.");
                    return pressure;
                }

                internal Complex NormalAt(Point receiver, Hare.Geometry.Vector normal)
                {
                    Complex value = Complex.Zero;
                    for (int j = 0; j < Sources.Length; j++) value += Strength[j] * NormalDerivative(receiver, Sources[j], normal, K);
                    return value;
                }
            }

            public Cabinet_Diffraction_Numerical(double width, double height, double depth, double sound_speed = 343.0, int samplesPerWavelength = 4, bool highFrequencyContinuation = true)
            {
                if (!(width > 0) || !(height > 0) || !(depth > 0) || !(sound_speed > 0) || double.IsInfinity(width + height + depth + sound_speed)) throw new ArgumentOutOfRangeException("width", "Numerical cabinet dimensions and sound speed must be finite and positive.");
                if (samplesPerWavelength < 3 || samplesPerWavelength > 8) throw new ArgumentOutOfRangeException("samplesPerWavelength", "Use 3 through 8 source samples per wavelength.");
                Width = width; Height = height; Depth = depth; SoundSpeed = sound_speed; SamplesPerWavelength = samplesPerWavelength;
                HighFrequencyContinuation = highFrequencyContinuation;
                DED = new Cabinet_Diffraction(width, height, depth, sound_speed);
                double maximumFrequency = highFrequencyContinuation ? 1000 : 8000;
                double spacing = Math.Min(0.7 * Math.Max(Width, Math.Max(Height, Depth)) / (1.5 * SamplesPerWavelength), SoundSpeed / (maximumFrequency * SamplesPerWavelength));
                double nx = Math.Max(2, Math.Ceiling(0.7 * Width / spacing)), ny = Math.Max(2, Math.Ceiling(0.7 * Depth / spacing)), nz = Math.Max(2, Math.Ceiling(0.7 * Height / spacing));
                double sourceCount = 2 * (nx * ny + nx * nz + ny * nz);
                double FaceCount(double a, double b, double inset)
                {
                    double step = Math.Min(SoundSpeed / (8 * maximumFrequency), inset * 0.5);
                    return Math.Max(8, Math.Ceiling(a / step)) * Math.Max(8, Math.Ceiling(b / step));
                }
                double integrationCount = 2 * (FaceCount(Width, Height, 0.15 * Depth) + FaceCount(Depth, Height, 0.15 * Width) + FaceCount(Width, Depth, 0.15 * Height));
                if (sourceCount * integrationCount > 16000000) throw new ArgumentException("Numerical cabinet exceeds the boundary-matrix memory limit. Use a smaller cabinet or Vanderkooy DED.");
                if (sourceCount > 2048) throw new ArgumentException("Numerical cabinet exceeds the 2048-source limit in its numerical bands. Use a smaller cabinet or Vanderkooy DED.");
            }

            //Log-frequency continuation at the available octave centers: 500 Hz -> 1 kHz -> 2 kHz.
            private double NumericalWeight(int octave)
            {
                return !HighFrequencyContinuation || octave <= 3 ? 1 : octave == 4 ? 0.5 : 0;
            }

            private static bool Finite(Complex value)
            {
                return !double.IsNaN(value.Real) && !double.IsNaN(value.Imaginary) && !double.IsInfinity(value.Real) && !double.IsInfinity(value.Imaginary);
            }

            private static Complex NormalDerivative(Point receiver, Point source, Hare.Geometry.Vector normal, double k)
            {
                Hare.Geometry.Vector path = receiver - source;
                double r = path.Length();
                double cosine = Hare.Geometry.Hare_math.Dot(path, normal) / r;
                return Complex.FromPolarCoordinates(1.0 / r, -k * r) * new Complex(-1.0 / r, -k) * cosine;
            }

            private void ValidateDriver(Point driver, double diameter)
            {
                double radius = diameter * 0.5;
                if (!(diameter > 0) || double.IsInfinity(diameter) || double.IsNaN(driver.x + driver.y + driver.z) || double.IsInfinity(driver.x + driver.y + driver.z) || Math.Abs(driver.y) > 1E-9 || Math.Abs(driver.x) + radius > Width * 0.5 + 1E-9 || Math.Abs(driver.z) + radius > Height * 0.5 + 1E-9) throw new ArgumentException("The finite piston must fit entirely on the front baffle at Y = 0.");
            }

            private Field Prepare(Point driver, int octave, double diameter)
            {
                ValidateDriver(driver, diameter);
                double radius = diameter * 0.5;
                double k = 2.0 * Math.PI * Cabinet_Diffraction.Frequencies[octave] / SoundSpeed;
                Band band;
                lock (Bands)
                {
                    band = Bands[octave];
                    if (band == null)
                    {
                        double spacing = Math.Min(0.7 * Math.Max(Width, Math.Max(Height, Depth)) / (1.5 * SamplesPerWavelength), SoundSpeed / (Cabinet_Diffraction.Frequencies[octave] * SamplesPerWavelength));
                        int nx = Math.Max(2, (int)Math.Ceiling(0.7 * Width / spacing)), ny = Math.Max(2, (int)Math.Ceiling(0.7 * Depth / spacing)), nz = Math.Max(2, (int)Math.Ceiling(0.7 * Height / spacing));
                        List<Point> integration = new List<Point>(), sources = new List<Point>();
                        List<Hare.Geometry.Vector> normals = new List<Hare.Geometry.Vector>();
                        List<double> weights = new List<double>();
                        void AddFace(Point origin, Hare.Geometry.Vector u, Hare.Geometry.Vector v, Hare.Geometry.Vector normal, int nu, int nv, double inset)
                        {
                            for (int i = 0; i < nu; i++) for (int j = 0; j < nv; j++)
                                {
                                    Point p = origin + u * ((i + 0.5) / nu) + v * ((j + 0.5) / nv);
                                    sources.Add(new Point(0.7 * p.x, -Depth * 0.5 + 0.7 * (p.y + Depth * 0.5), 0.7 * p.z));
                                }
                            double step = Math.Min(SoundSpeed / (8 * Cabinet_Diffraction.Frequencies[octave]), inset * 0.5);
                            var qu = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, Math.Max(8, (int)Math.Ceiling(u.Length() / step)));
                            var qv = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, Math.Max(8, (int)Math.Ceiling(v.Length() / step)));
                            for (int i = 0; i < qu.Order; i++) for (int j = 0; j < qv.Order; j++)
                                {
                                    integration.Add(origin + u * qu.GetAbscissa(i) + v * qv.GetAbscissa(j));
                                    normals.Add(normal);
                                    weights.Add(u.Length() * v.Length() * qu.GetWeight(i) * qv.GetWeight(j));
                                }
                        }
                        AddFace(new Point(-Width / 2, 0, -Height / 2), new Hare.Geometry.Vector(Width, 0, 0), new Hare.Geometry.Vector(0, 0, Height), new Hare.Geometry.Vector(0, 1, 0), nx, nz, 0.15 * Depth);
                        AddFace(new Point(-Width / 2, -Depth, -Height / 2), new Hare.Geometry.Vector(Width, 0, 0), new Hare.Geometry.Vector(0, 0, Height), new Hare.Geometry.Vector(0, -1, 0), nx, nz, 0.15 * Depth);
                        AddFace(new Point(-Width / 2, -Depth, -Height / 2), new Hare.Geometry.Vector(0, Depth, 0), new Hare.Geometry.Vector(0, 0, Height), new Hare.Geometry.Vector(-1, 0, 0), ny, nz, 0.15 * Width);
                        AddFace(new Point(Width / 2, -Depth, -Height / 2), new Hare.Geometry.Vector(0, Depth, 0), new Hare.Geometry.Vector(0, 0, Height), new Hare.Geometry.Vector(1, 0, 0), ny, nz, 0.15 * Width);
                        AddFace(new Point(-Width / 2, -Depth, -Height / 2), new Hare.Geometry.Vector(Width, 0, 0), new Hare.Geometry.Vector(0, Depth, 0), new Hare.Geometry.Vector(0, 0, -1), nx, ny, 0.15 * Height);
                        AddFace(new Point(-Width / 2, -Depth, Height / 2), new Hare.Geometry.Vector(Width, 0, 0), new Hare.Geometry.Vector(0, Depth, 0), new Hare.Geometry.Vector(0, 0, 1), nx, ny, 0.15 * Height);
                        //Test the prescribed normal derivative against conjugate fundamental solutions.
                        //The source-independent boundary matrix is reused for every driver in this cabinet.
                        var test = Matrix<Complex>.Build.Dense(sources.Count, integration.Count);
                        var gradient = Matrix<Complex>.Build.Dense(integration.Count, sources.Count);
                        System.Threading.Tasks.Parallel.For(0, integration.Count, i =>
                        {
                            for (int j = 0; j < sources.Count; j++)
                            {
                                double r = (integration[i] - sources[j]).Length();
                                test[j, i] = Complex.Conjugate(Complex.FromPolarCoordinates(1.0 / r, -k * r)) * weights[i];
                                gradient[i, j] = NormalDerivative(integration[i], sources[j], normals[i], k);
                            }
                        });
                        var matrix = test * gradient;
                        band = new Band { Sources = sources.ToArray(), Matrix = matrix, Factor = matrix.LU() };
                        Bands[octave] = band;
                    }
                }
                lock (band)
                {
                    if (band.LastField != null && band.LastField.Driver.x == driver.x && band.LastField.Driver.z == driver.z && band.LastField.Diameter == diameter) return band.LastField;
                    //Gauss-Legendre quadrature in squared radius; angular multiples of four preserve symmetry.
                    int radialCount = Math.Max(12, (int)Math.Ceiling(2 * k * radius));
                    int angularCount = 4 * Math.Max(12, (int)Math.Ceiling(Math.PI * k * radius));
                    Point[] piston = new Point[radialCount * angularCount];
                    double[] pistonWeights = new double[piston.Length];
                    var radialRule = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, radialCount);
                    for (int i = 0; i < radialCount; i++) for (int j = 0; j < angularCount; j++)
                        {
                            double r = radius * Math.Sqrt(radialRule.GetAbscissa(i)), angle = 2.0 * Math.PI * (j + 0.5) / angularCount;
                            piston[i * angularCount + j] = new Point(driver.x + r * Math.Cos(angle), 0, driver.z + r * Math.Sin(angle));
                            pistonWeights[i * angularCount + j] = radialRule.GetWeight(i) / angularCount;
                        }
                    var rhs = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(band.Sources.Length);
                    for (int i = 0; i < rhs.Count; i++)
                    {
                        //G = exp(-ikr)/r: unit forward infinite-baffle pressure gives piston flux -2*pi.
                        for (int j = 0; j < piston.Length; j++)
                        {
                            double r = (band.Sources[i] - piston[j]).Length();
                            rhs[i] -= 2 * Math.PI * pistonWeights[j] * Complex.Conjugate(Complex.FromPolarCoordinates(1.0 / r, -k * r));
                        }
                    }
                    var solution = band.Factor.Solve(rhs);
                    foreach (Complex value in solution) if (!Finite(value)) throw new InvalidOperationException("Numerical cabinet solve is singular or not finite.");
                    double residual = (band.Matrix * solution - rhs).L2Norm() / rhs.L2Norm();
                    if (!(residual <= 1E-7)) throw new InvalidOperationException("Numerical cabinet boundary solve failed its residual check.");
                    band.LastField = new Field { Driver = new Point(driver.x, 0, driver.z), Diameter = diameter, K = k, Sources = band.Sources, Strength = solution.ToArray() };
                    return band.LastField;
                }
            }

            public Complex Pressure(Point driver, Hare.Geometry.Vector direction, int octave, double driverDiameter, double distance)
            {
                double length = direction.Length();
                if (length <= 1E-12) return Complex.Zero;
                if (!(distance > 0) || double.IsInfinity(distance) || double.IsNaN(length) || double.IsInfinity(length)) throw new ArgumentOutOfRangeException("distance");
                ValidateDriver(driver, driverDiameter);
                Point receiver = driver + direction * (distance / length);
                if (Math.Abs(receiver.x) <= Width / 2 && receiver.y <= 0 && receiver.y >= -Depth && Math.Abs(receiver.z) <= Height / 2) throw new ArgumentException("Receiver must be outside the cabinet.");
                octave = Math.Max(0, Math.Min(7, octave));
                double weight = NumericalWeight(octave);
                if (weight == 0) return DED.Pressure(driver, direction, octave, driverDiameter, distance);
                Complex pressure = Prepare(driver, octave, driverDiameter).At(receiver);
                return weight == 1 ? pressure : weight * pressure + (1 - weight) * DED.Pressure(driver, direction, octave, driverDiameter, distance);
            }

            public string[] Driver_Balloon(Point driver, double driverDiameter)
            {
                return Cabinet_Diffraction.Balloon(Width, Height, Depth, SoundSpeed, (octave, distance) =>
                {
                    double weight = NumericalWeight(octave);
                    if (weight == 0) return DED.PrepareField(driver, octave, driverDiameter, distance);
                    Field field = Prepare(driver, octave, driverDiameter);
                    if (weight == 1) return direction => field.At(driver + direction * distance);
                    var edgeField = DED.PrepareField(driver, octave, driverDiameter, distance);
                    return direction => weight * field.At(driver + direction * distance) + (1 - weight) * edgeField(direction);
                });
            }
        }

        public class Cabinet_Diffraction : ICabinet_Diffraction
        {
            private struct EdgeElement
            {
                public Hare.Geometry.Point Point;
                public double Weight;
                public double SourceDistance;
                public double Length;
                public int Side;

                public EdgeElement(Hare.Geometry.Point point, double weight, double sourceDistance, double length, int side)
                {
                    Point = point;
                    Weight = weight;
                    SourceDistance = sourceDistance;
                    Length = length;
                    Side = side;
                }
            }

            private readonly double Width;
            private readonly double Height;
            private readonly double Depth;
            private readonly double SoundSpeed;

            internal static readonly double[] Frequencies = new double[]{ 62.5, 125, 250, 500, 1000, 2000, 4000, 8000 };

            /// <summary>
            /// Cabinet diffraction model based on the Distributed Edge Dipole
            /// formulation of Urban et al., following the distributed-edge
            /// interpretation introduced by Vanderkooy.
            ///
            /// Local coordinates:
            ///
            /// +Y = forward
            /// +Z = up
            /// +X = right
            ///
            /// The front baffle lies at Y = 0.
            ///
            /// The four front-baffle edges drive the corresponding rear edges and
            /// the four depth-running edges in a second distributed-edge integration.
            /// Depth edges retain separate dipoles for the two incident cabinet faces.
            /// Shared-face edge coupling continues through third and fourth order.
            /// Corner neighborhoods participate through the finite-edge integrals;
            /// a separate true point-vertex diffraction coefficient is not included.
            /// </summary>
            public Cabinet_Diffraction(double width, double height, double depth, double sound_speed = 343.0)
            {
                Width = width;
                Height = height;
                Depth = depth;
                SoundSpeed = sound_speed;
            }

            /// <summary>
            /// Signed circular piston aperture factor.
            ///
            /// Do not take Math.Abs() here. Above the first piston zero the sign
            /// represents a pi phase reversal and is important when the direct
            /// and diffracted fields are summed coherently.
            /// </summary>
            internal double PistonFactor(double frequency, double cosTheta, double diameter)
            {
                cosTheta = Math.Max(-1.0, Math.Min(1.0, cosTheta));

                double sinTheta = Math.Sqrt(Math.Max(0, 1.0 - cosTheta * cosTheta));
                double k = 2.0 * Math.PI * frequency / SoundSpeed;
                double x = k * diameter * 0.5 * sinTheta;

                if (Math.Abs(x) < 1E-10) return 1.0;
                return 2.0 * MathNet.Numerics.SpecialFunctions.BesselJ(1, x) / x;
            }

            private System.Numerics.Complex[] DrivenEdgeDrive(List<EdgeElement> edgeElements, List<EdgeElement> drivenElements, double frequency, double driverDiameter, int maxOrder = 4)
            {
                System.Numerics.Complex[] driven = new System.Numerics.Complex[drivenElements.Count];
                double k = 2.0 * Math.PI * frequency / SoundSpeed;
                double frontEdgeDrive = 0.5 * PistonFactor(frequency, 0, driverDiameter);

                for (int j = 0; j < drivenElements.Count; j++)
                {
                    EdgeElement element = drivenElements[j];
                    System.Numerics.Complex incident = System.Numerics.Complex.Zero;

                    //Rear sides 0..3 and depth dipoles 4..7 are driven only by their incident face.
                    //Each depth edge has two face-specific dipoles, rather than a point-corner source.
                    for (int i = 0; i < edgeElements.Count; i++)
                    {
                        EdgeElement frontElement = edgeElements[i];
                        if (frontElement.Side != element.Side % 4) continue;

                        Hare.Geometry.Vector path = element.Point - frontElement.Point;
                        double r = path.Length();
                        if (r <= 1E-9) continue;

                        //The first-order front dipole axis is +Y; preserve its polarity and phase.
                        double F = -path.dy / r;
                        System.Numerics.Complex frontStrength = frontEdgeDrive * frontElement.Weight * System.Numerics.Complex.FromPolarCoordinates(1, -k * frontElement.SourceDistance);
                        incident += frontStrength * F * System.Numerics.Complex.FromPolarCoordinates(1.0 / r, -k * r);
                    }

                    //Second DED integration uses physical dl, not the driver's projected dAlpha.
                    driven[j] = incident * element.Length / (2.0 * Math.PI);
                }

                //Keep each order separate: an order-n field drives order n+1, not the accumulated field.
                //Only distinct physical edges connected across an exterior cabinet face can couple.
                if (maxOrder > 2 && drivenElements.Count > 0)
                {
                    maxOrder = Math.Min(4, maxOrder);
                    List<KeyValuePair<int, System.Numerics.Complex>>[] transfer = new List<KeyValuePair<int, System.Numerics.Complex>>[driven.Length];
                    for (int j = 0; j < transfer.Length; j++) transfer[j] = new List<KeyValuePair<int, System.Numerics.Complex>>();
                    for (int i = 0; i < drivenElements.Count; i++)
                    {
                        EdgeElement source = drivenElements[i];
                        int outgoingSide = source.Side;
                        if (source.Side >= 4)
                        {
                            int xSide = source.Point.x >= 0 ? 1 : 3;
                            int zSide = source.Point.z >= 0 ? 0 : 2;
                            outgoingSide = source.Side % 4 == xSide ? zSide : xSide;
                        }

                        for (int j = 0; j < drivenElements.Count; j++)
                        {
                            EdgeElement target = drivenElements[j];
                            if (target.Side % 4 != outgoingSide) continue;
                            if (source.Side < 4 && target.Side < 4) continue;
                            if (source.Side >= 4 && target.Side >= 4 && source.Point.x == target.Point.x && source.Point.z == target.Point.z) continue;

                            //Midpoint quadrature at adjoining finite edges includes their corner neighborhoods.
                            //No arbitrary point-source strength or duplicate zero-length vertex path is added.
                            Hare.Geometry.Vector path = target.Point - source.Point;
                            double r = path.Length();
                            if (r <= 1E-9) continue;
                            double radiation = EdgeRadiation(source, path);
                            if (radiation == 0) continue;
                            transfer[j].Add(new KeyValuePair<int, System.Numerics.Complex>(i, radiation * target.Length / (2.0 * Math.PI) * System.Numerics.Complex.FromPolarCoordinates(1.0 / r, -k * r)));
                        }
                    }

                    System.Numerics.Complex[] previous = (System.Numerics.Complex[])driven.Clone();
                    for (int order = 3; order <= maxOrder; order++)
                    {
                        System.Numerics.Complex[] next = new System.Numerics.Complex[driven.Length];
                        for (int j = 0; j < next.Length; j++)
                        {
                            for (int i = 0; i < transfer[j].Count; i++) next[j] += transfer[j][i].Value * previous[transfer[j][i].Key];
                            driven[j] += next[j];
                        }
                        previous = next;
                    }
                }

                return driven;
            }
            /// <summary>
            /// Angle subtended by one elementary edge segment as viewed
            /// from the driver position in the baffle plane.
            ///
            /// Urban et al. use dAlpha / 2pi as the relative weighting of
            /// elementary edge sources for a rectangular baffle.
            /// </summary>
            private double SubtendedAngle(Hare.Geometry.Point driver, Hare.Geometry.Point a, Hare.Geometry.Point b)
            {
                double ax = a.x - driver.x;
                double az = a.z - driver.z;
                double bx = b.x - driver.x;
                double bz = b.z - driver.z;
                double la = Math.Sqrt(ax * ax + az * az);
                double lb = Math.Sqrt(bx * bx + bz * bz);

                if (la <= 1E-12 || lb <= 1E-12) return 0;

                double cross = ax * bz - az * bx;
                double dot = ax * bx + az * bz;

                return Math.Abs(Math.Atan2(cross, dot));
            }

            private void AddEdgeElements(List<EdgeElement> elements, Hare.Geometry.Point driver, Hare.Geometry.Point a, Hare.Geometry.Point b, double spacing, int side, bool driven = false)
            {
                Hare.Geometry.Vector span = b - a;
                double length = span.Length();

                if (length <= 1E-12) return;

                int count = Math.Max(1, (int)Math.Ceiling(length / spacing));
                double dl = length / count;

                for (int i = 0; i < count; i++)
                {
                    double t0 = (double)i / count;
                    double t1 = (double)(i + 1) / count;

                    Hare.Geometry.Point p0 = a + span * t0;
                    Hare.Geometry.Point p1 = a + span * t1;
                    Hare.Geometry.Point midpoint = (p0 + p1) / 2.0;

                    double dAlpha = driven ? 0 : SubtendedAngle(driver, p0, p1);

                    if (!driven && dAlpha <= 0) continue;
                    double sourceDistance = (midpoint - driver).Length();
                    if (sourceDistance <= 1E-12) continue;
                    elements.Add(new EdgeElement(midpoint, dAlpha / (2.0 * Math.PI), sourceDistance, dl, side));
                }
            }

            /// <summary>
            /// Constructs elementary dipoles along the four front-baffle edges.
            ///
            /// Edge spacing is frequency dependent. There is no need to use a
            /// single 8-kHz discretization for the lower octave bands.
            /// </summary>
            private List<EdgeElement> EdgeElements(Hare.Geometry.Point driver, int octave, bool driven = false)
            {
                octave = Math.Max(0, Math.Min( 7, octave));

                double wavelength = SoundSpeed / Frequencies[octave];

                //Urban correction to Vanderkooy method - specify elementary edge spacing much smaller
                //than wavelength.
                //Lambda / 8 is conservative enough for the coherent edge integration while the 20-mm cap keeps the low-frequencyrepresentation geometrically reasonable.
                
                double spacing = Math.Min(0.020, wavelength / 8.0);

                double x0 = -Width * 0.5;
                double x1 = Width * 0.5;
                double z0 = -Height * 0.5;
                double z1 = Height * 0.5;

                List<EdgeElement> elements = new List<EdgeElement>();

                // Top = 0
                AddEdgeElements(elements, driver, new Hare.Geometry.Point(x0, 0, z1), new Hare.Geometry.Point(x1, 0, z1), spacing, 0);
                // Right = 1
                AddEdgeElements(elements, driver, new Hare.Geometry.Point(x1, 0, z1), new Hare.Geometry.Point(x1, 0, z0), spacing, 1);
                // Bottom = 2
                AddEdgeElements(elements, driver, new Hare.Geometry.Point(x1, 0, z0), new Hare.Geometry.Point(x0, 0, z0), spacing, 2);
                // Left = 3
                AddEdgeElements(elements, driver, new Hare.Geometry.Point(x0, 0, z0), new Hare.Geometry.Point(x0, 0, z1), spacing, 3);

                double sum = 0;

                for (int i = 0; i < elements.Count; i++)
                {
                    sum += elements[i].Weight;
                }

                if (sum > 1E-12)
                {
                    for (int i = 0; i < elements.Count; i++)
                    {
                        EdgeElement e = elements[i];
                        e.Weight /= sum;
                        elements[i] = e;
                    }
                }

                if (driven)
                {
                    if (Depth <= 1E-9) return new List<EdgeElement>();

                    //Reuse the front discretization for the corresponding rear edges.
                    for (int i = 0; i < elements.Count; i++)
                    {
                        EdgeElement e = elements[i];
                        e.Point = new Hare.Geometry.Point(e.Point.x, -Depth, e.Point.z);
                        elements[i] = e;
                    }

                    Hare.Geometry.Point[] corners = new Hare.Geometry.Point[]{ new Hare.Geometry.Point(x1, 0, z1), new Hare.Geometry.Point(x1, 0, z0), new Hare.Geometry.Point(x0, 0, z0), new Hare.Geometry.Point(x0, 0, z1) };
                    for (int i = 0; i < corners.Length; i++)
                    {
                        Hare.Geometry.Point rear = new Hare.Geometry.Point(corners[i].x, -Depth, corners[i].z);
                        AddEdgeElements(elements, driver, corners[i], rear, spacing, 4 + i, true);
                        AddEdgeElements(elements, driver, corners[i], rear, spacing, 4 + (i + 1) % 4, true);
                    }
                }

                return elements;
            }

            //Shared by edge-to-edge propagation and the final receiver field.
            private double EdgeRadiation(EdgeElement element, Hare.Geometry.Vector direction)
            {
                double r = direction.Length();
                if (r <= 1E-12) return 0;
                double ux = direction.dx / r;
                double uy = direction.dy / r;
                double uz = direction.dz / r;
                double faceCosine;
                switch (element.Side % 4)
                {
                    case 0: faceCosine = uz; break;
                    case 1: faceCosine = ux; break;
                    case 2: faceCosine = -uz; break;
                    default: faceCosine = -ux; break;
                }

                //C1 face taper across direction cosines +/-0.1 (about +/-5.7 degrees).
                double Exterior(double cosine)
                {
                    double t = Math.Max(0, Math.Min(1, 0.5 + 5.0 * cosine));
                    return t * t * (3.0 - 2.0 * t);
                }

                if (element.Side < 4) return -uy * (1.0 - (1.0 - Exterior(-uy)) * (1.0 - Exterior(faceCosine)));
                double xExterior = Exterior(element.Point.x >= 0 ? ux : -ux);
                double zExterior = Exterior(element.Point.z >= 0 ? uz : -uz);
                return -faceCosine * (1.0 - (1.0 - xExterior) * (1.0 - zExterior));
            }
            private System.Numerics.Complex Pressure(Hare.Geometry.Point driver, Hare.Geometry.Vector direction, int octave, double driverDiameter, double distance, List<EdgeElement> edgeElements, List<EdgeElement> drivenElements, System.Numerics.Complex[] driven)
            {
                double length = direction.Length();

                if (length <= 1E-12) return System.Numerics.Complex.Zero;

                direction /= length;
                octave = Math.Max(0, Math.Min( 7, octave));
                double frequency = Frequencies[octave];
                double k = 2.0 * Math.PI * frequency / SoundSpeed;
                Hare.Geometry.Point receiver = driver + direction * distance;

                //FINITE - BAFFLE DRIVING FIELD
                //DED:
                //              1 + cos(theta)
                //  K(theta) = --------------
                //                  2
                //Unlike Vanderkooy's original hard front/shadow switch, this is continuous through 90 degrees.

                double cosTheta = Math.Max(-1.0, Math.Min(1.0, direction.dy));

                double driveFactor = 0.5 * (1.0 + cosTheta);
                double piston = PistonFactor(frequency, cosTheta, driverDiameter);
                System.Numerics.Complex pressure = piston * driveFactor / distance * System.Numerics.Complex.FromPolarCoordinates(1.0, -k * distance);
                double edgeDrive = 0.5 * PistonFactor(frequency, 0, driverDiameter);

                for (int i = 0; i < edgeElements.Count; i++)
                {
                    EdgeElement element = edgeElements[i];
                    Hare.Geometry.Vector toReceiver = receiver - element.Point;
                    double receiverDistance = toReceiver.Length();
                    if (receiverDistance <= 1E-12) continue;

                    double edgeCosTheta = toReceiver.dy / receiverDistance;
                    double edgeDirectivity = -edgeCosTheta;
                    double pathLength = element.SourceDistance + receiverDistance;
                    double amplitude = edgeDrive * edgeDirectivity * element.Weight / receiverDistance;
                    pressure += amplitude * System.Numerics.Complex.FromPolarCoordinates(1.0, -k * pathLength);
                }

                //SECOND THROUGH FOURTH ORDER DISTRIBUTED-EDGE DIFFRACTION
                for (int i = 0; i < drivenElements.Count; i++)
                {
                    EdgeElement element = drivenElements[i];
                    if (driven[i] == System.Numerics.Complex.Zero) continue;

                    Hare.Geometry.Vector toReceiver = receiver - element.Point;
                    double receiverDistance = toReceiver.Length();
                    if (receiverDistance <= 1E-12) continue;
                    pressure += driven[i] * EdgeRadiation(element, toReceiver) * System.Numerics.Complex.FromPolarCoordinates(1.0 / receiverDistance, -k * receiverDistance);
                }
                return pressure;
            }

            /// <summary>
            /// Returns the complex direct and distributed-edge cabinet field for one driver
            /// in one direction.
            ///
            /// This overload is useful for diagnostics and polar plotting.
            /// </summary>
            public System.Numerics.Complex Pressure(Hare.Geometry.Point driver, Hare.Geometry.Vector direction, int octave, double driverDiameter, double distance)
            {
                octave = Math.Max(0, Math.Min(7, octave));
                List<EdgeElement> edgeElements = EdgeElements(driver, octave);
                List<EdgeElement> drivenElements = EdgeElements(driver, octave, true);
                System.Numerics.Complex[] driven = DrivenEdgeDrive(edgeElements, drivenElements, Frequencies[octave], driverDiameter);

                return Pressure(driver, direction, octave, driverDiameter, distance, edgeElements, drivenElements, driven);
            }

            /// <summary>
            /// Generates Pachyderm full-sphere balloon strings for one driver
            /// mounted in this cabinet.
            ///
            /// The driver location is given in cabinet-local meters.
            /// </summary>
            public string[] Driver_Balloon(Hare.Geometry.Point driver, double driverDiameter)
            {
                return Balloon(Width, Height, Depth, SoundSpeed, (octave, distance) => PrepareField(driver, octave, driverDiameter, distance));
            }

            internal Func<Hare.Geometry.Vector, System.Numerics.Complex> PrepareField(Hare.Geometry.Point driver, int octave, double driverDiameter, double distance)
            {
                List<EdgeElement> edgeElements = EdgeElements(driver, octave);
                List<EdgeElement> drivenElements = EdgeElements(driver, octave, true);
                System.Numerics.Complex[] driven = DrivenEdgeDrive(edgeElements, drivenElements, Frequencies[octave], driverDiameter);
                return direction => Pressure(driver, direction, octave, driverDiameter, distance, edgeElements, drivenElements, driven);
            }
            internal static string[] Balloon(double width, double height, double depth, double sound_speed, Func<int, double, Func<Hare.Geometry.Vector, System.Numerics.Complex>> prepare)
            {
                const int umax = 37;
                const int vmax = 72;

                string[] result = new string[8];

                double largestDimension = Math.Max(width, Math.Max(height, depth));

                for (int octave = 0; octave < 8; octave++)
                {
                    double wavelength = sound_speed / Frequencies[octave];
                    double referenceDistance = Math.Max(1.0, Math.Max(10.0 * largestDimension, 2.0 * largestDimension * largestDimension / wavelength));

                    var field = prepare(octave, referenceDistance);
                    double[,] magnitude = new double[umax, vmax];

                    System.Threading.Tasks.Parallel.For(0, vmax, v =>
                        {
                            double phi = 2.0 * Math.PI * v / vmax + Math.PI / 2.0;

                            for (int u = 0; u < umax; u++)
                            {
                                double theta = Math.PI * u / (umax - 1);
                                Hare.Geometry.Vector direction = new Hare.Geometry.Vector(Math.Sin(theta) * Math.Cos(phi), Math.Cos(theta), Math.Sin(theta) * Math.Sin(phi));

                                System.Numerics.Complex p = field(direction);

                                double value = p.Magnitude;

                                if (double.IsNaN(value) || double.IsInfinity(value))
                                {
                                    value = 0;
                                }

                                magnitude[u, v] = value;
                            }
                        });

                    double max = 0;

                    for (int v = 0; v < vmax; v++)
                    {
                        for (int u = 0; u < umax; u++)
                        {
                            max = Math.Max(max, magnitude[u, v]);
                        }
                    }

                    if (max <= 1E-20)
                    {
                        max = 1;
                    }

                    System.Text.StringBuilder code = new System.Text.StringBuilder();

                    for (int v = 0; v < vmax; v++)
                    {
                        for (int u = 0; u < umax; u++)
                        {
                            double value = magnitude[u, v];
                            double attenuation;

                            if (value <= 1E-20)
                            {
                                attenuation = 60;
                            }
                            else
                            {
                                attenuation = 20.0 * Math.Log10(max / value);
                                attenuation = Math.Max(0, Math.Min(60, attenuation));
                            }

                            if (u > 0)
                            {
                                code.Append(" ");
                            }

                            code.Append(attenuation.ToString("0.000", System.Globalization.CultureInfo.InvariantCulture));
                        }

                        code.Append(";");
                    }

                    result[octave] = code.ToString();
                }

                return result;
            }
        }
    }
}
