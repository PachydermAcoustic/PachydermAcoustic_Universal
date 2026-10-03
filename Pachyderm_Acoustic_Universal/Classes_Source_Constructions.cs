using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pachyderm_Acoustic
{
    namespace Source_Constructions
    {
        public class Cabinet_Diffraction
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

            private readonly double[] Frequencies = new double[]{ 62.5, 125, 250, 500, 1000, 2000, 4000, 8000 };

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
            /// This implementation is first-order: the acoustic calculation uses
            /// the four front-baffle edges. Width, height and driver position are
            /// included directly. Depth is retained for the cabinet geometry and
            /// for later higher-order front-edge -> rear-edge diffraction.
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
            private double PistonFactor(double frequency, double cosTheta, double diameter)
            {
                cosTheta = Math.Max(-1.0, Math.Min(1.0, cosTheta));

                double sinTheta = Math.Sqrt(Math.Max(0, 1.0 - cosTheta * cosTheta));
                double k = 2.0 * Math.PI * frequency / SoundSpeed;
                double x = k * diameter * 0.5 * sinTheta;

                if (Math.Abs(x) < 1E-10) return 1.0;
                return 2.0 * MathNet.Numerics.SpecialFunctions.BesselJ(1, x) / x;
            }

            private System.Numerics.Complex[] RearEdgeDrive(List<EdgeElement> edgeElements, double frequency, double driverDiameter)
            {
                System.Numerics.Complex[] rearDrive = new System.Numerics.Complex[edgeElements.Count];
                if (Depth <= 1E-9) return rearDrive;

                double k = 2.0 * Math.PI * frequency / SoundSpeed;

                //Pressure that drives the first-order front - edge distribution.
                //theta = 90 degrees at the front edge.

                double frontEdgeDrive = 0.5 * PistonFactor(frequency, 0, driverDiameter);

                for (int j = 0; j < edgeElements.Count; j++)
                {
                    EdgeElement rearElement = edgeElements[j];
                    Hare.Geometry.Point rearPoint = new Hare.Geometry.Point(rearElement.Point.x, -Depth, rearElement.Point.z);
                    System.Numerics.Complex incident = System.Numerics.Complex.Zero;

                    //Couple only front and rear edges belonging to the same cabinet side.
                    //Top->top rear
                    //Right -> right rear
                    //Bottom->bottom rear
                    //Left -> left rear
                    //Cross - side paths would involve an additional corner event and should therefore not be included in this order.
                     
                    for (int i = 0; i < edgeElements.Count; i++)
                    {
                        EdgeElement frontElement = edgeElements[i];
                        if (frontElement.Side != rearElement.Side) continue;

                        Hare.Geometry.Vector path = rearPoint - frontElement.Point;
                        double r = path.Length();
                        if (r <= 1E-9) continue;

                        //DED edge dipole axis remains perpendicular to the original front baffle(+Y).
                        //Front edge->rear edge therefore lies in the rear hemisphere and has positive edge polarity.
                        double cosTheta = path.dy / r;
                        double F = -cosTheta;

                        //The source - independent part of the front edge element.
                        //Weight is dAlpha / 2pi, which is equivalent to dl/(2pi L) for a small element viewed from the driver.
                        System.Numerics.Complex frontStrength = frontEdgeDrive * frontElement.Weight * System.Numerics.Complex.FromPolarCoordinates(1, -k * frontElement.SourceDistance);

                        incident += frontStrength * F * System.Numerics.Complex.FromPolarCoordinates(1.0 / r, -k * r);
                    }

                    //Apply the second DED edge integration.
                    //Urban's elementary edge-source equation contains dl / (2 pi r)
                    //so the rear edge source strength contains its actual physical segment length here.

                    rearDrive[j] = incident * rearElement.Length / (2.0 * Math.PI);
                }

                return rearDrive;
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

            private void AddEdgeElements(List<EdgeElement> elements, Hare.Geometry.Point driver, Hare.Geometry.Point a, Hare.Geometry.Point b, double spacing, int side)
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

                    double dAlpha = SubtendedAngle(driver, p0, p1);

                    if (dAlpha <= 0) continue;
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
            private List<EdgeElement> EdgeElements(Hare.Geometry.Point driver, int octave)
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

                return elements;
            }

            private System.Numerics.Complex Pressure(Hare.Geometry.Point driver, Hare.Geometry.Vector direction, int octave, double driverDiameter, double distance, List<EdgeElement> edgeElements, System.Numerics.Complex[] rearDrive)
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

                //SECOND-ORDER FRONT EDGE -> REAR EDGE DIFFRACTION
                if (Depth > 1E-9 && rearDrive != null && rearDrive.Length == edgeElements.Count)
                {
                    for (int i = 0; i < edgeElements.Count; i++)
                    {
                        EdgeElement element = edgeElements[i];

                        if (rearDrive[i] == System.Numerics.Complex.Zero) continue;

                        Hare.Geometry.Point rearPoint = new Hare.Geometry.Point(element.Point.x, -Depth, element.Point.z);
                        Hare.Geometry.Vector toReceiver = receiver - rearPoint;
                        double receiverDistance = toReceiver.Length();

                        if (receiverDistance <= 1E-12) continue;

                        bool visible;
                        switch (element.Side)
                        {
                            // Top
                            case 0:
                                visible = toReceiver.dy <= 0 || toReceiver.dz >= 0;
                                break;

                            // Right
                            case 1:
                                visible = toReceiver.dy <= 0 || toReceiver.dx >= 0;
                                break;

                            // Bottom
                            case 2:
                                visible = toReceiver.dy <= 0 || toReceiver.dz <= 0;
                                break;
                     
                            // Left
                            case 3:
                                visible = toReceiver.dy <= 0 || toReceiver.dx <= 0;
                                break;

                            default:
                                visible = false;
                                break;
                        }

                        if (!visible)
                        {
                            continue;
                        }

                        double F = -toReceiver.dy / receiverDistance; ;
                        pressure += rearDrive[i] * F * System.Numerics.Complex.FromPolarCoordinates( 1.0 / receiverDistance, -k * receiverDistance);
                    }
                }

                return pressure;
            }

            /// <summary>
            /// Returns the complex first-order cabinet field for one driver
            /// in one direction.
            ///
            /// This overload is useful for diagnostics and polar plotting.
            /// </summary>
            public System.Numerics.Complex Pressure(Hare.Geometry.Point driver, Hare.Geometry.Vector direction, int octave, double driverDiameter, double distance)
            {
                List<EdgeElement> edgeElements = EdgeElements(driver, octave);
                System.Numerics.Complex[] rearDrive = RearEdgeDrive(edgeElements, Frequencies[octave], driverDiameter);

                return Pressure(driver, direction, octave, driverDiameter, distance, edgeElements, rearDrive);
}

            /// <summary>
            /// Generates Pachyderm full-sphere balloon strings for one driver
            /// mounted in this cabinet.
            ///
            /// The driver location is given in cabinet-local meters.
            /// </summary>
            public string[] Driver_Balloon(Hare.Geometry.Point driver, double driverDiameter)
            {
                const int umax = 37;
                const int vmax = 72;

                string[] result = new string[8];

                double largestDimension = Math.Max(Width, Math.Max(Height, Depth));

                for (int octave = 0; octave < 8; octave++)
                {
                    double wavelength = SoundSpeed / Frequencies[octave];
                    double referenceDistance = Math.Max(1.0, Math.Max(10.0 * largestDimension, 2.0 * largestDimension * largestDimension / wavelength));

                    List<EdgeElement> edgeElements = EdgeElements(driver, octave);
                    System.Numerics.Complex[] rearDrive = RearEdgeDrive(edgeElements, Frequencies[octave], driverDiameter);
                    double[,] magnitude = new double[umax, vmax];

                    System.Threading.Tasks.Parallel.For(0, vmax, v =>
                        {
                            double phi = 2.0 * Math.PI * v / vmax + Math.PI / 2.0;

                            for (int u = 0; u < umax; u++)
                            {
                                double theta = Math.PI * u / (umax - 1);
                                Hare.Geometry.Vector direction = new Hare.Geometry.Vector(Math.Sin(theta) * Math.Cos(phi), Math.Cos(theta), Math.Sin(theta) * Math.Sin(phi));

                                System.Numerics.Complex p = Pressure(driver, direction, octave, driverDiameter, referenceDistance, edgeElements, rearDrive);

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
