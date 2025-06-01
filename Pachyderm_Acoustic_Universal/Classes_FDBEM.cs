using System;
using System.Threading;
using System.Collections.Generic;
using Pachyderm_Acoustic.Environment;
using System.Numerics;
using System.Threading.Tasks;
using Eto.Forms;
using System.Configuration;

namespace Pachyderm_Acoustic.Simulation
{
    /// <summary>
    /// A simulation class implementing the Boundary Element Method (BEM) for acoustic analysis.
    /// </summary>
    public class BoundaryElementSimulation_FreqDom : Simulation_Type
    {
        private Polygon_Scene Room;
        private Source Source;
        private Receiver_Bank Receivers;
        private Thread SimulationThread;
        private string ProgressMessage;
        private AutoResetEvent SimulationResetEvent;
        public Complex[][] Results;//Frequency, Receiver
        private double[] frequency;
        public BoundaryElementSimulation_FreqDom(Scene room, Source source, Receiver_Bank receivers, double[] freq)
        {
            Room = (Polygon_Scene)room;
            Source = source;
            Receivers = receivers;
            ProgressMessage = "Simulation not started.";
            Results = new Complex[freq.Length][];
            for(int i = 0; i < freq.Length; i++) Results[i] = new Complex[Receivers.Count];
            SimulationResetEvent = new AutoResetEvent(false);
            frequency = freq;
        }

        /// <summary>
        /// Returns the simulation type as a string.
        /// </summary>
        /// <returns></returns>
        public override string Sim_Type()
        {
            return "Boundary Element Method (Frequency Domain) Simulation";
        }

        /// <summary>
        /// Begins the simulation process.
        /// </summary>
        public override void Begin()
        {
            if (SimulationThread != null && SimulationThread.ThreadState != System.Threading.ThreadState.Stopped)
            {
                ProgressMessage = "Simulation is already running.";
                return;
            }

            SimulationThread = new Thread(new ThreadStart(Simulate));
            SimulationThread.Start();
            ProgressMessage = "Simulation started.";
        }

        /// <summary>
        /// Returns a progress message.
        /// </summary>
        /// <returns></returns>
        public override string ProgressMsg()
        {
            return ProgressMessage;
        }

        /// <summary>
        /// Returns the current state of the simulation thread.
        /// </summary>
        /// <returns></returns>
        public override ThreadState ThreadState()
        {
            if (SimulationThread != null)
                return SimulationThread.ThreadState;
            else
                return System.Threading.ThreadState.Unstarted;
        }

        /// <summary>
        /// Combines results from different threads (if multithreading is used).
        /// </summary>
        public override void Combine_ThreadLocal_Results()
        {
        }

        /// <summary>
        /// Computes the maximum edge length of the polygon.
        /// </summary>
        /// <param name="vertices">Vertices of the polygon.</param>
        /// <returns>Maximum edge length.</returns>
        private double ComputePolygonSize(Hare.Geometry.Point[] vertices)
        {
            double maxLength = 0.0;
            int N = vertices.Length;
            for (int i = 0; i < N; i++)
            {
                Hare.Geometry.Point p1 = vertices[i];
                Hare.Geometry.Point p2 = vertices[(i + 1) % N]; // Wrap around to the first vertex
                double edgeLength = (p1 - p2).Length();
                if (edgeLength > maxLength)
                    maxLength = edgeLength;
            }
            return maxLength;
        }

        /// <summary>
        /// The main simulation method where BEM computations occur.
        /// </summary>
        private void Simulate()
        {
            try
            {
                ProgressMessage = "Initializing BEM simulation...";

                // Step 1: Define boundary elements and collocation points
                // Create boundary elements based on the room geometry.

                for (int f = 0; f < frequency.Length; f++)
                {
                    List<BoundaryElement> elements = new List<BoundaryElement>();

                    // Calculate wavelength and maximum element size
                    double wavelength = Room.Sound_speed(0) / frequency[f];
                    double maxElementSize = wavelength / 6.0;

                    double k = 2 * Math.PI * frequency[f] / Room.Sound_speed(0);
                    Complex i_k = Complex.ImaginaryOne / k;

                    for (int i = 0; i < Room.AbsorptionValue.Count; i++)
                    {
                        // Obtain polygon vertices
                        Hare.Geometry.Point[] vertices = Room.polygon(i);
                        Hare.Geometry.Vector normal = Room.Normal(i);
                        //Obtain material from room...
                        Environment.Material m = Room.Surface_Material(i);
                        Complex R = m.Reflection_Narrow(frequency[f], new Hare.Geometry.Vector(0, 0, 1), new Hare.Geometry.Vector(0, 0, 1));
                        Complex Y = (1 - R) / (Room.Rho_C(0) * (1 + R));

                        // Check the size of the polygon
                        double polygonSize = ComputePolygonSize(vertices);

                        if (polygonSize > maxElementSize)
                        {
                            // Subdivide the polygon
                            List<Hare.Geometry.Point[]> subdividedPolygons = SubdividePolygon(vertices, maxElementSize);

                            // Create boundary elements from subdivided polygons
                            foreach (var subVertices in subdividedPolygons)
                            {
                                BoundaryElement element = new BoundaryElement(subVertices, Y, normal, elements.Count);
                                elements.Add(element);
                            }
                        }
                        else
                        {
                            // Create boundary element without subdivision
                            BoundaryElement element = new BoundaryElement(vertices, Y, normal, elements.Count);
                            elements.Add(element);
                        }
                    }
                    ProgressMessage = $"Generated {elements.Count} boundary elements.";

                    // Step 2: Assemble the system of equations
                    ProgressMessage = "Assembling system of equations...";
                    int N = elements.Count;
                    Complex[,] matrix = new Complex[N, N];
                    Complex[] rhs = new Complex[N];

                    // Fill the matrix with the integral equations coefficients
                    for (int i = 0; i < N; i++)
                    {
                        Hare.Geometry.Point receiverPoint = elements[i].CollocationPoint;
                        Hare.Geometry.Vector Direct = (receiverPoint - Source.Origin);
                        double distance = Direct.Length();
                        Direct /= distance;

                        //Burton-Miller Correction
                        Complex BM_Derivative = ComputeNormalDerivative(receiverPoint, Source.Origin, elements[i].Normal, k);


                        // Check if the source element is on the same side as the receiver
                        double dot = Hare.Geometry.Hare_math.Dot(elements[i].Normal, Direct);
                        //Temporary check for the case of a sphere. More sophisticated check will be needed for more complex geometries

                        if (dot > 0 || distance == 0) rhs[i] = 0;
                        else rhs[i] = (-((1 / (4 * Math.PI * distance)) * Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * k * distance)) * elements[i].area) + i_k * BM_Derivative;


                        for (int j = 0; j < N; j++)
                        {
                            BoundaryElement sourceElement = elements[j];

                            // Check if the source element is on the same side as the receiver
                            double dotl = Hare.Geometry.Hare_math.Dot(sourceElement.CollocationPoint - receiverPoint, elements[j].Normal);
                            //Temporary check for the case of a sphere. More sophisticated check will be needed for more complex geometries

                            if (i == j)
                            {
                                // Diagonal term (singular value)
                                matrix[i, j] = (0.5 * elements[i].area + i_k * BM_Derivative) * elements[i].area;
                            }
                            else
                            {
                                // Compute the Green's function between collocation point i and source element j
                                // Integrate over the source element to compute the Green's function effect
                                Hare.Geometry.Point sourcePoint = sourceElement.CollocationPoint;
                                double re = (receiverPoint - sourcePoint).Length();

                                Complex G = 0;
                                if (re != 0) G = ((1 / (4 * Math.PI * re)) * Complex.Exp(-1.0 * Complex.ImaginaryOne * k * re));
                                matrix[i, j] = (G + i_k  * BM_Derivative) * sourceElement.area;
                            }

                            if (dotl < 0) matrix[i, j] = double.Epsilon;//*= -1;
                        }
                    }

                    // Step 3: Solve the system of equations
                    Complex[] boundaryPressures = SolveMatrixIteratively(matrix, rhs);

                    // Step 4: Compute pressures at receiver locations
                    ProgressMessage = "Computing pressures at receiver locations...";
                    ComputeReceiverPressures(elements, boundaryPressures, Receivers, f);
                }

                ProgressMessage = "Simulation completed successfully.";
            }
            catch (ThreadAbortException)
            {
                ProgressMessage = "Simulation aborted.";
            }
            catch (Exception ex)
            {
                ProgressMessage = $"Simulation failed: {ex.Message}";
            }
            finally
            {
                SimulationResetEvent.Set();
            }
        }

        private Complex[] SolveMatrixIteratively(Complex[,] matrix, Complex[] rhs)
        {
            int N = rhs.Length;
            Complex[] x = new Complex[N];
            Complex[] r = new Complex[N];
            Complex[] p = new Complex[N];
            Complex[] Ap = new Complex[N];

            // Initial guess x = 0
            Array.Copy(rhs, r, N);
            Array.Copy(r, p, N);

            //compute dot product of the vector r:
            Complex rsold = Complex.Zero;
            for (int i = 0; i < r.Length; i++) rsold += r[i] * r[i];

            for (int iter = 0; iter < N; iter++)
            {
                // Ap = A * p
                Parallel.For(0, N, i =>
                {
                    Ap[i] = Complex.Zero;
                    for (int j = 0; j < N; j++)
                    {
                        Ap[i] += matrix[i, j] * p[j];
                    }
                });

                Complex dp = Complex.Zero;
                for (int i = 0; i < p.Length; i++) dp += p[i] * Ap[i];
                Complex alpha = rsold / dp;

                Parallel.For(0, N, i =>
                {
                    x[i] += alpha * p[i];
                    r[i] -= alpha * Ap[i];
                });

                Complex rsnew = Complex.Zero;
                for (int i = 0; i < r.Length; i++) rsnew += r[i] * r[i];

                if (Math.Sqrt(rsnew.Real) < 1e-10)
                    break;

                Parallel.For(0, N, i =>
                {
                    p[i] = r[i] + (rsnew / rsold) * p[i];
                });

                rsold = rsnew;
            }

            return x;
        }

        /// <summary>
        /// Subdivides a polygon into smaller polygons based on the maximum element size.
        /// </summary>
        /// <param name="vertices">Vertices of the polygon.</param>
        /// <param name="maxSize">Maximum allowed size of an element.</param>
        /// <returns>List of subdivided polygons.</returns>
        private List<Hare.Geometry.Point[]> SubdividePolygon(Hare.Geometry.Point[] vertices, double maxSize)
        {
            List<Hare.Geometry.Point[]> subdividedPolygons = new List<Hare.Geometry.Point[]>();
            Queue<Hare.Geometry.Point[]> polygonsToSplit = new Queue<Hare.Geometry.Point[]>();
            polygonsToSplit.Enqueue(vertices);

            while (polygonsToSplit.Count > 0)
            {
                var currentPolygon = polygonsToSplit.Dequeue();

                if (ComputePolygonSize(currentPolygon) <= maxSize)
                {
                    subdividedPolygons.Add(currentPolygon);
                }
                else
                {
                    // Split the polygon
                    var splitPolygons = SplitPolygon(currentPolygon);
                    foreach (var poly in splitPolygons)
                    {
                        polygonsToSplit.Enqueue(poly);
                    }
                }
            }
            return subdividedPolygons;
        }

        /// <summary>
        /// Splits a polygon into smaller polygons by dividing along the longest edge.
        /// </summary>
        /// <param name="vertices">Vertices of the polygon to split.</param>
        /// <returns>List containing smaller polygons resulting from the split.</returns>
        private List<Hare.Geometry.Point[]> SplitPolygon(Hare.Geometry.Point[] vertices)
        {
            int N = vertices.Length;
            List<Hare.Geometry.Point[]> newPolygons = new List<Hare.Geometry.Point[]>();

            if (N == 3)
            {
                // Calculate midpoints of each edge
                Hare.Geometry.Point midPoint1 = new Hare.Geometry.Point(
                    (vertices[0].x + vertices[1].x) / 2.0,
                    (vertices[0].y + vertices[1].y) / 2.0,
                    (vertices[0].z + vertices[1].z) / 2.0
                );

                Hare.Geometry.Point midPoint2 = new Hare.Geometry.Point(
                    (vertices[1].x + vertices[2].x) / 2.0,
                    (vertices[1].y + vertices[2].y) / 2.0,
                    (vertices[1].z + vertices[2].z) / 2.0
                );

                Hare.Geometry.Point midPoint3 = new Hare.Geometry.Point(
                    (vertices[2].x + vertices[0].x) / 2.0,
                    (vertices[2].y + vertices[0].y) / 2.0,
                    (vertices[2].z + vertices[0].z) / 2.0
                );

                // Create four smaller triangles
                newPolygons.Add(new Hare.Geometry.Point[] { vertices[0], midPoint1, midPoint3 });
                newPolygons.Add(new Hare.Geometry.Point[] { midPoint1, vertices[1], midPoint2 });
                newPolygons.Add(new Hare.Geometry.Point[] { midPoint3, midPoint2, vertices[2] });
                newPolygons.Add(new Hare.Geometry.Point[] { midPoint1, midPoint2, midPoint3 });
            }
            else if (N == 4)
            {
                // Handle quadrilateral
                Hare.Geometry.Point midPoint1 = new Hare.Geometry.Point(
                    (vertices[0].x + vertices[1].x) / 2.0,
                    (vertices[0].y + vertices[1].y) / 2.0,
                    (vertices[0].z + vertices[1].z) / 2.0
                );

                Hare.Geometry.Point midPoint2 = new Hare.Geometry.Point(
                    (vertices[1].x + vertices[2].x) / 2.0,
                    (vertices[1].y + vertices[2].y) / 2.0,
                    (vertices[1].z + vertices[2].z) / 2.0
                );

                Hare.Geometry.Point midPoint3 = new Hare.Geometry.Point(
                    (vertices[2].x + vertices[3].x) / 2.0,
                    (vertices[2].y + vertices[3].y) / 2.0,
                    (vertices[2].z + vertices[3].z) / 2.0
                );

                Hare.Geometry.Point midPoint4 = new Hare.Geometry.Point(
                    (vertices[3].x + vertices[0].x) / 2.0,
                    (vertices[3].y + vertices[0].y) / 2.0,
                    (vertices[3].z + vertices[0].z) / 2.0
                );

                // Calculate the centroid of the quadrilateral
                Hare.Geometry.Point centroid = new Hare.Geometry.Point(
                    (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x) / 4.0,
                    (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y) / 4.0,
                    (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z) / 4.0
                );

                // Create four smaller quadrilaterals or triangles
                newPolygons.Add(new Hare.Geometry.Point[] { vertices[0], midPoint1, centroid, midPoint4 });
                newPolygons.Add(new Hare.Geometry.Point[] { midPoint1, vertices[1], midPoint2, centroid });
                newPolygons.Add(new Hare.Geometry.Point[] { centroid, midPoint2, vertices[2], midPoint3 });
                newPolygons.Add(new Hare.Geometry.Point[] { midPoint4, centroid, midPoint3, vertices[3] });
            }

            return newPolygons;
        }

        private Complex ComputeNormalDerivative(Hare.Geometry.Point observationPoint, Hare.Geometry.Point sourcePoint, Hare.Geometry.Vector boundaryNormal, double k)
        {
            // Vector from source to observation
            double Rx = observationPoint.x - sourcePoint.x;
            double Ry = observationPoint.y - sourcePoint.y;
            double Rz = observationPoint.z - sourcePoint.z;

            double re = Math.Sqrt(Rx * Rx + Ry * Ry + Rz * Rz);
            if (re < 1e-12) return Complex.Zero;

            // Dot product with boundary normal
            double dotRn = (Rx * boundaryNormal.dx + Ry * boundaryNormal.dy + Rz * boundaryNormal.dz) / re;

            // e^(i k r)/(4π r^3) * (i k r - 1) * [R ⋅ n]
            Complex expTerm = Complex.Exp(-Complex.ImaginaryOne * k * re);
            Complex factor = (expTerm / (4.0 * Math.PI * re * re * re)) * (Complex.ImaginaryOne * k * re - 1.0);
            return factor * dotRn;
        }

        /// <summary>
        /// Computes pressures at receiver locations using the solved boundary pressures.
        /// </summary>
        /// <param name="elements"></param>
        /// <param name="boundaryPressures"></param>
        /// <param name="receivers"></param>
        private void ComputeReceiverPressures(List<BoundaryElement> elements, Complex[] boundaryPressures, Receiver_Bank receivers, int f)
        {
            double k = 2 * Math.PI * frequency[f] / Room.Sound_speed(0);

            // Iterate over receivers and compute pressures
            for (int r = 0; r < receivers.Count; r++)
            {
                Hare.Geometry.Point receiverPosition = receivers.Origin(r);
                Complex pressure = 0.0;

                for (int i = 0; i < elements.Count; i++)
                {
                    double dot = Hare.Geometry.Hare_math.Dot(receiverPosition - elements[i].CollocationPoint, elements[i].Normal);
                    if (dot < 0) continue;
                    BoundaryElement element = elements[i];
                    // Integrate over the source element to compute the Green's function effect
                    Hare.Geometry.Point sourcePoint = element.CollocationPoint;
                    double re = (receivers.Origin(r) - sourcePoint).Length();
                    Complex G = re == 0 ? 0 : (1 / ((4 * Math.PI * re)) * Complex.Exp(-1.0 * Complex.ImaginaryOne * k * re));
                    pressure += boundaryPressures[i] * G;
                }

                // Store or process the pressure value as needed
                Results[f][r] = pressure;
            }
        }

        /// <summary>
        /// Represents a boundary element in the BEM simulation.
        /// </summary>
        public class BoundaryElement
        {
            public Hare.Geometry.Point[] Vertices { get; private set; }
            public Hare.Geometry.Point CollocationPoint { get; private set; }
            public int ElementID { get; private set; }
            public double area = 0;
            public Hare.Geometry.Vector Normal;
            public BoundaryElement(Hare.Geometry.Point[] vertices, Complex Admittance, Hare.Geometry.Vector Normal, int id)
            {
                Vertices = vertices;
                ElementID = id;
                this.Normal = Normal;
                CollocationPoint = ComputeCentroid(vertices);
                for (int i = 1, j = 2, k = 0; j < Vertices.Length; i++, j++)
                {
                    area += .5 * Hare.Geometry.Hare_math.Cross(Vertices[i] - Vertices[k], Vertices[i] - Vertices[j]).Length();
                }
            }

            private Hare.Geometry.Point ComputeCentroid(Hare.Geometry.Point[] vertices)
            {
                double x = 0, y = 0, z = 0;
                int N = vertices.Length;
                foreach (var vertex in vertices)
                {
                    x += vertex.x;
                    y += vertex.y;
                    z += vertex.z;
                }
                return new Hare.Geometry.Point(x / N, y / N, z / N);
            }
        }
    }
}