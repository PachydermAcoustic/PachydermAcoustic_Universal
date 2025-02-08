using System;
using System.Threading;
using System.Collections.Generic;
using Pachyderm_Acoustic.Environment;
using MathNet.Numerics.LinearAlgebra;
using System.Numerics;
using Hare.Geometry;
using System.Linq;
using System.Xml.Linq;
using System.Threading.Tasks;

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
        private bool IsRunning;
        private AutoResetEvent SimulationResetEvent;
        public double[][] Results;//Frequency, Receiver
        private double[] frequency;
        public BoundaryElementSimulation_FreqDom(Polygon_Scene room, Source source, Receiver_Bank receivers, double[] freq)
        {
            Room = room;
            Source = source;
            Receivers = receivers;
            ProgressMessage = "Simulation not started.";
            Results = new double[freq.Length][];
            for(int i = 0; i < freq.Length; i++) Results[i] = new double[Receivers.Count];
            IsRunning = false;
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
            if (IsRunning)
            {
                ProgressMessage = "Simulation is already running.";
                return;
            }

            IsRunning = true;
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
                //List<BoundaryElement> elements = GenerateBoundaryElements(Room);
                // Create boundary elements based on the room geometry.

                for (int f = 0; f < frequency.Length; f++)
                {
                    List<BoundaryElement> elements = new List<BoundaryElement>();

                    // Calculate wavelength and maximum element size
                    double wavelength = Room.Sound_speed(0) / frequency[f];
                    double maxElementSize = wavelength / 6.0;

                    double k = 2 * Math.PI * frequency[f] / Room.Sound_speed(0);

                    for (int i = 0; i < Room.AbsorptionValue.Count; i++)
                    {
                        // Obtain polygon vertices
                        Hare.Geometry.Point[] vertices = Room.polygon(i);
                        //Obtain material from room...
                        Environment.Material m = Room.Surface_Material(i);
                        Complex R = m.Reflection_Narrow(frequency[f], new Hare.Geometry.Vector(0,0,1), new Hare.Geometry.Vector(0,0,1));
                        Complex Y = (1 - R) / (413.3 * (1 + R));

                        // Check the size of the polygon
                        double polygonSize = ComputePolygonSize(vertices);

                        if (polygonSize > maxElementSize)
                        {
                            // Subdivide the polygon
                            List<Hare.Geometry.Point[]> subdividedPolygons = SubdividePolygon(vertices, maxElementSize);

                            // Create boundary elements from subdivided polygons
                            foreach (var subVertices in subdividedPolygons)
                            {
                                BoundaryElement element = new BoundaryElement(subVertices, Y, elements.Count);
                                elements.Add(element);
                            }
                        }
                        else
                        {
                            // Create boundary element without subdivision
                            BoundaryElement element = new BoundaryElement(vertices, Y, elements.Count);
                            elements.Add(element);
                        }
                    }
                    ProgressMessage = $"Generated {elements.Count} boundary elements.";

                    // Step 2: Assemble the system of equations
                    ProgressMessage = "Assembling system of equations...";
                    int N = elements.Count;
                    double[,] matrix = new double[N, N];
                    double[] rhs = new double[N];

                    // Fill the matrix with the integral equations coefficients
                    for (int i = 0; i < N; i++)
                    {
                        Hare.Geometry.Point receiverPoint = elements[i].CollocationPoint;
                        double distance = (receiverPoint - Source.Origin).Length();
                        
                        if (distance == 0) rhs[i] = 0;
                        else rhs[i] = -((1 / (4 * Math.PI * distance)) * Complex.Exp(-1.0 * System.Numerics.Complex.ImaginaryOne * k * distance)).Real * elements[i].area;

                        for (int j = 0; j < N; j++)
                        {
                            BoundaryElement sourceElement = elements[j];
                            
                            if (i == j)
                            {
                                // Diagonal term (singular value)
                                matrix[i, j] = 0.5 * elements[i].area;
                            }
                            else
                            {
                                // Compute the Green's function between collocation point i and source element j
                                // Integrate over the source element to compute the Green's function effect
                                Hare.Geometry.Point sourcePoint = sourceElement.CollocationPoint;
                                double re = (receiverPoint - sourcePoint).Length();

                                double G = 0;
                                if (re != 0) G = ((1 / (4 * Math.PI * re)) * Complex.Exp(-1.0 * Complex.ImaginaryOne * k * re)).Real;
                                matrix[i, j] = G * sourceElement.area;
                            }
                        }
                    }

                    // Step 3: Solve the system of equations
                    ProgressMessage = "Solving system of equations...";
                    // Use a linear solver (e.g., LU decomposition, Gaussian elimination)
                    // For simplicity, we'll use a basic solver (not optimized for large matrices)
                    // Convert to MathNet matrices and vectors
                    Matrix<double> mathNetMatrix = Matrix<double>.Build.DenseOfArray(matrix);
                    MathNet.Numerics.LinearAlgebra.Vector<double> mathNetVector = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(rhs);

                    // Solve the linear system
                    double[] boundaryPressures = mathNetMatrix.Solve(mathNetVector).ToArray();

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
                IsRunning = false;
                SimulationResetEvent.Set();
            }
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
        /// <returns>List containing two polygons resulting from the split.</returns>
        private List<Hare.Geometry.Point[]> SplitPolygon(Hare.Geometry.Point[] vertices)
        {
            int N = vertices.Length;
            int maxEdgeIndex = 0;
            double maxEdgeLength = 0.0;

            // Find the longest edge
            for (int i = 0; i < N; i++)
            {
                Hare.Geometry.Point p1 = vertices[i];
                Hare.Geometry.Point p2 = vertices[(i + 1) % N];
                double edgeLength = (p1 - p2).Length();
                if (edgeLength > maxEdgeLength)
                {
                    maxEdgeLength = edgeLength;
                    maxEdgeIndex = i;
                }
            }

            // Compute midpoint of the longest edge
            Hare.Geometry.Point midPoint = new Hare.Geometry.Point(
                (vertices[maxEdgeIndex].x + vertices[(maxEdgeIndex + 1) % N].x) / 2.0,
                (vertices[maxEdgeIndex].y + vertices[(maxEdgeIndex + 1) % N].y) / 2.0,
                (vertices[maxEdgeIndex].z + vertices[(maxEdgeIndex + 1) % N].z) / 2.0
            );

            // Create two new polygons by adding the midpoint
            List<Hare.Geometry.Point[]> newPolygons = new List<Hare.Geometry.Point[]>();

            // First new polygon
            List<Hare.Geometry.Point> poly1 = new List<Hare.Geometry.Point>();
            for (int i = 0; i <= maxEdgeIndex; i++)
            {
                poly1.Add(vertices[i]);
            }
            poly1.Add(midPoint);

            // Second new polygon
            List<Hare.Geometry.Point> poly2 = new List<Hare.Geometry.Point>();
            poly2.Add(midPoint);
            for (int i = maxEdgeIndex + 1; i < N; i++)
            {
                poly2.Add(vertices[i]);
            }

            newPolygons.Add(poly1.ToArray());
            newPolygons.Add(poly2.ToArray());

            return newPolygons;
        }

        /// <summary>
        /// Computes pressures at receiver locations using the solved boundary pressures.
        /// </summary>
        /// <param name="elements"></param>
        /// <param name="boundaryPressures"></param>
        /// <param name="receivers"></param>
        private void ComputeReceiverPressures(List<BoundaryElement> elements, double[] boundaryPressures, Receiver_Bank receivers, int f)
        {
            double k = 2 * Math.PI * frequency[f] / Room.Sound_speed(0);

            // Iterate over receivers and compute pressures
            for (int r = 0; r < receivers.Count; r++)
            {
                Hare.Geometry.Point receiverPosition = receivers.Origin(r);
                double pressure = 0.0;

                for (int i = 0; i < elements.Count; i++)
                {
                    BoundaryElement element = elements[i];
                    // Integrate over the source element to compute the Green's function effect
                    Hare.Geometry.Point sourcePoint = element.CollocationPoint;
                    double re = (receivers.Origin(r) - sourcePoint).Length();
                    double G = re == 0 ? 0 : (1 / ((4 * Math.PI * re)) * Complex.Exp(-1.0 * Complex.ImaginaryOne * k * re)).Real;
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
            public BoundaryElement(Hare.Geometry.Point[] vertices, Complex Admittance, int id)
            {
                Vertices = vertices;
                ElementID = id;
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