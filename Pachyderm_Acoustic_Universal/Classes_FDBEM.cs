using System;
using System.Threading;
using System.Collections.Generic;
using System.Numerics;
using Pachyderm_Acoustic.Environment;
using Hare.Geometry;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra.Complex;
using Vector = Hare.Geometry.Vector;

namespace Pachyderm_Acoustic.Simulation
{
    /// <summary>
    /// Constant-panel Helmholtz BEM with a surface-Laplacian admittance GIBC.
    /// Uses exp(i*omega*t), G = exp(-i*k*r)/(4*pi*r), and scene normals INTO the fluid.
    /// Absorbing velocity is opposite that normal: dp/dn = i*omega*rho*Y_Gamma*p.
    /// Geometry is in metres; Results are total pressures for a unit monopole Green field.
    /// </summary>
    public class BoundaryElementSimulation_FreqDom : Simulation_Type
    {
        private Polygon_Scene Room;
        private Source Source;
        private Receiver_Bank Receivers;
        private Thread SimulationThread;
        private string ProgressMessage;
        private AutoResetEvent SimulationResetEvent;
        private double[] frequency;
        private Dictionary<int, Complex>[] admittance;
        private static readonly GaussLegendreRule Quadrature = new GaussLegendreRule(0, 1, 8);
        public Complex[][] Results; // Frequency, receiver; NaN marks an unfinished/failed result.
        public Exception Failure { get; private set; }
        public double GIBC_Fit_Error { get; private set; }
        public int MaximumElements { get; set; } = 3000;

        public BoundaryElementSimulation_FreqDom(Scene room, Source source, Receiver_Bank receivers, double[] freq)
        {
            Room = room as Polygon_Scene ?? throw new ArgumentException("BEM requires a polygon scene.", nameof(room));
            Source = source ?? throw new ArgumentNullException(nameof(source));
            Receivers = receivers ?? throw new ArgumentNullException(nameof(receivers));
            frequency = (double[])(freq ?? throw new ArgumentNullException(nameof(freq))).Clone();
            foreach (double f in frequency) if (!(f > 0) || double.IsInfinity(f)) throw new ArgumentOutOfRangeException(nameof(freq), "BEM frequencies must be finite and positive.");
            Results = new Complex[frequency.Length][];
            for (int f = 0; f < frequency.Length; f++)
            {
                Results[f] = new Complex[Receivers.Count];
                for (int r = 0; r < Receivers.Count; r++) Results[f][r] = new Complex(double.NaN, double.NaN);
            }
            ProgressMessage = "Simulation not started.";
            SimulationResetEvent = new AutoResetEvent(false);
        }

        public override string Sim_Type() { return "Boundary Element Method (Frequency Domain) Simulation"; }
        public override string ProgressMsg() { return ProgressMessage; }
        public override ThreadState ThreadState() { return SimulationThread == null ? System.Threading.ThreadState.Unstarted : SimulationThread.ThreadState; }
        public override void Combine_ThreadLocal_Results() { }
        public override void Begin()
        {
            if (SimulationThread != null && SimulationThread.IsAlive) return;
            Failure = null;
            GIBC_Fit_Error = 0;
            for (int f = 0; f < Results.Length; f++) for (int r = 0; r < Results[f].Length; r++) Results[f][r] = new Complex(double.NaN, double.NaN);
            ProgressMessage = "Simulation started.";
            SimulationThread = new Thread(Simulate);
            SimulationThread.Start();
        }

        private double ComputePolygonSize(Point[] vertices)
        {
            double size = 0;
            for (int i = 0; i < vertices.Length; i++) size = Math.Max(size, (vertices[i] - vertices[(i + 1) % vertices.Length]).Length());
            return size;
        }

        private void Simulate()
        {
            try
            {
                double c = Room.Sound_speed(0), rho_c = Room.Rho_C(0);
                if (!(c > 0) || !(rho_c > 0) || double.IsInfinity(c) || double.IsInfinity(rho_c)) throw new InvalidOperationException("BEM requires a finite homogeneous fluid.");
                for (int f = 0; f < frequency.Length; f++)
                {
                    double k = Utilities.Numerics.PiX2 * frequency[f] / c;
                    Complex derivativeFactor = Complex.ImaginaryOne * k * rho_c;
                    double maxSize = c / (6 * frequency[f]), size = 0;
                    List<(Point[] Vertices, int Polygon)> triangles = new List<(Point[], int)>();
                    for (int p = 0; p < Room.Count(); p++)
                    {
                        Point[] v = Room.polygon(p);
                        if (v.Length < 3 || v.Length > 4) throw new NotSupportedException("BEM accepts triangles and planar convex quadrilaterals. Tessellate other faces first.");
                        Vector n = Room.Normal(p);
                        n = n / n.Length();
                        for (int a = 0; a < v.Length; a++)
                        {
                            if (Hare_math.Dot(Hare_math.Cross(v[(a + 1) % v.Length] - v[a], v[(a + 2) % v.Length] - v[(a + 1) % v.Length]), n) <= 0) throw new InvalidOperationException("BEM face is degenerate, concave, or has inconsistent winding.");
                            if (Math.Abs(Hare_math.Dot(v[a] - v[0], n)) > 1e-8 * Math.Max(1, ComputePolygonSize(v))) throw new InvalidOperationException("BEM quadrilaterals must be planar.");
                        }
                        for (int a = 1; a < v.Length - 1; a++)
                        {
                            Point[] triangle = new Point[] { v[0], v[a], v[a + 1] };
                            size = Math.Max(size, ComputePolygonSize(triangle));
                            triangles.Add((triangle, p));
                        }
                    }
                    // A common refinement depth preserves shared edges on a conforming input mesh.
                    int levels = 0;
                    double count = triangles.Count;
                    while (size > maxSize) { size *= .5; levels++; count *= 4; }
                    if (count > MaximumElements) throw new InvalidOperationException($"BEM needs {count} panels; the dense solver limit is {MaximumElements}. Reduce frequency or geometry size, or raise MaximumElements with sufficient memory.");
                    List<BoundaryElement> elements = new List<BoundaryElement>();
                    List<int> polygons = new List<int>();
                    foreach (var triangle in triangles)
                        foreach (Point[] v in SubdividePolygon(triangle.Vertices, levels))
                        {
                            elements.Add(new BoundaryElement(v, Complex.Zero, Room.Normal(triangle.Polygon), elements.Count));
                            polygons.Add(triangle.Polygon);
                        }
                    foreach (BoundaryElement element in elements) if (IsOnBoundary(Source.Origin, element)) throw new InvalidOperationException("BEM monopole sources must be off the boundary.");
                    int N = elements.Count;
                    ProgressMessage = $"Assembling {N} boundary elements at {frequency[f]} Hz...";
                    var fits = new Dictionary<Environment.Material, (Complex Y0, Complex Y2, double RelativeError)>();
                    admittance = new Dictionary<int, Complex>[N];
                    Complex[] slope = new Complex[N];
                    for (int i = 0; i < N; i++)
                    {
                        Environment.Material material = Room.Surface_Material(polygons[i]);
                        if (!fits.TryGetValue(material, out var fit))
                        {
                            fit = material.Surface_Admittance(frequency[f], rho_c);
                            fits.Add(material, fit);
                            GIBC_Fit_Error = Math.Max(GIBC_Fit_Error, fit.RelativeError);
                        }
                        admittance[i] = new Dictionary<int, Complex> { { i, fit.Y0 } };
                        slope[i] = fit.Y2 / (k * k);
                    }

                    // P1 cotangent stiffness K and lumped nodal mass M, projected onto panels:
                    // L = P M^-1 K M^-1 P^T A, where P averages the three triangle nodes.
                    // L approximates -Delta_Gamma, annihilates constants, and is A-self-adjoint.
                    // Object/plane keys impose zero tangential flux at planar patch edges.
                    // Curved scene objects share nodes across their tessellated facets.
                    // Exact coordinate keys require welded, conforming input; no proximity welding.
                    var nodeIDs = new Dictionary<(int Object, int Plane, double X, double Y, double Z), int>();
                    List<double> mass = new List<double>();
                    List<List<int>> incident = new List<List<int>>();
                    List<Dictionary<int, double>> stiffness = new List<Dictionary<int, double>>();
                    int[][] nodes = new int[N][];
                    void AddStiffness(int a, int b, double value)
                    {
                        stiffness[a].TryGetValue(b, out double previous);
                        stiffness[a][b] = previous + value;
                    }
                    for (int i = 0; i < N; i++)
                    {
                        nodes[i] = new int[3];
                        for (int a = 0; a < 3; a++)
                        {
                            Point v = elements[i].Vertices[a];
                            int objectID = Room.ObjectID(polygons[i]);
                            var key = (objectID, Room.IsPlanar(objectID) ? Room.PlaneID(polygons[i]) : -1, v.x, v.y, v.z);
                            if (!nodeIDs.TryGetValue(key, out int node))
                            {
                                node = mass.Count; nodeIDs.Add(key, node);
                                mass.Add(0); incident.Add(new List<int>()); stiffness.Add(new Dictionary<int, double>());
                            }
                            nodes[i][a] = node; mass[node] += elements[i].area / 3; incident[node].Add(i);
                        }
                        for (int a = 0; a < 3; a++)
                        {
                            int b = (a + 1) % 3, d = (a + 2) % 3;
                            double weight = Hare_math.Dot(elements[i].Vertices[b] - elements[i].Vertices[a], elements[i].Vertices[d] - elements[i].Vertices[a]) / (4 * elements[i].area);
                            int nb = nodes[i][b], nd = nodes[i][d];
                            AddStiffness(nb, nb, weight); AddStiffness(nd, nd, weight);
                            AddStiffness(nb, nd, -weight); AddStiffness(nd, nb, -weight);
                        }
                    }
                    for (int i = 0; i < N; i++)
                    {
                        if (slope[i] == Complex.Zero) continue;
                        foreach (int a in nodes[i]) foreach (var entry in stiffness[a]) foreach (int j in incident[entry.Key])
                        {
                            double laplacian = entry.Value * elements[j].area / (9 * mass[a] * mass[entry.Key]);
                            admittance[i].TryGetValue(j, out Complex previous);
                            admittance[i][j] = previous + slope[i] * laplacian;
                        }
                    }
                    // Interior trace with normals into the fluid: (1/2 I - D + S i*rho*omega*Y)p = p_inc.
                    Complex[,] matrix = new Complex[N, N];
                    Complex[] rhs = new Complex[N];
                    for (int i = 0; i < N; i++)
                    {
                        Point observation = elements[i].CollocationPoint;
                        rhs[i] = Green(observation, Source.Origin, k);
                        matrix[i, i] = .5;
                        for (int j = 0; j < N; j++)
                        {
                            var integral = IntegrateElement(observation, elements[j], k, i == j);
                            matrix[i, j] -= integral.D;
                            foreach (var entry in admittance[j]) matrix[i, entry.Key] += integral.G * derivativeFactor * entry.Value;
                        }
                    }
                    Complex[] pressure = N == 0 ? new Complex[0] : SolveMatrix(matrix, rhs);
                    ComputeReceiverPressures(elements, pressure, Receivers, f);
                }
                ProgressMessage = $"Simulation completed successfully. Maximum relative angular GIBC fit RMS: {GIBC_Fit_Error:P2}.";
            }
            catch (Exception ex)
            {
                Failure = ex;
                ProgressMessage = $"Simulation failed: {ex.Message}";
            }
            finally { SimulationResetEvent.Set(); }
        }

        private Complex[] SolveMatrix(Complex[,] matrix, Complex[] rhs)
        {
            // The Helmholtz/GIBC matrix is neither Hermitian nor positive definite; use pivoted LU.
            var A = DenseMatrix.OfArray(matrix);
            var b = DenseVector.OfArray(rhs);
            var solution = A.LU().Solve(b);
            double relativeResidual = (A * solution - b).L2Norm() / Math.Max(b.L2Norm(), 1e-30);
            if (double.IsNaN(relativeResidual) || double.IsInfinity(relativeResidual) || relativeResidual > 1e-8) throw new InvalidOperationException($"BEM solve failed its relative residual check ({relativeResidual:G3}).");
            return solution.ToArray();
        }

        private List<Point[]> SubdividePolygon(Point[] vertices, int levels)
        {
            List<Point[]> panels = new List<Point[]> { vertices };
            for (int level = 0; level < levels; level++)
            {
                List<Point[]> refined = new List<Point[]>();
                foreach (Point[] panel in panels) refined.AddRange(SplitPolygon(panel));
                panels = refined;
            }
            return panels;
        }

        private List<Point[]> SplitPolygon(Point[] v)
        {
            Point a = new Point((v[0].x + v[1].x) * .5, (v[0].y + v[1].y) * .5, (v[0].z + v[1].z) * .5), b = new Point((v[1].x + v[2].x) * .5, (v[1].y + v[2].y) * .5, (v[1].z + v[2].z) * .5), c = new Point((v[2].x + v[0].x) * .5, (v[2].y + v[0].y) * .5, (v[2].z + v[0].z) * .5);
            return new List<Point[]> { new Point[] { v[0], a, c }, new Point[] { a, v[1], b }, new Point[] { c, b, v[2] }, new Point[] { a, b, c } };
        }

        private bool IsOnBoundary(Point point, BoundaryElement element)
        {
            double tolerance = 1e-10 * Math.Max(1, ComputePolygonSize(element.Vertices));
            if (Math.Abs(Hare_math.Dot(point - element.Vertices[0], element.Normal)) > tolerance) return false;
            for (int a = 0; a < element.Vertices.Length; a++)
            {
                Vector edge = element.Vertices[(a + 1) % element.Vertices.Length] - element.Vertices[a];
                if (Hare_math.Dot(Hare_math.Cross(edge, point - element.Vertices[a]), element.Normal) < -tolerance * edge.Length()) return false;
            }
            return true;
        }
        private Complex Green(Point observation, Point source, double k)
        {
            double distance = (observation - source).Length();
            if (!(distance > 0)) throw new InvalidOperationException("BEM observation coincides with a monopole or quadrature point.");
            return Complex.Exp(-Complex.ImaginaryOne * k * distance) / (4 * Math.PI * distance);
        }

        private Complex ComputeNormalDerivative(Point observation, Point source, Vector normal, double k)
        {
            Vector delta = observation - source;
            double distance = delta.Length();
            return Green(observation, source, k) * (1 + Complex.ImaginaryOne * k * distance) * Hare_math.Dot(delta, normal) / (distance * distance);
        }

        private (Complex G, Complex D) IntegrateElement(Point observation, BoundaryElement element, double k, bool self)
        {
            Complex single = 0, dual = 0;
            // Duffy radial mapping about the panel centroid removes the 1/r self singularity.
            // This is regular Gauss quadrature for other panels; near-boundary receivers need refinement.
            Point center = element.CollocationPoint;
            for (int edge = 0; edge < element.Vertices.Length; edge++)
            {
                Vector a = element.Vertices[edge] - center, b = element.Vertices[(edge + 1) % element.Vertices.Length] - center;
                double jacobian = Hare_math.Cross(a, b).Length();
                for (int u = 0; u < Quadrature.Order; u++) for (int v = 0; v < Quadrature.Order; v++)
                {
                    double radial = Quadrature.GetAbscissa(u), transverse = Quadrature.GetAbscissa(v);
                    Point point = center + (a * (1 - transverse) + b * transverse) * radial;
                    double weight = jacobian * radial * Quadrature.GetWeight(u) * Quadrature.GetWeight(v);
                    single += weight * Green(observation, point, k);
                    if (!self) dual += weight * ComputeNormalDerivative(observation, point, element.Normal, k);
                }
            }
            return (single, dual);
        }

        private void ComputeReceiverPressures(List<BoundaryElement> elements, Complex[] pressure, Receiver_Bank receivers, int f)
        {
            double k = Utilities.Numerics.PiX2 * frequency[f] / Room.Sound_speed(0);
            Complex factor = Complex.ImaginaryOne * k * Room.Rho_C(0);
            Complex[] derivative = new Complex[elements.Count];
            for (int i = 0; i < elements.Count; i++) foreach (var entry in admittance[i]) derivative[i] += factor * entry.Value * pressure[entry.Key];
            for (int r = 0; r < receivers.Count; r++)
            {
                Point observation = receivers.Origin(r);
                Complex total = Green(observation, Source.Origin, k);
                for (int i = 0; i < elements.Count; i++)
                {
                    if (IsOnBoundary(observation, elements[i])) throw new InvalidOperationException("BEM receivers must be off the boundary.");
                    var integral = IntegrateElement(observation, elements[i], k, false);
                    total += integral.D * pressure[i] - integral.G * derivative[i];
                }
                Results[f][r] = total;
            }
        }

        public class BoundaryElement
        {
            public Point[] Vertices { get; private set; }
            public Point CollocationPoint { get; private set; }
            public int ElementID { get; private set; }
            public double area = 0;
            public Vector Normal;
            public Complex Admittance { get; private set; }
            public BoundaryElement(Point[] vertices, Complex Admittance, Vector Normal, int id)
            {
                if (vertices == null || vertices.Length < 3) throw new ArgumentException("A boundary element needs at least three vertices.", nameof(vertices));
                Vertices = vertices; ElementID = id; this.Admittance = Admittance;
                double length = Normal.Length();
                if (!(length > 0) || double.IsInfinity(length)) throw new ArgumentException("Boundary normal must be finite and nonzero.", nameof(Normal));
                this.Normal = Normal / length;
                double x = 0, y = 0, z = 0;
                foreach (Point v in vertices) { x += v.x; y += v.y; z += v.z; }
                CollocationPoint = new Point(x / vertices.Length, y / vertices.Length, z / vertices.Length);
                for (int a = 1; a < vertices.Length - 1; a++) area += .5 * Hare_math.Cross(vertices[a] - vertices[0], vertices[a + 1] - vertices[0]).Length();
                if (!(area > 0) || double.IsInfinity(area)) throw new ArgumentException("Boundary element has invalid area.", nameof(vertices));
            }
        }
    }
}