//using System;
//using System.Threading;
//using System.Collections.Generic;
//using Pachyderm_Acoustic.Environment;
//using Hare.Geometry;
//using MathNet.Numerics.LinearAlgebra;

//namespace Pachyderm_Acoustic.Simulation
//{
//    /// <summary>
//    /// A simulation class implementing the Time-Domain Boundary Element Method (TD-BEM) for acoustic analysis, including impedance boundary conditions.
//    /// </summary>
//    public class TimeDomainBoundaryElementSimulation : Simulation_Type
//    {
//        private Polygon_Scene Room;
//        private Source Source;
//        private Receiver_Bank Receivers;
//        private Thread SimulationThread;
//        private string ProgressMessage;
//        private List<double[]> Results; // Results at each time step
//        private bool IsRunning;
//        private AutoResetEvent SimulationResetEvent;
//        private double TotalTime;
//        private double TimeStep;
//        private int TimeSteps;
//        private double[][] Result;

//        public TimeDomainBoundaryElementSimulation(Polygon_Scene room, Source source, Receiver_Bank receivers, double totalTime, double timeStep)
//        {
//            Room = room;
//            Source = source;
//            Receivers = receivers;
//            TotalTime = totalTime;
//            TimeStep = timeStep;
//            TimeSteps = (int)(TotalTime / TimeStep);
//            ProgressMessage = "Simulation not started.";
//            Results = new List<double[]>();
//            IsRunning = false;
//            SimulationResetEvent = new AutoResetEvent(false);
//            Result = new double[Receivers.Count][];
//            for (int i = 0; i < Receivers.Count; i++) Result[i] = new double[TimeSteps];
//        }

//        /// <summary>
//        /// Returns the simulation type as a string.
//        /// </summary>
//        public override string Sim_Type()
//        {
//            return "Time-Domain Boundary Element Method Simulation with Impedance";
//        }

//        /// <summary>
//        /// Begins the simulation process.
//        /// </summary>
//        public override void Begin()
//        {
//            if (IsRunning)
//            {
//                ProgressMessage = "Simulation is already running.";
//                return;
//            }

//            IsRunning = true;
//            SimulationThread = new Thread(new ThreadStart(Simulate));
//            SimulationThread.Start();
//            ProgressMessage = "Simulation started.";
//        }

//        /// <summary>
//        /// Returns a progress message.
//        /// </summary>
//        public override string ProgressMsg()
//        {
//            return ProgressMessage;
//        }

//        /// <summary>
//        /// Returns the current state of the simulation thread.
//        /// </summary>
//        public override ThreadState ThreadState()
//        {
//            if (SimulationThread != null)
//                return SimulationThread.ThreadState;
//            else
//                return System.Threading.ThreadState.Unstarted;
//        }

//        /// <summary>
//        /// Combines results from different threads (if multithreading is used).
//        /// </summary>
//        public override void Combine_ThreadLocal_Results()
//        {
//            // Implementation for combining results from parallel computations if necessary
//        }

//        /// <summary>
//        /// The main simulation method where TD-BEM computations occur.
//        /// </summary>
//        private void Simulate()
//        {
//            double c = Room.Sound_speed(0);

//            try
//            {
//                ProgressMessage = "Initializing TD-BEM simulation...";

//                // Step 1: Define boundary elements and collocation points
//                //var elements = GenerateBoundaryElements(Room);
//                List<BoundaryElement> elements = new List<BoundaryElement>();

//                // Iterate over room polygons and create elements
//                for (int i = 0; i < Room.AbsorptionValue.Count; i++)
//                {
//                    Point[] poly = Room.polygon(i);
//                    Environment.Material material = Room.AbsorptionValue[i];

//                    // Compute admittance from absorption coefficient
//                    double absorptionCoefficient = material.Coefficient_A_Broad(0);
//                    double admittance = ComputeAdmittanceFromAbsorption(absorptionCoefficient);

//                    var element = new BoundaryElement(poly, i, admittance);
//                    elements.Add(element);
//                }
//                ProgressMessage = $"Generated {elements.Count} boundary elements.";

//                // Step 2: Initialize time-dependent variables
//                foreach (var element in elements) element.InitializeTimeVariables(TimeSteps);

//                // Step 3: Time-stepping loop
//                for (int timeIndex = 0; timeIndex < TimeSteps; timeIndex++)
//                {
//                    double currentTime = timeIndex * TimeStep;
//                    ProgressMessage = $"Time step {timeIndex + 1}/{TimeSteps} at time {currentTime:F4}s";

//                    // Step 3a: Assemble system at current time step
//                    //var system = AssembleSystem(elements, currentTime, timeIndex);
//                    int N = elements.Count;
//                    double[,] matrix = new double[N, N];
//                    double[] rhs = new double[N];

//                    // Fill the matrix and RHS using convolution integrals over past time steps
//                    for (int i = 0; i < N; i++)
//                    {
//                        var receiverElement = elements[i];
//                        for (int j = 0; j < N; j++)
//                        {
//                            var sourceElement = elements[j];

//                            // Compute the time-domain influence coefficient
//                            //double influence = ComputeTimeDomainInfluence(receiverElement, sourceElement, currentTime, timeIndex);
//                            double distance = (receiverElement.CollocationPoint - sourceElement.CollocationPoint).Length();
//                            double delay = distance / c;
//                            int delaySteps = (int)(delay / TimeStep);

//                            if (timeIndex >= delaySteps)
//                            {
//                                /// Time-domain Green's function for free-space propagation.
//                                // Dirac delta function represented as derivative of Heaviside function
//                                matrix[i, j] = (1.0 / (4 * Math.PI * distance)) * Math.Abs(delay - distance / c) < TimeStep ? 1.0 / TimeStep : 0.0;
//                                // Compute influence using the boundary condition
//                                matrix[i, j] *= receiverElement.Admittance * sourceElement.GetPressure(timeIndex - delaySteps);
//                            }
//                        }

//                        // Compute RHS from source and impedance boundary condition
//                        if (receiverElement.IsSourceElement(Source))
//                        {
//                            rhs[i] = SourceTimeFunction(currentTime) * receiverElement.Admittance;
//                        }
//                    }

//                    // Step 3b: Solve system
//                    // Convert to MathNet matrices and vectors
//                    Matrix<double> mathNetMatrix = Matrix<double>.Build.DenseOfArray(matrix);
//                    Vector<double> mathNetVector = Vector<double>.Build.Dense(rhs);

//                    // Solve the linear system
//                    double[] boundaryPressures = mathNetMatrix.Solve(mathNetVector).ToArray();

//                    // Step 3c: Store results
//                    Results.Add(boundaryPressures);

//                    // Step 3d: Update boundary elements with new pressures
//                    for (int i = 0; i < elements.Count; i++)
//                    {
//                        elements[i].RecordPressure(boundaryPressures[i]);//, timeIndex);
//                    }

//                    // Step 3e: Compute pressures at receiver locations
//                    for (int r = 0; r < Receivers.Count; r++)
//                    {
//                        Point receiverPosition = Receivers.Origin(r);
//                        double pressure = 0.0;

//                        for (int i = 0; i < elements.Count; i++)
//                        {
//                            double influence = 0.0;
//                            double distance = (receiverPosition - elements[i].CollocationPoint).Length();
//                            double tau = currentTime - distance / Room.Sound_speed(0);
//                            if (tau >= 0)
//                            {
//                                /// Time-domain Green's function for free-space propagation.
//                                // Dirac delta function represented as derivative of Heaviside function
//                                influence = (1.0 / (4 * Math.PI * distance)) * Math.Abs(tau - distance / Room.Sound_speed(0)) < TimeStep ? 1.0 / TimeStep : 0.0;
//                            }

//                            pressure += boundaryPressures[i] * influence;
//                        }

//                        // Store or process the pressure value as needed
//                        Result[r][(int)Math.Floor(currentTime / TimeStep)] = pressure;
//                    }

//                    // Step 3f: Apply attenuation using IIR filter
//                    // At the end of each time step, after solving for boundary pressures
//                    for (int i = 0; i < elements.Count; i++)
//                    {
//                        var element = elements[i];

//                        // Record the unfiltered pressure
//                        element.RecordPressure(boundaryPressures[i]);

//                        // Apply attenuation using the IIR filter
//                        double filteredPressure = element.ApplyAttenuation();

//                        // Update the boundary pressure with the filtered value
//                        boundaryPressures[i] = filteredPressure;

//                        // Optionally, if you need to update PastFilteredPressures for future computations
//                        element.RecordFilteredPressure(filteredPressure);
//                    }

//                }

//                ProgressMessage = "Simulation completed successfully.";
//            }
//            catch (ThreadAbortException)
//            {
//                ProgressMessage = "Simulation aborted.";
//            }
//            catch (Exception ex)
//            {
//                ProgressMessage = $"Simulation failed: {ex.Message}";
//            }
//            finally
//            {
//                IsRunning = false;
//                SimulationResetEvent.Set();
//            }
//        }

//        /// <summary>
//        /// Computes pressures at receiver locations using the solved boundary pressures.
//        /// </summary>
//        private void ComputeReceiverPressures(List<BoundaryElement> elements, double[] boundaryPressures, Receiver_Bank receivers, double currentTime)
//        {
//            for (int r = 0; r < Receivers.Count; r++)
//            {
//                Point receiverPosition = Receivers.Origin(r);
//                double pressure = 0.0;

//                for (int i = 0; i < elements.Count; i++)
//                {
//                    double influence = 0.0;
//                    double distance = (receiverPosition - elements[i].CollocationPoint).Length();
//                    double tau = currentTime - distance / Room.Sound_speed(0);
//                    if (tau >= 0)
//                    {
//                        /// Time-domain Green's function for free-space propagation.
//                        // Dirac delta function represented as derivative of Heaviside function
//                        influence = (1.0 / (4 * Math.PI * distance)) * Math.Abs(tau - distance / Room.Sound_speed(0)) < TimeStep ? 1.0 / TimeStep : 0.0;
//                    }

//                    pressure += boundaryPressures[i] * influence;
//                }

//                // Store or process the pressure value as needed
//                Result[r][(int)Math.Floor(currentTime / TimeStep)] = pressure;
//            }
//        }

//        /// <summary>
//        /// Time-dependent source function.
//        /// </summary>
//        private double SourceTimeFunction(double t)
//        {
//            // Example: Gaussian pulse
//            double amplitude = 1.0;
//            double pulseWidth = 0.01; // 10 ms pulse width
//            return amplitude * Math.Exp(-Math.Pow((t - pulseWidth) / (pulseWidth / 4), 2));
//        }

//        /// <summary>
//        /// Computes admittance from absorption coefficient.
//        /// </summary>
//        private double ComputeAdmittanceFromAbsorption(double absorptionCoefficient)
//        {
//            // Simple model: admittance Y = 2 * absorptionCoefficient
//            return 2.0 * absorptionCoefficient;
//        }
//    }

//    /// <summary>
//    /// Represents a boundary element in the TD-BEM simulation with IIR filtered absorption.
//    /// </summary>
//    /// 
//    public class BoundaryElement
//    {
//        public Point[] Vertices { get; private set; }
//        public Point CollocationPoint { get; private set; }
//        public int ElementID { get; private set; }
//        private double[] PastPressures;
//        private double[] PastFilteredPressures;
//        private int pressureIndex = 0;
//        public double Admittance { get; private set; }

//        private readonly double[] aCoefficients;
//        private readonly double[] bCoefficients;

//        public BoundaryElement(Point[] vertices, int id, double admittance, double[] bCoeffs, double[] aCoeffs)
//        {
//            Vertices = vertices;
//            ElementID = id;
//            Admittance = admittance;
//            CollocationPoint = ComputeCentroid(vertices);

//            aCoefficients = aCoeffs;
//            bCoefficients = bCoeffs;
//        }

//        public void InitializeTimeVariables(int filterOrder)
//        {
//            PastPressures = new double[filterOrder];
//            PastFilteredPressures = new double[filterOrder];
//        }

//        private Point ComputeCentroid(Point[] vertices)
//        {
//            double x = 0, y = 0, z = 0;
//            int N = vertices.Length;
//            foreach (var vertex in vertices)
//            {
//                x += vertex.x;
//                y += vertex.y;
//                z += vertex.z;
//            }
//            return new Point(x / N, y / N, z / N);
//        }

//        public void RecordPressure(double pressure)
//        {
//            PastPressures[pressureIndex] = pressure;
//            pressureIndex = (pressureIndex + 1) % PastPressures.Length;
//        }

//        public void RecordFilteredPressure(double pressure)
//        {
//            PastFilteredPressures[pressureIndex] = pressure;
//        }

//        public double ApplyAttenuation()
//        {
//            double output = 0.0;
//            int order = aCoefficients.Length;
//            int index;

//            // Numerator (b) coefficients
//            for (int i = 0; i < bCoefficients.Length; i++)
//            {
//                index = (pressureIndex - i + PastPressures.Length) % PastPressures.Length;
//                output += bCoefficients[i] * PastPressures[index];
//            }

//            // Denominator (a) coefficients (skip a[0], assumed to be 1)
//            for (int i = 1; i < aCoefficients.Length; i++)
//            {
//                index = (pressureIndex - i + PastPressures.Length) % PastPressures.Length;
//                output -= aCoefficients[i] * PastFilteredPressures[index];
//            }

//            // Store the filtered output
//            PastFilteredPressures[pressureIndex] = output;

//            return output;
//        }

//        public double GetPressure(int delaySteps)
//        {
//            int index = (pressureIndex - delaySteps + PastPressures.Length) % PastPressures.Length;
//            return PastFilteredPressures[index];
//        }
//    }

//    //public class BoundaryElement
//    //{
//    //    public Point[] Vertices { get; private set; }
//    //    public Point CollocationPoint { get; private set; }
//    //    public int ElementID { get; private set; }
//    //    private double[] PastPressures;
//    //    public double Admittance { get; private set; }

//    //    public BoundaryElement(Point[] vertices, int id, double admittance)
//    //    {
//    //        Vertices = vertices;
//    //        ElementID = id;
//    //        Admittance = admittance;
//    //        CollocationPoint = ComputeCentroid(vertices);
//    //    }

//    //    /// <summary>
//    //    /// Initializes time-dependent variables.
//    //    /// </summary>
//    //    public void InitializeTimeVariables(int timeSteps)
//    //    {
//    //        PastPressures = new double[timeSteps];
//    //    }

//    //    private Point ComputeCentroid(Point[] vertices)
//    //    {
//    //        double x = 0, y = 0, z = 0;
//    //        int N = vertices.Length;
//    //        foreach (var vertex in vertices)
//    //        {
//    //            x += vertex.x;
//    //            y += vertex.y;
//    //            z += vertex.z;
//    //        }
//    //        return new Point(x / N, y / N, z / N);
//    //    }

//    //    /// <summary>
//    //    /// Checks if this element is associated with the source.
//    //    /// </summary>
//    //    public bool IsSourceElement(Source source)
//    //    {
//    //        // Implement logic to check if this element is the source
//    //        return false;
//    //    }

//    //    /// <summary>
//    //    /// Records the pressure at a given time index.
//    //    /// </summary>
//    //    public void RecordPressure(double pressure, int timeIndex)
//    //    {
//    //        if (timeIndex < PastPressures.Length)
//    //        {
//    //            PastPressures[timeIndex] = pressure;
//    //        }
//    //    }

//    //    /// <summary>
//    //    /// Retrieves the pressure at a given time index.
//    //    /// </summary>
//    //    public double GetPressure(int timeIndex)
//    //    {
//    //        if (timeIndex < PastPressures.Length)
//    //        {
//    //            return PastPressures[timeIndex];
//    //        }
//    //        return 0.0;
//    //    }
//    //}
//}
