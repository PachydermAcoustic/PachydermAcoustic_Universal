using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Numerics;
using Vector = Hare.Geometry.Vector;

//General spectral contour with explicit Abel-Poisson damping, Bonner (2003) 6.1.1.
//Finite damping is an experimental approximation; uniform physical matching remains unfinished.
namespace Pachyderm_Acoustic.Source_Constructions
{
    internal sealed partial class Trihedral_Cone
    {
        static readonly MathNet.Numerics.Integration.GaussLegendreRule AngularRule = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, 64);

        readonly Vector[] Nodes, Normals;
        readonly Vector[][] Quadrature;
        readonly double[][] Weights;

        readonly bool Neumann;

        internal Trihedral_Cone(int panelsPerArc, bool neumann = true)
        {
            if (panelsPerArc < 4 || panelsPerArc > 128 || panelsPerArc % 2 != 0) throw new ArgumentOutOfRangeException(nameof(panelsPerArc));
            Neumann = neumann;
            int count = 3 * panelsPerArc;
            Nodes = new Vector[count]; Normals = new Vector[count];
            Quadrature = new Vector[count][]; Weights = new double[count][];
            var axes = new[] { new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(0, 0, 1) };
            //Cubic grading resolves the reentrant spherical corners of the Neumann domain.
            double Graded(double t) => Math.PI / 2 * (t <= 0.5 ? 0.5 * Math.Pow(2 * t, 3) : 1 - 0.5 * Math.Pow(2 * (1 - t), 3));
            for (int arc = 0; arc < 3; arc++) for (int j = 0; j < panelsPerArc; j++)
                {
                    int p = arc * panelsPerArc + j;
                    double a = Graded((double)j / panelsPerArc), b = Graded((double)(j + 1) / panelsPerArc);
                    Vector first = axes[arc], second = axes[(arc + 1) % 3];
                    Vector At(double t) => first * Math.Cos(t) + second * Math.Sin(t);
                    Nodes[p] = At((a + b) / 2); Normals[p] = axes[(arc + 2) % 3];

                    var rule = new MathNet.Numerics.Integration.GaussLegendreRule(a, b, 8);
                    Quadrature[p] = new Vector[rule.Order]; Weights[p] = new double[rule.Order];
                    for (int q = 0; q < rule.Order; q++) { Quadrature[p][q] = At(rule.GetAbscissa(q)); Weights[p][q] = rule.GetWeight(q); }
                }
        }

        internal static double Dot(Vector a, Vector b) => Hare.Geometry.Hare_math.Dot(a, b);

        internal static Vector Exterior(Vector direction)
        {
            double length = direction.Length();
            if (!(length > 0) || double.IsInfinity(length)) throw new ArgumentException("A finite nonzero direction is required.");
            Vector unit = direction / length;
            if (unit.dx >= 0 && unit.dy >= 0 && unit.dz >= 0) throw new ArgumentException("Direction is in the solid trihedral cone or on its boundary.");
            return unit;
        }

        //Scaled conical Legendre kernel and derivative with respect to dot(omega,omega').
        //The zero-balanced hypergeometric expansion near coincidence avoids a billion-term
        //ordinary series. Cos(pi*i*tau) and the gamma product cancel analytically there.
        internal static (double Value, double Derivative) Green(double dot, double tau)
        {
            if (!(tau >= 0) || tau > 24 || !(dot >= -1 - 1E-14 && dot < 1)) throw new ArgumentOutOfRangeException();
            double z = (1 + Math.Max(-1, dot)) / 2;
            if (z <= 0.995)
            {
                double term = 1, sum = 1, derivativeTerm = 0.25 + tau * tau, derivative = derivativeTerm;
                if (z == 0) return (-1 / (4 * Math.Cosh(Math.PI * tau)), -derivative / (8 * Math.Cosh(Math.PI * tau)));
                for (int n = 1; n < 16384; n++)
                {
                    term *= ((n - 0.5) * (n - 0.5) + tau * tau) * z / (n * n);
                    sum += term;
                    if (n > 1)
                    {
                        derivativeTerm *= ((n - 0.5) * (n - 0.5) + tau * tau) * z / (n * (n - 1));
                        derivative += derivativeTerm;
                    }
                    if (n > tau + 4 && Math.Abs(term) < 2E-15 * Math.Abs(sum) && Math.Abs(derivativeTerm) < 2E-15 * Math.Abs(derivative))
                        return (-sum / (4 * Math.Cosh(Math.PI * tau)), -derivative / (8 * Math.Cosh(Math.PI * tau)));
                }
                throw new InvalidOperationException("Conical Legendre series did not converge.");
            }
            //Complex digamma: recurrence into its asymptotic region, then Bernoulli series.
            //MathNet supplies only the real digamma; no complex helper exists in Universal.
            Complex w = new Complex(0.5, tau), correction = Complex.Zero;
            while (w.Real < 16) { correction -= 1 / w; w += 1; }
            Complex inv = 1 / w, inv2 = inv * inv;
            Complex psi = correction + Complex.Log(w) - inv / 2 - inv2 * (1.0 / 12 - inv2 * (1.0 / 120 - inv2 * (1.0 / 252 - inv2 * (1.0 / 240 - inv2 * (1.0 / 132)))));
            double h = -2 * 0.5772156649015328606 - 2 * psi.Real, qz = 1 - z, log = Math.Log(qz);
            double coefficient = 1, value = h - log, dz = 1 / qz;
            for (int n = 1; n < 16384; n++)
            {
                h += 2.0 / n - 2 * (n - 0.5) / ((n - 0.5) * (n - 0.5) + tau * tau);
                coefficient *= ((n - 0.5) * (n - 0.5) + tau * tau) * qz / (n * n);
                double term = coefficient * (h - log), dterm = coefficient / qz * (1 - n * (h - log));
                value += term; dz += dterm;
                if (n > tau + 4 && Math.Abs(term) < 2E-15 * Math.Max(1E-100, Math.Abs(value)) && Math.Abs(dterm) < 2E-15 * Math.Max(1E-100, Math.Abs(dz)))
                    return (-value / (4 * Math.PI), -dz / (8 * Math.PI));
            }
            throw new InvalidOperationException("Conical logarithmic kernel series did not converge.");
        }

        Matrix<double> BoundaryMatrix(double tau)
        {
            var matrix = Matrix<double>.Build.Dense(Nodes.Length, Nodes.Length);
            for (int i = 0; i < Nodes.Length; i++) for (int j = 0; j < Nodes.Length; j++)
                {
                    //Every arc is a great circle: its same-arc double-layer derivative is zero,
                    //including the principal-value diagonal. The jump term is exactly 1/2.
                    if (i / (Nodes.Length / 3) == j / (Nodes.Length / 3)) { if (i == j) matrix[i, j] = 0.5; continue; }
                    for (int q = 0; q < Quadrature[j].Length; q++)
                    {
                        Vector source = Quadrature[j][q];
                        double derivative = Green(Dot(Nodes[i], source), tau).Derivative;
                        matrix[i, j] += Weights[j][q] * derivative * (Neumann ? -Dot(Normals[i], source) : Dot(Nodes[i], Normals[j]));
                    }
                }
            return matrix;
        }

        double Potential(Vector observer, double tau, double[] density)
        {
            double value = 0;
            for (int j = 0; j < Nodes.Length; j++) for (int q = 0; q < Quadrature[j].Length; q++)
                {
                    Vector source = Quadrature[j][q];
                    var g = Green(Dot(observer, source), tau);
                    value += Weights[j][q] * density[j] * (Neumann ? g.Value : g.Derivative * Dot(observer, Normals[j]));
                }
            return value;
        }

        internal double Spectral(Vector observer, Vector incident, double tau)
        {
            observer = Exterior(observer); incident = Exterior(incident);
            var rhs = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(Nodes.Length);
            for (int i = 0; i < rhs.Count; i++)
            {
                var g = Green(Dot(Nodes[i], incident), tau);
                rhs[i] = Neumann ? g.Derivative * Dot(Normals[i], incident) : -g.Value;
            }
            var matrix = BoundaryMatrix(tau);
            var density = matrix.LU().Solve(rhs);
            double residual = (matrix * density - rhs).L2Norm() / rhs.L2Norm();
            if (!(residual < 1E-10)) throw new InvalidOperationException("Spherical boundary solve failed its residual check.");
            return Potential(observer, tau, density.ToArray());
        }

        //Check the boundary formulation against an independent, exactly harmonic field
        //whose pole is inside the excluded cone; neither the input nor target uses this BIE.
        internal double ManufacturedError(double tau, Vector observer)
        {
            observer = Exterior(observer);
            Vector pole = new Vector(1, 1, 1) / Math.Sqrt(3);
            var rhs = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(Nodes.Length);
            for (int i = 0; i < rhs.Count; i++)
            {
                var g = Green(Dot(Nodes[i], pole), tau);
                rhs[i] = Neumann ? -g.Derivative * Dot(Normals[i], pole) : g.Value;
            }
            double computed = Potential(observer, tau, BoundaryMatrix(tau).LU().Solve(rhs).ToArray());
            double exact = Green(Dot(observer, pole), tau).Value;
            return Math.Abs((computed - exact) / exact);
        }

        //Conservative M1 certificate: for every boundary point, both distances are at least
        //their respective minimum distances to the contour. It may reject valid M1 pairs.
        //This certificate must never be interpreted as a full-sphere visibility rule.
        internal static double M1Margin(Vector observer, Vector incident)
        {
            double Distance(Vector w)
            {
                w = Exterior(w);
                double best = -1;
                double[] c = { w.dx, w.dy, w.dz };
                for (int arc = 0; arc < 3; arc++)
                {
                    double a = c[arc], b = c[(arc + 1) % 3];
                    best = Math.Max(best, Math.Max(a, b));
                    if (a > 0 && b > 0) best = Math.Max(best, Math.Sqrt(a * a + b * b));
                }
                return Math.Acos(Math.Max(-1, Math.Min(1, best)));
            }
            return Distance(observer) + Distance(incident) - Math.PI;
        }

        internal Complex Coefficient(Vector observer, Vector incident, double cutoff = 12, int spectralOrder = 48)
        {
            observer = Exterior(observer); incident = Exterior(incident);
            if (!(M1Margin(observer, incident) > 0.05)) throw new ArgumentException("Pair is not certified in M1. M2/transition evaluation is not implemented; no coefficient is returned.");
            if (!(cutoff > 0 && cutoff <= 24) || spectralOrder < 8 || spectralOrder > 192) throw new ArgumentOutOfRangeException();
            var rule = new MathNet.Numerics.Integration.GaussLegendreRule(0, cutoff, spectralOrder);
            double sum = 0;
            for (int i = 0; i < rule.Order; i++)
            {
                double tau = rule.GetAbscissa(i);
                sum += rule.GetWeight(i) * 2 * Math.Sinh(Math.PI * tau) * tau * Spectral(observer, incident, tau);
            }
            return new Complex(0, -sum / Math.PI);
        }


        internal static (Complex Value, Complex Derivative) GeneralGreen(double dot, Complex nu, double separation = double.NaN)
        {
            //Optional separation = sin²(angle/2). Logarithmic trace quadrature supplies
            //it directly because the dot product can round to one at distinct points.
            double q = double.IsNaN(separation) ? (1 - dot) / 2 : separation;
            if (!(dot >= -1 - 1E-14 && dot <= 1 + 1E-14 && q > 0 && q <= 1 + 1E-14) || double.IsNaN(nu.Magnitude) || nu.Magnitude > 40) throw new ArgumentOutOfRangeException();
            dot = Math.Max(-1, Math.Min(1, dot));
            Complex cos = Complex.Cos(Math.PI * nu);
            if (cos.Magnitude < 1E-10) throw new ArgumentException("Contour must avoid full-sphere eigenvalues.");
            if (dot == -1) return (-1 / (4 * cos), -(0.25 - nu * nu) / (8 * cos));
            if (q < 0.005)
            {
                Complex Psi(Complex w)
                {
                    Complex correction = Complex.Zero;
                    while (w.Real < 16) { correction -= 1 / w; w += 1; }
                    Complex inv = 1 / w, square = inv * inv;
                    return correction + Complex.Log(w) - inv / 2 - square * (1.0 / 12 - square * (1.0 / 120 - square * (1.0 / 252 - square * (1.0 / 240 - square / 132))));
                }
                double log = Math.Log(q);
                Complex h = -2 * 0.5772156649015328606 - Psi(0.5 + nu) - Psi(0.5 - nu);
                Complex coefficient = 1, value = h - log, derivative = 1 / q;
                for (int n = 1; n < 1024; n++)
                {
                    Complex product = (n - 0.5) * (n - 0.5) - nu * nu;
                    h += 2.0 / n - 2 * (n - 0.5) / product;
                    coefficient *= product * q / (n * n);
                    Complex term = coefficient * (h - log), dterm = coefficient / q * (1 - n * (h - log));
                    value += term; derivative += dterm;
                    if (n > nu.Magnitude + 4 && term.Magnitude < 2E-14 * Math.Max(1E-100, value.Magnitude) && dterm.Magnitude < 2E-14 * Math.Max(1E-100, derivative.Magnitude))
                        return (-value / (4 * Math.PI), -derivative / (8 * Math.PI));
                }
                throw new InvalidOperationException("General conical logarithmic kernel did not converge.");
            }
            //Mehler-Dirichlet integral with a squared-endpoint substitution. Differentiation
            //under this regular integral avoids unstable differencing of complex Legendre P.
            double angle = Math.Acos(-dot);
            Complex sum = Complex.Zero, derivativeAngle = Complex.Zero;
            for (int i = 0; i < AngularRule.Order; i++)
            {
                double s = AngularRule.GetAbscissa(i), u = 1 - s * s, a = angle * (1 + u) / 2, b = angle * (1 - u) / 2;
                double factor = AngularRule.GetWeight(i) * 2 * angle * s / Math.Sqrt(2 * Math.Sin(a) * Math.Sin(b));
                Complex c = Complex.Cos(nu * angle * u);
                double logDerivative = 1 / angle - 0.25 * ((1 + u) / Math.Tan(a) + (1 - u) / Math.Tan(b));
                sum += factor * c;
                derivativeAngle += factor * (logDerivative * c - nu * u * Complex.Sin(nu * angle * u));
            }
            Complex scale = -Math.Sqrt(2) / (4 * Math.PI * cos);
            return (scale * sum, scale * derivativeAngle / Math.Sin(angle));
        }

        static Vector GeneralDirection(Vector direction)
        {
            double length = direction.Length();
            if (!(length > 0) || double.IsInfinity(length)) throw new ArgumentException("A finite nonzero direction is required.");
            Vector unit = direction / length;
            if (unit.dx > 0 && unit.dy > 0 && unit.dz > 0) throw new ArgumentException("Direction is inside the solid cone.");
            return unit;
        }

        internal Complex GeneralSpectral(Vector observer, Vector incident, Complex nu)
        {
            observer = GeneralDirection(observer); incident = GeneralDirection(incident);
            bool NearFace(Vector v) => (Math.Abs(v.dx) < 0.05 && v.dy >= 0 && v.dz >= 0) || (Math.Abs(v.dy) < 0.05 && v.dx >= 0 && v.dz >= 0) || (Math.Abs(v.dz) < 0.05 && v.dx >= 0 && v.dy >= 0);
            if (Neumann && NearFace(incident) && !NearFace(observer)) { Vector swap = observer; observer = incident; incident = swap; }
            if (Neumann && NearFace(incident)) throw new ArgumentException("Both directions approach the contour; the double boundary trace is not implemented.");
            var matrix = Matrix<Complex>.Build.Dense(Nodes.Length, Nodes.Length);
            var rhs = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(Nodes.Length);
            for (int i = 0; i < Nodes.Length; i++)
            {
                var g = GeneralGreen(Dot(Nodes[i], incident), nu);
                rhs[i] = Neumann ? g.Derivative * Dot(Normals[i], incident) : -g.Value;
                for (int j = 0; j < Nodes.Length; j++)
                {
                    if (i / (Nodes.Length / 3) == j / (Nodes.Length / 3)) { if (i == j) matrix[i, j] = 0.5; continue; }
                    for (int q = 0; q < Quadrature[j].Length; q++)
                    {
                        Vector source = Quadrature[j][q];
                        Complex derivative = GeneralGreen(Dot(Nodes[i], source), nu).Derivative;
                        matrix[i, j] += Weights[j][q] * derivative * (Neumann ? -Dot(Normals[i], source) : Dot(Nodes[i], Normals[j]));
                    }
                }
            }
            var density = matrix.LU().Solve(rhs);
            if (!((matrix * density - rhs).L2Norm() / rhs.L2Norm() < 1E-9)) throw new InvalidOperationException("General spherical boundary solve failed its residual check.");
            Complex value = Complex.Zero;
            var axes = new[] { new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(0, 0, 1) };
            double Graded(double t) => Math.PI / 2 * (t <= 0.5 ? 0.5 * Math.Pow(2 * t, 3) : 1 - 0.5 * Math.Pow(2 * (1 - t), 3));
            for (int j = 0; j < Nodes.Length; j++)
            {
                int count = Nodes.Length / 3, arc = j / count, panel = j % count;
                double alongA = Dot(observer, axes[arc]), alongB = Dot(observer, axes[(arc + 1) % 3]);
                if (Neumann && alongA >= 0 && alongB >= 0 && Math.Abs(Dot(observer, Normals[j])) < 0.05)
                {
                    //The Neumann single layer is continuous at a face. Resolve its logarithmic
                    //trace integral explicitly, including the exact source-on-baffle limit.
                    double a = Graded((double)panel / count), b = Graded((double)(panel + 1) / count);
                    double observerAngle = Math.Atan2(alongB, alongA), center = Math.Max(a, Math.Min(b, observerAngle));
                    var rule = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, 24);
                    void IntegrateTo(double end)
                    {
                        if (end == center) return;
                        for (int q = 0; q < rule.Order; q++)
                        {
                            double t = rule.GetAbscissa(q), angle = center + (end - center) * t * t;
                            Vector source = axes[arc] * Math.Cos(angle) + axes[(arc + 1) % 3] * Math.Sin(angle);
                            double difference = (center - observerAngle) + (end - center) * t * t;
                            double normal = Dot(observer, Normals[j]), projection = Math.Sqrt(alongA * alongA + alongB * alongB);
                            double separation = 0.5 * normal * normal / (1 + projection) + projection * Math.Pow(Math.Sin(difference / 2), 2);
                            value += density[j] * GeneralGreen(Dot(observer, source), nu, separation).Value * rule.GetWeight(q) * 2 * Math.Abs(end - center) * t;
                        }
                    }
                    IntegrateTo(a); IntegrateTo(b);
                }
                else for (int q = 0; q < Quadrature[j].Length; q++)
                    {
                        var g = GeneralGreen(Dot(observer, Quadrature[j][q]), nu);
                        value += Weights[j][q] * density[j] * (Neumann ? g.Value : g.Derivative * Dot(observer, Normals[j]));
                    }
            }
            return value;
        }

        internal Complex GeneralCoefficient(Vector observer, Vector incident, double damping, double cutoff = 16, double step = 0.5, double contourHeight = 0.5)
        {
            observer = GeneralDirection(observer); incident = GeneralDirection(incident);
            if (!(damping > 0 && damping <= 2 && cutoff >= 4 && cutoff <= 32 && step > 0 && step <= 1 && contourHeight >= 0.2 && contourHeight <= 1)) throw new ArgumentOutOfRangeException();
            Complex Integrand(Complex nu) => Complex.Exp(-Complex.ImaginaryOne * Math.PI * nu - damping * nu) * nu * GeneralSpectral(observer, incident, nu);
            Complex integral = Complex.Zero;
            //Contour runs from +infinity below the positive spectrum, through the left
            //half-plane, then returns above it. End tails must be checked by extending cutoff.
            int panels = (int)Math.Ceiling(cutoff / step);
            for (int p = 0; p < panels; p++)
            {
                var rule = new MathNet.Numerics.Integration.GaussLegendreRule(cutoff * p / panels, cutoff * (p + 1) / panels, 4);
                for (int q = 0; q < rule.Order; q++)
                {
                    double x = rule.GetAbscissa(q);
                    integral += rule.GetWeight(q) * (Integrand(new Complex(x, contourHeight)) - Integrand(new Complex(x, -contourHeight)));
                }
            }
            var vertical = new MathNet.Numerics.Integration.GaussLegendreRule(-contourHeight, contourHeight, 12);
            for (int q = 0; q < vertical.Order; q++) integral += Complex.ImaginaryOne * vertical.GetWeight(q) * Integrand(new Complex(0, vertical.GetAbscissa(q)));
            return Complex.ImaginaryOne * integral / Math.PI;
        }
        internal sealed class FaceField
        {
            readonly Trihedral_Cone Cone;
            readonly Vector Source;
            readonly Complex[][] Density;
            readonly List<ContourSample> Samples;
            internal FaceField(Trihedral_Cone cone, Vector source)
            {
                Cone = cone; source /= source.Length(); Source = source;
                if (Math.Abs(source.dy) > 1E-12 || !(source.dx > 0 && source.dz > 0)) throw new ArgumentException("Canonical source must be on the open Y=0 face.");
                Samples = cone.Contour(); Density = new Complex[Samples.Count][];
                for (int p = 0; p < Samples.Count; p++)
                {
                    //Reciprocity: prepare the single-layer trace on the source face,
                    //then transpose-solve once so each outgoing direction needs no solve.
                    var trace = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(cone.Nodes.Length);
                    var rule = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, 24);
                    int count = cone.Nodes.Length / 3;
                    double Graded(double t) => Math.PI / 2 * (t <= 0.5 ? 0.5 * Math.Pow(2 * t, 3) : 1 - 0.5 * Math.Pow(2 * (1 - t), 3));
                    for (int j = 0; j < trace.Count; j++)
                    {
                        if (j / count == 2)
                        {
                            double a = Graded((double)(j % count) / count), b = Graded((double)(j % count + 1) / count);
                            double sourceAngle = Math.Atan2(source.dx, source.dz), center = Math.Max(a, Math.Min(b, sourceAngle));
                            void IntegrateTo(double end)
                            {
                                if (end == center) return;
                                for (int q = 0; q < rule.Order; q++)
                                {
                                    double t = rule.GetAbscissa(q), angle = center + (end - center) * t * t;
                                    Vector node = new Vector(Math.Sin(angle), 0, Math.Cos(angle));
                                    double difference = (center - sourceAngle) + (end - center) * t * t;
                                    trace[j] += GeneralGreen(Dot(source, node), Samples[p].Nu, Math.Pow(Math.Sin(difference / 2), 2)).Value * rule.GetWeight(q) * 2 * Math.Abs(end - center) * t;
                                }
                            }
                            IntegrateTo(a); IntegrateTo(b);
                        }
                        else for (int q = 0; q < cone.Quadrature[j].Length; q++) trace[j] += cone.Weights[j][q] * GeneralGreen(Dot(source, cone.Quadrature[j][q]), Samples[p].Nu).Value;
                    }
                    var solution = Samples[p].Factor.Solve(trace);
                    if (!((Samples[p].Matrix.Transpose() * solution - trace).L2Norm() / Math.Max(1E-100, trace.L2Norm()) < 1E-9)) throw new InvalidOperationException("Trihedral face solve failed its residual check.");
                    Density[p] = solution.ToArray();
                }
            }
            internal Complex At(Vector direction)
            {
                direction /= direction.Length();
                Complex sum = Complex.Zero;
                for (int p = 0; p < Samples.Count; p++) for (int j = 0; j < Cone.Nodes.Length; j++)
                        sum += Samples[p].Weight * Density[p][j] * Samples[p].Kernel(Dot(direction, Cone.Nodes[j])) * Dot(Cone.Normals[j], direction);
                //Remove the known Abel-Poisson flat-face image before adding to DED.
                double decay = Math.Exp(-0.75);
                Complex flatFace = -Complex.ImaginaryOne * Math.Exp(-0.375) * (1 - decay * decay) / (4 * Math.PI * Math.Pow(1 + 2 * decay * Dot(direction, Source) + decay * decay, 1.5));
                return sum - flatFace;
            }
        }

        sealed class ContourSample
        {
            internal Complex Nu, Weight;
            internal Matrix<Complex> Matrix;
            internal MathNet.Numerics.LinearAlgebra.Factorization.LU<Complex> Factor;
            internal Complex[] Regular;
            internal Complex Kernel(double dot)
            {
                dot = Math.Max(-1, Math.Min(1 - 1E-15, dot));
                double position = Math.Acos(dot) * (Regular.Length - 1) / Math.PI;
                int i = Math.Min(Regular.Length - 2, (int)position);
                double t = position - i;
                Complex a = Regular[Math.Max(0, i - 1)], b = Regular[i], c = Regular[i + 1], d = Regular[Math.Min(Regular.Length - 1, i + 2)];
                return 0.5 * (2 * b + (-a + c) * t + (2 * a - 5 * b + 4 * c - d) * t * t + (-a + 3 * b - 3 * c + d) * t * t * t) - 1 / (4 * Math.PI * (1 - dot)) - (0.25 - Nu * Nu) * Math.Log((1 - dot) / 2) / (8 * Math.PI);
            }
        }

        List<ContourSample> ContourCache;
        List<ContourSample> Contour()
        {
            lock (this)
            {
                if (ContourCache != null) return ContourCache;
                var samples = new List<ContourSample>();
                void Add(Complex nu, Complex quadrature)
                {
                    var matrix = Matrix<Complex>.Build.Dense(Nodes.Length, Nodes.Length);
                    for (int i = 0; i < Nodes.Length; i++) for (int j = 0; j < Nodes.Length; j++)
                        {
                            if (i / (Nodes.Length / 3) == j / (Nodes.Length / 3)) { if (i == j) matrix[i, j] = 0.5; continue; }
                            for (int q = 0; q < Quadrature[j].Length; q++) matrix[i, j] -= Weights[j][q] * GeneralGreen(Dot(Nodes[i], Quadrature[j][q]), nu).Derivative * Dot(Normals[i], Quadrature[j][q]);
                        }
                    Complex[] regular = new Complex[1025];
                    for (int i = 0; i < regular.Length; i++)
                    {
                        double angle = i == 0 ? 1E-5 : Math.PI * i / (regular.Length - 1), dot = Math.Cos(angle);
                        regular[i] = GeneralGreen(dot, nu).Derivative + 1 / (4 * Math.PI * (1 - dot)) + (0.25 - nu * nu) * Math.Log((1 - dot) / 2) / (8 * Math.PI);
                    }
                    samples.Add(new ContourSample { Nu = nu, Weight = Complex.ImaginaryOne / Math.PI * quadrature * Complex.Exp(-Complex.ImaginaryOne * Math.PI * nu - 0.75 * nu) * nu, Matrix = matrix, Factor = matrix.Transpose().LU(), Regular = regular });
                }
                for (int p = 0; p < 32; p++)
                {
                    var rule = new MathNet.Numerics.Integration.GaussLegendreRule(p * 0.5, (p + 1) * 0.5, 4);
                    for (int q = 0; q < rule.Order; q++) { double x = rule.GetAbscissa(q), w = rule.GetWeight(q); Add(new Complex(x, 0.5), w); Add(new Complex(x, -0.5), -w); }
                }
                var vertical = new MathNet.Numerics.Integration.GaussLegendreRule(-0.5, 0.5, 12);
                for (int q = 0; q < vertical.Order; q++) Add(new Complex(0, vertical.GetAbscissa(q)), Complex.ImaginaryOne * vertical.GetWeight(q));
                ContourCache = samples;
                return samples;
            }
        }
        internal FaceField PrepareFace(Vector source) => new FaceField(this, source);
    }

    /// <summary>
    /// Experimental DED plus regularized Neumann trihedral residual at four front corners.
    /// Finite Abel damping=0.75, coarse spherical mesh, smooth engineering activation.
    /// Rear-vertex driving and uniform matching to finite DED endpoint terms remain unfinished.
    /// </summary>
    public sealed class Cabinet_Diffraction_Trihedral : ICabinet_Diffraction
    {
        readonly double Width, Height, Depth, SoundSpeed, VertexScale;
        readonly Cabinet_Diffraction DED;
        static readonly Trihedral_Cone Cone = new Trihedral_Cone(4);
        readonly Dictionary<double, Complex[,]> Patterns = new Dictionary<double, Complex[,]>();

        public Cabinet_Diffraction_Trihedral(double width, double height, double depth, double sound_speed = 343.0, double vertexScale = 1.0)
        {
            if (!(width > 0 && height > 0 && depth > 0 && sound_speed > 0) || double.IsInfinity(width + height + depth + sound_speed) || !(vertexScale >= 0 && vertexScale <= 2)) throw new ArgumentOutOfRangeException("width", "Finite positive cabinet dimensions and vertex scale from 0 through 2 are required.");
            Width = width; Height = height; Depth = depth; SoundSpeed = sound_speed; VertexScale = vertexScale;
            DED = new Cabinet_Diffraction(width, height, depth, sound_speed);
        }

        Complex[,] Pattern(double angle)
        {
            lock (Patterns)
            {
                if (Patterns.TryGetValue(angle, out Complex[,] pattern)) return pattern;
                var face = Cone.PrepareFace(new Vector(Math.Cos(angle), 0, Math.Sin(angle)));
                pattern = new Complex[25, 48];
                System.Threading.Tasks.Parallel.For(0, 25, u =>
                {
                    double t = Math.PI * u / 24;
                    for (int v = 0; v < 48; v++)
                    {
                        double phi = 2 * Math.PI * v / 48;
                        Vector direction = new Vector(Math.Sin(t) * Math.Cos(phi), Math.Cos(t), Math.Sin(t) * Math.Sin(phi));
                        pattern[u, v] = direction.dx >= -1E-12 && direction.dy >= -1E-12 && direction.dz >= -1E-12 ? Complex.Zero : face.At(direction);
                    }
                });
                Patterns.Add(angle, pattern);
                return pattern;
            }
        }

        internal Func<Vector, Complex> PrepareField(Hare.Geometry.Point driver, int octave, double diameter, double distance)
        {
            if (!(diameter > 0) || double.IsInfinity(diameter) || double.IsNaN(driver.x + driver.y + driver.z) || double.IsInfinity(driver.x + driver.y + driver.z) || Math.Abs(driver.y) > 1E-9 || Math.Abs(driver.x) + diameter / 2 > Width / 2 + 1E-9 || Math.Abs(driver.z) + diameter / 2 > Height / 2 + 1E-9) throw new ArgumentException("Driver must fit on the front baffle.");
            octave = Math.Max(0, Math.Min(7, octave));
            var baseline = DED.PrepareField(driver, octave, diameter, distance);
            if (VertexScale == 0) return baseline;
            double k = 2 * Math.PI * Cabinet_Diffraction.Frequencies[octave] / SoundSpeed;
            double ka = k / Math.Sqrt(1 / (Width * Width) + 1 / (Height * Height) + 1 / (Depth * Depth)), activation = Math.Pow(ka, 4) / (16 + Math.Pow(ka, 4));
            var corners = new Hare.Geometry.Point[4]; var tables = new Complex[4][,]; var drive = new Complex[4];
            for (int i = 0; i < 4; i++)
            {
                double sx = i < 2 ? 1 : -1, sz = i % 2 == 0 ? 1 : -1;
                corners[i] = new Hare.Geometry.Point(sx * Width / 2, 0, sz * Height / 2);
                double x = Width / 2 - sx * driver.x, z = Height / 2 - sz * driver.z, r = Math.Sqrt(x * x + z * z);
                tables[i] = Pattern(Math.Atan2(z, x));
                drive[i] = VertexScale * activation * (k * r) * (k * r) / (16 + (k * r) * (k * r)) * 0.5 * DED.PistonFactor(Cabinet_Diffraction.Frequencies[octave], 0, diameter) * Complex.FromPolarCoordinates(1 / r, -k * r);
            }
            return direction =>
            {
                double length = direction.Length(); if (!(length > 0)) return Complex.Zero;
                direction /= length; Hare.Geometry.Point receiver = driver + direction * distance;
                Complex pressure = baseline(direction);
                for (int i = 0; i < 4; i++)
                {
                    Vector path = receiver - corners[i]; double r = path.Length();
                    Vector canonical = new Vector((i < 2 ? -1 : 1) * path.dx / r, -path.dy / r, (i % 2 == 0 ? -1 : 1) * path.dz / r);
                    double t = Math.Max(0, Math.Min(1, Math.Max(-canonical.dx, Math.Max(-canonical.dy, -canonical.dz)) / 0.05)), exterior = t * t * (3 - 2 * t);
                    if (exterior == 0) continue;
                    double u = Math.Acos(Math.Max(-1, Math.Min(1, canonical.dy))) * 24 / Math.PI;
                    double v = Math.Atan2(canonical.dz, canonical.dx) * 48 / (2 * Math.PI); if (v < 0) v += 48;
                    Complex Cubic(Complex a, Complex b, Complex c, Complex d, double f) => 0.5 * (2 * b + (-a + c) * f + (2 * a - 5 * b + 4 * c - d) * f * f + (-a + 3 * b - 3 * c + d) * f * f * f);
                    int iu = Math.Min(23, (int)u), iv = (int)v;
                    Complex Row(int row) => Cubic(tables[i][Math.Max(0, Math.Min(24, row)), (iv + 47) % 48], tables[i][Math.Max(0, Math.Min(24, row)), iv % 48], tables[i][Math.Max(0, Math.Min(24, row)), (iv + 1) % 48], tables[i][Math.Max(0, Math.Min(24, row)), (iv + 2) % 48], v - iv);
                    Complex coefficient = Cubic(Row(iu - 1), Row(iu), Row(iu + 1), Row(iu + 2), u - iu);
                    pressure += drive[i] * exterior * 2 * Math.PI / (k * r) * Complex.Conjugate(coefficient) * Complex.FromPolarCoordinates(1, -k * r);
                }
                return pressure;
            };
        }

        public Complex Pressure(Hare.Geometry.Point driver, Vector direction, int octave, double driverDiameter, double distance)
        {
            if (!(distance > 0) || double.IsInfinity(distance) || double.IsNaN(direction.Length()) || double.IsInfinity(direction.Length())) throw new ArgumentOutOfRangeException("distance");
            return PrepareField(driver, octave, driverDiameter, distance)(direction);
        }

        public string[] Driver_Balloon(Hare.Geometry.Point driver, double driverDiameter) => Cabinet_Diffraction.Balloon(Width, Height, Depth, SoundSpeed, (octave, distance) => PrepareField(driver, octave, driverDiameter, distance));
        readonly Vector[] Nodes, Normals;
        readonly Vector[][] Quadrature;
        readonly double[][] Weights;

        readonly bool Neumann;


        internal static double Dot(Vector a, Vector b) => Hare.Geometry.Hare_math.Dot(a, b);

        internal static Vector Exterior(Vector direction)
        {
            double length = direction.Length();
            if (!(length > 0) || double.IsInfinity(length)) throw new ArgumentException("A finite nonzero direction is required.");
            Vector unit = direction / length;
            if (unit.dx >= 0 && unit.dy >= 0 && unit.dz >= 0) throw new ArgumentException("Direction is in the solid trihedral cone or on its boundary.");
            return unit;
        }

        //Scaled conical Legendre kernel and derivative with respect to dot(omega,omega').
        //The zero-balanced hypergeometric expansion near coincidence avoids a billion-term
        //ordinary series. Cos(pi*i*tau) and the gamma product cancel analytically there.
        internal static (double Value, double Derivative) Green(double dot, double tau)
        {
            if (!(tau >= 0) || tau > 24 || !(dot >= -1 - 1E-14 && dot < 1)) throw new ArgumentOutOfRangeException();
            double z = (1 + Math.Max(-1, dot)) / 2;
            if (z <= 0.995)
            {
                double term = 1, sum = 1, derivativeTerm = 0.25 + tau * tau, derivative = derivativeTerm;
                if (z == 0) return (-1 / (4 * Math.Cosh(Math.PI * tau)), -derivative / (8 * Math.Cosh(Math.PI * tau)));
                for (int n = 1; n < 16384; n++)
                {
                    term *= ((n - 0.5) * (n - 0.5) + tau * tau) * z / (n * n);
                    sum += term;
                    if (n > 1)
                    {
                        derivativeTerm *= ((n - 0.5) * (n - 0.5) + tau * tau) * z / (n * (n - 1));
                        derivative += derivativeTerm;
                    }
                    if (n > tau + 4 && Math.Abs(term) < 2E-15 * Math.Abs(sum) && Math.Abs(derivativeTerm) < 2E-15 * Math.Abs(derivative))
                        return (-sum / (4 * Math.Cosh(Math.PI * tau)), -derivative / (8 * Math.Cosh(Math.PI * tau)));
                }
                throw new InvalidOperationException("Conical Legendre series did not converge.");
            }
            //Complex digamma: recurrence into its asymptotic region, then Bernoulli series.
            //MathNet supplies only the real digamma; no complex helper exists in Universal.
            Complex w = new Complex(0.5, tau), correction = Complex.Zero;
            while (w.Real < 16) { correction -= 1 / w; w += 1; }
            Complex inv = 1 / w, inv2 = inv * inv;
            Complex psi = correction + Complex.Log(w) - inv / 2 - inv2 * (1.0 / 12 - inv2 * (1.0 / 120 - inv2 * (1.0 / 252 - inv2 * (1.0 / 240 - inv2 * (1.0 / 132)))));
            double h = -2 * 0.5772156649015328606 - 2 * psi.Real, qz = 1 - z, log = Math.Log(qz);
            double coefficient = 1, value = h - log, dz = 1 / qz;
            for (int n = 1; n < 16384; n++)
            {
                h += 2.0 / n - 2 * (n - 0.5) / ((n - 0.5) * (n - 0.5) + tau * tau);
                coefficient *= ((n - 0.5) * (n - 0.5) + tau * tau) * qz / (n * n);
                double term = coefficient * (h - log), dterm = coefficient / qz * (1 - n * (h - log));
                value += term; dz += dterm;
                if (n > tau + 4 && Math.Abs(term) < 2E-15 * Math.Max(1E-100, Math.Abs(value)) && Math.Abs(dterm) < 2E-15 * Math.Max(1E-100, Math.Abs(dz)))
                    return (-value / (4 * Math.PI), -dz / (8 * Math.PI));
            }
            throw new InvalidOperationException("Conical logarithmic kernel series did not converge.");
        }

        Matrix<double> BoundaryMatrix(double tau)
        {
            var matrix = Matrix<double>.Build.Dense(Nodes.Length, Nodes.Length);
            for (int i = 0; i < Nodes.Length; i++) for (int j = 0; j < Nodes.Length; j++)
                {
                    //Every arc is a great circle: its same-arc double-layer derivative is zero,
                    //including the principal-value diagonal. The jump term is exactly 1/2.
                    if (i / (Nodes.Length / 3) == j / (Nodes.Length / 3)) { if (i == j) matrix[i, j] = 0.5; continue; }
                    for (int q = 0; q < Quadrature[j].Length; q++)
                    {
                        Vector source = Quadrature[j][q];
                        double derivative = Green(Dot(Nodes[i], source), tau).Derivative;
                        matrix[i, j] += Weights[j][q] * derivative * (Neumann ? -Dot(Normals[i], source) : Dot(Nodes[i], Normals[j]));
                    }
                }
            return matrix;
        }

        double Potential(Vector observer, double tau, double[] density)
        {
            double value = 0;
            for (int j = 0; j < Nodes.Length; j++) for (int q = 0; q < Quadrature[j].Length; q++)
                {
                    Vector source = Quadrature[j][q];
                    var g = Green(Dot(observer, source), tau);
                    value += Weights[j][q] * density[j] * (Neumann ? g.Value : g.Derivative * Dot(observer, Normals[j]));
                }
            return value;
        }

        internal double Spectral(Vector observer, Vector incident, double tau)
        {
            observer = Exterior(observer); incident = Exterior(incident);
            var rhs = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(Nodes.Length);
            for (int i = 0; i < rhs.Count; i++)
            {
                var g = Green(Dot(Nodes[i], incident), tau);
                rhs[i] = Neumann ? g.Derivative * Dot(Normals[i], incident) : -g.Value;
            }
            var matrix = BoundaryMatrix(tau);
            var density = matrix.LU().Solve(rhs);
            double residual = (matrix * density - rhs).L2Norm() / rhs.L2Norm();
            if (!(residual < 1E-10)) throw new InvalidOperationException("Spherical boundary solve failed its residual check.");
            return Potential(observer, tau, density.ToArray());
        }

        //Check the boundary formulation against an independent, exactly harmonic field
        //whose pole is inside the excluded cone; neither the input nor target uses this BIE.
        internal double ManufacturedError(double tau, Vector observer)
        {
            observer = Exterior(observer);
            Vector pole = new Vector(1, 1, 1) / Math.Sqrt(3);
            var rhs = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(Nodes.Length);
            for (int i = 0; i < rhs.Count; i++)
            {
                var g = Green(Dot(Nodes[i], pole), tau);
                rhs[i] = Neumann ? -g.Derivative * Dot(Normals[i], pole) : g.Value;
            }
            double computed = Potential(observer, tau, BoundaryMatrix(tau).LU().Solve(rhs).ToArray());
            double exact = Green(Dot(observer, pole), tau).Value;
            return Math.Abs((computed - exact) / exact);
        }

        //Conservative M1 certificate: for every boundary point, both distances are at least
        //their respective minimum distances to the contour. It may reject valid M1 pairs.
        //This certificate must never be interpreted as a full-sphere visibility rule.
        internal static double M1Margin(Vector observer, Vector incident)
        {
            double Distance(Vector w)
            {
                w = Exterior(w);
                double best = -1;
                double[] c = { w.dx, w.dy, w.dz };
                for (int arc = 0; arc < 3; arc++)
                {
                    double a = c[arc], b = c[(arc + 1) % 3];
                    best = Math.Max(best, Math.Max(a, b));
                    if (a > 0 && b > 0) best = Math.Max(best, Math.Sqrt(a * a + b * b));
                }
                return Math.Acos(Math.Max(-1, Math.Min(1, best)));
            }
            return Distance(observer) + Distance(incident) - Math.PI;
        }

        internal Complex Coefficient(Vector observer, Vector incident, double cutoff = 12, int spectralOrder = 48)
        {
            observer = Exterior(observer); incident = Exterior(incident);
            if (!(M1Margin(observer, incident) > 0.05)) throw new ArgumentException("Pair is not certified in M1. M2/transition evaluation is not implemented; no coefficient is returned.");
            if (!(cutoff > 0 && cutoff <= 24) || spectralOrder < 8 || spectralOrder > 192) throw new ArgumentOutOfRangeException();
            var rule = new MathNet.Numerics.Integration.GaussLegendreRule(0, cutoff, spectralOrder);
            double sum = 0;
            for (int i = 0; i < rule.Order; i++)
            {
                double tau = rule.GetAbscissa(i);
                sum += rule.GetWeight(i) * 2 * Math.Sinh(Math.PI * tau) * tau * Spectral(observer, incident, tau);
            }
            return new Complex(0, -sum / Math.PI);
        }
    }
}
