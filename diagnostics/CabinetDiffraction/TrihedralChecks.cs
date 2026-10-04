using Pachyderm_Acoustic.Source_Constructions;
using System;
using System.Globalization;
using System.IO;
using System.Numerics;
using Vector = Hare.Geometry.Vector;

internal static class TrihedralChecks
{
    static void Require(bool condition, string message) { if (!condition) throw new Exception(message); }
    internal static void Run(string output)
    {
        Directory.CreateDirectory(output);
        //Independent Mehler-Dirichlet integral (NIST DLMF 14.12.1); no series/digamma.
        var kernelRule = new MathNet.Numerics.Integration.GaussLegendreRule(0, 1, 256);
        double KernelIntegral(double dot, double tau)
        {
            double angle = Math.Acos(-dot), sum = 0;
            for (int i = 0; i < kernelRule.Order; i++)
            {
                double s = kernelRule.GetAbscissa(i), t = angle * (1 - s * s);
                double denominator = Math.Sqrt(2 * Math.Sin((angle + t) / 2) * Math.Sin((angle - t) / 2));
                sum += kernelRule.GetWeight(i) * 2 * angle * s * Math.Cosh(tau * t) / (denominator * Math.Cosh(Math.PI * tau));
            }
            return -Math.Sqrt(2) * sum / (4 * Math.PI);
        }
        double kernelError = 0;
        foreach (double tau in new[] { 0.0, 1.0, 8.0, 16.0, 24.0 }) foreach (double dot in new[] { -0.9, 0.0, 0.98, 0.991, 0.999, 0.9999 })
        {
            var g = Trihedral_Cone.Green(dot, tau);
            double exact = KernelIntegral(dot, tau), error = Math.Abs((g.Value - exact) / exact);
            kernelError = Math.Max(kernelError, error);
            Require(error < 1E-7, "Independent conical Legendre kernel integral.");
            if (dot <= 0.991)
            {
                double step = Math.Min(1E-3, (1 - dot) * 0.01);
                double derivative = (-KernelIntegral(dot + 2 * step, tau) + 8 * KernelIntegral(dot + step, tau) - 8 * KernelIntegral(dot - step, tau) + KernelIntegral(dot - 2 * step, tau)) / (12 * step);
                Require(Math.Abs((g.Derivative - derivative) / derivative) < 2E-6, $"Independent kernel normal derivative: tau={tau}, dot={dot}, computed={g.Derivative:R}, integral={derivative:R}, relative={Math.Abs((g.Derivative - derivative) / derivative):R}");
            }
        }
        Console.WriteLine($"Maximum independent kernel relative error: {kernelError:G8}");
        Vector axis = new Vector(1, 1, 1) / Math.Sqrt(3);
        Vector lateral = new Vector(2, -1, -1) / Math.Sqrt(6);
        Vector incident = axis * -1;
        using var csv = new StreamWriter(Path.Combine(output, "trihedral-reference.csv"));
        csv.WriteLine("condition,panelsPerArc,thetaDegrees,cutoff,spectralOrder,real,imaginary");
        var previousErrors = new double[2, 4];
        foreach (bool neumann in new[] { false, true }) foreach (int panels in new[] { 8, 16, 32 })
        {
            var cone = new Trihedral_Cone(panels, neumann);
            foreach (double tau in new[] { 0.0, 1.0, 4.0, 8.0 })
            {
                double error = cone.ManufacturedError(tau, incident);
                Console.WriteLine($"{(neumann ? "N" : "D")} manufacture: n={panels}, tau={tau}, relative error={error:G8}");
                int boundary = neumann ? 1 : 0, spectral = tau == 0 ? 0 : tau == 1 ? 1 : tau == 4 ? 2 : 3;
                if (panels > 8) Require(error < 0.4 * previousErrors[boundary, spectral], "Manufactured field mesh convergence.");
                previousErrors[boundary, spectral] = error;
                if (panels == 32) Require(error < 0.01, "Manufactured spherical boundary solution.");
            }
            for (int theta = 0; theta <= 5; theta++)
            {
                double t = theta * Math.PI / 24;
                Vector observer = lateral * Math.Sin(t) - axis * Math.Cos(t);
                Complex f = cone.Coefficient(observer, incident);
                csv.WriteLine(FormattableString.Invariant($"{(neumann ? "Neumann" : "Dirichlet")},{panels},{theta * 7.5},12,48,{f.Real:R},{f.Imaginary:R}"));
                Console.WriteLine($"{(neumann ? "N" : "D")} n={panels}, theta={theta * 7.5}, f={f.Imaginary:G10}i");
                Require(!double.IsNaN(f.Imaginary) && !double.IsInfinity(f.Imaginary), "Finite coefficient.");
                //Bonner table 6.7 uses Dirichlet data, not the cabinet's Neumann condition.
                //Published 192-panel collocation values: useful independent amplitude/phase check.
                double[] reference = { -0.067187, -0.068032, -0.070720, -0.075706, -0.083874, -0.096960 };
                if (!neumann && panels == 32) Require(Math.Abs(f.Imaginary - reference[theta]) < 0.001, "Published Dirichlet trihedral benchmark.");
            }
        }
        Vector sample = lateral * Math.Sin(Math.PI / 12) - axis * Math.Cos(Math.PI / 12);
        var rigid = new Trihedral_Cone(32);
        Complex baseline = rigid.Coefficient(sample, incident);
        Complex quadrature = rigid.Coefficient(sample, incident, 12, 72);
        Complex tail = rigid.Coefficient(sample, incident, 16, 72);
        Complex reciprocal = rigid.Coefficient(incident, sample);
        Vector permutation = new Vector(sample.dz, sample.dx, sample.dy);
        Complex symmetric = rigid.Coefficient(permutation, incident);
        Console.WriteLine($"Neumann spectral refinement: {(quadrature - baseline).Magnitude:G8}; cutoff refinement: {(tail - quadrature).Magnitude:G8}; reciprocity: {(reciprocal - baseline).Magnitude:G8}; symmetry: {(symmetric - baseline).Magnitude:G8}");
        Require((quadrature - baseline).Magnitude < 1E-5, "Spectral quadrature convergence.");
        Require((tail - quadrature).Magnitude < 1E-4, "Spectral tail convergence.");
        Require((reciprocal - baseline).Magnitude < 0.001, "Canonical reciprocity.");
        Require((symmetric - baseline).Magnitude < 1E-10, "Trihedral axis permutation symmetry.");
        var finer = new Trihedral_Cone(64);
        Complex mesh = finer.Coefficient(sample, incident);
        double manufactured = finer.ManufacturedError(8, incident);
        Console.WriteLine($"Neumann coefficient mesh refinement: {(mesh - baseline).Magnitude:G8}; refined tau=8 manufactured error: {manufactured:G8}");
        Require((mesh - baseline).Magnitude < 1E-5, "Neumann coefficient mesh convergence.");
        Require(manufactured < 0.002, "Refined independent Neumann field.");
        //A front-baffle source lies on a tangent face. The face limit remains outside M1;
        //production must not quietly replace it with an axial incident plane wave.
        bool grazingRejected = false;
        try { rigid.Coefficient(incident, new Vector(1, -1E-6, 1)); } catch (ArgumentException) { grazingRejected = true; }
        Require(grazingRejected, "Cabinet face-grazing incidence requires the missing general-contour calculation.");
        bool rejected = false;
        try { rigid.Coefficient(new Vector(1, -0.1, 0), incident); } catch (ArgumentException) { rejected = true; }
        Require(rejected, "Unsupported transition/M2 directions must not return a fabricated coefficient.");
        Console.WriteLine("PASS: restricted canonical trihedral reference. Full cabinet matching is not implemented.");
    }
}

