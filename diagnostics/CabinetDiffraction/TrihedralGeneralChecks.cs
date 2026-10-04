using Pachyderm_Acoustic.Source_Constructions;
using System;
using System.Numerics;
using Vector = Hare.Geometry.Vector;

internal static class TrihedralGeneralChecks
{
    static void Require(bool condition, string message) { if (!condition) throw new Exception(message); }
    internal static void Run()
    {
        foreach (double tau in new[] { 0.0, 1.0, 8.0, 16.0, 24.0 }) foreach (double dot in new[] { -0.9, 0.0, 0.98, 0.991, 0.999 })
        {
            var expected = Trihedral_Cone.Green(dot, tau);
            var actual = Trihedral_Cone.GeneralGreen(dot, new Complex(0, tau));
            Require((actual.Value - expected.Value).Magnitude / Math.Abs(expected.Value) < 1E-7, $"General kernel value, tau={tau}, dot={dot}");
            Require((actual.Derivative - expected.Derivative).Magnitude / Math.Abs(expected.Derivative) < 1E-6, $"General kernel derivative, tau={tau}, dot={dot}");
        }
        foreach (double real in new[] { 0.0, 4.0, 12.0, 24.0 }) foreach (double dot in new[] { 0.0, 0.98, 0.991 })
        {
            Complex nu = new Complex(real, 0.5);
            double step = Math.Min(1E-3, 0.005 * (1 - dot));
            Complex derivative = (-Trihedral_Cone.GeneralGreen(dot + 2 * step, nu).Value + 8 * Trihedral_Cone.GeneralGreen(dot + step, nu).Value - 8 * Trihedral_Cone.GeneralGreen(dot - step, nu).Value + Trihedral_Cone.GeneralGreen(dot - 2 * step, nu).Value) / (12 * step);
            Complex actual = Trihedral_Cone.GeneralGreen(dot, nu).Derivative;
            Require((derivative - actual).Magnitude / actual.Magnitude < 1E-5, $"Complex kernel derivative, real={real}, dot={dot}");
        }
        Console.WriteLine("PASS: general complex spectral kernel agrees with independently validated imaginary kernel and derivative checks.");
        Vector incident = new Vector(-1, -1, -1) / Math.Sqrt(3);
        var cone = new Trihedral_Cone(4);
        const double damping = 0.75;
        var rule = new MathNet.Numerics.Integration.GaussLegendreRule(0, 16, 64);
        Complex reference = Complex.Zero;
        for (int i = 0; i < rule.Order; i++)
        {
            double tau = rule.GetAbscissa(i);
            reference -= Complex.ImaginaryOne / Math.PI * rule.GetWeight(i) * tau * cone.Spectral(incident, incident, tau) * (Complex.Exp(new Complex(Math.PI * tau, -damping * tau)) - Complex.Exp(new Complex(-Math.PI * tau, damping * tau)));
        }
        Complex general = cone.GeneralCoefficient(incident, incident, damping);
        Console.WriteLine($"Damped M1 imaginary contour={reference}; general contour={general}; difference={(general - reference).Magnitude:G8}");
        Require((general - reference).Magnitude < 2E-4, "General contour amplitude, orientation and Abel phase agree with deformed contour.");
        Vector side = new Vector(1, -0.1, 0);
        Complex m2 = cone.GeneralCoefficient(side, incident, damping);
        Complex extended = cone.GeneralCoefficient(side, incident, damping, 24);
        Console.WriteLine($"M2 damped coefficient={m2}; extended cutoff={extended}; difference={(extended - m2).Magnitude:G8}");
        Require((m2 - extended).Magnitude < 0.002, "M2 contour tail check.");
        Complex face = cone.GeneralCoefficient(incident, new Vector(1, 0, 1), damping);
        foreach (double offset in new[] { 0.01, 0.001, 0.0001 })
        {
            Complex approach = cone.GeneralCoefficient(incident, new Vector(1, -offset, 1), damping);
            Console.WriteLine($"Face approach offset={offset}, coefficient={approach}, distance from exact trace={(approach - face).Magnitude:G8}");
            if (offset == 0.0001) Require((approach - face).Magnitude < 0.001, "Continuous face-grazing Neumann trace.");
        }
        Console.WriteLine($"Exact face-grazing finite-damping trace: {face}");
        Console.WriteLine("PASS: finite-damping general contour. The zero-damping limit and physical transition matching remain to be checked.");
    }
}



