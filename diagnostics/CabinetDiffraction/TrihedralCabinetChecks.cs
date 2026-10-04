using System;
using System.IO;
using System.Numerics;
using Pachyderm_Acoustic.Source_Constructions;
using Point = Hare.Geometry.Point;
using Vector = Hare.Geometry.Vector;

internal static class TrihedralCabinetChecks
{
    static void Require(bool test, string message) { if (!test) throw new Exception(message); }
    internal static void DefaultColumn()
    {
        //Both separations are lost in a rounded dot product of one.
        Complex nu = new Complex(4, 0.5);
        var near = Trihedral_Cone.GeneralGreen(1, nu, 1E-20);
        var nearer = Trihedral_Cone.GeneralGreen(1, nu, 1E-24);
        Require((nearer.Value - near.Value - Math.Log(1E-4) / (4 * Math.PI)).Magnitude < 1E-12, "Coincidence logarithm retains sub-epsilon separation.");
        Require((nearer.Derivative * (8 * Math.PI * 1E-24) + 1).Magnitude < 1E-12, "Coincidence derivative retains analytic singular term.");
        bool rejected = false;
        try { Trihedral_Cone.GeneralGreen(1, nu, 0); } catch (ArgumentOutOfRangeException) { rejected = true; }
        Require(rejected, "An actual coincident kernel remains singular.");
        const int count = 8;
        const double spacing = 0.075, diameter = 0.05;
        var cabinet = new Cabinet_Diffraction_Trihedral(0.08, (count - 1) * spacing + diameter, 0.10);
        for (int i = 0; i < count; i++)
        {
            double offset = ((count - 1) * 0.5 - i) * spacing;
            Console.WriteLine($"Default column driver {i + 1}/{count}, z={offset:R}");
            string[] balloons = cabinet.Driver_Balloon(new Point(0, 0, offset), diameter);
            Require(balloons.Length == 8, "Default column returns all eight bands.");
            foreach (string balloon in balloons) Require(!string.IsNullOrEmpty(balloon) && !balloon.Contains("NaN") && !balloon.Contains("Infinity"), "Finite default-column balloons.");
        }
        Console.WriteLine("PASS: all eight default column drivers generate finite eight-band experimental balloons.");
    }
    internal static void Run(string output)
    {
        Directory.CreateDirectory(output);
        Vector incident = new Vector(1, 0, 1) / Math.Sqrt(2), observer = new Vector(-1, -1, -1) / Math.Sqrt(3);
        var cone = new Trihedral_Cone(4);
        var face = cone.PrepareFace(incident);
        Complex residual = face.At(observer), trace = cone.GeneralCoefficient(observer, incident, 0.75);
        double decay = Math.Exp(-0.75);
        Complex flatFace = -Complex.ImaginaryOne * Math.Exp(-0.375) * (1 - decay * decay) / (4 * Math.PI * Math.Pow(1 + 2 * decay * Trihedral_Cone.Dot(observer, incident) + decay * decay, 1.5));
        Console.WriteLine($"Prepared face residual={residual}; reference trace minus flat-face image={trace - flatFace}; difference={(residual - trace + flatFace).Magnitude:G8}");
        Require((residual - trace + flatFace).Magnitude < 1E-7, "Prepared Neumann face-source residual matches independent reciprocity trace.");
        var driver = new Point(0, 0, 0); const double diameter = 0.05, distance = 20;
        var cabinet = new Cabinet_Diffraction_Trihedral(0.08, 0.16, 0.10);
        var ded = new Cabinet_Diffraction(0.08, 0.16, 0.10);
        var disabled = new Cabinet_Diffraction_Trihedral(0.08, 0.16, 0.10, vertexScale: 0);
        using var csv = new StreamWriter(Path.Combine(output, "trihedral-cabinet.csv"));
        csv.WriteLine("octave,angleDegrees,dedReal,dedImaginary,experimentalReal,experimentalImaginary");
        double largestChange = 0;
        for (int octave = 0; octave < 8; octave++)
        {
            var field = cabinet.PrepareField(driver, octave, diameter, distance);
            for (int angle = 0; angle <= 360; angle++)
            {
                double a = Math.PI * angle / 180;
                Vector direction = new Vector(Math.Sin(a), Math.Cos(a), 0);
                Complex p = field(direction), reference = ded.Pressure(driver, direction, octave, diameter, distance);
                Require(!double.IsNaN(p.Magnitude) && !double.IsInfinity(p.Magnitude), "Finite experimental full-sphere pressure.");
                Require((p - field(new Vector(-direction.dx, direction.dy, direction.dz))).Magnitude < 1E-10, "Centered cabinet horizontal symmetry.");
                largestChange = Math.Max(largestChange, (p - reference).Magnitude * distance);
                csv.WriteLine(FormattableString.Invariant($"{octave},{angle},{reference.Real:R},{reference.Imaginary:R},{p.Real:R},{p.Imaginary:R}"));
            }
            Vector sample = new Vector(0.4, -0.6, 0.3);
            Require((disabled.Pressure(driver, sample, octave, diameter, distance) - ded.Pressure(driver, sample, octave, diameter, distance)).Magnitude < 1E-12, "Zero vertex scale exactly restores DED.");
            Complex left = field(new Vector(1, 1E-8, 0)), right = field(new Vector(1, -1E-8, 0));
            Require((left - right).Magnitude * distance < 1E-5, "Grazing continuity.");
        }
        Require(largestChange > 1E-4, "Selected experimental method actually changes the field.");
        //Depth only changes the smooth activation of the front-corner correction.
        //Subtract DED so its inherited subdivision seams do not obscure this check.
        double worstDepthSecondDifference = 0, coarseDepthSecondDifference = 0;
        var depthModels = new Cabinet_Diffraction_Trihedral[41];
        for (int step = 0; step <= 40; step++) depthModels[step] = new Cabinet_Diffraction_Trihedral(0.08, 0.16, 0.05 + step * 0.0025);
        using (var depthCsv = new StreamWriter(Path.Combine(output, "trihedral-depth.csv")))
        {
            depthCsv.WriteLine("octave,depth,angleDegrees,correctionReal,correctionImaginary");
            for (int octave = 0; octave < 8; octave++)
            {
                var previous = new Complex[3]; var beforePrevious = new Complex[3]; var coarsePrevious = new Complex[3]; var coarseBeforePrevious = new Complex[3];
                for (int step = 0; step <= 40; step++)
                {
                    double depth = 0.05 + step * 0.0025;
                    var trial = depthModels[step].PrepareField(driver, octave, diameter, distance);
                    var baseField = new Cabinet_Diffraction(0.08, 0.16, depth).PrepareField(driver, octave, diameter, distance);
                    for (int j = 0; j < 3; j++)
                    {
                        double angle = 90 + j * 45, a = angle * Math.PI / 180;
                        Vector direction = new Vector(Math.Sin(a), Math.Cos(a), 0);
                        Complex correction = (trial(direction) - baseField(direction)) * distance;
                        Require(!double.IsNaN(correction.Magnitude) && !double.IsInfinity(correction.Magnitude), "Finite side/rear depth sweep.");
                        if (step > 1) worstDepthSecondDifference = Math.Max(worstDepthSecondDifference, (correction - 2 * previous[j] + beforePrevious[j]).Magnitude);
                        beforePrevious[j] = previous[j]; previous[j] = correction;
                        if (step % 2 == 0)
                        {
                            if (step > 3) coarseDepthSecondDifference = Math.Max(coarseDepthSecondDifference, (correction - 2 * coarsePrevious[j] + coarseBeforePrevious[j]).Magnitude);
                            coarseBeforePrevious[j] = coarsePrevious[j]; coarsePrevious[j] = correction;
                        }
                        depthCsv.WriteLine(FormattableString.Invariant($"{octave},{depth:R},{angle},{correction.Real:R},{correction.Imaginary:R}"));
                    }
                }
            }
        }
        Require(worstDepthSecondDifference < 0.35 * coarseDepthSecondDifference, "Smooth experimental depth correction.");
        Console.WriteLine($"Depth sweep maximum correction second difference: {worstDepthSecondDifference:G8}; fine/coarse ratio={worstDepthSecondDifference / coarseDepthSecondDifference:G6}");
        var balloons = cabinet.Driver_Balloon(driver, diameter);
        Require(balloons.Length == 8, "Existing eight-octave balloon API.");
        foreach (string balloon in balloons) Require(!string.IsNullOrEmpty(balloon) && !balloon.Contains("NaN") && !balloon.Contains("Infinity"), "Finite balloon encoding.");
        Console.WriteLine($"PASS: selectable experimental cabinet field, symmetry, grazing continuity, zero-scale DED regression and eight balloon strings. Largest unit-pressure change: {largestChange:G6}");
    }
}
