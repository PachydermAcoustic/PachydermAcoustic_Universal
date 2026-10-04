using System;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Reflection;
using System.Globalization;
using Pachyderm_Acoustic.Source_Constructions;
using Point = Hare.Geometry.Point;
using Vector = Hare.Geometry.Vector;

static class NumericalChecks
{
    static void Require(bool condition, string message)
    {
        if (!condition) throw new Exception(message);
    }
    static Vector Polar(double degrees)
    {
        double a = degrees * Math.PI / 180;
        return new Vector(Math.Sin(a), Math.Cos(a), 0);
    }
    internal static void Run(string output)
    {
        Directory.CreateDirectory(output);
        var center = new Point(0, 0, 0);
        var compact = new Cabinet_Diffraction_Numerical(0.03, 0.03, 0.03, highFrequencyContinuation: false);
        var monopole = compact.Pressure(center, new Vector(0, 1, 0), 0, 0.01, 20);
        Console.WriteLine($"Compact-cabinet low-frequency distance-scaled magnitude: {20 * monopole.Magnitude:G8} (monopole limit 0.5)");
        Require(Math.Abs(20 * monopole.Magnitude - 0.5) < 0.02, "compact-cabinet monopole limit");
        Complex rear = compact.Pressure(center, new Vector(0, -1, 0), 0, 0.01, 20);
        Require(Math.Abs(rear.Magnitude / monopole.Magnitude - 1) < 0.01, "compact-cabinet near-omnidirectional low-frequency limit");
        var derivative = typeof(Cabinet_Diffraction_Numerical).GetMethod("NormalDerivative", BindingFlags.Static | BindingFlags.NonPublic);
        var receiver = new Point(0.2, 0.3, -0.4); var source = new Point(-0.1, -0.2, 0.05); var normal = new Vector(0, 1, 0);
        double k = 12, eps = 1E-6;
        Complex Green(Point p) { double r = (p - source).Length(); return Complex.FromPolarCoordinates(1 / r, -k * r); }
        Complex finiteDifference = (Green(receiver + normal * eps) - Green(receiver - normal * eps)) / (2 * eps);
        Complex analytical = (Complex)derivative.Invoke(null, new object[]{receiver, source, normal, k});
        Require((finiteDifference - analytical).Magnitude < 1E-7, "normal derivative sign and outgoing phase");
        var cabinet = new Cabinet_Diffraction_Numerical(0.08, 0.16, 0.10, highFrequencyContinuation: false);
        var coarse = new Cabinet_Diffraction_Numerical(0.08, 0.16, 0.10, 343, 3, highFrequencyContinuation: false);
        var refined = new Cabinet_Diffraction_Numerical(0.08, 0.16, 0.10, 343, 6, highFrequencyContinuation: false);
        var ded = new Cabinet_Diffraction(0.08, 0.16, 0.10);
        var hybrid = new Cabinet_Diffraction_Numerical(0.08, 0.16, 0.10);
        var prepare = typeof(Cabinet_Diffraction_Numerical).GetMethod("Prepare", BindingFlags.Instance | BindingFlags.NonPublic);
        using (var csv = new StreamWriter(Path.Combine(output, "numerical-polar.csv")))
        {
            csv.WriteLine("octave,angle,real,imaginary,magnitude,coarse,refined,ded,hybrid");
            for (int octave = 0; octave < 8; octave++)
            {
                double change = 0, reference = 0, coarseChange = 0, maxGrazing = 0;
                object field = prepare.Invoke(cabinet, new object[]{center, octave, 0.05});
                var normalAt = field.GetType().GetMethod("NormalAt", BindingFlags.Instance | BindingFlags.NonPublic);
                Complex rigidFlux = (Complex)normalAt.Invoke(field, new object[]{new Point(0.04, -0.043, 0.013), new Vector(1, 0, 0)});
                object refinedField = prepare.Invoke(refined, new object[]{center, octave, 0.05});
                Complex refinedFlux = (Complex)normalAt.Invoke(refinedField, new object[]{new Point(0.04, -0.043, 0.013), new Vector(1, 0, 0)});
                Console.WriteLine($"Band {octave}: off-grid rigid-wall flux / piston flux, default {rigidFlux.Magnitude / 3200:G4}, refined {refinedFlux.Magnitude / 3200:G4}");
                if (octave <= 5) Require(refinedFlux.Magnitude / 3200 < 0.05 && refinedFlux.Magnitude < rigidFlux.Magnitude, "low-band rigid-wall refinement reduces flux error below 5% of piston flux");
                //Report high-band surface errors separately; far-field convergence is not near-field validation.
                var oblique = new Vector(0.3, 0.5, 0.8);
                Require((cabinet.Pressure(center, oblique, octave, 0.05, 20) - cabinet.Pressure(center, new Vector(0.3, 0.5, -0.8), octave, 0.05, 20)).Magnitude < 1E-8, "up/down complex symmetry");
                for (int angle = 0; angle <= 360; angle++)
                {
                    Vector direction = Polar(angle);
                    Complex p = cabinet.Pressure(center, direction, octave, 0.05, 20);
                    Complex c = coarse.Pressure(center, direction, octave, 0.05, 20);
                    Complex q = refined.Pressure(center, direction, octave, 0.05, 20);
                    Require(double.IsFinite(p.Magnitude) && double.IsFinite(q.Magnitude), "finite pressure");
                    Complex mirror = cabinet.Pressure(center, new Vector(-direction.dx, direction.dy, 0), octave, 0.05, 20);
                    Require((p - mirror).Magnitude < 1E-8, "left/right complex symmetry");
                    change += (p - q).Magnitude * (p - q).Magnitude; reference += q.Magnitude * q.Magnitude;
                    coarseChange += (c - q).Magnitude * (c - q).Magnitude;
                    csv.WriteLine(FormattableString.Invariant($"{octave},{angle},{p.Real},{p.Imaginary},{p.Magnitude},{c.Magnitude},{q.Magnitude},{ded.Pressure(center, direction, octave, 0.05, 20).Magnitude},{hybrid.Pressure(center, direction, octave, 0.05, 20).Magnitude}"));
                    Complex edge = ded.Pressure(center, direction, octave, 0.05, 20);
                    Complex expected = octave <= 3 ? p : octave == 4 ? (p + edge) * 0.5 : edge;
                    Require((hybrid.Pressure(center, direction, octave, 0.05, 20) - expected).Magnitude < 1E-12, "hybrid selects and blends complex fields exactly");
                }
                foreach (double angle in new[]{90.0, 270.0})
                {
                    double jump = 20 * (cabinet.Pressure(center, Polar(angle - 1E-5), octave, 0.05, 20) - cabinet.Pressure(center, Polar(angle + 1E-5), octave, 0.05, 20)).Magnitude;
                    maxGrazing = Math.Max(maxGrazing, jump);
                    Require(jump < 1E-5, "grazing continuity");
                }
                Console.WriteLine($"Band {octave}: relative complex refinement change {Math.Sqrt(change / reference):G4}; coarse-to-refined {Math.Sqrt(coarseChange / reference):G4}; grazing change {maxGrazing:G4}");
                Require(Math.Sqrt(change / reference) < 0.10, "complex angular field convergence within 10% on diagnostic cabinet");
            }
        }
        //Compare depth changes, including where the source count changes.
        using (var csv = new StreamWriter(Path.Combine(output, "numerical-depth.csv")))
        {
            csv.WriteLine("octave,depth,angle,real,imaginary,magnitude,refined,hybrid");
            for (int octave = 0; octave < 8; octave++) for (int step = 0; step <= 16; step++)
            {
                double depth = 0.06 + step * 0.005;
                var model = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth, highFrequencyContinuation: false);
                var continued = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth);
                var fine = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth, 343, 6, highFrequencyContinuation: false);
                foreach (double angle in new[]{90.0, 135.0, 180.0})
                {
                    Complex pressure = model.Pressure(new Point(0.005, 0, 0.02), Polar(angle), octave, 0.05, 20);
                    csv.WriteLine(FormattableString.Invariant($"{octave},{depth},{angle},{pressure.Real},{pressure.Imaginary},{pressure.Magnitude},{fine.Pressure(new Point(0.005, 0, 0.02), Polar(angle), octave, 0.05, 20).Magnitude},{continued.Pressure(new Point(0.005, 0, 0.02), Polar(angle), octave, 0.05, 20).Magnitude}"));
                }
            }
        }
        double maxSeam = 0, maxHybridSeam = 0;
        for (int octave = 0; octave < 8; octave++)
        {
            double spacing = Math.Min(0.7 * 0.16 / 6, 343 / (62.5 * Math.Pow(2, octave) * 4));
            double depth = 4 * spacing / 0.7;
            if (depth < 0.04 || depth > 0.14) continue;
            var left = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth - 1E-8, highFrequencyContinuation: false);
            var right = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth + 1E-8, highFrequencyContinuation: false);
            double seam = 20 * (left.Pressure(center, Polar(135), octave, 0.05, 20) - right.Pressure(center, Polar(135), octave, 0.05, 20)).Magnitude;
            maxSeam = Math.Max(maxSeam, seam);
            var hybridLeft = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth - 1E-8);
            var hybridRight = new Cabinet_Diffraction_Numerical(0.08, 0.16, depth + 1E-8);
            double hybridSeam = 20 * (hybridLeft.Pressure(center, Polar(135), octave, 0.05, 20) - hybridRight.Pressure(center, Polar(135), octave, 0.05, 20)).Magnitude;
            maxHybridSeam = Math.Max(maxHybridSeam, hybridSeam);
            Require(hybridSeam < 0.005, "selected hybrid depth seam below 0.5% of unit pressure");
            Require(seam < 0.05, "numerical source-count seam below 5% of unit pressure");
        }
        Console.WriteLine($"Largest tested numerical depth subdivision seam: {maxSeam:G4} of unit pressure.");
        Console.WriteLine($"Largest tested selected-hybrid depth seam: {maxHybridSeam:G4} of unit pressure.");
        //The default column must complete for both a central and an end driver.
        var column = new Cabinet_Diffraction_Numerical(0.08, 0.575, 0.10, highFrequencyContinuation: false);
        foreach (var driver in new[]{new Point(0, 0, 0.0375), new Point(0, 0, 0.2625)})
        {
            for (int octave = 0; octave < 8; octave++)
            {
                var p = column.Pressure(driver, Polar(135), octave, 0.05, 20);
                Require(double.IsFinite(p.Magnitude), "finite default column field");
            }
        }
        Console.WriteLine("PASS: default-column central and end drivers across all eight bands.");
        Require(cabinet.Pressure(center, Polar(35), -1, 0.05, 20) == cabinet.Pressure(center, Polar(35), 0, 0.05, 20), "octave clamp");
        Require(cabinet.Pressure(center, new Vector(0, 0, 0), 0, 0.05, 20) == Complex.Zero, "zero direction");
        bool rejected = false;
        try { new Cabinet_Diffraction_Numerical(1, 1, 1, highFrequencyContinuation: false); } catch (ArgumentException) { rejected = true; }
        Require(rejected, "source budget rejects oversized cabinet");
        rejected = false;
        try { cabinet.Pressure(center, Polar(0), 0, 0.2, 20); } catch (ArgumentException) { rejected = true; }
        Require(rejected, "piston outside baffle rejected");
        string[] balloons = hybrid.Driver_Balloon(center, 0.05);
        Require(balloons.Length == 8, "eight numerical balloons");
        foreach (string balloon in balloons)
        {
            var rows = balloon.Split(';', StringSplitOptions.RemoveEmptyEntries);
            Require(rows.Length == 72, "72 meridians");
            foreach (string row in rows)
            {
                var values = row.Split(' ').Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
                Require(values.Length == 37 && values.All(v => double.IsFinite(v) && v >= 0 && v <= 60), "bounded numerical balloon format");
            }
        }
        Console.WriteLine("PASS: monopole limit, derivative, symmetry, convergence, grazing continuity, validation and balloon format.");
    }
}