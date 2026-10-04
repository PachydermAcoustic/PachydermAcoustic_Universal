using System;
using System.Collections;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Reflection;
using Pachyderm_Acoustic.Source_Constructions;
using Point = Hare.Geometry.Point;
using Vector = Hare.Geometry.Vector;

//Source-linked diagnostic: exercises the real cabinet code with Hare, without loading Rhino.
class Program
{
    const double Distance = 20;
    const double Diameter = 0.13;
    static readonly Point Center = new Point(0, 0, 0);
    static readonly Point Offset = new Point(0.03, 0, 0.17);
    static readonly Complex[] FrontGolden = new Complex[]{
        new Complex(-0.021897434086073773, 0.016329103460466825),
        new Complex(0.005568148389194235, -0.032363531250636235),
        new Complex(-0.045090219970135105, 0.0023207972032073566),
        new Complex(0.04324890382104409, -0.03776366960212772),
        new Complex(-0.029155998873850526, -0.05218433792085691),
        new Complex(-0.03855196021933678, 0.03338445374356917),
        new Complex(0.004828223093445394, -0.04958715503148994),
        new Complex(-0.049337401302435274, -0.00866863548657702),
        new Complex(-0.01274976487138348, 0.02073002140803053),
        new Complex(-0.009727497646946562, -0.020452322353195374),
        new Complex(-0.013932436190550608, 0.012291109609188134),
        new Complex(0.007745785549019751, -0.014153690531547835),
        new Complex(-0.007879406791487397, -0.014018410821543752),
        new Complex(-0.005723604847018169, 0.006679014770253815),
        new Complex(0.0001191711256908558, 0.0015911631601034527),
        new Complex(-0.0008923368090690576, -0.00011352826829130939)
    };

    //Cache the same per-band distributions that Driver_Balloon caches, for dense diagnostic sweeps.
    sealed class Field
    {
        readonly Cabinet_Diffraction Cabinet;
        readonly Point Driver;
        readonly int Octave;
        readonly object Front, Driven;
        readonly Complex[] Strength;
        readonly MethodInfo Evaluate;
        public readonly int FrontCount, DrivenCount;

        public Field(double depth, Point driver, int octave, bool frontOnly = false, int maxOrder = 4)
        {
            Cabinet = new Cabinet_Diffraction(0.32, 0.8, depth);
            Driver = driver;
            Octave = octave;
            const BindingFlags flags = BindingFlags.Instance | BindingFlags.NonPublic;
            var edges = typeof(Cabinet_Diffraction).GetMethod("EdgeElements", flags);
            Front = edges.Invoke(Cabinet, new object[]{ driver, octave, false });
            Driven = edges.Invoke(Cabinet, new object[]{ driver, octave, true });
            if (frontOnly) Driven = Activator.CreateInstance(Front.GetType());
            FrontCount = ((IList)Front).Count;
            DrivenCount = ((IList)Driven).Count;
            Strength = (Complex[])typeof(Cabinet_Diffraction).GetMethod("DrivenEdgeDrive", flags).Invoke(Cabinet, new object[]{ Front, Driven, 62.5 * Math.Pow(2, octave), Diameter, maxOrder });
            Evaluate = typeof(Cabinet_Diffraction).GetMethods(flags).Single(m => m.Name == "Pressure");
        }

        public Complex At(Vector direction, bool omitDepth = false)
        {
            Complex[] strength = Strength;
            if (omitDepth)
            {
                strength = (Complex[])Strength.Clone();
                Array.Clear(strength, FrontCount, DrivenCount - FrontCount);
            }
            return (Complex)Evaluate.Invoke(Cabinet, new object[]{ Driver, direction, Octave, Diameter, Distance, Front, Driven, strength });
        }

#if BASELINE
        public void CompareRearDrive(double depth)
        {
            var legacy = new Legacy_Cabinet_Diffraction(0.32, 0.8, depth);
            const BindingFlags flags = BindingFlags.Instance | BindingFlags.NonPublic;
            var front = typeof(Legacy_Cabinet_Diffraction).GetMethod("EdgeElements", flags).Invoke(legacy, new object[]{ Driver, Octave });
            var rear = (Complex[])typeof(Legacy_Cabinet_Diffraction).GetMethod("RearEdgeDrive", flags).Invoke(legacy, new object[]{ front, 62.5 * Math.Pow(2, Octave), Diameter });
            var second = (Complex[])typeof(Cabinet_Diffraction).GetMethod("DrivenEdgeDrive", flags).Invoke(Cabinet, new object[]{ Front, Driven, 62.5 * Math.Pow(2, Octave), Diameter, 2 });
            for (int i = 0; i < rear.Length; i++) Near(second[i], rear[i], "original rear drive");
        }
#endif
    }

    static void Require(bool condition, string name) { if (!condition) throw new Exception(name); }
    static void Near(Complex a, Complex b, string name) => Require((a - b).Magnitude < 1E-12, name);
    static Vector Polar(double degrees, double azimuth = 0)
    {
        double t = degrees * Math.PI / 180, p = azimuth * Math.PI / 180;
        return new Vector(Math.Sin(t) * Math.Cos(p), Math.Cos(t), Math.Sin(t) * Math.Sin(p));
    }

    static void Main(string[] args)
    {
        if (args.Length == 2 && args[0] == "--compare-balloon")
        {
            var originalType = Assembly.LoadFrom(Path.GetFullPath(args[1])).GetType("Pachyderm_Acoustic.Source_Constructions.Cabinet_Diffraction", true);
            var original = Activator.CreateInstance(originalType, new object[]{0.32, 0.8, 0.25, 343.0});
            var expected = (string[])originalType.GetMethod("Driver_Balloon").Invoke(original, new object[]{Offset, Diameter});
            var actual = new Cabinet_Diffraction(0.32, 0.8, 0.25).Driver_Balloon(Offset, Diameter);
            Require(actual.SequenceEqual(expected), "default DED balloon strings remain byte-for-byte identical");
            Console.WriteLine("PASS: all eight DED balloon strings are identical to the pre-selection implementation.");
            return;
        }
        if (args.Length > 0 && args[0] == "--numerical") { NumericalChecks.Run(args.Length > 1 ? args[1] : Path.Combine(Path.GetTempPath(), "pachyderm-numerical-cabinet")); return; }
        if (args.Length > 0 && args[0] == "--trihedral") { TrihedralChecks.Run(args.Length > 1 ? args[1] : Path.Combine(Path.GetTempPath(), "pachyderm-trihedral")); return; }
        if (args.Length > 0 && args[0] == "--trihedral-column") { TrihedralCabinetChecks.DefaultColumn(); return; }
        if (args.Length > 0 && args[0] == "--trihedral-general") { TrihedralGeneralChecks.Run(); return; }
        if (args.Length > 0 && args[0] == "--trihedral-cabinet") { TrihedralCabinetChecks.Run(args.Length > 1 ? args[1] : Path.Combine(Path.GetTempPath(), "pachyderm-trihedral-cabinet")); return; }
        string output = args.Length > 0 ? args[0] : Path.Combine(Path.GetTempPath(), "pachyderm-cabinet-diagnostics");
        Directory.CreateDirectory(output);
        var directions = new[]{ new Vector(0, 1, 0), new Vector(1, -0.4, 0.3) };
        for (int octave = 0; octave < 8; octave++)
        {
            var zero = new Field(0, Offset, octave);
            var front = new Field(0.25, Offset, octave, true);
            Require(zero.DrivenCount == 0, "zero depth must omit driven edges");
            for (int i = 0; i < directions.Length; i++)
            {
                Near(zero.At(directions[i]), FrontGolden[i * 8 + octave], "direct/front golden");
                Near(front.At(directions[i]), zero.At(directions[i]), "front field independent of depth");
            }
            var centered = new Field(0.25, Center, octave);
            double spacing = Math.Min(0.020, 343 / (62.5 * Math.Pow(2, octave)) / 8);
            Require(centered.DrivenCount == centered.FrontCount + 8 * (int)Math.Ceiling(0.25 / spacing), "four depth edges with two incident faces each");
            var offsetField = new Field(0.25, Offset, octave);
            var secondOrder = new Field(0.25, Offset, octave, false, 2);
            var thirdOrder = new Field(0.25, Offset, octave, false, 3);
            Complex thirdContribution = thirdOrder.At(Polar(135, 35)) - secondOrder.At(Polar(135, 35));
            Complex fourthContribution = offsetField.At(Polar(135, 35)) - thirdOrder.At(Polar(135, 35));
            Require(thirdContribution.Magnitude > 1E-12 && fourthContribution.Magnitude > 1E-12, "nonzero third and fourth orders");
            Console.WriteLine($"Band {octave}: scaled |order3|={Distance * thirdContribution.Magnitude:G4}, |order4|={Distance * fourthContribution.Magnitude:G4}");
            Complex depthContribution = offsetField.At(Polar(135, 35)) - offsetField.At(Polar(135, 35), true);
            Require(depthContribution.Magnitude > 1E-10 && Math.Abs(depthContribution.Imaginary) > 1E-12, "depth field must contribute with complex phase");
            for (int angle = 0; angle < 360; angle += 13)
            {
                var d = Polar(angle, 37);
                var p = centered.At(d);
                Require(double.IsFinite(p.Real) && double.IsFinite(p.Imaginary), "finite pressure");
                Near(p, centered.At(new Vector(-d.dx, d.dy, d.dz)), "left/right symmetry");
                Near(p, centered.At(new Vector(d.dx, d.dy, -d.dz)), "top/bottom symmetry");
            }
#if BASELINE
            centered.CompareRearDrive(0.25);
#endif
        }
        Console.WriteLine("PASS: golden direct/front field, zero depth, depth-edge topology, finite complex pressure, and reflection symmetry (all octaves).");
#if BASELINE
        Console.WriteLine("PASS: complex rear-edge drive unchanged from the saved working baseline (all octaves).");
#endif
        double maxAngularChange = 0, maxDepthJump = 0;
        using (var polar = new StreamWriter(Path.Combine(output, "polar.csv")))
        using (var depth = new StreamWriter(Path.Combine(output, "depth.csv")))
        {
            polar.WriteLine("octave,depth_m,azimuth_deg,angle_deg,real,imaginary,scaled_magnitude");
            depth.WriteLine("octave,depth_m,angle_deg,real,imaginary,scaled_magnitude");
            for (int octave = 0; octave < 8; octave++)
            {
                foreach (double cabinetDepth in new[]{0.10, 0.25, 0.50})
                {
                    var field = new Field(cabinetDepth, Offset, octave);
                    foreach (double azimuth in new[]{0.0, 45.0, 90.0, 135.0})
                        for (int angle = 0; angle <= 720; angle++)
                        {
                            var p = field.At(Polar(angle * 0.5, azimuth));
                            polar.WriteLine(FormattableString.Invariant($"{octave},{cabinetDepth},{azimuth},{angle * 0.5},{p.Real:R},{p.Imaginary:R},{Distance * p.Magnitude:R}"));
                        }
                    //Convergence at the old side-face visibility boundary in the forward hemisphere.
                    double grazing = Math.Asin((0.16 - Offset.x) / Distance) * 180 / Math.PI;
                    double broad = (field.At(Polar(grazing + 1E-4)) - field.At(Polar(grazing - 1E-4))).Magnitude;
                    double narrow = (field.At(Polar(grazing + 1E-5)) - field.At(Polar(grazing - 1E-5))).Magnitude;
                    Require(narrow <= 0.12 * broad + 1E-14, "angular grazing continuity");
                    maxAngularChange = Math.Max(maxAngularChange, narrow * Distance);
                }
                for (int step = 0; step <= 100; step++)
                {
                    double cabinetDepth = 0.05 + step * 0.005;
                    var field = new Field(cabinetDepth, Offset, octave);
                    foreach (double angle in new[]{90.0, 135.0, 180.0})
                    {
                        var p = field.At(Polar(angle, 35));
                        depth.WriteLine(FormattableString.Invariant($"{octave},{cabinetDepth},{angle},{p.Real:R},{p.Imaginary:R},{Distance * p.Magnitude:R}"));
                    }
                }
                //Check re-discretization seams, rather than merely sampling between them.
                double spacing = Math.Min(0.020, 343 / (62.5 * Math.Pow(2, octave)) / 8);
                foreach (double cabinetDepth in new[]{10 * spacing, 20 * spacing})
                {
                    var before = new Field(cabinetDepth - 1E-8, Offset, octave);
                    var after = new Field(cabinetDepth + 1E-8, Offset, octave);
                    double jump = Distance * (before.At(Polar(135, 35)) - after.At(Polar(135, 35))).Magnitude;
                    maxDepthJump = Math.Max(maxDepthJump, jump);
                    Require(jump < 0.005, "depth discretization jump below 0.5% of unit reference pressure");
                }
            }
        }
        Console.WriteLine($"PASS: grazing angular convergence; largest scaled change across 0.00002 degrees = {maxAngularChange:G4}.");
        Console.WriteLine($"PASS: depth sweep 0.05..0.55 m; largest tested subdivision jump = {maxDepthJump:G4} of unit reference pressure.");
        var cabinet = new Cabinet_Diffraction(0.32, 0.8, 0.25);
        Near(cabinet.Pressure(Offset, Polar(135), -1, Diameter, Distance), cabinet.Pressure(Offset, Polar(135), 0, Diameter, Distance), "lower octave clamp");
        Near(cabinet.Pressure(Offset, Polar(135), 8, Diameter, Distance), cabinet.Pressure(Offset, Polar(135), 7, Diameter, Distance), "upper octave clamp");
        Near(cabinet.Pressure(Offset, new Vector(0, 0, 0), 0, Diameter, Distance), Complex.Zero, "zero direction");
        string[] balloons = cabinet.Driver_Balloon(Offset, Diameter);
        Require(balloons.Length == 8, "eight balloon strings");
        foreach (string balloon in balloons)
        {
            string[] rows = balloon.Split(';', StringSplitOptions.RemoveEmptyEntries);
            Require(rows.Length == 72, "72 balloon meridians");
            foreach (string row in rows)
            {
                double[] values = row.Split(' ').Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
                Require(values.Length == 37 && values.All(v => double.IsFinite(v) && v >= 0 && v <= 60), "37 bounded balloon samples");
            }
        }
        Console.WriteLine("PASS: public octave clamps, zero direction, and 8 x 72 x 37 balloon string format.");
        Console.WriteLine("Diagnostic CSVs: " + output);
    }
}
