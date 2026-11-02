// Ticket 09 (Deliverable A): T6 fixed-count BEFORE/AFTER regression — shared scenario.
// This file is compiled BOTH by tests/T6Baseline (standalone runner, also copied into the
// baseline worktree rt-09-before to produce the BEFORE fixture from commit b6490b5) and by
// tests/AnalysisCancellation (T6RegressionSuite). The scenario, worker split, seeds and
// fixture format are therefore byte-identical on both sides by construction.
//
// Fixed per spec T6: scene, worker count (Environment.ProcessorCount), per-worker seeds
// (T6Shared.BASE_SEED + p), disjoint ray ranges, source initial state, total ray count.
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;


using System.Text;
using Hare.Geometry;
using Material = Pachyderm_Acoustic.Environment.Material;

namespace RayTraceTests;

internal static class T6Scenario
{
    public const int SOURCES = 3, RECEIVERS = 4, SAMPLE_RATE = 44100, CUTOFF_MS = 500, BANDS = 8;

    /// <summary>Runs the fixed scenario once and returns the fixture text.</summary>
    public static async Task<string> Run()
    {
        int workers = Environment.ProcessorCount;
        Pachyderm_Acoustic.Audio.Pach_SP.Initialize_FFTW();
        Pachyderm_Acoustic.Audio.Pach_SP.Filter = new Pachyderm_Acoustic.Audio.Pach_SP.Minimum_Phase_System();

        Point[] receiverPts = new Point[RECEIVERS]
        {
            new Point(2, 4, 1.2), new Point(4, 3, 1.2), new Point(6, 4, 1.2), new Point(4, 1.5, 1.2),
        };
        double[] swl = new double[8] { 100, 101, 102, 103, 104, 105, 106, 107 };
        RandomSource[] sources = new RandomSource[SOURCES]
        {
            new RandomSource(swl, new Point(2, 2, 1.7), 0, false),
            new RandomSource(swl, new Point(6, 2, 1.7), 1, false),
            new RandomSource(swl, new Point(4, 5, 1.7), 2, false),
        };

        Polygon_Scene scene = new Polygon_Scene(20, 50, 101325, 0, false, true);
        scene.Construct(
            new Point[][][] { new Point[][] { T6Shared.Face(0,0,0, 8,0,0, 8,6,0, 0,6,0), T6Shared.Face(0,0,4, 0,6,4, 8,6,4, 8,0,4),
                                              T6Shared.Face(0,0,0, 0,0,4, 8,0,4, 8,0,0), T6Shared.Face(0,6,0, 8,6,0, 8,6,4, 0,6,4),
                                              T6Shared.Face(0,0,0, 0,6,0, 0,6,4, 0,0,4), T6Shared.Face(8,0,0, 8,0,4, 8,6,4, 8,6,0) } },
            new List<Material> { new Basic_Material(new double[] { 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26 }) },
            new List<Scattering> { new Lambert_Scattering(new double[] { 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12 }) },
            new List<double[]> { new double[8] });
        List<Point> partitionPts = new List<Point>(receiverPts);
        foreach (RandomSource s in sources) s.AppendPts(ref partitionPts);
        scene.partition(partitionPts);

        var fixture = new StringBuilder();
        fixture.AppendLine($"workers={workers} rays={T6Shared.RAY_COUNT} sources={SOURCES} receivers={RECEIVERS} bands={BANDS}");
        for (int s = 0; s < SOURCES; s++)
        {
            Receiver_Bank rec = new Receiver_Bank(receiverPts, sources[s], scene, SAMPLE_RATE, CUTOFF_MS, Receiver_Bank.Type.Stationary, false);
            Direct_Sound DS = new Direct_Sound(sources[s], rec, scene, Enumerable.Range(0, BANDS).ToArray());
            await T6Shared.Run_To_Completion(DS);
            SplitRayTracer RT = new SplitRayTracer(sources[s], rec, scene, (CUTOFF_MS / 1000.0) * scene.Sound_speed(0), new int[] { 0, BANDS - 1 }, 0, T6Shared.RAY_COUNT, null);
            List<Task> tasks = new List<Task>();
            for (int p = 0; p < workers; p++)
            {
                int start = (int)Math.Floor((double)p * T6Shared.RAY_COUNT / workers);
                int end = p == workers - 1 ? T6Shared.RAY_COUNT : (int)Math.Floor((double)(p + 1) * T6Shared.RAY_COUNT / workers);
                Calc_Params cp = new Calc_Params(scene, start, end, p, T6Shared.BASE_SEED + p);
                tasks.Add(Task.Run(async () => await RT.Calculate(cp)));
            }
            await Task.WhenAll(tasks);
            RT.Combine_ThreadLocal_Results();
            DS.Create_Filter();
            rec.Create_Filter(null);

            T6Shared.Assert(RT.PercentLost < 100, $"source {s}: all energy lost");
            T6Shared.Assert(RT._currentRay[0] == T6Shared.RAY_COUNT, $"source {s}: traced {RT._currentRay[0]} rays, expected {T6Shared.RAY_COUNT}");
            fixture.AppendLine($"S{s} SUBRAYS {RT._rayTotal[0]:R}");
            fixture.AppendLine($"S{s} SRCPOWER {string.Join(" ", sources[s].SWL().Select(v => v.ToString("R")))}");
            fixture.AppendLine($"S{s} BANDS {string.Join(" ", rec.SWL.Select(v => v.ToString("R")))}");
            for (int r = 0; r < RECEIVERS; r++)
            {
                Spherical_Receiver R = rec.Rec_List[r];
                for (int oct = 0; oct < BANDS; oct++)
                {
                    double sum = 0;
                    for (int t = 0; t < rec.SampleCT; t++) sum += R.Energy(t, oct);
                    T6Shared.Assert(!double.IsNaN(sum) && !double.IsInfinity(sum), $"S{s} R{r} B{oct}: NaN/Inf energy");
                    fixture.AppendLine($"S{s} RT R{r} B{oct} E {sum:R}");
                }
                for (int oct = 0; oct < BANDS; oct++) fixture.AppendLine($"S{s} DS R{r} B{oct} E {DS.EnergyValue(oct, r).Sum():R}");
                double[] filter = R.Create_Filter(swl, 0);
                for (int oct = 0; oct < BANDS; oct++)
                {
                    double sum = 0;
                    for (int t = 0; t < filter.Length; t++) sum += filter[t] * filter[t];
                    T6Shared.Assert(!double.IsNaN(sum) && !double.IsInfinity(sum), $"S{s} R{r} B{oct}: NaN/Inf filter");
                    fixture.AppendLine($"S{s} RT R{r} B{oct} P {sum:R}");
                }
                for (int cat = 0; cat < 3; cat++)
                    for (int oct = 0; oct < BANDS; oct++)
                    {
                        double sum = 0;
                        for (int t = 0; t < R.Reflection_Analysis_Data[cat][oct].Length; t++) sum += R.Reflection_Analysis_Data[cat][oct][t];
                        T6Shared.Assert(!double.IsNaN(sum) && !double.IsInfinity(sum), $"S{s} R{r} B{oct} RAD{cat}: NaN/Inf");
                        fixture.AppendLine($"S{s} RT R{r} B{oct} RAD{cat} {sum:R}");
                    }
                // Direction histograms (T6 "方向"): per band, sum of squared components of the
                // + and - axis directions over all samples; identical formatting on both sides.
                for (int oct = 0; oct < BANDS; oct++)
                {
                    double pos = 0, neg = 0;
                    for (int t = 0; t < rec.SampleCT; t++)
                    {
                        Hare.Geometry.Vector p = R.Directions_Pos(oct, t), n = R.Directions_Neg(oct, t);
                        pos += p.dx * p.dx + p.dy * p.dy + p.dz * p.dz;
                        neg += n.dx * n.dx + n.dy * n.dy + n.dz * n.dz;
                    }
                    fixture.AppendLine($"S{s} RT R{r} B{oct} DIR {pos:R} {neg:R}");
                }
            }
        }
        return fixture.ToString();
    }
}

/// <summary>Element-wise fixture comparison shared by the runner and the harness suite.</summary>
internal static class T6Compare
{
    public static void Files(string referencePath, string resultPath, string label) => Files(referencePath, resultPath, label, radReportOnly: false);

    public static void Files(string referencePath, string resultPath, string label, bool radReportOnly)
    {
        string[] refLines = File.ReadAllLines(referencePath);
        string[] newLines = File.ReadAllLines(resultPath);
        T6Shared.Assert(refLines.Length == newLines.Length, $"{label}: line count {refLines.Length} -> {newLines.Length}");
        int compared = 0, radReported = 0;
        for (int i = 0; i < refLines.Length; i++)
        {
            string o = refLines[i].TrimEnd('\r'), n = newLines[i].TrimEnd('\r');
            string[] op = o.Split(' '), np = n.Split(' ');
            T6Shared.Assert(op.Length == np.Length && string.Equals(string.Join(" ", op.Take(Math.Min(4, op.Length))), string.Join(" ", np.Take(Math.Min(4, np.Length))), StringComparison.Ordinal),
                $"{label}: fixture key changed: {o} -> {n}");
            for (int f = Math.Min(4, op.Length); f < op.Length; f++)
            {
                // "sub-rays" counts are integers and must match exactly; floating tokens
                // are compared with tolerance below.
                if (!double.TryParse(op[f], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double a)
                    || !double.TryParse(np[f], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double b)
                    || (op[f].IndexOfAny(new[] { '.', 'E', 'e' }) < 0 && np[f].IndexOfAny(new[] { '.', 'E', 'e' }) < 0))
                {
                    T6Shared.Assert(string.Equals(op[f], np[f], StringComparison.Ordinal), $"{label}: fixture text mismatch: {o} -> {n}");
                    continue;
                }
                bool radLine = radReportOnly && np.Length >= 2 && np[1] == "RT" && np.Length > 4 && np[4].StartsWith("RAD", StringComparison.Ordinal);
                if (!double.IsNaN(a)) T6Shared.Assert(!double.IsNaN(b) && !double.IsInfinity(b), $"{label}: newly appearing NaN/Infinity: {o} -> {n}");
                if (radLine)
                {
                    // Baseline RAD is the racy pre-D5 output (spec D5): report the delta
                    // instead of failing against a nondeterministic reference.
                    if (a != b) radReported++;
                    continue;
                }
                T6Shared.Assert(Math.Abs(a - b) <= 1e-12 + 1e-9 * Math.Max(Math.Abs(a), Math.Abs(b)),
                    $"{label}: T6 regression at {o} -> {n}");
                compared++;
            }
        }
        Console.WriteLine(radReported > 0
            ? $"  T6 compare ({label}): {compared} elements within tolerance; {radReported} RAD element(s) differ from the racy baseline (expected, spec D5)"
            : $"  T6 compare ({label}): {compared} numeric elements identical within tolerance (1e-12 + 1e-9*max)");
    }
}
