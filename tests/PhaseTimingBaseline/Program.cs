// Ticket 02 (phase timing): headless fixed scenario for the BEFORE/AFTER baseline.
// Stages mirror the Hybrid/Mapping per-source summary: prepare/construct (receivers vs
// tracer+pools), trace work, filter, end-to-end wall. Fixed-count ray tracing is driven
// through SplitRayTracer.Calculate with fixed worker seeds so results are comparable
// across builds (Begin() seeds from the wall clock).
// Run: dotnet run --project tests/PhaseTimingBaseline [-- <fixture-file>]
// Exit code 0 = pass. First run records the fixture; later runs compare (tolerance
// 1e-12 + 1e-9 * max, per spec) and fail on drift.

using System.Diagnostics;
using System.Text;
using Hare.Geometry;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;
using Material = Pachyderm_Acoustic.Environment.Material;

const int SOURCES = 3, RECEIVERS = 4, SAMPLE_RATE = 44100, CUTOFF_MS = 500, RAY_COUNT = 2000, BANDS = 8;
int workers = Environment.ProcessorCount;
string fixturePath = args.Length > 0 ? args[0] : "phase-timing-fixture.txt";
// The Rhino UI assigns the DSP system before calculating; default to minimum phase here.
Pachyderm_Acoustic.Audio.Pach_SP.Initialize_FFTW();
Pachyderm_Acoustic.Audio.Pach_SP.Filter = new Pachyderm_Acoustic.Audio.Pach_SP.Minimum_Phase_System();

Console.WriteLine($"scenario: sources={SOURCES} receivers={RECEIVERS} sample_rate={SAMPLE_RATE} cutoff_ms={CUTOFF_MS} bands=octave({BANDS}) ray_count={RAY_COUNT} mode=fixed-count workers={workers} seeds=0x20251026+p runtime={Environment.Version} os={Environment.OSVersion}");

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

Stopwatch batch = Stopwatch.StartNew();
Polygon_Scene scene = new Polygon_Scene(20, 50, 101325, 0, false, true);
scene.Construct(
    new Point[][][] { new Point[][] { Face(0,0,0, 8,0,0, 8,6,0, 0,6,0), Face(0,0,4, 0,6,4, 8,6,4, 8,0,4),
                                      Face(0,0,0, 0,0,4, 8,0,4, 8,0,0), Face(0,6,0, 8,6,0, 8,6,4, 0,6,4),
                                      Face(0,0,0, 0,6,0, 0,6,4, 0,0,4), Face(8,0,0, 8,0,4, 8,6,4, 8,6,0) } },
    new List<Material> { new Basic_Material(new double[] { 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26 }) },
    new List<Scattering> { new Lambert_Scattering(new double[] { 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12 }) },
    new List<double[]> { new double[8] });
List<Point> partitionPts = new List<Point>(receiverPts);
foreach (RandomSource s in sources) s.AppendPts(ref partitionPts);
scene.partition(partitionPts);
batch.Stop();
Console.WriteLine($"batch prep (scene construct + partition): {batch.Elapsed.TotalSeconds:0.###} s");

var fixture = new StringBuilder();
fixture.AppendLine($"workers={workers} rays={RAY_COUNT}");
for (int s = 0; s < SOURCES; s++)
{
    Stopwatch e2e = Stopwatch.StartNew();
    Stopwatch sw = Stopwatch.StartNew();
    Receiver_Bank rec = new Receiver_Bank(receiverPts, sources[s], scene, SAMPLE_RATE, CUTOFF_MS, Receiver_Bank.Type.Stationary, false);
    sw.Stop();
    TimeSpan prepareRecv = sw.Elapsed;

    sw = Stopwatch.StartNew();
    Direct_Sound DS = new Direct_Sound(sources[s], rec, scene, Enumerable.Range(0, BANDS).ToArray());
    sw.Stop();
    TimeSpan dsCtor = sw.Elapsed;
    sw = Stopwatch.StartNew();
    await Run_To_Completion(DS);
    sw.Stop();
    TimeSpan dsWork = sw.Elapsed;

    sw = Stopwatch.StartNew();
    SplitRayTracer RT = new SplitRayTracer(sources[s], rec, scene, (CUTOFF_MS / 1000.0) * scene.Sound_speed(0), new int[] { 0, BANDS - 1 }, 0, RAY_COUNT, null);
    sw.Stop();
    TimeSpan rtCtor = sw.Elapsed;
    sw = Stopwatch.StartNew();
    List<Task> tasks = new List<Task>();
    for (int p = 0; p < workers; p++)
    {
        int start = (int)Math.Floor((double)p * RAY_COUNT / workers);
        int end = p == workers - 1 ? RAY_COUNT : (int)Math.Floor((double)(p + 1) * RAY_COUNT / workers);
        Calc_Params cp = new Calc_Params(scene, start, end, p, 0x20251026 + p);
        tasks.Add(Task.Run(async () => await RT.Calculate(cp)));
    }
    await Task.WhenAll(tasks);
    RT.Combine_ThreadLocal_Results();
    sw.Stop();
    TimeSpan rtWork = sw.Elapsed;

    sw = Stopwatch.StartNew();
    DS.Create_Filter();
    rec.Create_Filter(null);
    sw.Stop();
    TimeSpan filter = sw.Elapsed;
    e2e.Stop();

    Console.WriteLine($"source {s}/{SOURCES - 1} end-to-end wall: {e2e.Elapsed.TotalSeconds:0.###} s (includes stages below; do not add to them)");
    Console.WriteLine($"  prepare receiver bank: {prepareRecv.TotalSeconds:0.###} s");
    Console.WriteLine($"  direct sound: construct {dsCtor.TotalSeconds:0.###} s; work {dsWork.TotalSeconds:0.###} s");
    Console.WriteLine($"  ray trace: construct (tracer + pools) {rtCtor.TotalSeconds:0.###} s; work {rtWork.TotalSeconds:0.###} s");
    Console.WriteLine($"  filter (direct + receiver IRs): {filter.TotalSeconds:0.###} s");
    Console.WriteLine($"  image source / IS-trace: not run (disabled in this scenario)");
    Console.WriteLine($"  rays={RT._currentRay[0]} sub-rays={RT._rayTotal[0]} lost={RT.PercentLost:0.##}%"); // library _ts invalid here: Begin() is bypassed for fixed seeds

    Assert(RT.PercentLost < 100, $"source {s}: all energy lost - scene winding is wrong");
    for (int r = 0; r < RECEIVERS; r++)
    {
        for (int oct = 0; oct < BANDS; oct++)
        {
            double sum = 0;
            for (int t = 0; t < rec.SampleCT; t++) sum += rec.Rec_List[r].Energy(t, oct);
            Assert(!double.IsNaN(sum) && !double.IsInfinity(sum), $"source {s} rec {r} band {oct}: NaN/Inf energy");
            fixture.AppendLine($"S{s} RT R{r} B{oct} {sum:R}");
        }
        for (int oct = 0; oct < BANDS; oct++) fixture.AppendLine($"S{s} DS R{r} B{oct} {DS.EnergyValue(oct, r).Sum():R}");
    }
    double total = 0;
    for (int r = 0; r < RECEIVERS; r++) for (int oct = 0; oct < BANDS; oct++) total += rec.Rec_List[r].Energy(rec.SampleCT / 2, oct);
    Assert(total > 0, $"source {s}: no energy accumulated");
}

Compare_Or_Record_Fixture(fixturePath, fixture.ToString());
Console.WriteLine("PASS: stage timings above; acoustic fixture matches (or was recorded).");
return 0;

static Point[] Face(params double[] xyz)
{
    // 4 corners given as x,y,z triplets: c0,c1,c2 (cross product of c1-c0 and c2-c0 is the interior normal), c3.
    return new Point[] { new Point(xyz[0], xyz[1], xyz[2]), new Point(xyz[3], xyz[4], xyz[5]), new Point(xyz[6], xyz[7], xyz[8]), new Point(xyz[9], xyz[10], xyz[11]) };
}

static async Task Run_To_Completion(Simulation_Type sim)
{
    // Replicates RCPachTools.RunSimulation coordination without Rhino.
    sim.Begin();
    while (sim.ThreadState() != System.Threading.ThreadState.Stopped) await Task.Delay(25);
    sim.Combine_ThreadLocal_Results();
    while (sim.ThreadState() == System.Threading.ThreadState.Running) await Task.Delay(25);
}

static void Compare_Or_Record_Fixture(string path, string current)
{
    if (!File.Exists(path))
    {
        File.WriteAllText(path, current);
        Console.WriteLine($"fixture recorded to {Path.GetFullPath(path)} (first run).");
        return;
    }
    string[] oldLines = File.ReadAllLines(path);
    string[] newLines = current.Split('\n', StringSplitOptions.RemoveEmptyEntries);
    Assert(oldLines.Length == newLines.Length, $"fixture line count changed: {oldLines.Length} -> {newLines.Length}");
    for (int i = 0; i < oldLines.Length; i++)
    {
        string o = oldLines[i].TrimEnd('\r'), n = newLines[i].TrimEnd('\r');
        if (o.StartsWith("workers=")) { Assert(o == n, $"worker count changed: {o} -> {n}"); continue; }
        string[] op = o.Split(' '), np = n.Split(' ');
        Assert(op.Length == np.Length && op[0] == np[0] && op[1] == np[1] && op[2] == np[2] && op[3] == np[3], $"fixture key changed: {o} -> {n}");
        double a = double.Parse(op[4]), b = double.Parse(np[4]);
        Assert(Math.Abs(a - b) <= 1e-12 + 1e-9 * Math.Max(Math.Abs(a), Math.Abs(b)), $"acoustic result drifted: {o} -> {n}");
    }
    Console.WriteLine("fixture compared: identical within tolerance.");
}

static void Assert(bool condition, string message)
{
    if (!condition) throw new Exception(message);
}
