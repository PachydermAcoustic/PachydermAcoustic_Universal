using Hare.Geometry;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Deterministic fixtures shared by every suite (spec § Testing Decisions, matrix T1–T8).
/// Fixed scene, explicit-seed ray-tracing entry via <see cref="Calc_Params"/>, direct
/// receiver construction, pool identity access, and a fixed worker-count / per-worker-seed /
/// ray-range split. Only public production surface is used — no Rhino, no private reflection.
/// </summary>
internal static class Fixtures
{
    public const int WorkerCount = 4;
    public const int RaysPerWorker = 500;
    public const int TotalRays = WorkerCount * RaysPerWorker;
    public const int BaseSeed = 11;          // worker w gets seed BaseSeed + w
    public const int SampleRate = 1000;      // receiver histogram sample rate (Hz)
    public const double CutOffTime_ms = 1000;// receiver cutoff time (ms)

    /// <summary>Fixed scene with no surfaces: every ray misses, so each ray's fate is fixed by its seed alone.</summary>
    public static Empty_Scene Scene() => new(20, 50, 1013.25, 0, false, false,
        new Point(-10, -10, -10), new Point(10, 10, 10));

    /// <summary>8-band octave spectrum in dB; offset shifts every band so T2 can use distinct spectra.</summary>
    public static double[] OctaveSpectrum(double dbOffset = 0) =>
        Enumerable.Range(0, 8).Select(b => 100.0 + b + dbOffset).ToArray();

    /// <summary>24-band third-octave spectrum in dB.</summary>
    public static double[] ThirdOctaveSpectrum(double dbOffset = 0) =>
        Enumerable.Range(0, 24).Select(b => 100.0 + b + dbOffset).ToArray();

    /// <summary>Pseudo-random-direction source at a fixed point; directions come only from the seeded Random.</summary>
    public static RandomSource Source(double dbOffset = 0, int id = 0) =>
        new(OctaveSpectrum(dbOffset), new Point(0, 0, 2), id, false);

    /// <summary>Third-octave variant of <see cref="Source(double, int)"/> for 24-band paths.</summary>
    public static RandomSource ThirdOctaveSource(double dbOffset = 0, int id = 0) =>
        new(ThirdOctaveSpectrum(dbOffset), new Point(0, 0, 2), id, true);

    public static readonly Point[] ReceiverPoints = { new(3, 1, 1.5), new(-2, 4, 1) };

    /// <summary>Direct receiver-bank construction over fixed points (no UI, no Rhino doc).</summary>
    public static Receiver_Bank ReceiverBank(Source src, Scene room, bool thirdOctave = false) =>
        new(ReceiverPoints, src, room, SampleRate, CutOffTime_ms, Receiver_Bank.Type.Stationary, thirdOctave);

    /// <summary>Direct single-receiver construction (bank-free path).</summary>
    public static Spherical_Receiver Receiver(Scene room, Point at, Point sourceAt) =>
        new(at, sourceAt, room.Attenuation(at, false), room.Sound_speed(at), room.Rho(at), SampleRate, CutOffTime_ms);

    /// <summary>
    /// Fixed worker split (T6): worker count, per-worker seed (BaseSeed + w) and disjoint
    /// [start, end) ray ranges over TotalRays, exactly like Begin() would divide them.
    /// </summary>
    public static Calc_Params[] WorkerParams(Scene room, int workerCount = WorkerCount,
        int totalRays = TotalRays, int baseSeed = BaseSeed)
    {
        var ps = new Calc_Params[workerCount];
        for (int w = 0; w < workerCount; w++)
            ps[w] = new Calc_Params(room, w * totalRays / workerCount, (w + 1) * totalRays / workerCount, w, baseSeed + w);
        return ps;
    }

    /// <summary>
    /// Fresh tracer for the fixed scene. imageSourceOrder = -1 so ray order 0 is already logged
    /// (this harness runs no image-source pass). Raycount &gt;= 1 keeps it a fixed-count calculation.
    /// </summary>
    public static SplitRayTracer Tracer(Source src, Receiver_Bank rec, Scene room, int rayCount = TotalRays,
        int bandLow = 0, int? bandHigh = null, int imageSourceOrder = -1) =>
        new(src, rec, room, CutOffTime_ms, new[] { bandLow, bandHigh ?? src.BandCount - 1 }, imageSourceOrder, rayCount, null);

    /// <summary>
    /// One fully deterministic fixed-count run: fresh scene/source/receiver/tracer, then
    /// <see cref="SplitRayTracer.Calculate(object)"/> awaited per worker sequentially
    /// (single-threaded accumulation, so repeat runs are bit-identical).
    /// </summary>
    public static async Task<SplitRayTracer> RunSeeded()
    {
        var room = Scene();
        var src = Source();
        var tracer = Tracer(src, ReceiverBank(src, room), room);
        foreach (var cp in WorkerParams(room)) await tracer.Calculate(cp);
        return tracer;
    }
}

/// <summary>
/// Pool identity accessors (T1). Suites compare object identities across these —
/// this file asserts nothing and changes no pool behavior.
/// </summary>
internal static class PoolFixtures
{
    public static BroadRayPool Broad => BroadRayPool.Instance;
    public static OctaveRayPool Octave => OctaveRayPool.Instance;
    /// <summary>Captured at GeometryPool's static initialization — always the live instance since ticket 05 made Initialize() idempotent.</summary>
    public static BroadRayPool SharedBroad => GeometryPool.MainPool;
    public static OctaveRayPool SharedOctave => GeometryPool.OctPool;
    public static void Reinitialize()
    {
        BroadRayPool.Initialize();
        OctaveRayPool.Initialize();
    }
}
