using System.Collections.Concurrent;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Ticket 05 (spec D2, matrix T1/T2): pool lifecycle. T1 — one thread-safe lazy
/// creation path per pool: concurrent first access yields one complete instance,
/// repeated <see cref="BroadRayPool.Initialize"/> only ensures existence, a borrowed
/// ray survives re-initialization, borrows pair with releases, and an in-use ray is
/// never re-rented (that last failure would be an independent correctness blocker —
/// the ring-rent/anonymous-release algorithm is intentionally NOT changed here).
/// T2 — two successive sources with different spectra, cold and warm pool, 8 and 24
/// bands; the later source never overwrites the earlier source's recorded results.
/// </summary>
internal static class PoolInitSuite
{
    public static void Register()
    {
        Harness.Suite("pool lazy init: one idempotent instance (T1)", () => { Lifecycle(); return Task.CompletedTask; });
        Harness.Suite("successive sources, cold/warm pool, 8+24 bands (T2)", SuccessiveSources);
    }

    // ---------------- T1 ----------------

    private static void Lifecycle()
    {
        ConcurrentFirstAccess();
        RepeatedInitializeKeepsInstanceAndBorrowedRays();
        GeometryPoolCaptureStaysLive();
        ConcurrentBorrowReleaseBalanced();
    }

    /// <summary>Many threads hitting Initialize()/Instance simultaneously must all see the same, complete instance.</summary>
    private static void ConcurrentFirstAccess()
    {
        const int threads = 8;
        var start = new Barrier(threads);
        var broad = new BroadRayPool[threads];
        var oct = new OctaveRayPool[threads];
        var workers = new Task[threads];
        for (int t = 0; t < threads; t++)
        {
            int i = t;
            workers[i] = Task.Run(() =>
            {
                start.SignalAndWait();
                if ((i & 1) == 0) { BroadRayPool.Initialize(); OctaveRayPool.Initialize(); }
                broad[i] = BroadRayPool.Instance;
                oct[i] = OctaveRayPool.Instance;
            });
        }
        Task.WaitAll(workers);

        for (int i = 1; i < threads; i++)
        {
            Harness.Assert(ReferenceEquals(broad[0], broad[i]),
                $"concurrent BroadRayPool access returned different instances to threads 0 and {i}");
            Harness.Assert(ReferenceEquals(oct[0], oct[i]),
                $"concurrent OctaveRayPool access returned different instances to threads 0 and {i}");
        }

        // "Complete" = the borrow entry points work on the instance every thread saw.
        BroadRay br = broad[0].pull;
        broad[0].release();
        OctaveRay or = oct[0].pull;
        oct[0].release(or);
    }

    /// <summary>Initialize() may not replace an existing pool or disturb rays still borrowed from it.</summary>
    private static void RepeatedInitializeKeepsInstanceAndBorrowedRays()
    {
        BroadRayPool bp = BroadRayPool.Instance;
        OctaveRayPool op = OctaveRayPool.Instance;

        BroadRay heldBroad = bp.pull;
        heldBroad.Ray_ID = 4242;
        heldBroad.t_sum = 7.5;
        OctaveRay heldOct = op.pull;
        heldOct.Intensity = 3.25;
        heldOct.Octave = 5;

        BroadRayPool.Initialize();
        OctaveRayPool.Initialize();
        BroadRayPool.Initialize();
        OctaveRayPool.Initialize();

        Harness.Assert(ReferenceEquals(bp, BroadRayPool.Instance),
            "BroadRayPool.Initialize() replaced an existing pool (must only ensure existence)");
        Harness.Assert(ReferenceEquals(op, OctaveRayPool.Instance),
            "OctaveRayPool.Initialize() replaced an existing pool (must only ensure existence)");
        Harness.Assert(heldBroad.Ray_ID == 4242 && heldBroad.t_sum == 7.5,
            "Initialize() mutated a BroadRay still borrowed from the pool");
        Harness.Assert(heldOct.Intensity == 3.25 && heldOct.Octave == 5,
            "Initialize() mutated an OctaveRay still borrowed from the pool");

        // The held rays were not silently re-rented: releasing them keeps the pools responsive.
        bp.release();
        op.release(heldOct);
        BroadRay r2 = bp.pull;
        bp.release();
        OctaveRay o2 = op.pull;
        op.release(o2);
    }

    /// <summary>GeometryPool captures Instance at type-init; after re-initialization the capture must still be the live pool.</summary>
    private static void GeometryPoolCaptureStaysLive()
    {
        BroadRayPool captured = PoolFixtures.SharedBroad;
        OctaveRayPool capturedOct = PoolFixtures.SharedOctave;
        BroadRayPool.Initialize();
        OctaveRayPool.Initialize();
        Harness.Assert(ReferenceEquals(captured, BroadRayPool.Instance),
            "GeometryPool.MainPool no longer refers to the live BroadRayPool after Initialize()");
        Harness.Assert(ReferenceEquals(capturedOct, OctaveRayPool.Instance),
            "GeometryPool.OctPool no longer refers to the live OctaveRayPool after Initialize()");
    }

    /// <summary>
    /// Existing-concurrency borrow/release churn with identity tracking: every pull returns a
    /// ray nobody else holds, and all borrows are paired with releases so both pools stay usable.
    /// A failure here means an in-use ray was re-rented before release — an independent
    /// correctness blocker (ticket 05 explicitly leaves the lease algorithm untouched).
    /// </summary>
    private static void ConcurrentBorrowReleaseBalanced()
    {
        const int workers = 8, iterations = 250;
        var activeBroad = new ConcurrentDictionary<BroadRay, int>();
        var activeOct = new ConcurrentDictionary<OctaveRay, int>();
        var tasks = new Task[workers];
        for (int w = 0; w < workers; w++)
        {
            tasks[w] = Task.Run(() =>
            {
                for (int i = 0; i < iterations; i++)
                {
                    BroadRay b = BroadRayPool.Instance.pull;
                    Harness.Assert(activeBroad.TryAdd(b, i),
                        "BLOCKER (report separately, do not fix the lease): BroadRay re-rented while still borrowed");
                    b.t_sum = i;
                    Harness.Assert(activeBroad.TryRemove(b, out _), "broad identity bookkeeping failure");
                    BroadRayPool.Instance.release();

                    OctaveRay o = OctaveRayPool.Instance.pull;
                    Harness.Assert(activeOct.TryAdd(o, i),
                        "BLOCKER (report separately, do not fix the lease): OctaveRay re-rented while still borrowed");
                    o.Intensity = i;
                    Harness.Assert(activeOct.TryRemove(o, out _), "octave identity bookkeeping failure");
                    OctaveRayPool.Instance.release(o);
                }
            });
        }
        Task.WaitAll(tasks);

        // Pools still answer borrow/release after all churn (no leaked semaphore count).
        BroadRay fin = BroadRayPool.Instance.pull;
        BroadRayPool.Instance.release();
        OctaveRay fin2 = OctaveRayPool.Instance.pull;
        OctaveRayPool.Instance.release(fin2);
    }

    // ---------------- T2 ----------------

    private static async Task SuccessiveSources()
    {
        await RunSuccessive(thirdOctave: false); // 8-band octave
        await RunSuccessive(thirdOctave: true);  // 24-band third-octave
    }

    private static async Task RunSuccessive(bool thirdOctave)
    {
        int bands = thirdOctave ? 24 : 8;
        string mode = thirdOctave ? "24-band" : "8-band";

        SplitRayTracer cold = await RunSource(0, thirdOctave);      // first (cold-pool) source A
        SplitRayTracer warm = await RunSource(0, thirdOctave);      // same source again, warm pool
        SplitRayTracer next = await RunSource(10, thirdOctave);     // source B: +10 dB in every band

        double[][] coldSnap = Snapshot(cold, bands);
        double[][] warmSnap = Snapshot(warm, bands);
        double[][] nextSnap = Snapshot(next, bands);

        double total = coldSnap.Sum(h => h.Sum());
        Harness.Assert(total > 0, $"{mode}: no receiver energy recorded; comparison would be vacuous");

        // Warm pool reproduces the cold-pool run exactly (no cross-run pool contamination).
        AssertEqual(coldSnap, warmSnap, $"{mode}: warm-pool rerun differs from cold-pool run");

        // B's spectrum is +10 dB per band -> ~10x the energy per recorded sample.
        for (int row = 0; row < coldSnap.Length; row++)
            for (int s = 0; s < coldSnap[row].Length; s++)
            {
                double a = coldSnap[row][s], b = nextSnap[row][s];
                if (a == 0) Harness.Assert(b == 0, $"{mode}: row {row} sample {s}: source B recorded {b} where A recorded 0");
                else Harness.Assert(Math.Abs(b - 10 * a) <= 1e-9 * 10 * a,
                    $"{mode}: row {row} sample {s}: +10 dB source gave {b}, expected 10x A's {a}");
            }

        // Source A's already-recorded results are untouched by source B's run.
        AssertEqual(coldSnap, Snapshot(cold, bands), $"{mode}: source A's results changed after source B ran");
    }

    /// <summary>Fresh scene/source/receiver/tracer + seeded fixed-count Calculate, mirroring <see cref="Fixtures.RunSeeded"/> with a spectrum offset.</summary>
    private static async Task<SplitRayTracer> RunSource(double dbOffset, bool thirdOctave)
    {
        var room = Fixtures.Scene();
        var src = thirdOctave ? Fixtures.ThirdOctaveSource(dbOffset) : Fixtures.Source(dbOffset);
        var tracer = Fixtures.Tracer(src, Fixtures.ReceiverBank(src, room, thirdOctave), room);
        foreach (var cp in Fixtures.WorkerParams(room)) await tracer.Calculate(cp);
        return tracer;
    }

    /// <summary>Deep copy of every receiver's per-band energy histograms: receiver-major, [r * bands + band].</summary>
    private static double[][] Snapshot(SplitRayTracer tracer, int bands)
    {
        var recs = tracer.GetReceiver.Rec_List;
        var snap = new double[recs.Length * bands][];
        for (int r = 0; r < recs.Length; r++)
            for (int b = 0; b < bands; b++)
                snap[r * bands + b] = (double[])recs[r].GetEnergyHistogram(b).Clone();
        return snap;
    }

    private static void AssertEqual(double[][] expected, double[][] actual, string what)
    {
        Harness.Assert(expected.Length == actual.Length, $"{what}: row counts {expected.Length} vs {actual.Length}");
        for (int row = 0; row < expected.Length; row++)
        {
            Harness.Assert(expected[row].Length == actual[row].Length, $"{what}: row {row} lengths differ");
            for (int s = 0; s < expected[row].Length; s++)
                Harness.Assert(expected[row][s] == actual[row][s],
                    $"{what}: row {row} sample {s}: {expected[row][s]} vs {actual[row][s]}");
        }
    }
}
