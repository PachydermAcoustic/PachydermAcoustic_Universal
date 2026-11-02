using Hare.Geometry;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Ticket 07 (spec D6, acceptance T7 + T8): a convergence run's task set describes only the
/// current round. T7 runs both convergence modes (minimum: Raycount &lt; 0, detailed:
/// Raycount == 0) over a real-energy fixture and asserts each round's check operates on the
/// bank accumulated by all previous rounds, the final round's contributions are counted, and
/// the run ends through the original convergence condition (no cancellation, no failure).
/// Convergence through the original condition is itself the accumulation proof: the minimum
/// check needs six consecutive rounds within 2%/2%/10% of the CUMULATIVE window sums, which
/// per-round-only data (relative drift O(1)) can never satisfy.
///
/// T8 asserts the lifecycle: no per-round growth of retained work after many rounds (live
/// bytes sampled right after each round's explicit GC.Collect stay flat — the pre-fix code
/// retained and re-awaited every historical round's tasks), cancellation and a controlled
/// round-task failure never declare success early, settle all started receiver writers, and
/// a failure is observable through SplitRayTracer.Failure once ThreadState() is Stopped.
/// Zero-energy fixtures keep the existing rules: the loop never self-converges and runs
/// until cancelled; the bank stays exactly empty.
/// </summary>
internal static class ConvergenceLifecycleSuite
{
    private const int ConvSeed = 4242;      // fixed master seed for Calculate_Conv (round/task seeds derive from it)
    private const int RoundsPerCheck = 80;  // Source.Revolution_Period default; asserted, not assumed
    private const double BoxAbsorption = 0.3;
    private const double BoxScattering = 0.15;
    private const int RetentionRounds = 25;
    private const int RetentionGrowthLimitBytes = 20_000; // pre-fix retention grew ~3 KB/round beyond this window
    private const string FailureMessage = "T8 injected shoot failure";

    public static void Register()
    {
        Harness.Suite("round task lifecycle: retention, cancellation, failure (T8)", RunT8);
        Harness.Suite("convergence rounds over accumulated data (T7: minimum + detailed)", RunT7);
    }

    private static async Task RunT7()
    {
        await RunConvergingScenario("minimum", rayCount: -1, minRounds: 30, checkPayload: true);
        await RunConvergingScenario("detailed", rayCount: 0, minRounds: 63, checkPayload: false);
    }

    private static async Task RunT8()
    {
        await RunRetentionScenario();      // also the zero-energy minimum-mode case
        await RunZeroEnergyScenario(-1);   // zero-energy minimum mode: no auto-converge
        await RunZeroEnergyScenario(0);    // detailed mode: no auto-converge rule added
        await RunCancellationScenario();
        await RunFailureScenario();
    }

    /// <summary>A real-energy shoebox run that must end through the original convergence condition.</summary>
    private static async Task RunConvergingScenario(string mode, int rayCount, int minRounds, bool checkPayload)
    {
        (SplitRayTracer tracer, ConvFeedback feedback, RandomSource src, Receiver_Bank bank) = ConvFixture(rayCount);

        await tracer.Calculate_Conv(4, ConvSeed);

        // Rounds come from the ray counter (each round traces exactly Revolution_Period
        // rays): the detailed check may return early without a Fill, so Payloads undercounts
        // rounds in that mode.
        int rounds = tracer._currentRay[0] / src.Revolution_Period();
        Harness.Assert(tracer._currentRay[0] % src.Revolution_Period() == 0,
            $"{mode}: ray count {tracer._currentRay[0]} is not a whole number of rounds");
        Harness.Assert(tracer.Failure == null, $"{mode}: run ended by task failure: {tracer.Failure?.Message}");
        Harness.Assert(tracer.ThreadState() == System.Threading.ThreadState.Stopped, $"{mode}: simulation not stopped");
        Harness.Assert(!feedback.Canceler.IsCancellationRequested, $"{mode}: run ended by cancellation, not convergence");
        Harness.Assert(rounds >= minRounds, $"{mode}: converged after only {rounds} rounds; accumulated-data convergence impossible that early");
        Harness.Assert(src.Revolution_Period() == RoundsPerCheck, $"{mode}: unexpected revolution period {src.Revolution_Period()}");
        Harness.Assert(tracer._currentRay[0] == rounds * RoundsPerCheck,
            $"{mode}: {tracer._currentRay[0]} rays traced over {rounds} rounds; expected {rounds * RoundsPerCheck}");

        // Fill-based evidence exists only where the check actually reports (minimum mode
        // and the non-final rounds that pass its early exits); detailed mode converges
        // through the same original condition without necessarily reporting.
        if (feedback.Totals.Count > 0)
        {
            for (int r = 1; r < feedback.Totals.Count; r++)
                Harness.Assert(feedback.Totals[r] > feedback.Totals[r - 1],
                    $"{mode}: round {r + 1} did not add to the bank accumulated by earlier rounds");
            double finalTotal = BankTotal(bank);
            Harness.Assert(finalTotal > feedback.Totals[feedback.Totals.Count - 1], $"{mode}: final round's contributions not counted into the bank");
        }

        if (checkPayload)
        {
            // The minimum check's Fill payload is {conv1, conv2, convinf}. The last six rounds
            // before convergence all passed the 2%/2%/10% thresholds, so the last three recorded
            // payloads must each be within the largest threshold.
            for (int p = feedback.Payloads.Count - 3; p < feedback.Payloads.Count; p++)
            {
                double drift = MaxFinite(feedback.Payloads[p]);
                Harness.Assert(drift <= 0.1, $"{mode}: round {p + 1} drift {drift} exceeds the 0.1 convergence threshold");
            }
        }
    }

    /// <summary>
    /// Zero-energy minimum-mode run across many rounds: the bank must stay exactly empty, the
    /// loop must keep running (no auto-converge) until cancelled, and the live heap sampled
    /// right after each round's explicit GC.Collect must not grow with the round count.
    /// </summary>
    private static async Task RunRetentionScenario()
    {
        (SplitRayTracer tracer, ConvFeedback feedback, _, Receiver_Bank bank) = ZeroEnergyFixture(-1);
        feedback.OnRound = r => { if (r == RetentionRounds) feedback.Canceler.Cancel(); };

        await tracer.Calculate_Conv(4, ConvSeed);

        Harness.Assert(feedback.Payloads.Count == RetentionRounds,
            $"zero-energy run did {feedback.Payloads.Count} rounds before cancellation; expected exactly {RetentionRounds} (no auto-converge)");
        Harness.Assert(tracer.Failure == null, $"zero-energy run failed: {tracer.Failure?.Message}");
        Harness.Assert(tracer.ThreadState() == System.Threading.ThreadState.Stopped, "zero-energy run not stopped after cancellation");
        Harness.Assert(BankTotal(bank) == 0, "zero-energy fixture recorded energy");

        // Measured after the fix: flat within ~2 KB from round 8 on (bounded IR queue in
        // the checks plus jitter). The pre-fix code retained every historical round's tasks
        // (~3 KB/round → ~50 KB by round 25 on top of this), so this limit discriminates.
        long lateGrowth = feedback.Heap[feedback.Heap.Count - 1] - feedback.Heap[7];
        Harness.Assert(lateGrowth < RetentionGrowthLimitBytes,
            $"retained work grew by {lateGrowth} bytes between round 8 and round {RetentionRounds}; round task sets must not accumulate");
    }

    /// <summary>
    /// Zero-energy run (both modes): the existing no-convergence rule is preserved — no
    /// auto-converge is added. Minimum mode reports each round via Fill (its check's drift
    /// values include finite numbers), so cancellation lands after an exact round count.
    /// Detailed mode's SpecMin check returns early on zero data (before its Fill), so a
    /// wall-clock timer cancels instead; the run must still stop promptly without declaring
    /// success and leave the bank exactly empty.
    /// </summary>
    private static async Task RunZeroEnergyScenario(int rayCount)
    {
        (SplitRayTracer tracer, ConvFeedback feedback, _, Receiver_Bank bank) = ZeroEnergyFixture(rayCount);
        int cancelAt = 5;
        if (rayCount < 0)
        {
            feedback.OnRound = r => { if (r == cancelAt) feedback.Canceler.Cancel(); };
        }
        else
        {
            // Wall-clock cancel: detailed mode's SpecMin check returns before its Fill on
            // zero data (conv == 0), so no round-count hook exists for this mode.
            using System.Threading.Timer timer = new System.Threading.Timer(
                _ => feedback.Canceler.Cancel(), null, 3000, Timeout.Infinite);
            await tracer.Calculate_Conv(4, ConvSeed);
            Harness.Assert(tracer.Failure == null, $"zero-energy detailed run failed: {tracer.Failure?.Message}");
            Harness.Assert(tracer.ThreadState() == System.Threading.ThreadState.Stopped,
                "zero-energy detailed run did not stop after cancellation");
            Harness.Assert(BankTotal(bank) == 0, "zero-energy detailed fixture recorded energy");
            return;
        }

        await tracer.Calculate_Conv(4, ConvSeed);

        Harness.Assert(feedback.Payloads.Count == cancelAt,
            $"zero-energy run did {feedback.Payloads.Count} rounds; expected exactly {cancelAt} (no auto-converge)");
        Harness.Assert(tracer.Failure == null, $"zero-energy run failed: {tracer.Failure?.Message}");
        Harness.Assert(tracer.ThreadState() == System.Threading.ThreadState.Stopped, "zero-energy run not stopped");
        Harness.Assert(BankTotal(bank) == 0, "zero-energy fixture recorded energy");
    }

    /// <summary>Real-energy run cancelled mid-flight: the current round still settles, success is not declared.</summary>
    private static async Task RunCancellationScenario()
    {
        (SplitRayTracer tracer, ConvFeedback feedback, _, Receiver_Bank bank) = ConvFixture(-1);
        int cancelAt = 3;
        feedback.OnRound = r => { if (r == cancelAt) feedback.Canceler.Cancel(); };

        await tracer.Calculate_Conv(4, ConvSeed);

        Harness.Assert(feedback.Payloads.Count == cancelAt,
            $"cancelled run did {feedback.Payloads.Count} rounds; expected exactly {cancelAt} (no early convergence)");
        Harness.Assert(tracer.Failure == null, $"cancelled run reported a failure: {tracer.Failure?.Message}");
        Harness.Assert(tracer.ThreadState() == System.Threading.ThreadState.Stopped, "cancelled run not stopped");
        Harness.Assert(tracer._currentRay[0] == cancelAt * RoundsPerCheck,
            $"cancelled run traced {tracer._currentRay[0]} rays; the in-flight round must settle all {cancelAt * RoundsPerCheck}");

        double settled = BankTotal(bank);
        await Task.Delay(150);
        Harness.Assert(BankTotal(bank) == settled, "receiver writes continued after the cancelled run stopped");
    }

    /// <summary>
    /// Controlled round-task failure: the injected shoot exception surfaces through
    /// SplitRayTracer.Failure, the simulation reaches Stopped instead of hanging, and no
    /// scene/ray activity continues after the run returns.
    /// </summary>
    private static async Task RunFailureScenario()
    {
        FlakyScene room = new FlakyScene();
        RandomSource src = Fixtures.Source();
        Receiver_Bank bank = Fixtures.ReceiverBank(src, room);
        ConvFeedback feedback = new ConvFeedback { Bank = bank };
        SplitRayTracer tracer = new SplitRayTracer(src, bank, room, Fixtures.CutOffTime_ms,
            new[] { 0, src.BandCount - 1 }, -1, -1, feedback);

        await tracer.Calculate_Conv(4, ConvSeed);

        Harness.Assert(tracer.Failure is InvalidOperationException && tracer.Failure.Message == FailureMessage,
            $"injected failure not surfaced; Failure = {tracer.Failure}");
        // A poller must not read a faulted run as a finished one: ThreadState() rethrows the
        // stored failure, the same contract ImageSourceData already uses.
        bool rethrown = false;
        try { tracer.ThreadState(); }
        catch (InvalidOperationException error) { rethrown = error.Message == FailureMessage; }
        Harness.Assert(rethrown, "failed run reported Stopped instead of rethrowing its failure");
        Harness.Assert(feedback.Payloads.Count == 0, "a failed round must not run convergence checks");

        double settled = BankTotal(bank);
        int shots = room.Shots;
        await Task.Delay(150);
        Harness.Assert(room.Shots == shots, "ray work continued after the failed run returned");
        Harness.Assert(BankTotal(bank) == settled, "receiver writes continued after the failed run returned");
    }

    /// <summary>Real-energy fixture: 5 x 4 x 3 m shoebox, reverberant enough to fill all convergence windows.</summary>
    private static (SplitRayTracer tracer, ConvFeedback feedback, RandomSource src, Receiver_Bank bank) ConvFixture(int rayCount)
    {
        Polygon_Scene room = BoxRoom();
        RandomSource src = new RandomSource(Fixtures.OctaveSpectrum(), new Point(0, 0, 0.8), 0, false);
        Receiver_Bank bank = new Receiver_Bank(new[] { new Point(1.2, 0.7, 0.3), new Point(-1.2, -0.7, -0.3) },
            src, room, Fixtures.SampleRate, Fixtures.CutOffTime_ms, Receiver_Bank.Type.Stationary, false);
        ConvFeedback feedback = new ConvFeedback { Bank = bank };
        SplitRayTracer tracer = new SplitRayTracer(src, bank, room, Fixtures.CutOffTime_ms,
            new[] { 0, src.BandCount - 1 }, -1, rayCount, feedback);
        return (tracer, feedback, src, bank);
    }

    /// <summary>Zero-energy fixture: silent source (power underflows to exactly 0) in an empty scene.</summary>
    private static (SplitRayTracer tracer, ConvFeedback feedback, RandomSource src, Receiver_Bank bank) ZeroEnergyFixture(int rayCount)
    {
        Scene room = Fixtures.Scene();
        RandomSource src = new RandomSource(Enumerable.Repeat(-4000.0, 8).ToArray(), new Point(0, 0, 2), 0, false);
        Receiver_Bank bank = Fixtures.ReceiverBank(src, room);
        ConvFeedback feedback = new ConvFeedback { Bank = bank };
        SplitRayTracer tracer = new SplitRayTracer(src, bank, room, Fixtures.CutOffTime_ms,
            new[] { 0, src.BandCount - 1 }, -1, rayCount, feedback);
        return (tracer, feedback, src, bank);
    }

    private static Polygon_Scene BoxRoom()
    {
        // One object, six inward-facing quads. Receivers (radius 1 m) and source clear of all walls.
        Point[][] faces =
        {
            new[] { new Point(-2.5, -2, -1.5), new Point(-2.5, 2, -1.5), new Point(2.5, 2, -1.5), new Point(2.5, -2, -1.5) },
            new[] { new Point(-2.5, -2, 1.5), new Point(2.5, -2, 1.5), new Point(2.5, 2, 1.5), new Point(-2.5, 2, 1.5) },
            new[] { new Point(-2.5, -2, -1.5), new Point(-2.5, -2, 1.5), new Point(-2.5, 2, 1.5), new Point(-2.5, 2, -1.5) },
            new[] { new Point(2.5, -2, -1.5), new Point(2.5, 2, -1.5), new Point(2.5, 2, 1.5), new Point(2.5, -2, 1.5) },
            new[] { new Point(-2.5, -2, -1.5), new Point(2.5, -2, -1.5), new Point(2.5, -2, 1.5), new Point(-2.5, -2, 1.5) },
            new[] { new Point(-2.5, 2, -1.5), new Point(-2.5, 2, 1.5), new Point(2.5, 2, 1.5), new Point(2.5, 2, -1.5) },
        };
        Polygon_Scene room = new Polygon_Scene(20, 50, 1013.25, 0, false, true);
        room.Construct(new Point[][][] { faces },
            new List<Pachyderm_Acoustic.Environment.Material> { new Basic_Material(Enumerable.Repeat(BoxAbsorption, 8).ToArray()) },
            new List<Scattering> { new Lambert_Scattering(Enumerable.Repeat(BoxScattering, 8).ToArray()) },
            new List<double[]> { new double[8] });
        room.partition();
        return room;
    }

    private static double BankTotal(Receiver_Bank bank)
    {
        double total = 0;
        foreach (Spherical_Receiver rec in bank.Rec_List)
            for (int oct = 0; oct < 8; oct++)
                total += rec.GetEnergyHistogram(oct).Sum();
        return total;
    }

    private static double MaxFinite(double[] values)
    {
        double max = double.MinValue;
        foreach (double v in values) if (!double.IsNaN(v) && !double.IsInfinity(v) && v > max) max = v;
        return max;
    }

    /// <summary>I_Conv_Progress stub: records each non-final round's Fill payload, the bank
    /// total at that moment, and the live heap right after the round's explicit GC.Collect.</summary>
    private sealed class ConvFeedback : I_Conv_Progress
    {
        public readonly CancellationTokenSource Canceler = new CancellationTokenSource();
        public readonly List<double[]> Payloads = new List<double[]>();
        public readonly List<double> Totals = new List<double>();
        public readonly List<long> Heap = new List<long>();
        public Receiver_Bank Bank;
        public Action<int> OnRound;

        public CancellationToken CancellationToken => Canceler.Token;
        CancellationTokenSource I_Conv_Progress.Canceler => Canceler;

        public void Update() { }
        public void Set_Targets(double[] X, double[] Y) { }
        public void Update_Targets() { }

        public void Fill(string title, double[] YArray, int ID, Queue<double[]> IR_Clear = null, int Ray_Req = 0, List<double[]> Add_Functions = null)
        {
            Payloads.Add(YArray == null ? Array.Empty<double>() : (double[])YArray.Clone());
            Heap.Add(GC.GetTotalMemory(true)); // production GC.Collect ran immediately before this check
            Totals.Add(BankTotal(Bank));
            OnRound?.Invoke(Payloads.Count);
        }
    }

    /// <summary>Empty scene whose shoot throws once the fault point is past the tracer
    /// constructor's validity probes, failing one round's ray task deterministically.</summary>
    private sealed class FlakyScene : Empty_Scene
    {
        private int shots;
        public int Shots => shots;

        public FlakyScene() : base(20, 50, 1013.25, 0, false, false, new Point(-10, -10, -10), new Point(10, 10, 10)) { }

        public override bool shoot(Ray R, out double u, out double v, out int Poly_ID, out List<Point> X_PT, out List<double> t, out List<int> code)
        {
            if (Interlocked.Increment(ref shots) > 50) throw new InvalidOperationException(FailureMessage);
            return base.shoot(R, out u, out v, out Poly_ID, out X_PT, out t, out code);
        }
    }
}
