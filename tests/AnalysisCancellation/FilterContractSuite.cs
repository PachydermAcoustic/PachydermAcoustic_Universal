using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Ticket 03 (spec D4, matrix T4): the three <see cref="Receiver_Bank.Create_Filter"/>
/// entries complete synchronously — normal return means every requested filter is
/// available, controllable failures surface to the caller (a background batch call's
/// failure is observable via await), and <c>CreatingFilter</c> is reset after both
/// normal and failing exits, so no half-built filter state is published as complete.
///
/// Failure seam: the <see cref="Spherical_Receiver"/> filter trio is public virtual, so a
/// subclass whose Create_Filter throws, installed through the public Rec_List field of a
/// normally constructed bank, makes filter creation fail controllably.
/// </summary>
internal static class FilterContractSuite
{
    private const string FailureMessage = "T4 injected receiver filter failure";

    private sealed class ThrowingReceiver : Spherical_Receiver
    {
        public ThrowingReceiver(double rhoC)
            : base(Fixtures.SampleRate, Fixtures.CutOffTime_ms, rhoC, 8) { }

        public override void Create_Filter() => throw new InvalidOperationException(FailureMessage);
        public override double[] Create_Filter(double[] SWL, int dim) => throw new InvalidOperationException(FailureMessage);
        public override double[][] Create_Filter(double[] SWL) => throw new InvalidOperationException(FailureMessage);
    }

    public static void Register() => Harness.Suite("filter completion and failure contract (three entries)", Run);

    private static async Task Run()
    {
        Pachyderm_Acoustic.Audio.Pach_SP.Initialize_FFTW();
        double[] swl = Fixtures.OctaveSpectrum();

        // Batch entry, normal return: a seeded run's bank must have every filter installed.
        SplitRayTracer tracer = await Fixtures.RunSeeded();
        Receiver_Bank traced = tracer.GetReceiver;
        traced.Create_Filter();
        for (int r = 0; r < traced.Rec_List.Length; r++)
            Harness.Assert(ReturnsPromptly(() => traced.HasFilter(r), $"HasFilter({r}) after normal batch Create_Filter"),
                $"receiver {r} has no filter after a normally returned batch Create_Filter");

        // Array-returning overloads, normal return: keep their result semantics and reset the flag.
        double[] oneAxis = traced.Create_Filter(swl, 0, 0);
        Harness.Assert(oneAxis != null && oneAxis.Length > 0, "double[] Create_Filter returned no filter array");
        double[][] threeAxis = traced.Create_Filter(swl, 0);
        Harness.Assert(threeAxis != null && threeAxis.Length > 0 && threeAxis[0].Length > 0,
            "double[][] Create_Filter returned no filter arrays");
        ReturnsPromptly(() => traced.HasFilter(0), "HasFilter(0) after overload Create_Filter");

        // Batch entry, synchronous failure observed directly by the caller.
        Receiver_Bank failing = BankWithThrowerAt(0);
        ExpectInjectedFailure(() => failing.Create_Filter(), "batch Create_Filter");
        Harness.Assert(!ReturnsPromptly(() => failing.HasFilter(0), "HasFilter(0) after failing batch Create_Filter"),
            "failing receiver reports a filter after a failed batch Create_Filter");

        // Batch entry executed as background work: the failure is observable through await.
        Receiver_Bank background = BankWithThrowerAt(0);
        Exception observed = null;
        try { await Task.Run(() => background.Create_Filter()); }
        catch (Exception error) { observed = error; }
        Harness.Assert(observed is InvalidOperationException && observed.Message == FailureMessage,
            $"awaited background batch Create_Filter did not surface the injected failure (got: {observed?.GetType().Name ?? "no exception"})");
        Harness.Assert(!ReturnsPromptly(() => background.HasFilter(0), "HasFilter(0) after awaited failing batch Create_Filter"),
            "failing receiver reports a filter after an awaited failed batch Create_Filter");

        // double[] overload, failure: propagates synchronously and resets the flag.
        Receiver_Bank oneD = BankWithThrowerAt(0);
        ExpectInjectedFailure(() => oneD.Create_Filter(swl, 0, 0), "double[] Create_Filter");
        ReturnsPromptly(() => oneD.HasFilter(0), "HasFilter(0) after failing double[] Create_Filter");

        // double[][] overload, failure: propagates synchronously and resets the flag.
        Receiver_Bank twoD = BankWithThrowerAt(1);
        ExpectInjectedFailure(() => twoD.Create_Filter(swl, 1), "double[][] Create_Filter");
        ReturnsPromptly(() => twoD.HasFilter(1), "HasFilter(1) after failing double[][] Create_Filter");
    }

    private static Receiver_Bank BankWithThrowerAt(int index)
    {
        Scene room = Fixtures.Scene();
        Receiver_Bank bank = Fixtures.ReceiverBank(Fixtures.Source(), room);
        bank.Rec_List[index] = new ThrowingReceiver(room.Rho(Fixtures.ReceiverPoints[index]));
        return bank;
    }

    private static void ExpectInjectedFailure(Action call, string entry)
    {
        Exception seen = null;
        try { call(); }
        catch (Exception error) { seen = error; }
        Harness.Assert(seen is InvalidOperationException && seen.Message == FailureMessage,
            $"{entry} did not surface the injected failure to its caller (got: {seen?.GetType().Name ?? "no exception"})");
    }

    /// <summary>
    /// CreatingFilter is private; its reset is observable because Receiver_Bank.HasFilter
    /// spin-waits on it. A stuck-true flag would hang the call, hence the timeout probe.
    /// </summary>
    private static T ReturnsPromptly<T>(Func<T> operation, string what)
    {
        Task<T> probe = Task.Run(operation);
        Harness.Assert(probe.Wait(TimeSpan.FromSeconds(10)), $"{what} did not return within 10 s (CreatingFilter not reset?)");
        return probe.Result;
    }
}
