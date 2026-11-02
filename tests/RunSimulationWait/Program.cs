using System;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Utilities;

// Ticket 08 (spec D7) headless proof for RCPachTools.RunSimulation, linked from the Rhino clone.
// Scenarios: completed work pays no fixed delay and combines exactly once; still-running work
// (first phase, and follow-up work started by the combine) waits on the 3,000 ms poll; failure
// and cancellation propagate before any combine; autosave is restored on every exit path.
// Exit code 0 = all scenarios passed.

int failed = 0;
void Check(bool ok, string name)
{
    Console.WriteLine($"{(ok ? "PASS" : "FAIL")}  {name}");
    if (!ok) failed++;
}

// 1/2: already-completed simulation => no first or post-combine fixed delay, combine exactly once.
foreach (bool feedback in new[] { true, false })
{
    Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled = true;
    Rhino.RhinoApp.PromptCalls = 0;
    FakeSim sim = new FakeSim();
    Stopwatch sw = Stopwatch.StartNew();
    Simulation_Type result = await RCPachTools.RunSimulation(sim, feedback);
    sw.Stop();
    Check(sw.ElapsedMilliseconds < 1000, $"completed sim pays no fixed delay, either phase (needfeedback={feedback}) [{sw.ElapsedMilliseconds} ms]");
    Check(ReferenceEquals(result, sim) && sim.CombineCalls == 1, $"completed sim returns after exactly one combine (needfeedback={feedback})");
    // One prompt announces the start; a stopped simulation must add no polling prompts.
    Check(Rhino.RhinoApp.PromptCalls == 1, $"no polling prompt churn when already stopped (needfeedback={feedback})");
    Check(Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled, $"autosave restored after success (needfeedback={feedback})");
}

// 3: still-running first phase waits on the 3,000 ms poll (a 500 ms job must cost one full poll).
{
    Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled = true;
    Rhino.RhinoApp.PromptCalls = 0;
    FakeSim sim = new FakeSim { RunPhase1 = TimeSpan.FromMilliseconds(500) };
    Stopwatch sw = Stopwatch.StartNew();
    await RCPachTools.RunSimulation(sim, true);
    sw.Stop();
    Check(sw.ElapsedMilliseconds >= 2900 && sw.ElapsedMilliseconds <= 9000, $"running first-phase work polls at 3,000 ms [{sw.ElapsedMilliseconds} ms]");
    Check(sim.CombineCalls == 1, "running first-phase work still combines exactly once");
    Check(Rhino.RhinoApp.PromptCalls == 2, "progress prompt still updates once per poll while running");
    Check(Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled, "autosave restored after polling wait");
}

// 4: the combine starts follow-up work (500 ms); the second phase waits for it, still one combine.
{
    Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled = true;
    FakeSim sim = new FakeSim { RunPhase2 = TimeSpan.FromMilliseconds(500) };
    Stopwatch sw = Stopwatch.StartNew();
    await RCPachTools.RunSimulation(sim, true);
    sw.Stop();
    Check(sw.ElapsedMilliseconds >= 2900 && sw.ElapsedMilliseconds <= 9000, $"combine-started work is waited for on the 3,000 ms poll [{sw.ElapsedMilliseconds} ms]");
    Check(sim.CombineCalls == 1, "combine-started work path combines exactly once");
    Check(Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled, "autosave restored after combine-started work");
}

// 5: a failure reported through the state protocol propagates; no combine while a writer failed.
{
    Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled = true;
    FakeSim sim = new FakeSim { FailWith = new InvalidOperationException("simulated filter failure") };
    Exception caught = null;
    try { await RCPachTools.RunSimulation(sim, true); }
    catch (Exception e) { caught = e; }
    Check(ReferenceEquals(caught, sim.FailWith), "failure propagates out of RunSimulation");
    Check(sim.CombineCalls == 0, "no combine forced after failure");
    Check(Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled, "autosave restored after failure");
}

// 6: workers observe cancellation and stop; the combine is skipped (no combine with writers alive).
{
    Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled = true;
    CancellationTokenSource cts = new CancellationTokenSource();
    FakeSim sim = new FakeSim { OnBegin = () => cts.Cancel() };
    OperationCanceledException cancelled = null;
    try { await RCPachTools.RunSimulation(sim, true, cts.Token); }
    catch (OperationCanceledException e) { cancelled = e; }
    Check(cancelled != null && cancelled.CancellationToken == cts.Token, "cancellation surfaces before the combine");
    Check(sim.CombineCalls == 0, "no combine forced after cancellation");
    Check(Rhino.ApplicationSettings.FileSettings.AutoSaveEnabled, "autosave restored after cancellation");
}

Console.WriteLine(failed == 0 ? "RunSimulationWait: all scenarios passed" : $"RunSimulationWait: {failed} check(s) failed");
return failed;

/// <summary> Simulation_Type double whose ThreadState follows the same protocol as the real sims:
/// Running while the phase window is open, Stopped after, and exceptions surfaced from the state
/// check (the ticket-03 failure contract). </summary>
internal sealed class FakeSim : Simulation_Type
{
    private readonly Stopwatch clock = new Stopwatch();
    private TimeSpan combineAt;

    internal TimeSpan RunPhase1;  // Running while clock < this after Begin
    internal TimeSpan RunPhase2;  // after Combine: Running while clock < combineAt + this
    internal Exception FailWith;  // if set, ThreadState() throws (failure contract)
    internal Action OnBegin;

    public int CombineCalls;

    public override void Begin()
    {
        clock.Restart();
        OnBegin?.Invoke();
    }

    public override string Sim_Type() => "Fake";
    public override string ProgressMsg() => "fake progress";

    public override System.Threading.ThreadState ThreadState()
    {
        if (FailWith != null) throw FailWith;
        if (clock.Elapsed < RunPhase1) return System.Threading.ThreadState.Running;
        if (CombineCalls > 0 && clock.Elapsed - combineAt < RunPhase2) return System.Threading.ThreadState.Running;
        return System.Threading.ThreadState.Stopped;
    }

    public override void Combine_ThreadLocal_Results()
    {
        CombineCalls++;
        combineAt = clock.Elapsed;
    }
}
