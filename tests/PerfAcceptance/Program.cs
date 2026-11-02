// Ticket 09 (Deliverable B): performance acceptance measurement — AFTER side.
// Fixed scenario identical to tests/PhaseTimingBaseline (3 sources × 4 receivers,
// 44100 Hz, 500 ms, octave, 2000 rays fixed count, fixed worker seeds), plus a
// 1-source variant, run cold (fresh process per run) and warm (in-process repeats).
using RayTraceTests;
// Reports per-source stage wall clock, GC/alloc deltas, and process peak memory.
//
// Usage:
//   dotnet run -- cold            # one cold 3-source run (pool built lazily inside)
//   dotnet run -- cold1           # one cold 1-source run
//   dotnet run -- warm [n]        # n (default 4) warm in-process 3-source runs
//   dotnet run -- warm1 [n]       # n (default 4) warm in-process 1-source runs
// Exit code 0 = pass (scenario completed; numbers are the output).

using System.Diagnostics;

string mode = args.Length > 0 ? args[0] : "warm";
int repeats = args.Length > 1 && int.TryParse(args[1], out int n) ? n : 4;
int sources = mode.EndsWith('1') ? 1 : 3;

Pachyderm_Acoustic.Audio.Pach_SP.Initialize_FFTW();
Pachyderm_Acoustic.Audio.Pach_SP.Filter = new Pachyderm_Acoustic.Audio.Pach_SP.Minimum_Phase_System();

Console.WriteLine($"perf-after: mode={mode} sources={sources} workers={Environment.ProcessorCount} rays/source={2000} runtime={Environment.Version}");

for (int run = 1; run <= (mode.StartsWith("cold") ? 1 : repeats); run++)
{
    GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();
    long alloc0 = GC.GetTotalAllocatedBytes(precise: true);
    int gc0 = GC.CollectionCount(0) + GC.CollectionCount(1) + GC.CollectionCount(2);
    var swTotal = Stopwatch.StartNew();

    var summary = await Perf.Scenario.Run(sources);

    swTotal.Stop();
    long alloc1 = GC.GetTotalAllocatedBytes(precise: true);
    int gc1 = GC.CollectionCount(0) + GC.CollectionCount(1) + GC.CollectionCount(2);
    Console.WriteLine($"run {run}: total {swTotal.Elapsed.TotalSeconds:0.###} s | alloc {(alloc1 - alloc0) / 1048576.0:0.0} MB | GC {gc1 - gc0} | peak WS {Process.GetCurrentProcess().PeakWorkingSet64 / 1048576.0:0.0} MB");
    foreach (var line in summary) Console.WriteLine(line);
}
return 0;
