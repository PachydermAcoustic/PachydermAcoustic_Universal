// Ticket 09 (Deliverable A): T6 fixed-count BEFORE/AFTER regression runner.
// Standalone entry: dumps the fixed-count fixture to a file (record mode) or compares
// against an existing fixture (compare mode). The scenario and comparison code live in
// T6Scenario.cs / T6Shared.cs, which are shared verbatim with tests/AnalysisCancellation
// (T6RegressionSuite) and with the baseline worktree copy that produces the BEFORE side.
//
// Usage: dotnet run [-- <fixture-file> [strict]]   (default t6-fixture.txt; first run
// records, later runs compare; "strict" fails on RAD deltas too — use it for same-build
// comparisons. Cross-build BEFORE compares default to reporting RAD deltas because the
// baseline RAD accumulation is racy, spec D5).
// Exit code 0 = pass.

using RayTraceTests;

string fixturePath = args.Length > 0 ? args[0] : "t6-fixture.txt";
bool strict = args.Contains("strict") || fixturePath == "t6-fixture.txt";
Console.WriteLine($"T6 runner: workers={Environment.ProcessorCount} rays={T6Shared.RAY_COUNT} sources={T6Scenario.SOURCES} receivers={T6Scenario.RECEIVERS} bands={T6Scenario.BANDS} seeds={T6Shared.BASE_SEED:x}+p runtime={Environment.Version} os={Environment.OSVersion}");
string current = await T6Scenario.Run();

if (!File.Exists(fixturePath))
{
    File.WriteAllText(fixturePath, current);
    Console.WriteLine($"fixture recorded to {Path.GetFullPath(fixturePath)} (first run).");
}
else
{
    string fresh = Path.Combine(Path.GetTempPath(), $"t6-run-{Guid.NewGuid():N}.txt");
    File.WriteAllText(fresh, current);
    T6Compare.Files(fixturePath, fresh, $"compare vs {Path.GetFileName(fixturePath)}", radReportOnly: !strict);
}
Console.WriteLine("PASS: T6 fixture matches (or was recorded).");
return 0;
