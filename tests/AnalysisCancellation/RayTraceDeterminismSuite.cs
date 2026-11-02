using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Ticket 01 smoke assertion: the same explicit seed and fixed ray count produce a
/// repeatable ray-tracing result (receiver energy histograms and first RAD category
/// identical, element for element), with progress exactly at the fixed ray count.
/// </summary>
internal static class RayTraceDeterminismSuite
{
    public static void Register() => Harness.Suite("ray trace determinism (explicit seed, fixed count)", Run);

    private static async Task Run()
    {
        SplitRayTracer first = await Fixtures.RunSeeded();
        SplitRayTracer second = await Fixtures.RunSeeded();

        Harness.Assert(first._currentRay[0] == Fixtures.TotalRays,
            $"first run traced {first._currentRay[0]} rays, expected {Fixtures.TotalRays}");
        Harness.Assert(second._currentRay[0] == Fixtures.TotalRays,
            $"second run traced {second._currentRay[0]} rays, expected {Fixtures.TotalRays}");

        var a = first.GetReceiver.Rec_List;
        var b = second.GetReceiver.Rec_List;
        Harness.Assert(a.Length == b.Length, "receiver counts differ between runs");
        double total = 0;
        for (int r = 0; r < a.Length; r++)
        {
            for (int oct = 0; oct < 8; oct++)
            {
                string where = $"receiver {r} octave {oct}";
                double[] ha = a[r].GetEnergyHistogram(oct);
                double[] hb = b[r].GetEnergyHistogram(oct);
                Harness.Assert(ha.Length == hb.Length, $"{where}: histogram lengths {ha.Length} vs {hb.Length}");
                for (int s = 0; s < ha.Length; s++)
                    Harness.Assert(ha[s] == hb[s], $"{where} sample {s}: {ha[s]} vs {hb[s]}");
                total += ha.Sum();
            }

            double[][] radA = a[r].Reflection_Analysis_Data[0];
            double[][] radB = b[r].Reflection_Analysis_Data[0];
            for (int oct = 0; oct < radA.Length; oct++)
            {
                for (int s = 0; s < radA[oct].Length; s++)
                    Harness.Assert(radA[oct][s] == radB[oct][s],
                        $"receiver {r} RAD octave {oct} sample {s}: {radA[oct][s]} vs {radB[oct][s]}");
            }
        }
        Harness.Assert(total > 0, "no receiver energy recorded; determinism comparison would be vacuous");
    }
}
