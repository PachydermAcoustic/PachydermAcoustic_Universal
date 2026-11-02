using System.Collections.Concurrent;
using Hare.Geometry;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Ticket 04 (spec D5, acceptance T5): single-band OctaveRay hits must accumulate all
/// three Reflection_Analysis_Data categories under the same per-band locking as the
/// multi-band path, and out-of-window samples must be neither written nor thrown on.
/// Controlled inputs make every per-cell contribution an exactly representable dyadic
/// double (1, 0.5 or 0.25), so multi-worker totals are compared to the single-worker
/// reference and to the closed-form total with exact equality — a lost read-modify-write
/// changes the integer-dyadic total; no float tolerance is used.
///
/// Fixture design: receiver of radius 1 at the origin, ray starting 2 m away aimed dead
/// center (chord half-length tsphere is exactly 1), zero attenuation spectrum
/// (air term Math.Pow(10, ±0) is exactly 1), so the recorded energy is
/// fl(Intensity * SizeMod) = 1 for Intensity = π; third-octave fractions
/// {0.5, 0.25, 0.25} then contribute exactly 0.5 / 0.25 / 0.25 per band.
/// </summary>
internal static class RadAccumulationSuite
{
    private const int Workers = 8;
    private const int HitsPerWorker = 5000;               // per-cell add count 40,000 (exact in double)
    private const int HitsPerCell = Workers * HitsPerWorker;
    private const int SampleRate = 1000;
    private const int CutOffMs = 1000;                    // RAD window = 1000 samples (0..999)
    private const int ThirdOctaveStart = 21;              // bands 21..23 of a 24-band receiver

    private static readonly Point Wall = new(5, 0, 0);    // ray end, far beyond the receiver sphere

    public static void Register() =>
        Harness.Suite("RAD concurrent accumulation (T5: single-band + third-octave, window edges)", Run);

    private static async Task Run()
    {
        // (name, speed of sound, ray t_sum, expected sample index; -1 means outside the window)
        (string name, double cSound, double tSum, int sample)[] scenarios =
        {
            ("sample 0",       1e6, 0.0,      0),
            ("mid sample",     1e3, 0.0,      2),
            ("last sample",    1e6, 0.9995, 999),
            ("past window",    1e6, 5.0,     -1),
            ("before window",  1e6, -0.5,    -1),
        };

        foreach (var s in scenarios)
        {
            SingleBandScenario(s.name, s.cSound, s.tSum, s.sample);
            ThirdOctaveScenario(s.name, s.cSound, s.tSum, s.sample);
        }
        await Task.CompletedTask;
    }

    /// <summary>BandCount == 1 rays: all eight octaves × all three RAD categories.</summary>
    private static void SingleBandScenario(string name, double cSound, double tSum, int sample)
    {
        Spherical_Receiver reference = Receiver(8, cSound);
        Spherical_Receiver shared = Receiver(8, cSound);
        OctaveRay[] rays = new OctaveRay[8 * 3];
        for (int oct = 0; oct < 8; oct++)
            for (int kind = 0; kind < 3; kind++)
                rays[3 * oct + kind] = Ray(oct, (short)kind, tSum);

        RunScenario($"single-band, {name}", rays, reference, shared, sample,
            (kind, band) => HitsPerCell);                 // every octave: N adds of exactly 1.0
    }

    /// <summary>BandCount == 3 ray on a 24-band receiver: bands 21..23 × all three categories.</summary>
    private static void ThirdOctaveScenario(string name, double cSound, double tSum, int sample)
    {
        double[] fractions = { 0.5, 0.25, 0.25 };
        Spherical_Receiver reference = Receiver(24, cSound);
        Spherical_Receiver shared = Receiver(24, cSound);
        OctaveRay[] rays = new OctaveRay[3];
        for (int kind = 0; kind < 3; kind++)
        {
            OctaveRay r = Ray(ThirdOctaveStart, (short)kind, tSum);
            r.BandStart = ThirdOctaveStart;
            r.BandCount = 3;
            r.Fractions = fractions;
            rays[kind] = r;
        }

        RunScenario($"third-octave, {name}", rays, reference, shared, sample,
            (kind, band) => band >= ThirdOctaveStart ? HitsPerCell * fractions[band - ThirdOctaveStart] : 0);
    }

    private static void RunScenario(string where, OctaveRay[] rays, Spherical_Receiver reference,
        Spherical_Receiver shared, int sample, Func<int, int, double> expectedCell)
    {
        // Single-worker reference: each touched cell receives HitsPerCell identical adds.
        for (int h = 0; h < HitsPerCell; h++)
            foreach (OctaveRay r in rays) reference.CheckRay(r, Wall);

        HammerConcurrently(shared, rays, where);

        AssertRadMatchesReference(shared, reference, where);
        if (sample < 0)
        {
            AssertAllZero(reference.Reflection_Analysis_Data, where);
            return;
        }
        for (int kind = 0; kind < 3; kind++)
            for (int band = 0; band < reference.Reflection_Analysis_Data[kind].Length; band++)
                AssertOnlyCell(reference.Reflection_Analysis_Data, kind, band, sample, expectedCell(kind, band), where);
    }

    /// <summary>All workers hammer the shared receiver at once; any worker exception fails the scenario.</summary>
    private static void HammerConcurrently(Spherical_Receiver rec, OctaveRay[] rays, string where)
    {
        var errors = new ConcurrentQueue<Exception>();
        var start = new Barrier(Workers);
        var threads = new Thread[Workers];
        for (int w = 0; w < Workers; w++)
        {
            threads[w] = new Thread(() =>
            {
                start.SignalAndWait();
                try
                {
                    for (int h = 0; h < HitsPerWorker; h++)
                        foreach (OctaveRay r in rays) rec.CheckRay(r, Wall);
                }
                catch (Exception e) { errors.Enqueue(e); }
            });
            threads[w].Start();
        }
        foreach (Thread t in threads) t.Join();
        Harness.Assert(errors.IsEmpty, $"{where}: worker exception(s): {string.Join("; ", errors)}");
    }

    private static void AssertRadMatchesReference(Spherical_Receiver shared, Spherical_Receiver reference, string where)
    {
        double[][][] a = shared.Reflection_Analysis_Data, b = reference.Reflection_Analysis_Data;
        for (int kind = 0; kind < 3; kind++)
            for (int band = 0; band < a[kind].Length; band++)
                for (int s = 0; s < a[kind][band].Length; s++)
                    Harness.Assert(a[kind][band][s] == b[kind][band][s],
                        $"{where}: RAD[{kind}][{band}][{s}] multi-worker {a[kind][band][s]} vs single-worker {b[kind][band][s]}");
    }

    /// <summary>The expected sample holds exactly the closed-form dyadic total; every other sample stays zero.</summary>
    private static void AssertOnlyCell(double[][][] rad, int kind, int band, int sample, double expected, string where)
    {
        for (int s = 0; s < rad[kind][band].Length; s++)
        {
            double want = s == sample ? expected : 0;
            Harness.Assert(rad[kind][band][s] == want,
                $"{where}: RAD[{kind}][{band}][{s}] is {rad[kind][band][s]}, expected exactly {want}");
        }
    }

    private static void AssertAllZero(double[][][] rad, string where)
    {
        for (int kind = 0; kind < 3; kind++)
            AssertOnlyCell(rad, kind, 0, -1, 0, where);   // sample -1 never matches: asserts every sample is 0
    }

    private static Spherical_Receiver Receiver(int bands, double cSound) =>
        new(new Point(0, 0, 0), new Point(0, 0, 0), new double[bands], cSound, 1.2, SampleRate, CutOffMs);

    private static OctaveRay Ray(int octave, short mask, double tSum) =>
        new(-2, 0, 0, 1, 0, 0, 0, 0, Math.PI, tSum, octave, 0, 0) { Scatter_Mask = mask };
}
