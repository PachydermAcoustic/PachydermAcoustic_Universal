using Hare.Geometry;
using Pachyderm_Acoustic;
using Pachyderm_Acoustic.Audio;
using Pachyderm_Acoustic.Environment;

namespace RayTraceTests;

/// <summary>
/// Ticket 06 (spec T3): Directional_Histogram keeps exactly one Energy allocation (the base
/// class one). Both constructors × 8/24 bands must produce the base semantics — dimensions,
/// all-zero content, writable last sample — and the receiver file format must round-trip.
/// The allocation part proves by differential measurement (GC.GetAllocatedBytesForCurrentThread,
/// warm constructors, fixed input sizes) that the second Energy array is gone: with one Energy
/// set per receiver, Spherical_Receiver allocates 11 band·sample double arrays (1 Energy +
/// 1 Pressure + 6 direction + 3 RAD) and Map_Receiver 8 (no RAD); a surviving duplicate
/// Energy would read 12 / 9. No byte counts are hardcoded — only array counts per band·sample.
/// </summary>
internal static class HistogramEnergySuite
{
    public static void Register() => Harness.Suite("directional histogram single energy allocation (T3)", Run);

    private static Task Run()
    {
        // --- double-time constructor (Spherical_Receiver public surface), 8 and 24 bands ---
        CheckFreshHistogram(new Spherical_Receiver(Fixtures.SampleRate, Fixtures.CutOffTime_ms, 1.2, 8),
            8, Fixtures.SampleRate, "double-time ctor, 8 bands");
        CheckFreshHistogram(new Spherical_Receiver(Fixtures.SampleRate, Fixtures.CutOffTime_ms, 1.2, 24),
            24, Fixtures.SampleRate, "double-time ctor, 24 bands");
        // fractional cutoff truncates to whole samples: 333 Hz * 1.5 s = 499.5 -> 499 (exact FP: 1.5 and 333 are dyadic)
        CheckFreshHistogram(new Spherical_Receiver(333, 1500.0, 1.2, 8),
            8, 499, "double-time ctor, fractional cutoff");

        // bank construction path (source + scene) for both band resolutions
        var scene = Fixtures.Scene();
        var bank8 = Fixtures.ReceiverBank(Fixtures.Source(), scene);
        foreach (var r in bank8.Rec_List) CheckFreshHistogram(r, 8, Fixtures.SampleRate, "bank ctor, 8 bands");
        var bank24 = Fixtures.ReceiverBank(Fixtures.ThirdOctaveSource(), scene, thirdOctave: true);
        foreach (var r in bank24.Rec_List) CheckFreshHistogram(r, 24, Fixtures.SampleRate, "bank ctor, 24 bands");

        // --- int-sample constructor (Map_Receiver public surface), 8 and 24 bands ---
        var origin = new Point(1, 2, 3);
        var srcAt = new Point(0, 0, 2);
        CheckFreshHistogram(new PachMapReceiver.Map_Receiver(origin, srcAt, 1.5, 1.2, 0, Fixtures.SampleRate, 777, true, false),
            8, 777, "int-sample ctor, 8 bands");
        CheckFreshHistogram(new PachMapReceiver.Map_Receiver(origin, srcAt, 1.5, 1.2, 0, Fixtures.SampleRate, 777, true, true),
            24, 777, "int-sample ctor, 24 bands");

        CheckReadWriteRoundTrip(bank8);
        CheckReadWriteRoundTrip(bank24);

        AssertSingleEnergyAllocation();
        return Task.CompletedTask;
    }

    /// <summary>Base-semantics structure of a freshly constructed directional histogram.</summary>
    private static void CheckFreshHistogram(Spherical_Receiver rec, int bands, int expectedSamples, string label)
    {
        for (int oct = 0; oct < bands; oct++)
        {
            double[] e = rec.GetEnergyHistogram(oct);
            double[] p = rec.GetRMSPressureHistogram(oct);
            Harness.Assert(ReferenceEquals(e, rec.GetEnergyHistogram(oct)) && ReferenceEquals(p, rec.GetRMSPressureHistogram(oct)),
                $"{label}: band {oct} histogram accessors must return the live storage, not a copy");
            Harness.Assert(e.Length == expectedSamples && p.Length == expectedSamples,
                $"{label}: band {oct} lengths {e.Length}/{p.Length}, expected {expectedSamples}");
            Harness.Assert(e.All(x => x == 0) && p.All(x => x == 0),
                $"{label}: band {oct} energy/pressure not all zero after construction");
            for (int c = 0; c < 3; c++)
            {
                Harness.Assert(rec.Directions_Pos(oct, 0, c) == 0 && rec.Directions_Pos(oct, expectedSamples - 1, c) == 0 &&
                               rec.Directions_Neg(oct, 0, c) == 0 && rec.Directions_Neg(oct, expectedSamples - 1, c) == 0,
                    $"{label}: band {oct} axis {c} directions not zero after construction");
            }
        }

        if (rec.Reflection_Analysis_Data != null)
        {
            for (int type = 0; type < 3; type++)
                for (int oct = 0; oct < bands; oct++)
                {
                    double[] rad = rec.Reflection_Analysis_Data[type][oct];
                    Harness.Assert(rad.Length == expectedSamples && rad.All(x => x == 0),
                        $"{label}: RAD {type} band {oct} length {rad.Length}, expected {expectedSamples}, or not all zero");
                }
        }

        // last sample is writable through the public surface and shows up in the returned arrays
        int last = expectedSamples - 1;
        Vector pos = new Vector(0.25, 0.5, 1.0);
        Vector neg = new Vector(-0.125, -0.75, -2.0);
        for (int oct = 0; oct < bands; oct++)
        {
            rec.Combine_Sample(last, (oct + 1) * 0.5, (oct + 1) * 0.25, pos, neg, oct);
            Harness.Assert(rec.GetEnergyHistogram(oct)[last] == (oct + 1) * 0.5, $"{label}: band {oct} last energy sample not writable");
            Harness.Assert(rec.GetRMSPressureHistogram(oct)[last] == (oct + 1) * 0.25, $"{label}: band {oct} last pressure sample not writable");
            Harness.Assert(rec.Directions_Pos(oct, last, 0) == 0.25 && rec.Directions_Pos(oct, last, 1) == 0.5 && rec.Directions_Pos(oct, last, 2) == 1.0 &&
                           rec.Directions_Neg(oct, last, 0) == -0.125 && rec.Directions_Neg(oct, last, 1) == -0.75 && rec.Directions_Neg(oct, last, 2) == -2.0,
                $"{label}: band {oct} last direction samples not writable");
        }
    }

    /// <summary>Result read/write format: Write_Data then Read_Data reproduces every energy, pressure and direction sample.</summary>
    private static void CheckReadWriteRoundTrip(Receiver_Bank bank)
    {
        Pach_SP.Initialize_FFTW(); // Read_Data finishes with filter creation, which needs the FFT cache
        int bands = bank.SWL.Length;
        int last = bank.SampleCT - 1;
        Vector pos = new Vector(0.25, 0.5, 1.0);
        Vector neg = new Vector(-0.125, -0.75, -2.0);
        foreach (var r in bank.Rec_List)
            for (int oct = 0; oct < bands; oct++)
                r.Combine_Sample(last, (oct + 1) * 0.5, (oct + 1) * 0.25, pos, neg, oct);

        using var stream = new MemoryStream();
        var writer = new BinaryWriter(stream);
        bank.Write_Data(ref writer);
        writer.Flush();
        stream.Position = 0;
        var reader = new BinaryReader(stream);
        Harness.Assert(reader.ReadString() == "Ray-Traced_Data", "round-trip: stream does not start with the ray-traced data indicator");

        int sampleRate = 0;
        var read = Receiver_Bank.Read_Data(ref reader, bank.SWL.ToArray(), bank.Rec_List.Length,
            Fixtures.ReceiverPoints, Enumerable.Repeat(1.2, bank.Rec_List.Length).ToArray(), 0.0,
            ref sampleRate, "2.6.0.0", null);
        Harness.Assert(sampleRate == bank.SampleRate && read.SampleCT == bank.SampleCT,
            $"round-trip: sample rate/count changed ({sampleRate}/{read.SampleCT} vs {bank.SampleRate})");

        for (int r = 0; r < bank.Rec_List.Length; r++)
            for (int oct = 0; oct < bands; oct++)
            {
                double[] e0 = bank.Rec_List[r].GetEnergyHistogram(oct);
                double[] e1 = read.Rec_List[r].GetEnergyHistogram(oct);
                double[] p0 = bank.Rec_List[r].GetRMSPressureHistogram(oct);
                double[] p1 = read.Rec_List[r].GetRMSPressureHistogram(oct);
                Harness.Assert(e1.Length == e0.Length && p1.Length == p0.Length,
                    $"round-trip: receiver {r} band {oct} lengths changed ({e1.Length}/{p1.Length} vs {e0.Length}/{p0.Length})");
                for (int s = 0; s < e0.Length; s++)
                {
                    string where = $"round-trip: receiver {r} band {oct} sample {s}";
                    Harness.Assert(e1[s] == e0[s], $"{where}: energy {e1[s]} vs {e0[s]}");
                    Harness.Assert(p1[s] == p0[s], $"{where}: pressure {p1[s]} vs {p0[s]}");
                    for (int c = 0; c < 3; c++)
                    {
                        Harness.Assert(read.Rec_List[r].Directions_Pos(oct, s, c) == bank.Rec_List[r].Directions_Pos(oct, s, c),
                            $"{where}: pos axis {c} changed");
                        Harness.Assert(read.Rec_List[r].Directions_Neg(oct, s, c) == bank.Rec_List[r].Directions_Neg(oct, s, c),
                            $"{where}: neg axis {c} changed");
                    }
                }
            }
    }

    /// <summary>
    /// Differential allocation check: doubling the histogram length must add exactly the
    /// per-band·sample arrays the (single-Energy) structure accounts for — 11 for
    /// Spherical_Receiver (E + P + 6 dir + 3 RAD), 8 for Map_Receiver (no RAD). A surviving
    /// duplicate Energy set would measure 12 / 9. Fixed-size overhead (headers, receiver
    /// shells, points) cancels in the difference; constructors are warmed first.
    /// </summary>
    private static void AssertSingleEnergyAllocation()
    {
        const int n = 64, bands = 8, rate = 1000, s1 = 2000, s2 = 4000;
        var sink = new Spherical_Receiver[n];
        for (int i = 0; i < n; i++) sink[i] = new Spherical_Receiver(rate, 4000.0, 1.2, bands);

        long b0 = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < n; i++) sink[i] = new Spherical_Receiver(rate, 2000.0, 1.2, bands);
        long b1 = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < n; i++) sink[i] = new Spherical_Receiver(rate, 4000.0, 1.2, bands);
        long b2 = GC.GetAllocatedBytesForCurrentThread();
        double coef = ((b2 - b1) - (b1 - b0)) / ((double)n * bands * (s2 - s1) * sizeof(double));
        Harness.Assert(Math.Abs(coef - 11) < 0.1,
            $"double-time ctor: expected 11 band·sample arrays (1 energy + 1 pressure + 6 direction + 3 RAD), measured {coef:0.###} — duplicate Energy allocation still present?");

        var origin = new Point(1, 2, 3);
        var srcAt = new Point(0, 0, 2);
        var mapSink = new PachMapReceiver.Map_Receiver[n];
        for (int i = 0; i < n; i++) mapSink[i] = new PachMapReceiver.Map_Receiver(origin, srcAt, 1.5, 1.2, 0, rate, s2, true, false);

        b0 = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < n; i++) mapSink[i] = new PachMapReceiver.Map_Receiver(origin, srcAt, 1.5, 1.2, 0, rate, s1, true, false);
        b1 = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < n; i++) mapSink[i] = new PachMapReceiver.Map_Receiver(origin, srcAt, 1.5, 1.2, 0, rate, s2, true, false);
        b2 = GC.GetAllocatedBytesForCurrentThread();
        coef = ((b2 - b1) - (b1 - b0)) / ((double)n * bands * (s2 - s1) * sizeof(double));
        Harness.Assert(Math.Abs(coef - 8) < 0.1,
            $"int-sample ctor: expected 8 band·sample arrays (1 energy + 1 pressure + 6 direction), measured {coef:0.###} — duplicate Energy allocation still present?");
    }
}
