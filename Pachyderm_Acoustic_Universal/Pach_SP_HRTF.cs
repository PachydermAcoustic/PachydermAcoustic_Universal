using System;
using System.IO;
using System.Collections.Generic;
using System.Numerics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Buffers;
using System.Collections.Concurrent;
using Pachyderm_Acoustic.Utilities;
using static Pachyderm_Acoustic.Audio.SystemResponseCompensation;
using Hare.Geometry;
using Vector = Hare.Geometry.Vector;
using NAudio.Wave;
using NAudio.Wave.SampleProviders;
using Eto.Forms;

namespace Pachyderm_Acoustic
{
    namespace Audio
    {
        public partial class Pach_SP_HRTF
        {
            public static double ComputeTotalEnergy(double[][][] HRIRs)
            {
                double totalEnergy = 0.0;

                foreach (var direction in HRIRs)
                {
                    foreach (var channel in direction)
                    {
                        foreach (var sample in channel)
                        {
                            totalEnergy += sample * sample;
                        }
                    }
                }

                return totalEnergy;
            }

            /// <summary>
            /// Resamples a HRIR using the WDL Resampler from NAudio.
            /// </summary>
            /// <param name="hrirs"></param>
            /// <param name="srcFs"></param>
            /// <param name="targetFs"></param>
            /// <returns></returns>
            public static double[][] ResampleHRIRWDL(float[][] hrir, int srcFs, int targetFs)
            {
                int channels = hrir.Length;
                int samples = hrir[0].Length;

                // No resampling needed
                if (srcFs == targetFs)
                {
                    return hrir
                        .Select(channel => channel.Select(x => (double)x).ToArray())
                            .ToArray();
                }

                // Flatten channels into [samples, channels] float[,] for resampler
                float[,] multiChannel = new float[samples, channels];
                for (int ch = 0; ch < channels; ch++)
                    for (int i = 0; i < samples; i++)
                        multiChannel[i, ch] = hrir[ch][i];

                // Set up resampler
                var sourceProvider = new MutableArraySampleProvider(srcFs, channels);
                sourceProvider.SetSamples(multiChannel);
                var resampler = new WdlResamplingSampleProvider(sourceProvider, targetFs);

                int outputFrames = (int)(samples * (double)targetFs / srcFs + 0.5);
                float[] tempBuffer = new float[1024 * channels];
                float[] outputBuffer = new float[outputFrames * channels];

                int totalRead = 0;
                while (totalRead < outputBuffer.Length)
                {
                    int read = resampler.Read(tempBuffer, 0, Math.Min(tempBuffer.Length, outputBuffer.Length - totalRead));
                    if (read == 0) break;
                    Array.Copy(tempBuffer, 0, outputBuffer, totalRead, read);
                    totalRead += read;
                }

                // If fewer samples were read, clear the rest
                if (totalRead < outputBuffer.Length)
                    Array.Clear(outputBuffer, totalRead, outputBuffer.Length - totalRead);

                // Convert to double[][] per channel
                double[][] resampled = new double[channels][];
                for (int ch = 0; ch < channels; ch++)
                {
                    resampled[ch] = new double[outputFrames];
                    for (int i = 0; i < outputFrames; i++)
                        resampled[ch][i] = outputBuffer[i * channels + ch];
                }

                return resampled;

            }

            /// <summary>
            /// Converts time-domain signal to minimum-phase using the cepstral method.
            /// </summary>
            /// <param name="signal"></param>
            /// <param name="threadId"></param>
            /// <returns></returns>
            /// <exception cref="ArgumentException"></exception>
            public static double[] ComputeMinimumPhaseSignal(double[] signal, int threadId)
            {
                if (signal == null || signal.Length == 0)
                {
                    throw new ArgumentException("HRIR is null or empty. Could not continue to minimum phase conversion.");
                }

                int n = signal.Length;

                // FFT
                Complex[] S = Pach_SP.FFT_General(signal, threadId);

                // Log magnitude
                for (int i = 0; i < S.Length; i++)
                {
                    double mag = Math.Max(S[i].Magnitude, 1e-12);
                    S[i] = new Complex(Math.Log(mag), 0);
                }

                // IFFT
                double[] cepstrum = Pach_SP.IFFT_Real_General(S, threadId);

                // Normalise IFFT
                for (int i = 0; i < n; i++)
                    cepstrum[i] /= n;

                // Keep only causal part
                int half = (cepstrum.Length + 1) / 2;
                for (int i = half; i < cepstrum.Length; i++)
                    cepstrum[i] = 0.0;
                for (int i = 1; i < half - 1; i++)
                    cepstrum[i] *= 2.0;

                // FFT of modified cepstrum
                Complex[] S_min = Pach_SP.FFT_General(cepstrum, threadId);

                // Exponentiate
                for (int i = 0; i < S_min.Length; i++)
                    S_min[i] = Complex.Exp(S_min[i]);

                // IFFT to get minimum-phase signal
                double[] minPhaseSignal = Pach_SP.IFFT_Real_General(S_min, threadId);

                // Normalise IFFT
                for (int i = 0; i < n; i++)
                    minPhaseSignal[i] /= n;

                // Scale to match original energy
                double energyOriginal = signal.Select(x => x * x).Sum();
                double energyMinPhase = minPhaseSignal.Select(x => x * x).Sum();
                double scale = Math.Sqrt(energyOriginal / (energyMinPhase + 1e-20));
                for (int i = 0; i < n; i++)
                    minPhaseSignal[i] *= scale;

                return minPhaseSignal;
            }

            /// <summary>
            /// Applies a fractional delay via windowed sinc interpolation.
            /// </summary>
            /// <param name="signal"></param>
            /// <param name="delaySeconds"></param>
            /// <param name="Fs"></param>
            /// <param name="kernelSize"></param>
            /// <returns></returns>
            /// <exception cref="ArgumentNullException"></exception>
            /// <exception cref="ArgumentException"></exception>
            public static double[] ApplyFractionalDelaySinc(double[] signal, double delaySeconds, int Fs, int kernelSize = 15)
            {
                if (signal == null) throw new ArgumentNullException(nameof(signal));
                if (Fs <= 0) throw new ArgumentException("Fs must be > 0", nameof(Fs));
                if (kernelSize < 3) kernelSize = 3;
                if (kernelSize % 2 == 0) kernelSize++;

                int len = signal.Length;
                double delaySamples = delaySeconds * Fs;
                int N = (int)Math.Floor(delaySamples);
                double frac = delaySamples - N;

                int R = kernelSize / 2;

                double[] taps = new double[2 * R + 1];
                for (int k = -R; k <= R; k++)
                {
                    double x = k - frac;
                    double sinc = (Math.Abs(x) < 1e-16) ? 1.0 : Math.Sin(Math.PI * x) / (Math.PI * x);

                    // Hann window centered on 0 with radius R
                    double w = 0.5 * (1.0 + Math.Cos(Math.PI * k / (double)R));

                    double val = sinc * w;
                    taps[k + R] = val;
                }

                // Normalise taps to unity DC
                double norm = 0.0;
                for (int i = 0; i < taps.Length; i++) norm += taps[i];
                if (Math.Abs(norm) > 1e-16)
                {
                    for (int i = 0; i < taps.Length; i++) taps[i] /= norm;
                }

                // Output
                double[] y = new double[len];

                for (int n = 0; n < len; n++)
                {
                    double sum = 0.0;
                    for (int k = -R; k <= R; k++)
                    {
                        int idx = n - N + k;
                        int safeIdx = MirrorIndex(idx, len);
                        sum += taps[k + R] * signal[safeIdx];
                    }
                    y[n] = sum;
                }

                return y;
            }

            /// <summary>
            /// Applies a fractional delay via Lagrange interpolation.
            /// </summary>
            /// <param name="h"></param>
            /// <param name="delay"></param>
            /// <param name="Sample_Freq"></param>
            /// <param name="thread"></param>
            /// <returns></returns>
            public static double[] ApplyFractionalDelayLagrange(double[] signal, double delaySeconds, int sampleFreq, int lagrangeOrder = 3)
            {
                double delaySamples = delaySeconds * sampleFreq;
                int N = (int)Math.Floor(delaySamples);
                double delta = delaySamples - N;

                int M = 3;
                int halfM = M / 2;
                int tapsCount = M + 1;
                double[] w = new double[tapsCount];
                int kStart = -halfM;
                int kEnd = kStart + tapsCount - 1;

                for (int k = kStart; k <= kEnd; k++)
                {
                    double prod = 1.0;
                    for (int m = kStart; m <= kEnd; m++)
                    {
                        if (m == k) continue;
                        prod *= (delta - m) / (k - m);
                    }
                    w[k - kStart] = prod;
                }

                int len = signal.Length;
                double[] y = new double[len];

                for (int n = 0; n < len; n++)
                {
                    double sum = 0.0;
                    for (int k = kStart; k <= kEnd; k++)
                    {
                        int idx = n - N - k;
                        int safeIdx = MirrorIndex(idx, len);
                        sum += w[k - kStart] * signal[safeIdx];
                    }
                    y[n] = sum;
                }

                return y;
            }

            public static double[] ApplyFractionalDelay(double[] signal, double delaySeconds, int Fs, bool useSinc = true, int kernelSize = 31, int lagrangeOrder = 3)
            {
                if (useSinc)
                    return ApplyFractionalDelaySinc(signal, delaySeconds, Fs, kernelSize);
                else
                    return ApplyFractionalDelayLagrange(signal, delaySeconds, Fs, lagrangeOrder);
            }

            public static double[] BuildDirectionalComponent(int i, double[][] filt, Vector[] Directions, double[][] Translation)
            {
                double[] s = new double[filt.Length];

                for (int t = 0; t < filt.Length; t++)
                {
                    s[t] += (Directions[i].dx > 0 ? filt[t][0] : filt[t][1]) * Translation[i][0];
                    s[t] += (Directions[i].dy > 0 ? filt[t][2] : filt[t][3]) * Translation[i][1];
                    s[t] += (Directions[i].dz > 0 ? filt[t][4] : filt[t][5]) * Translation[i][2];
                }

                return s;
            }

            public static double[] BuildDrySignal(double[][] directionalSignals)
            {
                int length = directionalSignals.First(s => s != null).Length;
                double[] drySignal = new double[length];

                for (int i = 0; i < directionalSignals.Length; i++)
                {
                    if (directionalSignals[i] == null) continue;

                    for (int t = 0; t < length; t++)
                        drySignal[t] += directionalSignals[i][t];
                }

                return drySignal;
            }

            public static double ComputeRMS(double[] signal)
            {
                if (signal == null || signal.Length == 0)
                    return 0.0;

                double sumSquares = 0.0;
                for (int t = 0; t < signal.Length; t++)
                    sumSquares += signal[t] * signal[t];

                return Math.Sqrt(sumSquares / signal.Length);
            }

            public static double[][] SumAcrossDirections(double[][][] signals)
            {
                if (signals == null || signals.Length == 0)
                    return new double[2][] { new double[0], new double[0] };

                int numDirs = signals.Length;

                // Determine max length across all directions/channels
                int maxLength = 0;
                for (int i = 0; i < numDirs; i++)
                {
                    if (signals[i] == null) continue;
                    maxLength = Math.Max(maxLength, signals[i][0]?.Length ?? 0);
                    maxLength = Math.Max(maxLength, signals[i][1]?.Length ?? 0);
                }

                double[] leftSum = new double[maxLength];
                double[] rightSum = new double[maxLength];

                for (int i = 0; i < numDirs; i++)
                {
                    if (signals[i] == null) continue;

                    var left = signals[i][0];
                    var right = signals[i][1];

                    if (left != null)
                    {
                        for (int n = 0; n < left.Length; n++)
                            leftSum[n] += left[n];
                    }

                    if (right != null)
                    {
                        for (int n = 0; n < right.Length; n++)
                            rightSum[n] += right[n];
                    }
                }

                return new double[][] { leftSum, rightSum };
            }

            public static void NormaliseStereoByDryRMS(double[][] Signal, double dryRMS)
            {
                if (Signal.Length != 2)
                    throw new ArgumentException("Signal must have 2 channels (stereo).");

                int length = Signal[0].Length;
                double outSumSquares = 0;

                for (int t = 0; t < length; t++)
                {
                    double avgSample = 0.5 * (Signal[0][t] + Signal[1][t]);
                    outSumSquares += avgSample * avgSample;
                }

                double outRMS = Math.Sqrt(outSumSquares / length);
                double gainFactor = dryRMS / (outRMS + 1e-12);

                for (int ch = 0; ch < 2; ch++)
                    for (int t = 0; t < length; t++)
                        Signal[ch][t] *= gainFactor;
            }

            private static int MirrorIndex(int idx, int len)
            {
                while (idx < 0 || idx >= len)
                {
                    if (idx < 0)
                        idx = -idx - 1;
                    else
                        idx = 2 * len - idx - 1;
                }
                return idx;
            }

            public static void ShiftHRIRPairs(double[][][] hrirData, double frac)
            {
                int numDirs = hrirData.Length;
                int numCh = 2; // left/right
                int maxSamples = 0;

                // Find the maximum length across all HRIRs
                for (int d = 0; d < numDirs; d++)
                {
                    if (hrirData[d] == null) continue;
                    for (int ch = 0; ch < numCh; ch++)
                    {
                        if (hrirData[d][ch] != null)
                            maxSamples = Math.Max(maxSamples, hrirData[d][ch].Length);
                    }
                }

                // Process each direction
                for (int d = 0; d < numDirs; d++)
                {
                    if (hrirData[d] == null) continue;

                    double[] hL = hrirData[d][0];
                    double[] hR = hrirData[d][1];
                    if (hL == null || hR == null) continue;

                    // Detect first significant peaks
                    int idxPeakL = DetectFirstSignificantPeak(hL, frac);
                    int idxPeakR = DetectFirstSignificantPeak(hR, frac);

                    // Identify earliest ear
                    int earliestCh = idxPeakL <= idxPeakR ? 0 : 1;
                    int earliestPeakIdx = Math.Min(idxPeakL, idxPeakR);

                    double[] earliestH = earliestCh == 0 ? hL : hR;
                    double[] laterH = earliestCh == 0 ? hR : hL;

                    // Find last zero crossing before peak
                    int idxZero = -1;
                    for (int i = earliestPeakIdx - 1; i > 0; i--)
                    {
                        if (earliestH[i - 1] * earliestH[i] <= 0.0)
                        {
                            idxZero = i;
                            break;
                        }
                    }
                    if (idxZero < 0) idxZero = Math.Max(earliestPeakIdx - 1, 0);

                    int shiftSamples = Math.Max(idxZero - 1, 0);

                    // Crop earliest ear
                    double[] earliestShifted = new double[earliestH.Length - shiftSamples];
                    Array.Copy(earliestH, shiftSamples, earliestShifted, 0, earliestShifted.Length);

                    // Fractional shift
                    double y0 = earliestShifted[Math.Max(earliestPeakIdx - shiftSamples - 1, 0)];
                    double y1 = earliestShifted[earliestPeakIdx - shiftSamples];
                    double denom = y0 - y1;
                    double fracShift = (Math.Abs(denom) > 1e-12) ? y0 / denom : 0.0;
                    fracShift = Clamp(fracShift, -2.0, 2.0);

                    earliestShifted = ApplyFractionalDelay(earliestShifted, -fracShift, 1, useSinc: true);

                    // Shift later ear by same integer amount
                    double[] laterShifted = new double[laterH.Length - shiftSamples];
                    Array.Copy(laterH, shiftSamples, laterShifted, 0, laterShifted.Length);

                    // Pad both arrays to maxSamples and overwrite original HRIRs
                    double[] earliestPadded = new double[maxSamples];
                    double[] laterPadded = new double[maxSamples];

                    Array.Copy(earliestShifted, earliestPadded, earliestShifted.Length);
                    Array.Copy(laterShifted, laterPadded, laterShifted.Length);

                    hrirData[d][earliestCh] = earliestPadded;
                    hrirData[d][1 - earliestCh] = laterPadded;
                }
            }

            public static double Clamp(double value, double min, double max)
            {
                if (value < min) return min;
                if (value > max) return max;
                return value;
            }

            public static int DetectFirstSignificantPeak(double[] h, double frac)
            {
                int N = h.Length;
                double[] env = h.Select(Math.Abs).ToArray();
                double thresh = frac * env.Max();

                var localMaxIdxs = new System.Collections.Generic.List<int>();
                for (int i = 1; i < N - 1; i++)
                {
                    if (env[i] > env[i - 1] && env[i] >= env[i + 1])
                    {
                        localMaxIdxs.Add(i);
                    }
                }

                var significantPeaks = localMaxIdxs.Where(idx => env[idx] >= thresh);
                return significantPeaks.Any() ? significantPeaks.First() : -1;
            }

            /// <summary>
            /// 
            /// Equalises a set of HRIRs to compensate for the measurement system.
            /// 
            /// References:
            /// Gardner, W.G., 1997. 3-D Audio Using Loudspeakers. PhD thesis. Massachusetts Institute of Technology, Program in Media Arts and Sciences. Available at: https://dspace.mit.edu/handle/1721.1/29134
            /// Blauert, J., 1983. Spatial Hearing.Cambridge, MA: MIT Press.
            /// Moller, H., 1992. ‘Fundamentals of binaural technology’, Applied Acoustics, 36, pp. 171–218.
            /// Jot, J. M., Larcher, V.and Warusfel, O., 1995. ‘Digital signal processing issues in the context of binaural and transaural stereophony’, Proceedings of the Audio Engineering Society Convention, Preprint 3980.
            ///
            /// </summary>
            /// <param name="hrirs"></param>
            /// <param name="Directions"></param>
            /// <param name="Fs"></param>
            /// <param name="threadId"></param>
            /// <param name="settings"></param>
            /// <param name="auto"></param>
            public static void ApplySystemCompensation(double[][][] hrirs, Vector[] Directions, int Fs, int threadId, SystemResponseCompensation.SystemCompensationSettings settings, bool auto)
            {
                if (auto) // Automatic routine used for calculation of IACC
                {
                    SystemResponseCompensation.DiffuseFieldEqualisation(
                        hrirs,
                        threadId,
                        Fs,
                        1.0 / 3.0,
                        12.0,
                        true,
                        50.0
                    );

                    settings.ConvertToDTF = false;
                }
                else
                {
                    switch (settings.SelectedEQ)
                    {
                        case 1:
                            SystemResponseCompensation.MeasurementEqualisation(
                                hrirs,
                                settings.TFReference,
                                Fs,
                                threadId,
                                settings.IsCalibrated,
                                settings.SmoothingOct
                            );
                            break;

                        case 2:
                            SystemResponseCompensation.FreeFieldEqualisation(
                                hrirs,
                                Directions,
                                Fs,
                                threadId,
                                settings.FreeFieldIncidence,
                                settings.MaxBoostDb,
                                settings.LowFreqHz,
                                1e-6,
                                settings.SmoothingOct
                            );
                            break;

                        case 3:
                            SystemResponseCompensation.DiffuseFieldEqualisation(
                                hrirs,
                                threadId,
                                Fs,
                                settings.SmoothingOct,
                                settings.MaxBoostDb,
                                settings.MinPhase,
                                settings.LowFreqHz
                            );
                            break;

                        case 0:
                        default:
                            break;
                    }
                }

                if (settings.ConvertToDTF)
                {
                    SystemResponseCompensation.ConvertToDTF(hrirs, threadId, Fs);
                }
            }
        }

        public class SystemResponseCompensation
        {
            public struct SystemCompensationSettings
            {
                public int SelectedEQ;
                public double? SmoothingOct;
                public double? MaxBoostDb;
                public double? LowFreqHz;
                public bool MinPhase;
                public bool IsCalibrated;
                public double[] TFReference; // For measurement EQ
                public bool ConvertToDTF;
                public bool InputIsDTF;
                public int FreeFieldIncidence;
            }

            public static void MeasurementEqualisation(
                double[][][] hrirs,
                double[] tfReference,
                int Fs,
                int threadId,
                bool isCalibrated,
                double? smoothingOct = null,
                double regFactor = 1e-6)
            {
                int numDirs = hrirs.Length;

                // Precompute reference spectrum once
                int maxN = 0;
                for (int dir = 0; dir < numDirs; dir++)
                {
                    int len = hrirs[dir][0].Length;
                    if (len > maxN) maxN = len;
                }

                int fftLen = NextPowerOfTwo(maxN + tfReference.Length - 1);

                double[] paddedTF = new double[fftLen];
                Array.Copy(tfReference, paddedTF, tfReference.Length);
                Complex[] tfRefFreq = Pach_SP.FFT_General(paddedTF, threadId);

                int nFreqs = fftLen / 2 + 1;

                double[] tfMag = new double[nFreqs];
                for (int i = 0; i < nFreqs; i++)
                    tfMag[i] = tfRefFreq[i].Magnitude;

                if (smoothingOct.HasValue)
                    tfMag = SmoothSpectrumOctave(tfMag, Fs, fftLen, smoothingOct.Value);

                for (int i = 0; i < nFreqs; i++)
                {
                    double phase = tfRefFreq[i].Phase;
                    tfRefFreq[i] = Complex.FromPolarCoordinates(tfMag[i], phase);
                }

                // Regularise / floor the reference magnitude
                double maxRefPower = 0.0;
                for (int i = 0; i < fftLen; i++)
                {
                    double mag = tfRefFreq[i].Magnitude;
                    double pow = mag * mag;
                    if (pow > maxRefPower) maxRefPower = pow;
                }

                double minReg = 1e-20;
                double reg = Math.Max(maxRefPower * regFactor, minReg);

                for (int dir = 0; dir < numDirs; dir++)
                {
                    if (hrirs[dir] == null) continue;
                    int channels = hrirs[dir].Length;
                    if (channels == 0) continue;

                    int N = hrirs[dir][0].Length;
                    for (int ch = 0; ch < channels; ch++)
                    {
                        double[] orig = hrirs[dir][ch];
                        if (orig == null || orig.Length == 0) continue;

                        double origRms = !isCalibrated ? Pach_SP_HRTF.ComputeRMS(orig) : 0.0;

                        double[] paddedHRIR = new double[fftLen];
                        Array.Copy(orig, paddedHRIR, N);
                        Complex[] hrirFreq = Pach_SP.FFT_General(paddedHRIR, threadId);

                        // Regularised inverse (Wiener/Tikhonov style):
                        Complex[] eqFreq = new Complex[fftLen];
                        for (int i = 0; i < fftLen; i++)
                        {
                            Complex H = hrirFreq[i];
                            Complex G = tfRefFreq[i];

                            double denom = G.Real * G.Real + G.Imaginary * G.Imaginary + reg;
                            eqFreq[i] = (H * Complex.Conjugate(G)) / denom;
                        }

                        // IFFT back to time domain (real)
                        double[] eqTime = Pach_SP.IFFT_Real_General(eqFreq, threadId);

                        for (int n = 0; n < fftLen; n++)
                            eqTime[n] /= fftLen;

                        // If requested, preserve original RMS
                        if (!isCalibrated)
                        {
                            double newRms = 0.0;
                            double sumSquares = 0.0;
                            for (int n = 0; n < N; n++) // Only loop through the first N samples
                            {
                                sumSquares += eqTime[n] * eqTime[n];
                            }

                            newRms = Math.Sqrt(sumSquares / N);

                            if (newRms > 0.0 && origRms > 0.0)
                            {
                                double scale = origRms / newRms;
                                for (int n = 0; n < N; n++) eqTime[n] *= scale;
                            }
                        }

                        Array.Copy(eqTime, 0, hrirs[dir][ch], 0, N);
                    }
                }
            }

            /// <summary>
            /// Free-field equalization: normalize HRTFs with respect to a particular incident direction. Usually frontal to 30 degrees azimuth.
            /// 
            /// Things to be improved:
            /// Use a frequency-dependent regularisation factor.
            /// Consider also non-uniform spatial sampling density, similar to diffuse-field EQ.
            /// 
            /// </summary>
            /// <param name="hrirs">[direction][ear][samples]</param>
            /// <param name="directionIndex">Index of the reference direction (e.g., frontal)</param>
            /// <param name="threadId">Thread for FFT calls</param>
            public static void FreeFieldEqualisation(
                double[][][] hrirs,
                Vector[] Directions,
                int Fs,
                int threadId,
                int FreeFieldIncidence,
                double? maxBoostdB = null,
                double? lowFreqHz = null,
                double regFactor = 1e-6,
                double? smoothingOct = null
                )
            {
                if (hrirs == null || hrirs.Length == 0 || Directions == null || Directions.Length == 0)
                    return;

                int nDirs = hrirs.Length;
                int nEars = 2;
                int fftSize = hrirs[0][0].Length;
                double fRes = Fs / (double)fftSize;
                int lowFreqBin = (lowFreqHz.HasValue)
                    ? (int)Math.Round(lowFreqHz.Value / fRes)
                    : 0;

                double azRad = FreeFieldIncidence * Math.PI / 180.0;
                double elRad = 0.0; // zero elevation
                Vector referenceDir = new Vector(
                    Math.Cos(elRad) * Math.Cos(azRad),
                    Math.Cos(elRad) * Math.Sin(azRad),
                    Math.Sin(elRad)
                );

                // Find the HRIR index closest to incident direction (generally frontal)
                int referenceIndex = 0;
                double minAngle = double.MaxValue;
                for (int i = 0; i < Directions.Length; i++)
                {
                    double angle = HRTF.AngularDistanceDegrees(referenceDir, Directions[i]);
                    if (angle < minAngle)
                    {
                        minAngle = angle;
                        referenceIndex = i;
                    }
                }

                // FFT of reference HRIR (closest to incident direction)
                Complex[] refLeft = Pach_SP.FFT_General(hrirs[referenceIndex][0], threadId);
                Complex[] refRight = Pach_SP.FFT_General(hrirs[referenceIndex][1], threadId);

                double[] refLeftMag = new double[fftSize];
                double[] refRightMag = new double[fftSize];
                for (int f = 0; f < fftSize; f++)
                {
                    refLeftMag[f] = refLeft[f].Magnitude;
                    refRightMag[f] = refRight[f].Magnitude;
                }

                if (smoothingOct.HasValue)
                {
                    refLeftMag = SmoothSpectrumOctave(refLeftMag, Fs, fftSize, smoothingOct.Value);
                    refRightMag = SmoothSpectrumOctave(refRightMag, Fs, fftSize, smoothingOct.Value);
                }

                // Compute max reference power (|G|^2) for relative regularisation
                double maxRefPower = 0.0;
                for (int f = 0; f < fftSize; f++)
                {
                    maxRefPower = Math.Max(maxRefPower, refLeftMag[f] * refLeftMag[f]);
                    maxRefPower = Math.Max(maxRefPower, refRightMag[f] * refRightMag[f]);
                }

                double minReg = 1e-20;
                double reg = Math.Max(maxRefPower * regFactor, minReg);

                // Apply equalisation to all HRIRs
                for (int dir = 0; dir < nDirs; dir++)
                {
                    if (hrirs[dir] == null) continue;
                    for (int ear = 0; ear < nEars; ear++)
                    {
                        Complex[] H = Pach_SP.FFT_General(hrirs[dir][ear], threadId);
                        Complex[] refFFT = (ear == 0) ? refLeft : refRight;

                        for (int f = 0; f < fftSize; f++)
                        {
                            double refMag = refFFT[f].Magnitude;
                            double regulatedMag = Math.Sqrt(refMag * refMag + reg);
                            double mag = H[f].Magnitude / regulatedMag;

                            if (maxBoostdB.HasValue)
                            {
                                double gaindB = 20.0 * Math.Log10(mag + 1e-300);
                                gaindB = Math.Max(-maxBoostdB.Value, Math.Min(maxBoostdB.Value, gaindB));
                                mag = Math.Pow(10.0, gaindB / 20.0);
                            }

                            if (lowFreqHz.HasValue && f < lowFreqBin)
                                mag = 1.0;

                            H[f] = Complex.FromPolarCoordinates(mag, H[f].Phase);
                        }

                        double[] equalisedHRIR = Pach_SP.IFFT_Real_General(H, threadId);
                        for (int n = 0; n < equalisedHRIR.Length; n++)
                            equalisedHRIR[n] /= fftSize;

                        double origRms = Pach_SP_HRTF.ComputeRMS(hrirs[dir][ear]);
                        double newRms = Pach_SP_HRTF.ComputeRMS(equalisedHRIR);

                        if (newRms > 0.0 && origRms > 0.0)
                        {
                            double scale = origRms / newRms;
                            for (int n = 0; n < hrirs[dir][ear].Length; n++)
                                equalisedHRIR[n] *= scale;
                        }

                        Array.Copy(equalisedHRIR, 0, hrirs[dir][ear], 0, hrirs[dir][ear].Length);
                    }
                }
            }

            /// <summary>
            /// Equalises the dataset with respect to the diffuse-field average of HRTFs across all incident directions for both ears.
            /// 
            /// For greater robustness, we could consider using a regularised Kirkeby inverse filter. See: https://github.com/Sinnerboy89/DFE/tree/master
            /// Diffuse-field equalisation really requires a uniform spherical distribution of incident directions.
            /// In practice, if the measured directions do not sample the sphere uniformly, each HRTF can be given a weight 
            /// in the averaging process (proportional to a solid angle associated to the corresponding direction).
            /// Plan for later work: Use MIConvexHull to compute Voronoi-based solid angle weights. Additionally, add a low-frequency 
            /// roll-off instead of hard cut-off.
            /// 
            /// References:
            /// M. Mein, "Perception de l'information binaurale liee aux reflexions precoces dans une salle. Application a la simulation de la qualite acoustique", 
            /// memoire de DEA, Univ. du Maine, Le Mans, Sept. 1993.
            /// 
            /// </summary>
            /// <param name="hrirs"></param>
            /// <param name="threadId"></param>
            /// <param name="Fs"></param>
            /// <param name="smoothingOct"></param>
            /// <param name="maxBoostdB"></param>
            /// <param name="minPhase"></param>
            /// <param name="lowFreqHz"></param>
            public static void DiffuseFieldEqualisation(double[][][] hrirs, int threadId, int Fs,
                                double? smoothingOct = 1.0 / 3.0,
                                double? maxBoostdB = 12.0,
                                bool minPhase = false,
                                double? lowFreqHz = 50.0)
            {
                if (hrirs == null || hrirs.Length == 0) return;
                int nDirs = hrirs.Length;
                int nEars = 2;
                int fftSize = hrirs[0][0].Length;
                int nFreqs = fftSize / 2 + 1;
                double fRes = Fs / (double)fftSize;
                int lowFreqBin = (lowFreqHz.HasValue)
                    ? (int)Math.Round(lowFreqHz.Value / fRes)
                    : 0;

                // Compute FFTs
                Complex[][][] H = new Complex[nDirs][][];
                for (int i = 0; i < nDirs; i++)
                {
                    if (hrirs[i] == null) continue;
                    H[i] = new Complex[nEars][];
                    for (int ch = 0; ch < nEars; ch++)
                        H[i][ch] = Pach_SP.FFT_General(hrirs[i][ch], threadId);
                }

                // Compute diffuse-field average per ear
                double[][] avgMag = new double[nEars][];
                for (int ch = 0; ch < nEars; ch++)
                {
                    avgMag[ch] = new double[nFreqs];
                    for (int f = 0; f < nFreqs; f++)
                    {
                        double sumPower = 0.0;
                        int count = 0;
                        for (int i = 0; i < nDirs; i++)
                        {
                            if (H[i] == null) continue;
                            double mag = H[i][ch][f].Magnitude;
                            sumPower += mag * mag;
                            count++;
                        }
                        avgMag[ch][f] = (count > 0) ? Math.Sqrt(sumPower / count) : 1e-12;
                    }

                    if (smoothingOct.HasValue)
                        avgMag[ch] = SmoothSpectrumOctave(avgMag[ch], Fs, fftSize, smoothingOct.Value);
                }

                // Build and apply EQ filter for each ear separately
                for (int ch = 0; ch < nEars; ch++)
                {
                    Complex[] eqFilterSpec = new Complex[fftSize];
                    for (int f = 0; f < nFreqs; f++)
                    {
                        double gaindB = -20.0 * Math.Log10(avgMag[ch][f] + 1e-12);

                        if (maxBoostdB.HasValue)
                            gaindB = Math.Max(-maxBoostdB.Value, Math.Min(maxBoostdB.Value, gaindB));

                        if (lowFreqHz.HasValue && f < lowFreqBin)
                            gaindB = 0.0;

                        double gainLin = Math.Pow(10.0, gaindB / 20.0);
                        eqFilterSpec[f] = new Complex(gainLin, 0.0);
                    }

                    for (int f = nFreqs; f < fftSize; f++)
                        eqFilterSpec[f] = Complex.Conjugate(eqFilterSpec[fftSize - f]);

                    // Minimum-phase or linear-phase EQ filter
                    double[] eqFilter;
                    if (minPhase)
                        eqFilter = Pach_SP.Minimum_Phase_Complex(eqFilterSpec, threadId);
                    else
                        eqFilter = Pach_SP.IFFT_Real_General(eqFilterSpec, threadId);

                    for (int n = 0; n < eqFilter.Length; n++)
                        eqFilter[n] /= fftSize;

                    // Apply EQ filter to all HRIRs for this ear
                    for (int i = 0; i < nDirs; i++)
                    {
                        if (H[i] == null) continue;
                        double[] filtered = Pach_SP.FFT_Convolution_double(hrirs[i][ch], eqFilter, 0);

                        // RMS normalization
                        double origRms = Pach_SP_HRTF.ComputeRMS(hrirs[i][ch]);
                        double newRms = Pach_SP_HRTF.ComputeRMS(filtered);

                        if (newRms > 0.0 && origRms > 0.0)
                        {
                            double scale = origRms / newRms;
                            for (int n = 0; n < filtered.Length; n++)
                                filtered[n] *= scale;
                        }

                        hrirs[i][ch] = filtered;
                    }
                }
            }

            private static double[] SmoothSpectrumOctave(double[] mag, int sampleRate, int fftSize, double oct)
            {
                int nFreqs = mag.Length;
                double[] smoothed = new double[nFreqs];
                double fRes = sampleRate / (double)fftSize;

                for (int f = 0; f < nFreqs; f++)
                {
                    double freq = f * fRes;
                    if (freq <= 0) { smoothed[f] = mag[f]; continue; }
                    double sum = 0.0, wSum = 0.0;
                    for (int k = 0; k < nFreqs; k++)
                    {
                        double freqK = k * fRes;
                        if (freqK <= 0) continue;
                        double logDist = Math.Log(freq / freqK) / Math.Log(2.0);
                        double sigma = oct / 2.0;
                        double w = Math.Exp(-0.5 * (logDist / sigma) * (logDist / sigma));
                        sum += mag[k] * w;
                        wSum += w;
                    }
                    smoothed[f] = (wSum > 0) ? sum / wSum : mag[f];
                }
                return smoothed;
            }

            /// <summary>
            /// Resamples a mono time-domain signal using the WDL ResamplingSampleProvider from NAudio. This is used for resampling the
            /// measurement system compensation transfer function, which is typically a mono signal.
            /// </summary>
            /// <param name="tfSignal"></param>
            /// <param name="srcFs"></param>
            /// <param name="targetFs"></param>
            /// <returns></returns>
            public static double[] ResampleMonoTFWdl(double[] tfSignal, int srcFs, int targetFs)
            {
                if (srcFs == targetFs)
                    return (double[])tfSignal.Clone(); // no resampling needed

                int channels = 1;
                int samples = tfSignal.Length;

                // Convert double -> float
                float[] inputFloat = tfSignal.Select(x => (float)x).ToArray();

                // Wrap in a simple mono sample provider
                var sourceProvider = new MutableArraySampleProvider(srcFs, channels);
                sourceProvider.SetSamples(new float[,] { { 0, } }); // placeholder
                                                                    // Since it's mono, reshape to [samples, channels]
                float[,] mono2D = new float[samples, 1];
                for (int i = 0; i < samples; i++)
                    mono2D[i, 0] = inputFloat[i];
                sourceProvider.SetSamples(mono2D);

                // Create the resampler
                var resampler = new WdlResamplingSampleProvider(sourceProvider, targetFs);

                int outputFrames = (int)(samples * (double)targetFs / srcFs + 0.5);
                float[] outputBuffer = new float[outputFrames]; // mono ? no interleaving
                int totalRead = 0;
                int chunkSize = 1024;
                float[] tempBuffer = new float[chunkSize];

                while (totalRead < outputFrames)
                {
                    int read = resampler.Read(tempBuffer, 0, Math.Min(chunkSize, outputFrames - totalRead));
                    if (read == 0) break;
                    Array.Copy(tempBuffer, 0, outputBuffer, totalRead, read);
                    totalRead += read;
                }

                // Convert back to double
                double[] resampledTF = outputBuffer.Take(totalRead).Select(x => (double)x).ToArray();

                return resampledTF;
            }

            /// <summary>
            /// Converts a set of HRTFs to Directional Transfer Functions (DTFs). Middlebrooks explicitly removes the measurement system response before this step,
            /// so the user is made to apply some form of system equalisation before calling this function.
            ///
            /// References:
            /// 
            /// The Institut fur Schallforschung, Osterreichische Akademie der Wissenschaften give an explanation of the HRTF to DTF conversion process on their website, 
            /// as supporting material for the ARI HRTF database:
            /// https://www.oeaw.ac.at/isf/outreach/software/hrtf-database/detailed-description
            /// 
            /// Middlebrooks, J.C., 1999. Individual differences in external-ear transfer functions reduced by scaling in frequency. 
            /// The Journal of the Acoustical Society of America, 106(3), pp.1480–1491. https://doi.org/10.1121/1.427176
            /// 
            /// </summary>
            /// <param name="hrirs"></param>
            /// <param name="threadId"></param>
            /// <param name="Fs"></param>
            /// <param name="smoothingOct"></param>
            public static void ConvertToDTF(double[][][] hrirs, int threadId, int Fs, double? smoothingOct = null)
            {
                if (hrirs == null || hrirs.Length == 0) return;

                int nDirs = hrirs.Length;
                int nEars = hrirs[0].Length;
                int fftLen = hrirs[0][0].Length;
                int nFreqs = fftLen / 2 + 1;

                // Compute FFTs
                Complex[][][] H = new Complex[nDirs][][];
                for (int dir = 0; dir < nDirs; dir++)
                {
                    H[dir] = new Complex[nEars][];
                    for (int ear = 0; ear < nEars; ear++)
                        H[dir][ear] = Pach_SP.FFT_General(hrirs[dir][ear], threadId);
                }

                // Compute average log amplitude spectrum
                double[][] avgLogMag = new double[nEars][];
                for (int ear = 0; ear < nEars; ear++)
                {
                    avgLogMag[ear] = new double[nFreqs];
                    for (int f = 0; f < nFreqs; f++)
                    {
                        double sumLogMag = 0.0;
                        for (int dir = 0; dir < nDirs; dir++)
                        {
                            // Convert magnitude to log scale
                            double mag = H[dir][ear][f].Magnitude;
                            sumLogMag += Math.Log(mag + 1e-12);
                        }
                        avgLogMag[ear][f] = sumLogMag / nDirs;
                    }

                    if (smoothingOct.HasValue)
                        avgLogMag[ear] = SmoothSpectrumOctave(avgLogMag[ear], Fs, fftLen, smoothingOct.Value);
                }

                // Convert average log amplitude spectrum to minimum-phase complex spectrum
                double[][] H_minPhase = new double[nEars][];
                for (int ear = 0; ear < nEars; ear++)
                {
                    H_minPhase[ear] = Pach_SP_HRTF.ComputeMinimumPhaseSignal(avgLogMag[ear], threadId);
                }

                // Filter HRTFs with inverse CTF to obtain DTFs
                for (int dir = 0; dir < nDirs; dir++)
                {
                    double rmsLeft = Pach_SP_HRTF.ComputeRMS(hrirs[dir][0]);
                    double rmsRight = Pach_SP_HRTF.ComputeRMS(hrirs[dir][1]);

                    for (int ear = 0; ear < nEars; ear++)
                    {
                        Complex[] H_DTF = new Complex[fftLen];
                        for (int f = 0; f < nFreqs; f++)
                        {
                            double epsilon = 1e-12;
                            H_DTF[f] = H[dir][ear][f] / (new Complex(H_minPhase[ear][f] + epsilon, 0.0));
                        }

                        // Mirror for negative frequencies
                        for (int f = nFreqs; f < fftLen; f++)
                            H_DTF[f] = Complex.Conjugate(H_DTF[fftLen - f]);

                        // Back to time domain
                        double[] timeDomain = Pach_SP.IFFT_Real_General(H_DTF, threadId);

                        for (int i = 0; i < timeDomain.Length; i++)
                            timeDomain[i] /= fftLen;

                        hrirs[dir][ear] = timeDomain;
                    }

                    double newRmsLeft = Pach_SP_HRTF.ComputeRMS(hrirs[dir][0]);
                    double newRmsRight = Pach_SP_HRTF.ComputeRMS(hrirs[dir][1]);

                    double scale = Math.Sqrt((rmsLeft * rmsLeft + rmsRight * rmsRight) / (newRmsLeft * newRmsLeft + newRmsRight * newRmsRight));

                    // Apply scaling to both ears to preserve overall energy
                    for (int i = 0; i < hrirs[dir][0].Length; i++)
                    {
                        hrirs[dir][0][i] *= scale;
                        hrirs[dir][1][i] *= scale;
                    }
                }
            }

            public static int NextPowerOfTwo(int n)
            {
                if (n < 1)
                    return 1;

                n--;
                n |= n >> 1;
                n |= n >> 2;
                n |= n >> 4;
                n |= n >> 8;
                n |= n >> 16;
                n++;
                return n;
            }
        }
    }
}