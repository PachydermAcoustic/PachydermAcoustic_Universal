//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2026, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
//'Pachyderm-Acoustic is free software; you can redistribute it and/or modify 
//'it under the terms of the GNU General Public License as published 
//'by the Free Software Foundation; either version 3 of the License, or 
//'(at your option) any later version. 
//'Pachyderm-Acoustic is distributed in the hope that it will be useful, 
//'but WITHOUT ANY WARRANTY; without even the implied warranty of 
//'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
//'GNU General Public License for more details. 
//' 
//'You should have received a copy of the GNU General Public 
//'License along with Pachyderm-Acoustic; if not, write to the Free Software 
//'Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 

using System;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using Pachyderm_Acoustic.Environment;

namespace Pachyderm_Acoustic
{
    /// <summary>Shared source phase references, prepared before tracing coherent speaker paths.</summary>
    public sealed class Coherent_Source_Waveforms : IDisposable
    {
        private const int Length = 16384;
        private const int SampleRate = 44100;
        private readonly double[][] spectra;
        private readonly long[] offsets;
        private readonly Audio.Pach_SP.DSP_System phaseSystem;
        private FileStream file;
        private MemoryMappedFile mapping;
        private MemoryMappedViewAccessor view;
        public int BandCount { get; private set; }
        internal double TimeOriginSeconds { get { return phaseSystem is Audio.Pach_SP.Linear_Phase_System ? Length * .5 / SampleRate : 0; } }

        public Coherent_Source_Waveforms(SourceCluster source, bool thirdOctave)
            : this(source.Sources.Select(s => s.IsSpeaker ? Utilities.AcousticalMath.ConvertBandEnergy(s.SoundPower, thirdOctave).Select(Math.Sqrt).ToArray() : null).ToArray(), Audio.Pach_SP.Filter)
        { }

        private Coherent_Source_Waveforms(double[][] spectra, Audio.Pach_SP.DSP_System phaseSystem)
        {
            this.spectra = spectra;
            this.phaseSystem = phaseSystem;
            BandCount = spectra.First(s => s != null).Length;
            offsets = new long[spectra.Length];
            long capacity = 0;
            var references = new System.Collections.Generic.Dictionary<string, long>();
            for (int source = 0; source < spectra.Length; source++)
            {
                if (spectra[source] == null) continue;
                double maximum = spectra[source].Max();
                string key = string.Join(",", spectra[source].Select(p => (maximum > 0 ? p / maximum : 0).ToString("R", System.Globalization.CultureInfo.InvariantCulture)));
                if (!references.TryGetValue(key, out offsets[source]))
                {
                    offsets[source] = capacity;
                    references.Add(key, capacity);
                    capacity += (long)BandCount * Length * sizeof(double);
                }
            }
            try
            {
                string path = Path.Combine(Path.GetTempPath(), "Pachyderm-source-" + Guid.NewGuid().ToString("N") + ".tmp");
                file = new FileStream(path, FileMode.CreateNew, FileAccess.ReadWrite, FileShare.Read | FileShare.Delete, 4096, FileOptions.DeleteOnClose);
                file.SetLength(capacity);
                mapping = MemoryMappedFile.CreateFromFile(file, null, capacity, MemoryMappedFileAccess.ReadWrite, HandleInheritability.None, true);
                view = mapping.CreateViewAccessor(0, capacity, MemoryMappedFileAccess.ReadWrite);
                var prepared = new System.Collections.Generic.HashSet<long>();
                for (int source = 0; source < spectra.Length; source++)
                {
                    if (spectra[source] == null || !prepared.Add(offsets[source])) continue;
                    double maximum = spectra[source].Max();
                    double[] response = phaseSystem.Transfer_Function(spectra[source].Select(p => maximum > 0 ? p / maximum : 0).ToArray(), SampleRate, Length, 0);
                    for (int band = 0; band < BandCount; band++)
                    {
                        double[] kernel = BandCount == 24 ? Audio.Pach_SP.FIR_Third_Bandpass(response, band, SampleRate, 0) : Audio.Pach_SP.FIR_Bandpass(response, band, SampleRate, 0);
                        double norm = Math.Sqrt(kernel.Sum(p => p * p));
                        if (norm > 0) for (int i = 0; i < kernel.Length; i++) kernel[i] /= norm;
                        view.WriteArray(offsets[source] + (long)band * Length * sizeof(double), kernel, 0, Length);
                    }
                }
            }
            catch { Dispose(); throw; }
        }

        /// <summary>Apply received amplitudes and source/path delay to the fixed source phase.</summary>
        internal double[] Response(int source, double[] pressure, int band, double[] delays, double fraction, int sampleRate, int thread)
        {
            if (source < 0 || source >= spectra.Length || spectra[source] == null) throw new ArgumentOutOfRangeException(nameof(source));
            if (pressure.Length != BandCount) throw new ArgumentException("Source and path band counts must match.", nameof(pressure));
            double[] result = new double[1];
            int first = band < BandCount ? band : 0, last = band < BandCount ? band + 1 : BandCount;
            for (int b = first; b < last; b++)
            {
                if (pressure[b] == 0) continue;
                double[] kernel = new double[Length];
                if (view != null) view.ReadArray(offsets[source] + (long)b * Length * sizeof(double), kernel, 0, Length);
                else
                {
                    // Completed simulations retain only source descriptors; later edits can reconstruct the reference.
                    double maximum = spectra[source].Max();
                    kernel = phaseSystem.Transfer_Function(spectra[source].Select(p => maximum > 0 ? p / maximum : 0).ToArray(), SampleRate, Length, thread);
                    kernel = BandCount == 24 ? Audio.Pach_SP.FIR_Third_Bandpass(kernel, b, SampleRate, thread) : Audio.Pach_SP.FIR_Bandpass(kernel, b, SampleRate, thread);
                    double norm = Math.Sqrt(kernel.Sum(p => p * p));
                    if (norm > 0) for (int i = 0; i < kernel.Length; i++) kernel[i] /= norm;
                }
                if (sampleRate != SampleRate)
                {
                    kernel = Audio.Pach_SP.Resample(kernel, SampleRate, sampleRate, thread, false);
                    double norm = Math.Sqrt(kernel.Sum(p => p * p));
                    if (norm > 0) for (int i = 0; i < kernel.Length; i++) kernel[i] /= norm;
                }
                kernel = Audio.Pach_SP.Fractional_Delay(kernel, delays, fraction, sampleRate, thread);
                if (result.Length < kernel.Length) Array.Resize(ref result, kernel.Length);
                for (int i = 0; i < kernel.Length; i++) result[i] += kernel[i] * pressure[b];
            }
            return result;
        }

        internal void Write_Data(BinaryWriter writer)
        {
            writer.Write(phaseSystem is Audio.Pach_SP.Linear_Phase_System);
            writer.Write(BandCount);
            writer.Write(spectra.Length);
            foreach (double[] spectrum in spectra)
            {
                writer.Write(spectrum != null);
                if (spectrum != null) foreach (double pressure in spectrum) writer.Write(pressure);
            }
        }

        internal static Coherent_Source_Waveforms Read_Data(BinaryReader reader)
        {
            bool linear = reader.ReadBoolean();
            int bands = reader.ReadInt32(), count = reader.ReadInt32();
            if ((bands != 8 && bands != 24) || count < 1 || count > 100000) throw new InvalidDataException("Invalid coherent source reference dimensions.");
            double[][] spectra = new double[count][];
            for (int i = 0; i < count; i++)
                if (reader.ReadBoolean())
                {
                    spectra[i] = new double[bands];
                    for (int band = 0; band < bands; band++) spectra[i][band] = reader.ReadDouble();
                    if (spectra[i].Any(p => p < 0 || double.IsNaN(p) || double.IsInfinity(p))) throw new InvalidDataException("Invalid coherent source reference spectrum.");
                }
            if (!spectra.Any(s => s != null)) throw new InvalidDataException("No coherent speaker references.");
            return new Coherent_Source_Waveforms(spectra, linear ? (Audio.Pach_SP.DSP_System)new Audio.Pach_SP.Linear_Phase_System() : new Audio.Pach_SP.Minimum_Phase_System());
        }

        public void Dispose()
        {
            view?.Dispose(); view = null;
            mapping?.Dispose(); mapping = null;
            file?.Dispose(); file = null;
        }
    }
}
