using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Pachyderm_Acoustic
{
    namespace Audio
    {
        public partial class Pach_SP
        {
            public static class Measurement
            {
                public class IO_Tester
                {
                    NAudio.Wave.WaveOutEvent WO = new NAudio.Wave.WaveOutEvent();
                    NAudio.Wave.WaveInEvent WI;
                    List<Int16>[] Response;
                    double[] TD_Signal;
                    int Channels_in;
                    int CT_Averages = 1;
                    int block;
                    int SampleFreq = 44100;
                    double Signal_Length = 1.1;
                    public double[][] IR;
                    double Gain_out = 1;
                    public double[][] SNR;
                    public int[] Direct_Sample;
                    public bool Running = false;

                    public IO_Tester()
                    {
                    }

                    bool StopCalled = false;

                    void WI_RecordingStopped(object sender, EventArgs e)
                    {
                        if (StopCalled) return;
                        StopCalled = true;
                        IR = new double[Channels_in][];
                        Direct_Sample = new int[Channels_in];
                        SNR = new double[Channels_in][];
                        for (int c = 0; c < Channels_in; c++)
                        {
                            double[] Recording = new double[Response[c].Count];

                            for (int n = 0; n < Response[c].Count; n++) Recording[n] = (float)Response[c][n] / (float)Int16.MaxValue;

                            double[] DCR = Pach_SP.FFT_Deconvolution(Recording, TD_Signal);

                            int length = (int)(SampleFreq * Signal_Length);

                            IR[c] = new double[length];

                            for (int a = 0; a < CT_Averages; a++)
                            {
                                //int Choice = ;
                                for (int s = 0; s < length; s++)
                                {
                                    //if (s > IR.Length) 
                                    IR[c][s] += DCR[(a + 1) * length + s];
                                }
                            }

                            for (int n = 0; n < length; n++) IR[c][n] /= (double)CT_Averages;
                            double[] ETC = new double[IR[c].Length];
                            SNR[c] = new double[8];

                            for (int oct = 0; oct < 8; oct++)
                            {
                                double[] octIR = Pach_SP.FIR_Bandpass(IR[c], oct, SampleFreq, 0);
                                for (int i = 0; i < IR[c].Length; i++) ETC[i] = octIR[i] * octIR[i];
                                Direct_Sample[c] = FindDirect(ETC);
                                SNR[c][oct] = Signal_Noise_Ratio(ETC, Direct_Sample[c]);
                            }
                        }
                        StopCalled = false;
                        Running = false;
                    }

                    void WI_DataAvailable(object sender, NAudio.Wave.WaveInEventArgs e)
                    {
                        float sum = 0;
                        int g = (int)(Math.Pow(2, Gain));

                        for (int n = 0; n < e.BytesRecorded; n += block)
                        {
                            for (int c = 0; c < Channels_in; c++)
                            {
                                short sample = (short)(BitConverter.ToInt16(e.Buffer, n + 2 * c) * g);
                                Response[c].Add(sample);
                                sum += sample;
                            }
                        }
                    }

                    public static double Signal_Noise_Ratio(double[] etc, int d_sound)
                    {
                        int noise_blk = 128;
                        double sum = 0;
                        for (int i = 0; i < noise_blk; i++)
                        {
                            sum += etc[etc.Length - noise_blk + i];
                        }
                        double mean = sum / noise_blk;
                        return 10 * Math.Log10(etc[d_sound] / mean);
                    }

                    public static int FindDirect(double[] etc)
                    {
                        double deltaEMax = 0;
                        int D_Sound = 0;
                        for (int i = 1; i < etc.Length; i++)
                        {
                            double deltaE = (etc[i] - etc[i - 1]);
                            if (deltaE > deltaEMax)
                            {
                                deltaEMax = deltaE;
                                D_Sound = i;
                            }
                        }
                        return D_Sound;
                    }

                    public void Acquire(int input_device, Signal_Type ST, int output_device)
                    {
                        Running = true;
                        Channels_in = NAudio.Wave.WaveIn.GetCapabilities(input_device).Channels;
                        Response = new List<short>[Channels_in];
                        block = 2 * Channels_in;
                        WI = new NAudio.Wave.WaveInEvent();
                        WI.WaveFormat = new NAudio.Wave.WaveFormat(SampleFreq, 16, Channels_in);
                        WI.DeviceNumber = input_device;

                        WI.BufferMilliseconds = 100;
                        WI.NumberOfBuffers = 3;
                        WI.RecordingStopped += WI_RecordingStopped;
                        WI.DataAvailable += WI_DataAvailable;
                        WO.DeviceNumber = output_device;
                        for (int c = 0; c < Channels_in; c++) Response[c] = new List<short>();

                        SignalProvider Signal;

                        switch (ST)
                        {
                            case Signal_Type.Pink_Noise:
                                Signal = new Noise_Provider(1, (int)CT_Averages, SampleFreq);
                                break;
                            case Signal_Type.MLS:
                                throw new NotImplementedException();
                            case Signal_Type.Swept_Sine:
                                Signal = new Sweep_Provider((float)Signal_Length, CT_Averages, 63, 20000, SampleFreq);
                                break;
                            case Signal_Type.Logarithmic_Swept_Sine:
                                throw new NotImplementedException();
                            default:
                                System.Windows.Forms.MessageBox.Show("Select a Signal...");
                                return;
                        }

                        TD_Signal = Signal.Signal;
                        
                        WO.NumberOfBuffers = 1;
                        WO.DesiredLatency = 3000 * CT_Averages;
                        WO.Volume = 1.0f;
                        WO.Init(Signal);
                        WI.StartRecording();
                        WO.Play();
                        System.Threading.Thread.Sleep((int)(Signal_Time_s * (3 + CT_Averages) * 1000));
                        WO.Stop();
                        WI.StopRecording();
                        System.Threading.Thread.Sleep(100);
                        WI_RecordingStopped(this, null);
                    }

                    public int Sample_Frequency
                    {
                        get { return this.SampleFreq; }
                        set { SampleFreq = value; }
                    }

                    public double Signal_Time_s
                    {
                        get { return this.Signal_Length; }
                        set { Signal_Length = value; }
                    }

                    public int No_of_Channels
                    {
                        get { return Channels_in; }
                        set { Channels_in = value; }
                    }

                    public int No_of_Averages
                    {
                        get { return CT_Averages; }
                        set { CT_Averages = value; }
                    }

                    public double Gain
                    {
                        get { return Gain_out; }
                        set { Gain_out = value; }
                    }
                }

                public enum Signal_Type
                {
                    Pink_Noise,
                    MLS,
                    Swept_Sine,
                    Logarithmic_Swept_Sine
                }

                public static List<String> Get_Input_Devices()
                {
                    List<String> InputSelection = new List<string>();
                    int DC = NAudio.Wave.WaveIn.DeviceCount;
                    for (int i = 0; i < DC; i++)
                    {
                        NAudio.Wave.WaveInCapabilities DCap = NAudio.Wave.WaveIn.GetCapabilities(i);
                        InputSelection.Add(String.Format("{0}, {1} Channels", DCap.ProductName, DCap.Channels));
                    }
                    return InputSelection;
                }

                public static List<String> Get_Output_Devices()
                {
                    List<String> OutputSelection = new List<string>();
                    int DC = NAudio.Wave.WaveOut.DeviceCount;
                    for (int i = 0; i < DC; i++)
                    {
                        NAudio.Wave.WaveOutCapabilities DCap = NAudio.Wave.WaveOut.GetCapabilities(i);
                        OutputSelection.Add(String.Format("{0}, {1} Channels", DCap.ProductName, DCap.Channels));
                    }
                    return OutputSelection;
                }

                private class Sweep_Provider : SignalProvider
                {
                    protected int sample_rate;

                    public Sweep_Provider(float seconds, int no_of_avgs, float f_start, float f_finish, int _sample_rate)
                        : base(_sample_rate, 1)
                    {
                        sample_rate = _sample_rate;
                        float amplitude = 1.0f;
                        float initial_phase = (float)Math.PI;

                        int length = (int)Math.Round(seconds * sample_rate);
                        float[] frequency = new float[length];
                        signal = new float[length * (1 + no_of_avgs)];
                        signal1 = new float[length * no_of_avgs];

                        for (int t = 0; t < length; t++) frequency[t] = (float)(Math.Exp(Math.Log(f_start) * (1 - (t / (float)sample_rate) / seconds) + Math.Log(f_finish) * (t / (float)sample_rate) / seconds));

                        double[] phase = new double[length];
                        phase[0] = initial_phase;

                        for (int i = 1; i < length; i++)
                        {
                            phase[i] = (phase[i - 1] + 2.0f * Math.PI * frequency[i] / sample_rate) % (2 * Math.PI);
                            signal[i + length] = (float)(amplitude * Math.Sin(phase[i]));
                            signal1[i] = signal[i + length];
                        }

                        float taperCount = (float)sample_rate / 200;
                        for (int n = 1; n <= taperCount; n++) signal[2 * length - n] *= (n / taperCount);

                        for (int a = 1; a < no_of_avgs; a++) for (int i = 0; i < length; i++) signal[(a + 1) * length + i] = signal[length + i];
                    }

                    int sample;

                    public override int Read(float[] buffer, int offset, int sampleCount)
                    {
                        if (sample + offset >= signal.Length)
                            return 0;

                        for (int n = 0; n < sampleCount; n++)
                        {
                            int nsample = n + sample;
                            if (nsample >= signal.Length)
                            {
                                buffer[n + offset] = 0;
                                sample += n;
                                return n;
                            }
                            else
                            {
                                buffer[n + offset] = signal[nsample];
                            }
                        }
                        sample += sampleCount;
                        return sampleCount;
                    }
                }

                private class Noise_Provider : SignalProvider
                {
                    protected int sample_rate;

                    public Noise_Provider(float seconds, int no_of_avgs, int _sample_rate)
                        : base(_sample_rate, 1)
                    {

                        sample_rate = _sample_rate;
                        System.Random R = new Random();
                        int length = (int)Math.Round(seconds * sample_rate);
                        signal = new float[length * (1 + no_of_avgs)];
                        signal1 = new float[length * no_of_avgs];
                        for (int i = 0; i < length; i++)
                        {
                            signal[i + length] = (float)R.NextDouble();// *short.MaxValue;
                            signal1[i] = signal[i + length];
                        }
                        for (int a = 1; a < no_of_avgs; a++) for (int i = 0; i < length; i++) signal[(a + 1) * length + i] = signal[length + i];
                    }

                    int sample;

                    public override int Read(float[] buffer, int offset, int sampleCount)
                    {
                        if (sample + offset >= signal.Length)
                            return 0;

                        for (int n = 0; n < sampleCount; n++)
                        {
                            int nsample = n + sample;
                            if (nsample >= signal.Length)
                            {
                                buffer[n + offset] = 0;
                                sample += n;
                                return n;
                            }
                            else
                            {
                                buffer[n + offset] = signal[nsample];
                            }
                        }
                        sample += sampleCount;
                        return sampleCount;
                    }
                }

                private abstract class SignalProvider : NAudio.Wave.WaveProvider32
                {
                    protected internal float[] signal;
                    protected internal float[] signal1;

                    public SignalProvider(int samplerate, int channels)
                        : base(samplerate, channels) { }

                    public double[] Signal
                    {
                        get
                        {
                            double[] s_out = new double[signal1.Length];
                            for (int i = 0; i < signal1.Length; i++) s_out[i] = signal1[i];
                            return s_out;
                        }
                    }
                }
            }
        }
    }
}