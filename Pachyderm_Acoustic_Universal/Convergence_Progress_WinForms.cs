using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using ZedGraph;

namespace Pachyderm_Acoustic
{
    public partial class Convergence_Progress_WinForms : Form
    {
        private static Convergence_Progress_WinForms instance = null;
        private Queue<double> history1 = new Queue<double>();
        private Queue<double> history2 = new Queue<double>();
        private List<double> CT = new List<double>();
        bool conclude = false;

        public Convergence_Progress_WinForms()
        {
            InitializeComponent();

            IR_View.GraphPane.Title.Text = "Impulse Response Status";
            IR_View.GraphPane.XAxis.Title.Text = "Time (s)";
            IR_View.GraphPane.YAxis.Title.Text = "Ratio of Change";
            Conv_View.GraphPane.Title.Text = "Convergence History (Last 20 records)";
            Conv_View.GraphPane.XAxis.Title.Text = "Iteration";
            Conv_View.GraphPane.YAxis.Title.Text = "Maximum Change";
        }
        
        public static Convergence_Progress_WinForms Instance
        {
            get
            {
                if (instance == null)
                {
                    instance = new Convergence_Progress_WinForms();
                }
                return instance;
            }
        }

        public bool Populate(double Conv1, double Conv2, double ConvInf, int ID, double corr = 0)
        {
            if (this.Visible == false) return false;

            if (corr != 0) IR_View.GraphPane.Title.Text = "Impulse Response Status - Schroeder correlation = " + Math.Round(corr, 3);

            double[] t = new double[6] {0, 0.05, 0.05, 0.08, 0.08, 3};
            double[] conv = new double[6] { Conv1, Conv1, Conv2, Conv2, ConvInf, ConvInf };

            double convmax = conv.Max();
            if (ID == 0)
            {
                history1.Enqueue(convmax);
                if (history1.Count > 20) history1.Dequeue();
                IR_View.GraphPane.CurveList.Clear();
                Conv_View.GraphPane.CurveList.Clear();
                IR_View.GraphPane.AddCurve("Data", t, conv, System.Drawing.Color.Red, SymbolType.None);
                IR_View.GraphPane.XAxis.Scale.Max = t.Last();
                IR_View.GraphPane.YAxis.Scale.Max = Math.Max(convmax, 0.3);
                if (history1.Count < 20) CT.Add(CT.Count -10);
                Conv_View.GraphPane.AddCurve("Clear", CT.ToArray(), history1.ToArray(), System.Drawing.Color.Red, SymbolType.None);
                Conv_View.GraphPane.XAxis.Scale.Min = -10;
                Conv_View.GraphPane.XAxis.Scale.Max = CT.Count() - 10;
                Conv_View.GraphPane.YAxis.Scale.Max = history1.Max();
            }
            if (ID == 1)
            {
                history2.Enqueue(convmax);
                if (history2.Count > 20) history2.Dequeue();
                IR_View.GraphPane.AddCurve("Data", t, conv, System.Drawing.Color.Blue, SymbolType.None);
                IR_View.GraphPane.XAxis.Scale.Max = Math.Max(IR_View.GraphPane.XAxis.Scale.Max, 1);
                IR_View.GraphPane.YAxis.Scale.Max = Math.Max(IR_View.GraphPane.YAxis.Scale.Max, Math.Max(convmax, 0.3));
                Conv_View.GraphPane.YAxis.Scale.Max = Math.Max(Conv_View.GraphPane.YAxis.Scale.Max, history2.Max());
                Conv_View.GraphPane.AddCurve("Clear", CT.ToArray(), history2.ToArray(), System.Drawing.Color.Blue, SymbolType.None);
            }
            IR_View.GraphPane.XAxis.Title = new AxisLabel("Time (ms)", "Arial", 10, System.Drawing.Color.Black, false, false, false);
            IR_View.GraphPane.YAxis.Title = new AxisLabel("Magnitude", "Arial", 10, System.Drawing.Color.Black, false, false, false);
            
            if (ID == 0)
            {
                if (conv.Max() < 100) IR_View.GraphPane.YAxis.Type = AxisType.Linear;
                else IR_View.GraphPane.YAxis.Type = AxisType.Log;
                if (history1.Max() < 100) Conv_View.GraphPane.YAxis.Type = AxisType.Linear;
                else Conv_View.GraphPane.YAxis.Type = AxisType.Log;
            }
            else
            {
                if (conv.Max() > 100) IR_View.GraphPane.YAxis.Type = AxisType.Log;
                if (history2.Max() > 100) Conv_View.GraphPane.YAxis.Type = AxisType.Log;
            }

            Color Conv_Color = Color.Gray;
            switch (conv_count)
            {
                case 1:
                    Conv_Color = Color.Red;
                    break;
                case 2:
                    Conv_Color = Color.OrangeRed;
                    break;
                case 3:
                    Conv_Color = Color.Orange;
                    break;
                case 4:
                    Conv_Color = Color.Yellow;
                    break;
                case 5:
                    Conv_Color = Color.YellowGreen;
                    break;
                case 6:
                    Conv_Color = Color.Green;
                    break;
                case 7:
                    Conv_Color = Color.ForestGreen;
                    break;
                case 8:
                    Conv_Color = Color.DarkSeaGreen;
                    break;
                case 9:
                    Conv_Color = Color.Blue;
                    break;
                case 10:
                    Conv_Color = Color.DarkBlue;
                    break;
                default:
                    break;
            }

            Conv_View.GraphPane.AddCurve("Convergence", new double[2] { -5, 25 }, new double[2] { 0.2, 0.2 }, Conv_Color);
            Conv_View.GraphPane.AddCurve("Convergence", new double[2] { -5, 25 }, new double[2] { 0.02, 0.02 }, Conv_Color);
            IR_View.GraphPane.AddCurve("Convergence", t, new double[6] { 0.02, 0.02, 0.02, 0.02, 0.1, 0.1 }, Conv_Color);

            if (conclude)
            {
                conclude = false;
                return true;
            }
            else return false;
        }

        int conv_count = 0;

        public bool Populate(double[] IR_p, double Convergence, double sample_freq, int count, int ID)
        {
            if (this.Visible == false) return false;

            double[] t = new double[IR_p.Length];
            for (int i = 0; i < IR_p.Length; i++) t[i] = (i / sample_freq);
            if (ID == 0)
            {
                conv_count = count;
                history1.Enqueue(Convergence);
                if (history1.Count > 20) history1.Dequeue();
                IR_View.GraphPane.CurveList.Clear();
                Conv_View.GraphPane.CurveList.Clear();
                IR_View.GraphPane.AddCurve("Data", t, IR_p, System.Drawing.Color.Red, SymbolType.None);
                IR_View.GraphPane.XAxis.Scale.Max = t.Last();
                IR_View.GraphPane.YAxis.Scale.Max = Math.Max(IR_p.Max(), 0.3);
                if (history1.Count < 20)CT.Add(CT.Count - 10);
                Conv_View.GraphPane.AddCurve("Clear", CT.ToArray(), history1.ToArray(), System.Drawing.Color.Red, SymbolType.None);
                IR_View.ZoomOutAll(IR_View.GraphPane);
                Conv_View.GraphPane.XAxis.Scale.Min = -10;
                Conv_View.GraphPane.XAxis.Scale.Max = CT.Count() - 10;
                Conv_View.GraphPane.YAxis.Scale.Max = history1.Max();
            }
            if (ID == 1)
            {
                conv_count = Math.Min(conv_count, count);
                history2.Enqueue(Convergence);
                if (history2.Count > 20) history2.Dequeue();
                IR_View.GraphPane.AddCurve("Data", t, IR_p, System.Drawing.Color.Blue, SymbolType.None);
                IR_View.GraphPane.XAxis.Scale.Max = Math.Max(IR_View.GraphPane.XAxis.Scale.Max, t.Last());
                IR_View.GraphPane.YAxis.Scale.Max = Math.Max(IR_View.GraphPane.YAxis.Scale.Max, Math.Max(IR_p.Max(), 0.3));
                Conv_View.GraphPane.YAxis.Scale.Max = Math.Max(Conv_View.GraphPane.YAxis.Scale.Max, history2.Max());
                Conv_View.GraphPane.AddCurve("Clear", CT.ToArray(), history2.ToArray(), System.Drawing.Color.Blue, SymbolType.None);
            }
            IR_View.GraphPane.XAxis.Title = new AxisLabel("Time (ms)", "Arial", 10, System.Drawing.Color.Black, false, false, false);
            IR_View.GraphPane.YAxis.Title = new AxisLabel("Magnitude", "Arial", 10, System.Drawing.Color.Black, false, false, false);

            if (ID == 0)
            {
                if (Convergence < 100) IR_View.GraphPane.YAxis.Type = AxisType.Linear;
                else IR_View.GraphPane.YAxis.Type = AxisType.Log;
                if (history1.Max() < 100) Conv_View.GraphPane.YAxis.Type = AxisType.Linear;
                else Conv_View.GraphPane.YAxis.Type = AxisType.Log;
            }
            else
            {
                if (Convergence > 100) IR_View.GraphPane.YAxis.Type = AxisType.Log;
                if (history2.Max() > 100) Conv_View.GraphPane.YAxis.Type = AxisType.Log;
            }

            Color Conv_Color = Color.Gray;
            switch (conv_count)
            {
                case 1:
                    Conv_Color = Color.Red;
                    break;
                case 2:
                    Conv_Color = Color.OrangeRed;
                    break;
                case 3:
                    Conv_Color = Color.Orange;
                    break;
                case 4:
                    Conv_Color = Color.Yellow;
                    break;
                case 5:
                    Conv_Color = Color.YellowGreen;
                    break;
                case 6:
                    Conv_Color = Color.Green;
                    break;
                case 7:
                    Conv_Color = Color.ForestGreen;
                    break;
                case 8:
                    Conv_Color = Color.DarkSeaGreen;
                    break;
                case 9:
                    Conv_Color = Color.Blue;
                    break;
                case 10:
                    Conv_Color = Color.DarkBlue;
                    break;
                default:
                    break;
            }

            IR_View.GraphPane.AddCurve("Convergence", new double[2] { -10000, 10000 }, new double[2] { 0.1, 0.1 }, Conv_Color);
            Conv_View.GraphPane.AddCurve("Convergence", new double[2] { -5, 25 }, new double[2] { 0.1, 0.1 }, Conv_Color);

            if (conclude)
            {
                conclude = false;
                return true;
            }
            else return false;
        }

        private void Conclude_Click(object sender, EventArgs e)
        {
            conclude = true;
        }
    }
}
