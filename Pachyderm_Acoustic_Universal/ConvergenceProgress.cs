//using System;
//using System.Collections.Generic;
//using Eto.Forms;
//using System.Linq;
//using OxyPlot;
//using Eto.OxyPlot;

//namespace Pachyderm_Acoustic
//{
//	public partial class ConvergenceProgress : Eto.Forms.Form
//    { 
//        PlotModel PM;
//        PlotModel History;
//        List<OxyPlot.Series.LineSeries> Hist = new List<OxyPlot.Series.LineSeries>();
//        List<OxyPlot.Series.LineSeries> ConvergenceLines = new List<OxyPlot.Series.LineSeries>();
//        bool Context = false;

//        private static ConvergenceProgress instance = null;

//        private ConvergenceProgress()
//        {
//            PM = new PlotModel() { Title = "Echogram" };
//            PM.Axes.Add(new OxyPlot.Axes.LinearAxis() { Position = OxyPlot.Axes.AxisPosition.Left, Minimum = 0, Maximum = double.NaN, MajorStep = 10, MinorStep = 5, Title = "SPL (dB)" });
//            PM.Axes.Add(new OxyPlot.Axes.LinearAxis() { Position = OxyPlot.Axes.AxisPosition.Bottom, Minimum = 0, Maximum = 1, MajorStep = .1, MinorStep = .05, Title = "Time (s)" });
//            History = new PlotModel() { Title = "Convergence Plot" };
//            PM.Axes.Add(new OxyPlot.Axes.LinearAxis() { Position = OxyPlot.Axes.AxisPosition.Left, MajorStep = 10, MinorStep = 5, Title = "Convergence (%)" });
//            PM.Axes.Add(new OxyPlot.Axes.LinearAxis() { Position = OxyPlot.Axes.AxisPosition.Bottom, Minimum = 0, Maximum = 1, Title = "Iterations"});
//            Hist.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Red, StrokeThickness = 0.1, Title = "Convergence with Farthest Visible Receiver" });
//            Hist.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Red, StrokeThickness = 0.1, Title = "Convergence with Farthest Occluded Receiver" });
//            ConvergenceLines.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Gray, StrokeThickness = 0.1, Title = "Early Convergence" , LineStyle = LineStyle.Dash});
//            ConvergenceLines.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Gray, StrokeThickness = 0.1, Title = "Early Convergence", LineStyle = LineStyle.Dash });
//            ConvergenceLines.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Gray, StrokeThickness = 0.1, Title = "Convergence", LineStyle = LineStyle.Dash });
//            ConvergenceLines[0].Points.Add(new OxyPlot.DataPoint(0, 0.02));
//            ConvergenceLines[0].Points.Add(new OxyPlot.DataPoint(10000, 0.02));
//            ConvergenceLines[0].Points.Add(new OxyPlot.DataPoint(0, 0.1));
//            ConvergenceLines[0].Points.Add(new OxyPlot.DataPoint(10000, 0.1));
//            ConvergenceLines[0].Points.Add(new OxyPlot.DataPoint(0, 0.1));
//            ConvergenceLines[0].Points.Add(new OxyPlot.DataPoint(10000, 0.1));
//            History.Series.Add(Hist[0]);
//            History.Series.Add(Hist[1]);
//            History.Series.Add(ConvergenceLines[0]);
//            PM.Series.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Blue, StrokeThickness = 0.1, Title = "Impulse Response" });
//            InitializeComponent();
//        }

//        public static ConvergenceProgress Instance
//        {
//            get
//            {
//                if (instance == null)
//                {
//                    instance = new ConvergenceProgress();
//                }
//                return instance;
//            }
//        }

//        public void reset_IR()
//        {
//            //this.ShowModalAsync();
//            //PM.Series.Clear();
//            Context = false;
//        }

//        public void reset()
//        {
//            //Context = false;
//            //PM.Series.Clear();
//            //History.Series.Clear();
//            //Hist.Clear();
//            //Hist.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Red, StrokeThickness = 0.1, Title = "Convergence with Farthest Visible Receiver" });
//            //Hist.Add(new OxyPlot.Series.LineSeries() { Color = OxyColors.Red, StrokeThickness = 0.1, Title = "Convergence with Farthest Occluded Receiver" });
//        }

//        public void Populate(double[] IR_p, double Convergence, double sample_freq, int ID)
//        {
//            //if (!Context)
//            //{
//            //    Context = true;
//            //    History.Series.Add(ConvergenceLines[2]);
//            //}
//            PM.Axes[1].Maximum = IR_p.Last();
//            //OxyPlot.Series.LineSeries pts = new OxyPlot.Series.LineSeries() { Color = OxyColors.Blue, StrokeThickness = 0.1, Title = "Impulse Response" };
//            //Hist[ID].Points.Add(new OxyPlot.DataPoint(Hist[ID].Points.Count(), Convergence));
//            (PM.Series[0] as OxyPlot.Series.LineSeries).Points.Clear();
//            for (int i = 0; i < IR_p.Length; i++) (PM.Series[0] as OxyPlot.Series.LineSeries).Points.Add(new OxyPlot.DataPoint(i / sample_freq, IR_p[i]));
//            //PM.Series.Add(pts);
//            (History.Series[ID] as OxyPlot.Series.LineSeries).Points.Add(new OxyPlot.DataPoint(Hist[ID].Points.Count(), Convergence));
//        }
//    }
//}
