//using Eto.Drawing;
//using Eto.Forms;
//using Eto;
//using Eto.Wpf;
//using Eto.Wpf.Forms;
//using OxyPlot;
//using OxyPlot.Wpf;

//namespace Pachyderm_Acoustic
//{
//	partial class ConvergenceProgress : Form
//	{
//        void InitializeComponent()
//        {
//            Title = "Convergence Progress";
//            ClientSize = new Eto.Drawing.Size(400, 350);

//            Maximizable = false;
//            Minimizable = false;
//            Padding = 5;
//            Resizable = false;
//            ShowInTaskbar = false;
//            Title = GetType().Name;
//            WindowStyle = WindowStyle.Default;

//            var pf = Eto.Platform.Detect;

//            if (pf.IsWpf)
//            {
//                this.Platform.Add(typeof(Eto.OxyPlot.Plot.IHandler), () => new Eto.OxyPlot.Wpf.PlotHandler());
//                this.Platform.Add(typeof(Eto.Forms.StackLayout.IHandler), () => new Eto.Forms.StackLayout());
//            }
//            else if (pf.IsMac)
//            {
//                //this.Platform.Add(typeof(Plot.IHandler), () => new Eto.OxyPlot.Mac.PlotHandler());
//            }

//            PlotView P1 = new PlotView();
//            Eto.OxyPlot.Plot P2 = new Eto.OxyPlot.Plot();

//            P1.Model = PM;
//            P2.Model = History;

//            Content = new StackLayout
//            {
//                Padding = 5,
//                //Spacing = new Size(5, 5),
//                Items = { P1, P2 }
//            };
//        }
//	}
//}