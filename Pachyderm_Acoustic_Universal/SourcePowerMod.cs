using System;
using System.Windows.Forms;

namespace Pachyderm_Acoustic
{
    public partial class SourcePowerMod : Form
    {
        public double[] Power;
        public bool accept = false;

        public SourcePowerMod(double[] initpower)
        {
            //KeyPreview = true;
            //KeyDown += keypressed;
            AcceptButton = OK;
            CancelButton = Cancel;
            InitializeComponent();
            if (initpower == null || initpower.Length != 8) initpower = new double[8] { 120, 120, 120, 120, 120, 120, 120, 120 };
            SWL0.Value = (decimal)initpower[0];
            SWL1.Value = (decimal)initpower[1];
            SWL2.Value = (decimal)initpower[2];
            SWL3.Value = (decimal)initpower[3];
            SWL4.Value = (decimal)initpower[4];
            SWL5.Value = (decimal)initpower[5];
            SWL6.Value = (decimal)initpower[6];
            SWL7.Value = (decimal)initpower[7];
        }

        //public void keypressed(object sender, KeyEventArgs e)
        //{
        //    switch (e.KeyCode)
        //    {
        //        case Keys.Enter:
        //            e.Handled = true;
        //            OK.PerformClick();
        //            break;
        //        case Keys.NumPad8:
        //            e.Handled = true;
        //            Cancel.PerformClick();
        //            break;
        //    }
        //}

    protected override void OnFormClosed(FormClosedEventArgs e)
        {
            Power = new double[8];
            Power[0] = (double)SWL0.Value;
            Power[1] = (double)SWL1.Value;
            Power[2] = (double)SWL2.Value;
            Power[3] = (double)SWL3.Value;
            Power[4] = (double)SWL4.Value;
            Power[5] = (double)SWL5.Value;
            Power[6] = (double)SWL6.Value;
            Power[7] = (double)SWL7.Value;
            base.OnFormClosed(e);
        }

        private void swl63_Scroll(object sender, EventArgs e)
        {
            SWL0.Value = swl63.Value;
        }
        private void swl125_Scroll(object sender, EventArgs e)
        {
            SWL1.Value = swl125.Value;
        }
        private void swl250_Scroll(object sender, EventArgs e)
        {
            SWL2.Value = swl250.Value;
        }
        private void swl500_Scroll(object sender, EventArgs e)
        {
            SWL3.Value = swl500.Value;
        }
        private void swl1k_Scroll(object sender, EventArgs e)
        {
            SWL4.Value = swl1k.Value;
        }
        private void swl2k_Scroll(object sender, EventArgs e)
        {
            SWL5.Value = swl2k.Value;
        }
        private void swl4k_Scroll(object sender, EventArgs e)
        {
            SWL6.Value = swl4k.Value;
        }
        private void swl8k_Scroll(object sender, EventArgs e)
        {
            SWL7.Value = swl8k.Value;
        }

        private void swl0_updown(object sender, EventArgs e)
        {
            swl63.Value = (int)SWL0.Value;
            SPL0.Text = (Math.Round(SWL0.Value, 2) - 11).ToString();
        }
        private void swl1_updown(object sender, EventArgs e)
        {
            swl125.Value = (int)SWL1.Value;
            SPL1.Text = (Math.Round(SWL1.Value, 2) - 11).ToString();
        }
        private void swl2updown(object sender, EventArgs e)
        {
            swl250.Value = (int)SWL2.Value;
            SPL2.Text = (Math.Round(SWL2.Value, 2) - 11).ToString();
        }
        private void swl3updown(object sender, EventArgs e)
        {
            swl500.Value = (int)SWL3.Value;
            SPL3.Text = (Math.Round(SWL3.Value, 2) - 11).ToString();
        }
        private void swl4updown(object sender, EventArgs e)
        {
            swl1k.Value = (int)SWL4.Value;
            SPL4.Text = (Math.Round(SWL4.Value, 2) - 11).ToString();
        }
        private void swl5updown(object sender, EventArgs e)
        {
            swl2k.Value = (int)SWL5.Value;
            SPL5.Text = (Math.Round(SWL5.Value, 2) - 11).ToString();
        }
        private void swl6updown(object sender, EventArgs e)
        {
            swl4k.Value = (int)SWL6.Value;
            SPL6.Text = (Math.Round(SWL6.Value, 2) - 11).ToString();
        }
        private void swl7updown(object sender, EventArgs e)
        {
            swl8k.Value = (int)SWL7.Value;
            SPL7.Text = (Math.Round(SWL7.Value, 2) - 11).ToString();
        }

        private void OK_Click(object sender, EventArgs e)
        {
            accept = true;
            Close();
        }

        private void Cancel_Click(object sender, EventArgs e)
        {
            accept = false;
            Close();
        }
    }
}
