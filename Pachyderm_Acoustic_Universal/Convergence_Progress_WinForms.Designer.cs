namespace Pachyderm_Acoustic
{
    partial class Convergence_Progress_WinForms
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.tableLayoutPanel1 = new System.Windows.Forms.TableLayoutPanel();
            this.Conv_View = new ZedGraph.ZedGraphControl();
            this.IR_View = new ZedGraph.ZedGraphControl();
            this.Conclude = new System.Windows.Forms.Button();
            this.tableLayoutPanel1.SuspendLayout();
            this.SuspendLayout();
            // 
            // tableLayoutPanel1
            // 
            this.tableLayoutPanel1.ColumnCount = 1;
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Absolute, 20F));
            this.tableLayoutPanel1.Controls.Add(this.Conv_View, 0, 1);
            this.tableLayoutPanel1.Controls.Add(this.IR_View, 0, 0);
            this.tableLayoutPanel1.Controls.Add(this.Conclude, 0, 2);
            this.tableLayoutPanel1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tableLayoutPanel1.Location = new System.Drawing.Point(0, 0);
            this.tableLayoutPanel1.Name = "tableLayoutPanel1";
            this.tableLayoutPanel1.RowCount = 3;
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 50F));
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 50F));
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Absolute, 43F));
            this.tableLayoutPanel1.Size = new System.Drawing.Size(546, 724);
            this.tableLayoutPanel1.TabIndex = 0;
            // 
            // Conv_View
            // 
            this.Conv_View.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.Conv_View.AutoSize = true;
            this.Conv_View.AutoSizeMode = System.Windows.Forms.AutoSizeMode.GrowAndShrink;
            this.Conv_View.EditButtons = System.Windows.Forms.MouseButtons.Left;
            this.Conv_View.Location = new System.Drawing.Point(8, 347);
            this.Conv_View.Margin = new System.Windows.Forms.Padding(8, 7, 8, 7);
            this.Conv_View.Name = "Conv_View";
            this.Conv_View.PanModifierKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Shift | System.Windows.Forms.Keys.None)));
            this.Conv_View.ScrollGrace = 0D;
            this.Conv_View.ScrollMaxX = 0D;
            this.Conv_View.ScrollMaxY = 120D;
            this.Conv_View.ScrollMaxY2 = 0D;
            this.Conv_View.ScrollMinX = 0D;
            this.Conv_View.ScrollMinY = 0D;
            this.Conv_View.ScrollMinY2 = 0D;
            this.Conv_View.Size = new System.Drawing.Size(530, 326);
            this.Conv_View.TabIndex = 44;
            this.Conv_View.UseExtendedPrintDialog = true;
            this.Conv_View.UseWaitCursor = true;
            // 
            // IR_View
            // 
            this.IR_View.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.IR_View.AutoSize = true;
            this.IR_View.AutoSizeMode = System.Windows.Forms.AutoSizeMode.GrowAndShrink;
            this.IR_View.EditButtons = System.Windows.Forms.MouseButtons.Left;
            this.IR_View.Location = new System.Drawing.Point(8, 7);
            this.IR_View.Margin = new System.Windows.Forms.Padding(8, 7, 8, 7);
            this.IR_View.Name = "IR_View";
            this.IR_View.PanModifierKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Shift | System.Windows.Forms.Keys.None)));
            this.IR_View.ScrollGrace = 0D;
            this.IR_View.ScrollMaxX = 0D;
            this.IR_View.ScrollMaxY = 120D;
            this.IR_View.ScrollMaxY2 = 0D;
            this.IR_View.ScrollMinX = 0D;
            this.IR_View.ScrollMinY = 0D;
            this.IR_View.ScrollMinY2 = 0D;
            this.IR_View.Size = new System.Drawing.Size(530, 326);
            this.IR_View.TabIndex = 43;
            this.IR_View.UseExtendedPrintDialog = true;
            this.IR_View.UseWaitCursor = true;
            // 
            // Conclude
            // 
            this.Conclude.Dock = System.Windows.Forms.DockStyle.Fill;
            this.Conclude.Location = new System.Drawing.Point(3, 683);
            this.Conclude.Name = "Conclude";
            this.Conclude.Size = new System.Drawing.Size(540, 38);
            this.Conclude.TabIndex = 45;
            this.Conclude.Text = "Conclude_Simulation";
            this.Conclude.UseVisualStyleBackColor = true;
            this.Conclude.Click += new System.EventHandler(this.Conclude_Click);
            // 
            // Convergence_Progress_WinForms
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(546, 724);
            this.Controls.Add(this.tableLayoutPanel1);
            this.Name = "Convergence_Progress_WinForms";
            this.Text = "Convergence Progress";
            this.tableLayoutPanel1.ResumeLayout(false);
            this.tableLayoutPanel1.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.TableLayoutPanel tableLayoutPanel1;
        private ZedGraph.ZedGraphControl Conv_View;
        private ZedGraph.ZedGraphControl IR_View;
        private System.Windows.Forms.Button Conclude;
    }
}