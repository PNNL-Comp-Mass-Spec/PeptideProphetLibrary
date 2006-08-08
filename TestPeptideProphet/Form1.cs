using System;
using System.Drawing;
using System.Collections;
using System.ComponentModel;
using System.Windows.Forms;
using System.Data;
using System.Threading;

namespace TestPeptideProphet
{
	/// <summary>
	/// Summary description for Form1.
	/// </summary>
	public class Form1 : System.Windows.Forms.Form
	{
		private System.Windows.Forms.Button button1;
		private System.Windows.Forms.Button button2;
		/// <summary>
		/// Required designer variable.
		/// </summary>
		private System.ComponentModel.Container components = null;

		public Form1()
		{
			//
			// Required for Windows Form Designer support
			//
			InitializeComponent();

			//
			// TODO: Add any constructor code after InitializeComponent call
			//
		}

		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		protected override void Dispose( bool disposing )
		{
			if( disposing )
			{
				if (components != null) 
				{
					components.Dispose();
				}
			}
			base.Dispose( disposing );
		}

		#region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		private void InitializeComponent()
		{
			this.button1 = new System.Windows.Forms.Button();
			this.button2 = new System.Windows.Forms.Button();
			this.SuspendLayout();
			// 
			// button1
			// 
			this.button1.Location = new System.Drawing.Point(36, 24);
			this.button1.Name = "button1";
			this.button1.TabIndex = 0;
			this.button1.Text = "button1";
			this.button1.Click += new System.EventHandler(this.button1_Click);
			// 
			// button2
			// 
			this.button2.Location = new System.Drawing.Point(164, 28);
			this.button2.Name = "button2";
			this.button2.TabIndex = 1;
			this.button2.Text = "Abort";
			this.button2.Click += new System.EventHandler(this.button2_Click);
			// 
			// Form1
			// 
			this.AutoScaleBaseSize = new System.Drawing.Size(5, 13);
			this.ClientSize = new System.Drawing.Size(344, 106);
			this.Controls.Add(this.button2);
			this.Controls.Add(this.button1);
			this.Name = "Form1";
			this.Text = "Form1";
			this.ResumeLayout(false);

		}
		#endregion

		/// <summary>
		/// The main entry point for the application.
		/// </summary>
		[STAThread]
		static void Main() 
		{
			Application.Run(new Form1());
		}

		private PeptideProphetLibrary.PeptideProphet p;
		private Thread t;

		private void button1_Click(object sender, System.EventArgs e)
		{
			try
			{
				p = new PeptideProphetLibrary.PeptideProphet();
				PeptideProphetLibrary.InitializationParams parms = new PeptideProphetLibrary.InitializationParams();
				parms.InputFileName = @"C:\d3p019\DuXiuxia\AID_MPH_002_071906_23_24July06_Griffin_06-03-06_syn.txt";
				parms.OutputFilePath = @"C:\d3p019\DuXiuxia";
				parms.Enzyme = "tryptic";
				p.Setup(parms);
				t = new Thread(new ThreadStart(Start));				
				t.Start();
				
				//p.Start();
			}
			catch(Exception ex)
			{
				MessageBox.Show(ex.ToString());
			}
		}

		private void Start()
		{
			p.Start();
			if (p.FileStatus == PeptideProphetLibrary.IPeptideProphet.ProcessCheckFile.PP_IOFileNonexistent)
			{
				MessageBox.Show(p.FileStatus.ToString());
				MessageBox.Show(p.Results.ToString());
			}
			else
			{
				MessageBox.Show(p.Results.ToString());
			}
		}	

		private void button2_Click(object sender, System.EventArgs e)
		{
			try
			{
				p.Abort();
			}
			catch(Exception ex)
			{
				MessageBox.Show(ex.ToString());
			}
		}
	}
}
