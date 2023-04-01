//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL) by Arthur van der Harten 
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2023, Arthur van der Harten 
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

using System.Collections.Generic;

namespace Pachyderm_Acoustic
{
    /// <summary>
    /// A class to contain the acoustic materials library.
    /// </summary>
    public class Acoustics_Library
    {
        public List<Material> Abs_List = new List<Material>();
        public List<Material> TL_List = new List<Material>();

        public Acoustics_Library()
            : base()
        {
            //Read Absorption values in from library file designated on the Options Page...
            Load_Library();
        }

        /// <summary>
        /// Loads the user defined materials library.
        /// </summary>
        public void Load_Library()
        {
            Abs_List.Clear();
            TL_List.Clear();
            System.IO.StreamReader ML_Reader;
            try
            {
                string MLPath = Pach_Properties.Instance.Lib_Path();
                MLPath += "\\Pach_Materials_Library.txt";
                ML_Reader = new System.IO.StreamReader(MLPath);
                do
                {
                    try
                    {
                        string Material = ML_Reader.ReadLine();
                        string[] D_Mat = Material.Split(new char[] { ':' });
                        string Name = D_Mat[0].Trim();
                        string Abs_Code = D_Mat[1].Trim();
                        Abs_Code += "0000000000000000";
                        double[] Abs = new double[8];
                        double[] Sct = new double[8];
                        double[] Trns = new double[8];
                        Utilities.PachTools.DecodeAcoustics(Abs_Code, ref Abs, ref Sct, ref Trns);
                        this.Add_Unique_Abs(Name, Abs);
                    }
                    catch (System.Exception)
                    { continue; }
                } while (!ML_Reader.EndOfStream);
                ML_Reader.Close();
            }
            catch (System.Exception)
            {
            }

            System.IO.StreamReader IL_Reader;
            try
            {
                string ILPath = Pach_Properties.Instance.Lib_Path();
                ILPath += "\\Pach_Isolation_Library.txt";
                IL_Reader = new System.IO.StreamReader(ILPath);
                do
                {
                    try
                    {
                        string Material = IL_Reader.ReadLine();
                        string[] D_Mat = Material.Split(new char[] { ':' });
                        string Name = D_Mat[0].Trim();
                        string TL_Code = D_Mat[1].Trim();
                        //TL_Code += "0;0;0;0;0;0;0;0";
                        double[] TL = Utilities.PachTools.DecodeTransmissionLoss(TL_Code);
                        this.Add_Unique_TL(Name, TL);
                    }
                    catch (System.Exception)
                    { continue; }
                } while (!IL_Reader.EndOfStream);
                IL_Reader.Close();
            }
            catch (System.Exception)
            {
            }
        }

        public void Add_Unique_Abs(string Name, double[] Abs)
        {
            for (int i = 0; i < Abs_List.Count; i++)
            {
                if (!string.Equals(Abs_List[i].Name, Name, System.StringComparison.OrdinalIgnoreCase)) continue;
                Abs_List.RemoveAt(i);
                break;
            }

            Abs_List.Add(new Material(Name, Abs));
        }

        public void Delete_Abs(string Name)
        {
            for (int i = 0; i < Abs_List.Count; i++)
            {
                if (!string.Equals(Abs_List[i].Name, Name, System.StringComparison.OrdinalIgnoreCase)) continue;
                Abs_List.RemoveAt(i);
                break;
            }
        }

        public void Delete_Abs(int index)
        {
            Abs_List.RemoveAt(index);
        }

        public Material Abs_byKey(string Selection)
        {
            foreach (Material Mat in Abs_List) if (Mat.Name == Selection) return Mat;
            throw new System.Exception();
        }

        /// <summary>
        /// Saves the user defined materials library.
        /// </summary>
        public void Save_Abs_Library()
        {
            //Enter an external file saver here... 
            string MLPath = Pach_Properties.Instance.Lib_Path();
            MLPath += "\\Pach_Materials_Library.txt";

            System.IO.StreamWriter Writer;
            Writer = new System.IO.StreamWriter(MLPath);

            for (int i = 0; i < Abs_List.Count; i++)
            {
                string Entry = Abs_List[i].Name + ':';
                int[] sct = new int[8];
                int[] abs = new int[8];
                int[] trns = new int[1];
                for (int oct = 0; oct < 8; oct++)
                {
                    abs[oct] = (int)(Abs_List[i].Values[oct] * 100.0);
                }
                string Abs_Code = Utilities.PachTools.EncodeAcoustics(abs, sct, trns);
                Entry += Abs_Code.Substring(0, 16);
                Writer.WriteLine(Entry);
            }
            Writer.Close();
        }

        /// <summary>
        /// Returns the names of all user defined materials for display in an interface component.
        /// </summary>
        /// <returns></returns>
        public List<string> Names_Abs()
        {
            List<string> Catalog = new List<string>();
            foreach (Material obj in Abs_List)
            {
                Catalog.Add(obj.Name);
            }
            return Catalog;
        }

        public void Add_Unique_TL(string Name, double[] TL)
        {
            for (int i = 0; i < TL_List.Count; i++)
            {
                if (!string.Equals(TL_List[i].Name, Name, System.StringComparison.OrdinalIgnoreCase)) continue;
                TL_List.RemoveAt(i);
                break;
            }

            TL_List.Add(new Material(Name, TL));
        }

        public void Delete_TL(string Name)
        {
            for (int i = 0; i < TL_List.Count; i++)
            {
                if (!string.Equals(TL_List[i].Name, Name, System.StringComparison.OrdinalIgnoreCase)) continue;
                TL_List.RemoveAt(i);
                break;
            }
        }

        public void Delete_TL(int index)
        {
            TL_List.RemoveAt(index);
        }

        public Material TL_byKey(string Selection)
        {
            foreach (Material Mat in TL_List) if (Mat.Name == Selection) return Mat;
            throw new System.Exception();
        }

        /// <summary>
        /// Saves the user defined materials library.
        /// </summary>
        public void Save_TL_Library()
        {
            //Enter an external file saver here... 
            string MLPath = Pach_Properties.Instance.Lib_Path();
            MLPath += "\\Pach_Isolation_Library.txt";

            System.IO.StreamWriter Writer;
            Writer = new System.IO.StreamWriter(MLPath);

            for (int i = 0; i < TL_List.Count; i++)
            {
                string Entry = TL_List[i].Name + ':';
                double[] TL = new double[8];
                for (int oct = 0; oct < 8; oct++)
                {
                    TL[oct] = (double)(TL_List[i].Values[oct]);
                }
                string TL_Code = Utilities.PachTools.EncodeTransmissionLoss(TL);
                Entry += TL_Code;
                Writer.WriteLine(Entry);
            }
            Writer.Close();
        }

        /// <summary>
        /// Returns the names of all user defined materials for display in an interface component.
        /// </summary>
        /// <returns></returns>
        public List<string> Names_TL()
        {
            List<string> Catalog = new List<string>();
            foreach (Material obj in TL_List)
            {
                Catalog.Add(obj.Name);
            }
            return Catalog;
        }
    }
    /// <summary>
    /// a structure defining a material.
    /// </summary>
    public struct Material
    {
        public double[] Values;
        public string Name;

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="nym">the name of the material</param>
        /// <param name="Abs">the absorption coefficients, from 62.5 Hz.[0] to 8 khz.[7]</param>
        public Material(string nym, double[] vals)
        {
            if (vals.Length != 8) throw new System.Exception("Material properties must be specified for all octave bands form 63 to 8k.");
            Name = nym;
            Values = vals;
        }
    }
}