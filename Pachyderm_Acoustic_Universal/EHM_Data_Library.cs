//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2023, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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

using Pachyderm_Acoustic.AbsorptionModels;
using System.Collections.Generic;

namespace Pachyderm_Acoustic
{
    /// <summary>
    /// A class to contain the acoustic materials library.
    /// </summary>
    public class EHM_Data_Library
    {
        public readonly List<string> Names = new List<string>();
        public readonly List<double> Thickness = new List<double>();
        public readonly List<double> Density_kgm = new List<double>();
        public readonly List<string> Category = new List<string>();
        public readonly List<string> Application = new List<string>();
        public readonly List<double> Flow_Resistivity = new List<double>();
        public readonly List<double> FR_Dev = new List<double>();
        public readonly List<double[]> Abs_Coef = new List<double[]>(); //Currently 13 coef - 3rd oct 125 to 2k
        public readonly List<int> Abs_Method = new List<int>(); //0 for ISO 10534 (tube), 1 for ISO 354 (chamber)
        public readonly List<double> ECC_Spec = new List<double>(); //kgCO2/m2
        public readonly List<double> ECC_Abs = new List<double>(); //kgCO2/m3
        public readonly List<string> Source = new List<string>(); //Where did acoustic data come from?
        public readonly List<ABS_Layer[]> Substrates = new List<ABS_Layer[]>(); //List of substrate wall assemblies 
        public readonly List<double> substrate_ECC = new List<double>(); //ECC of substrate wall assemblies

        public EHM_Data_Library()
            : base()
        {
            //Read Absorption values in from library file designated on the Options Page...
            System.IO.StreamReader ML_Reader;
            try
            {
                string MLPath = Pach_Properties.Instance.Lib_Path;
                MLPath += "\\Acoustic_Materials_and_ECCs.csv";
                ML_Reader = new System.IO.StreamReader(MLPath);
                ML_Reader.ReadLine();
                do
                {
                    try
                    {
                        string Material = ML_Reader.ReadLine();
                        string[] D_Mat = Material.Split(new char[] { ',' });
                        string Name = D_Mat[0].Trim();
                        Names.Add(Name);
                        Category.Add(D_Mat[1].Trim());
                        Application.Add(D_Mat[2].Trim());
                        Density_kgm.Add(double.Parse(D_Mat[3].Trim()));
                        double T = double.Parse(D_Mat[4].Trim());
                        Thickness.Add(T);
                        Flow_Resistivity.Add(double.Parse(D_Mat[5].Trim()));
                        FR_Dev.Add(double.Parse(D_Mat[6].Trim('(', ')', '±')));
                        Abs_Coef.Add(new double[13] { double.Parse(D_Mat[7].Trim()), double.Parse(D_Mat[8].Trim()), double.Parse(D_Mat[9].Trim()), double.Parse(D_Mat[10].Trim()), double.Parse(D_Mat[11].Trim()), double.Parse(D_Mat[12].Trim()), double.Parse(D_Mat[13].Trim()), double.Parse(D_Mat[14].Trim()), double.Parse(D_Mat[15].Trim()), double.Parse(D_Mat[16].Trim()), double.Parse(D_Mat[17].Trim()), double.Parse(D_Mat[18].Trim()), double.Parse(D_Mat[19].Trim()) });
                        Abs_Method.Add(D_Mat[20].Contains("534") ? 0 : D_Mat[29].Contains("354") ? 1 : -1);
                        double Simp_ECC_GWP = double.Parse(D_Mat[23]);
                        ECC_Spec.Add(Simp_ECC_GWP);
                        ECC_Abs.Add(Simp_ECC_GWP / T);
                        Source.Add(D_Mat[22]);
                    }
                    catch (System.Exception)
                    { continue; }
                } while (!ML_Reader.EndOfStream);
                ML_Reader.Close();
            }
            catch (System.Exception)
            {
            }

            try
            {
                string MLPath = Pach_Properties.Instance.Lib_Path;
                MLPath += "\\Wall_Assemblies_EC.csv";
                ML_Reader = new System.IO.StreamReader(MLPath);
                ML_Reader.ReadLine();
                List<string> lines = new List<string>();
                do
                {
                    lines.Add(ML_Reader.ReadLine());
                } while (!ML_Reader.EndOfStream);
                ML_Reader.Close();

                int i = 0;
                int id = 1;
                List<ABS_Layer> layers = new List<ABS_Layer>();
                do
                {
                    do
                    {
                        string line = lines[i];
                        try
                        {
                            string[] D_Mat = line.Split(new char[] { ',' });
                            if (double.Parse(D_Mat[0]) == id)
                            {
                                string Name = D_Mat[1].Trim();
                                //Assess type of layer
                                switch (D_Mat[2].Trim())
                                {
                                    case "N/A":
                                        substrate_ECC.Add(double.Parse(D_Mat[9]));
                                        Names.Add(Name);
                                        break;
                                    case "Normal Weight CMU":
                                        layers.Add(ABS_Layer.CreateSolid(double.Parse(D_Mat[6]), double.Parse(D_Mat[3]), double.Parse(D_Mat[11]), double.Parse(D_Mat[12]), double.Parse(D_Mat[7]), "Normal Weight CMU"));
                                        break;
                                    case "Gypsum Wall Board":
                                        layers.Add(ABS_Layer.CreateSolid(double.Parse(D_Mat[6]), double.Parse(D_Mat[3]), double.Parse(D_Mat[11]), double.Parse(D_Mat[12]), double.Parse(D_Mat[7]), "Gypsum Wall Board"));
                                        break;
                                    case "Batt Insulation":
                                        layers.Add(ABS_Layer.Create_Miki(double.Parse(D_Mat[6]), 25000, .99, 0, 0, 0, 0, 0, 0, 0, 22, "Generic Batt"));
                                        break;
                                    case "airspace":
                                        layers.Add(ABS_Layer.Airspace(double.Parse(D_Mat[6])));
                                        break;
                                }
                            }
                        }
                        catch (System.Exception)
                        {
                        }
                        i++;
                    } while (i < lines.Count && lines[i].Substring(0, 2).Contains(id.ToString()));
                    Substrates.Add(layers.ToArray());
                    layers = new List<ABS_Layer>();
                    id++;
                } while (i < lines.Count);
            }
            catch (System.Exception)
            {
            }
        }

        /// <summary>
        /// Returns the names of all user defined materials for display in an interface component.
        /// </summary>
        /// <returns></returns>
        public List<string> Names_Abs()
        {
            return Names;
        }
    }
}