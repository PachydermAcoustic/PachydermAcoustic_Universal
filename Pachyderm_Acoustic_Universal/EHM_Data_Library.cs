//'Pachyderm-Acoustic: Geometrical Acoustics for Rhinoceros (GPL)   
//' 
//'This file is part of Pachyderm-Acoustic. 
//' 
//'Copyright (c) 2008-2025, Open Research in Acoustical Science and Education, Inc. - a 501(c)3 nonprofit 
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
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

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

        // Online file URLs
        private const string ACOUSTIC_MATERIALS_URL = "https://data.mendeley.com/public-files/datasets/dj3cs2psm9/files/ef58fa43-af34-4378-94f4-42e30d07f9e5/file_downloaded";
        private const string WALL_ASSEMBLIES_URL = "https://data.mendeley.com/public-files/datasets/dj3cs2psm9/files/d95e4f01-4c5e-4117-b3ce-4e887e465aae/file_downloaded";

        public EHM_Data_Library()
            : base()
        {
            // Check for updated files before loading
            CheckForUpdatedFiles();

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
                                        layers.Add(ABS_Layer.CreateSolid(double.Parse(D_Mat[6]), double.Parse(D_Mat[3]), double.Parse(D_Mat[11])*1E9, double.Parse(D_Mat[12]), double.Parse(D_Mat[7]), "Normal Weight CMU"));
                                        break;
                                    case "Gypsum Wall Board":
                                        layers.Add(ABS_Layer.CreateSolid(double.Parse(D_Mat[6]), double.Parse(D_Mat[3]), double.Parse(D_Mat[11])*1E9, double.Parse(D_Mat[12]), double.Parse(D_Mat[7]), "Gypsum Wall Board"));
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
        /// Checks for updated versions of the CSV files online and downloads them if newer versions are available.
        /// Uses shell execution to avoid obsolete HTTP libraries.
        /// </summary>
        private void CheckForUpdatedFiles()
        {
            try
            {
                string libPath = Pach_Properties.Instance.Lib_Path;

                // Check and update Acoustic Materials file
                CheckAndUpdateFile(ACOUSTIC_MATERIALS_URL,
                    Path.Combine(libPath, "Acoustic_Materials_and_ECCs.csv"));

                // Check and update Wall Assemblies file
                CheckAndUpdateFile(WALL_ASSEMBLIES_URL,
                    Path.Combine(libPath, "Wall_Assemblies_EC.csv"));
            }
            catch (System.Exception ex)
            {
                // Log error but don't prevent loading local files
                System.Diagnostics.Debug.WriteLine($"Failed to check for file updates: {ex.Message}");
            }
        }

        /// <summary>
        /// Checks if an online file is newer than the local version and downloads it if so.
        /// Uses shell execution with curl or PowerShell to avoid obsolete HTTP libraries.
        /// </summary>
        private void CheckAndUpdateFile(string url, string localFilePath)
        {
            try
            {
                // Get the last modified time of the local file
                DateTime localFileTime = DateTime.MinValue;
                if (File.Exists(localFilePath))
                {
                    localFileTime = File.GetLastWriteTimeUtc(localFilePath);
                }

                // Get remote file last modified date using HEAD request
                DateTime? remoteFileTime = GetRemoteFileLastModified(url);

                // Download if remote file is newer, couldn't get timestamp, or local file doesn't exist
                if (!File.Exists(localFilePath) || 
                    (remoteFileTime.HasValue && remoteFileTime.Value > localFileTime) ||
                    !remoteFileTime.HasValue)
                {
                    DownloadFile(url, localFilePath);
                }
            }
            catch (System.Exception ex)
            {
                System.Diagnostics.Debug.WriteLine($"Failed to check file {url}: {ex.Message}");
            }
        }

        /// <summary>
        /// Gets the last modified date of a remote file using shell execution.
        /// </summary>
        private DateTime? GetRemoteFileLastModified(string url)
        {
            try
            {
                // Try curl first (available on Windows 10+ and most Unix systems)
                string curlOutput = ExecuteShellCommand("curl", $"-sI \"{url}\"");
                if (!string.IsNullOrEmpty(curlOutput))
                {
                    return ParseLastModifiedFromHeaders(curlOutput);
                }

                // Fallback to PowerShell on Windows
                if (System.Environment.OSVersion.Platform == PlatformID.Win32NT)
                {
                    string psScript = $"try {{ $response = Invoke-WebRequest -Uri '{url}' -Method Head -TimeoutSec 30; $response.Headers['Last-Modified'] }} catch {{ $null }}";
                    string psOutput = ExecuteShellCommand("powershell", $"-Command \"{psScript}\"");
                    if (!string.IsNullOrEmpty(psOutput) && DateTime.TryParse(psOutput.Trim(), out DateTime result))
                    {
                        return result.ToUniversalTime();
                    }
                }
            }
            catch (System.Exception ex)
            {
                System.Diagnostics.Debug.WriteLine($"Failed to get remote file timestamp: {ex.Message}");
            }

            return null; // Could not determine remote file timestamp
        }

        /// <summary>
        /// Downloads a file using shell execution.
        /// </summary>
        private void DownloadFile(string url, string localFilePath)
        {
            try
            {
                // Ensure directory exists
                string directory = Path.GetDirectoryName(localFilePath);
                if (!Directory.Exists(directory))
                {
                    Directory.CreateDirectory(directory);
                }

                // Try curl first (available on Windows 10+ and most Unix systems)
                string curlResult = ExecuteShellCommand("curl", $"-s -L -o \"{localFilePath}\" \"{url}\"");
                if (File.Exists(localFilePath) && new FileInfo(localFilePath).Length > 0)
                {
                    System.Diagnostics.Debug.WriteLine($"Successfully downloaded updated file using curl: {localFilePath}");
                    return;
                }

                // Fallback to PowerShell on Windows
                if (System.Environment.OSVersion.Platform == PlatformID.Win32NT)
                {
                    string psScript = $"try {{ Invoke-WebRequest -Uri '{url}' -OutFile '{localFilePath}' -TimeoutSec 60 }} catch {{ Write-Error $_.Exception.Message }}";
                    string psResult = ExecuteShellCommand("powershell", $"-Command \"{psScript}\"");
                    
                    if (File.Exists(localFilePath) && new FileInfo(localFilePath).Length > 0)
                    {
                        System.Diagnostics.Debug.WriteLine($"Successfully downloaded updated file using PowerShell: {localFilePath}");
                        return;
                    }
                }

                // Fallback to wget on Unix systems
                string wgetResult = ExecuteShellCommand("wget", $"-q -O \"{localFilePath}\" \"{url}\"");
                if (File.Exists(localFilePath) && new FileInfo(localFilePath).Length > 0)
                {
                    System.Diagnostics.Debug.WriteLine($"Successfully downloaded updated file using wget: {localFilePath}");
                    return;
                }

                System.Diagnostics.Debug.WriteLine($"Failed to download file using any available method: {url}");
            }
            catch (System.Exception ex)
            {
                System.Diagnostics.Debug.WriteLine($"Failed to download file {url}: {ex.Message}");
            }
        }

        /// <summary>
        /// Executes a shell command and returns the output.
        /// </summary>
        private string ExecuteShellCommand(string command, string arguments)
        {
            try
            {
                using (Process process = new Process())
                {
                    process.StartInfo.FileName = command;
                    process.StartInfo.Arguments = arguments;
                    process.StartInfo.UseShellExecute = false;
                    process.StartInfo.RedirectStandardOutput = true;
                    process.StartInfo.RedirectStandardError = true;
                    process.StartInfo.CreateNoWindow = true;
                    process.StartInfo.WindowStyle = ProcessWindowStyle.Hidden;

                    process.Start();
                    string output = process.StandardOutput.ReadToEnd();
                    process.WaitForExit(30000); // 30 second timeout

                    return output;
                }
            }
            catch (System.Exception)
            {
                return string.Empty;
            }
        }

        /// <summary>
        /// Parses the Last-Modified header from HTTP response headers.
        /// </summary>
        private DateTime? ParseLastModifiedFromHeaders(string headers)
        {
            try
            {
                foreach (string line in headers.Split('\n'))
                {
                    if (line.StartsWith("Last-Modified:", StringComparison.OrdinalIgnoreCase))
                    {
                        string dateString = line.Substring("Last-Modified:".Length).Trim();
                        if (DateTime.TryParse(dateString, out DateTime result))
                        {
                            return result.ToUniversalTime();
                        }
                    }
                }
            }
            catch (System.Exception)
            {
                // Ignore parsing errors
            }

            return null;
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