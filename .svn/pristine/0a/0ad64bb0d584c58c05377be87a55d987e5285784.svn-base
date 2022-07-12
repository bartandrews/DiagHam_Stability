#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ArrayTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/FQHEOnCylinderFileTools.h"

#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// extract the spectrum from a file
//
// spectrumFile = spectrum file name
// nbrLzSectors = reference on the number of momentum sectors
// minLzSector = reference on the minimal momentum
// nbrEnergie =  reference on the array where the number of levels per momentum sector
// energies =  reference on the array where the emergies will be stored
// return value = true if no error occured
bool FQHECylinderQuasiholesWithSpinTimeReversalSymmetryExtractSpectrum(char* spectrumFile, int& nbrLzSectors, int& minLzSector,
								       int*& nbrEnergie, double**& energies);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* OutputGroup  = new OptionGroup ("ouput options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "singlelayer-spectra", "name of the file containing the list of single layer spectra");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-states", "maximum number of states to generate per momentum sector", 100);  
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the two layer system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-ky", "twice the total momentum of the two layer system", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "charging-energy", "factor in front of the charging energy (i.e 1/(2C))", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "average-nbrparticles", "average number of particles", 0.0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "optimize-chargingenergy", "find the optimal value of the charging energy such that it is minimal at the indicated number of particles (only valid if positive)", -1);
  (*SystemGroup) += new BooleanOption ('\n', "only-optimize", "only compute the optimal value of the charging energy without generating the effective basis");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "fix-nbrparticles", "fix the number of particles for the two layer system (no restriction if negative)", -1);  
  (*SystemGroup) += new SingleDoubleOption ('\n', "degeneracy-error", "difference below which two energies are considered to be degenerate", 1.0e-12);
  (*SystemGroup) += new SingleStringOption ('\n', "directory", "use a specific directory for the input data instead of the current one (only useful when building the eigenstates in the full quasihole basis)");
  
  (*OutputGroup) += new BooleanOption ('\n', "build-eigenstates", "build the eigenstates in the full quasihole basis");
  (*OutputGroup) += new BooleanOption ('\n', "write-eigenstatebasis", "write a file that describes the effective eigenstate basis");
  (*OutputGroup) += new BooleanOption ('\n', "build-effectivehamiltonian", "compute and store all the  building blocks to generate the Hamiltonian in the effective basis");
  (*OutputGroup) += new BooleanOption ('\n', "write-fullspectrum", "write the full spectrum of the two decoupled layer problem (including the charging energy)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("singlelayer-spectra") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type FQHECylinderQuasiholesWithSpinTimeReversalSymmetryGenerateBasis -h" << endl;
      return -1;
    }

  int TwoLayerSzSector = Manager.GetInteger("total-sz");
  int TwoLayerLzSector = Manager.GetInteger("total-ky");
  int FixTotalNbrParticles = Manager.GetInteger("fix-nbrparticles");

  int MaxNbrParticlesPerLayer = -1;
  int NbrFluxQuanta = 0;
  bool Statistics = true;
  double Ratio = 0.0;
  double Perimeter = 0.0;
  int KValue = 1;
  int RValue = 2;

  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);

  MultiColumnASCIIFile SingleLayerSpectrumFile;
  if (SingleLayerSpectrumFile.Parse(Manager.GetString("singlelayer-spectra")) == false)
    {
      SingleLayerSpectrumFile.DumpErrors(cout);
      return -1;
    }
  for (int i = 0; i < SingleLayerSpectrumFile.GetNbrLines(); ++i)
    {
      int TmpNbrParticles = 0;
      if (FQHEOnCylinderFindSystemInfoFromFileName(SingleLayerSpectrumFile(0, i), TmpNbrParticles, NbrFluxQuanta,
						   Statistics, Ratio, Perimeter) == false)
	{
	  cout << "can't extract system information from file name " << SingleLayerSpectrumFile(0, i) << endl;
	  return -1;
	}
      if (TmpNbrParticles > MaxNbrParticlesPerLayer)
	{
	  MaxNbrParticlesPerLayer = TmpNbrParticles;
	}
    }
  char** SpectrumFileNames = new char* [MaxNbrParticlesPerLayer + 1];
  char* TmpOutputFileName;
  int* NbrLzSectors = new int [MaxNbrParticlesPerLayer + 1];
  int* MinLzValues = new int [MaxNbrParticlesPerLayer + 1];
  int** NbrEnergies = new int* [MaxNbrParticlesPerLayer + 1];
  double*** Energies = new double** [MaxNbrParticlesPerLayer + 1];
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      NbrLzSectors[i] = 0;
    }
  
  for (int i = 0; i < SingleLayerSpectrumFile.GetNbrLines(); ++i)
    {
      int TmpNbrParticles = 0;
      if (FQHEOnCylinderFindSystemInfoFromFileName(SingleLayerSpectrumFile(0, i), TmpNbrParticles, NbrFluxQuanta,
						   Statistics, Ratio, Perimeter) == false)
	{
	  cout << "can't extract system information from file name " << SingleLayerSpectrumFile(0, i) << endl;
	  return -1;
	}
      if (FQHECylinderQuasiholesWithSpinTimeReversalSymmetryExtractSpectrum(SingleLayerSpectrumFile(0, i), NbrLzSectors[TmpNbrParticles], MinLzValues[TmpNbrParticles],
									    NbrEnergies[TmpNbrParticles], Energies[TmpNbrParticles]) == false)
	{
	  return -1;
	}
      SpectrumFileNames[TmpNbrParticles] = new char[strlen(SingleLayerSpectrumFile(0, i)) + 1];
      strcpy (SpectrumFileNames[TmpNbrParticles], SingleLayerSpectrumFile(0, i));
      if (TmpNbrParticles == 0)
	{
	  char* TmpString = new char[128];
	  sprintf(TmpString, "sz_%d_lz_%d", TwoLayerSzSector, TwoLayerLzSector);
	  TmpOutputFileName = ReplaceString(SingleLayerSpectrumFile(0, i), "sz_0_lz", TmpString);
	  delete[] TmpString;
	}
    }

  int TotalNbrLevels = 0;

  for (int UpLayerNbrParticles = 0; UpLayerNbrParticles <= MaxNbrParticlesPerLayer; ++UpLayerNbrParticles)
    {
      if (NbrLzSectors[UpLayerNbrParticles] != 0)
	{
	  int DownLayerNbrParticles = UpLayerNbrParticles - TwoLayerSzSector;
	  if ((DownLayerNbrParticles >= 0) && (DownLayerNbrParticles <= MaxNbrParticlesPerLayer) && (NbrLzSectors[DownLayerNbrParticles] != 0) && ((FixTotalNbrParticles < 0) || (FixTotalNbrParticles == (DownLayerNbrParticles + UpLayerNbrParticles))))
	    {
	      for (int UpLayerLzSector = 0; UpLayerLzSector < NbrLzSectors[UpLayerNbrParticles]; ++UpLayerLzSector)
		{
		  int DownLayerLzSector = ((((UpLayerLzSector * 2) + MinLzValues[UpLayerNbrParticles]) - TwoLayerLzSector) -  MinLzValues[DownLayerNbrParticles]) / 2;
		  if ((DownLayerLzSector >= 0) && (DownLayerLzSector < NbrLzSectors[DownLayerNbrParticles]))
		    {
		      TotalNbrLevels += NbrEnergies[UpLayerNbrParticles][UpLayerLzSector] *  NbrEnergies[DownLayerNbrParticles][DownLayerLzSector];
		    }
		}
	    }
	} 
    }
  int* TwoLayerIndices = new int [TotalNbrLevels];
  int* TwoLayerNbrParticles = new int [TotalNbrLevels];
  int* TwoLayerUpIndices = new int [TotalNbrLevels];
  int* TwoLayerUpLzValues = new int [TotalNbrLevels];
  int* TwoLayerDownIndices = new int [TotalNbrLevels];
  double* TwoLayerEnergies = new double [TotalNbrLevels];
  double* TwoLayerSortedEnergies = new double [TotalNbrLevels];
  double* TwoLayerBareEnergies = new double [TotalNbrLevels];
  TotalNbrLevels = 0;
  double ChargingEnergy = Manager.GetDouble("charging-energy");
  for (int UpLayerNbrParticles = 0; UpLayerNbrParticles <= MaxNbrParticlesPerLayer; ++UpLayerNbrParticles)
    {
      if (NbrLzSectors[UpLayerNbrParticles] != 0)
	{
	  int DownLayerNbrParticles = UpLayerNbrParticles - TwoLayerSzSector;
	  if ((DownLayerNbrParticles >= 0) && (DownLayerNbrParticles <= MaxNbrParticlesPerLayer) && (NbrLzSectors[DownLayerNbrParticles] != 0) && ((FixTotalNbrParticles < 0) || (FixTotalNbrParticles == (DownLayerNbrParticles + UpLayerNbrParticles))))
	    {
	      int TmpNbrParticles = UpLayerNbrParticles + DownLayerNbrParticles;
	      double TmpEnergyShift = 0.0;
	      if (Manager.GetInteger("optimize-chargingenergy") < 0)
		{
		  TmpEnergyShift = (Manager.GetDouble("average-nbrparticles") - ((double) TmpNbrParticles));
		  TmpEnergyShift *= TmpEnergyShift;	      
		  TmpEnergyShift *= ChargingEnergy;
		}
	      for (int UpLayerLzSector = 0; UpLayerLzSector < NbrLzSectors[UpLayerNbrParticles]; ++UpLayerLzSector)
		{
		  int ShiftedUpLayerLzSector = ((UpLayerLzSector * 2) + MinLzValues[UpLayerNbrParticles]);
		  int DownLayerLzSector = ((ShiftedUpLayerLzSector - TwoLayerLzSector) - MinLzValues[DownLayerNbrParticles]) / 2;
		  if ((DownLayerLzSector >= 0) && (DownLayerLzSector < NbrLzSectors[DownLayerNbrParticles]))
		    {
		      for (int i = 0; i < NbrEnergies[UpLayerNbrParticles][UpLayerLzSector]; ++i)
			{
			  for (int j = 0; j < NbrEnergies[DownLayerNbrParticles][DownLayerLzSector]; ++j)
			    {
			      TwoLayerEnergies[TotalNbrLevels] = Energies[UpLayerNbrParticles][UpLayerLzSector][i] + Energies[DownLayerNbrParticles][DownLayerLzSector][j] + TmpEnergyShift;
			      TwoLayerSortedEnergies[TotalNbrLevels] = TwoLayerEnergies[TotalNbrLevels];
			      TwoLayerBareEnergies[TotalNbrLevels] = Energies[UpLayerNbrParticles][UpLayerLzSector][i] + Energies[DownLayerNbrParticles][DownLayerLzSector][j];
			      TwoLayerNbrParticles[TotalNbrLevels] = TmpNbrParticles;
			      TwoLayerIndices[TotalNbrLevels] = TotalNbrLevels;
			      TwoLayerUpLzValues[TotalNbrLevels] = ShiftedUpLayerLzSector;
			      TwoLayerUpIndices[TotalNbrLevels] = i;
			      TwoLayerDownIndices[TotalNbrLevels] = j;
			      ++TotalNbrLevels;
			    }
			}
		    }
		}
	    }
	} 
    }

  cout << "total number of levels in the Sz=" << TwoLayerSzSector << " Ky=" << TwoLayerLzSector << " : " << TotalNbrLevels << endl; 

  if (Manager.GetInteger("optimize-chargingenergy") >= 0)
    {
      double MinEnergyNbrParticlesMinus2 = 1.0e200;
      double MinEnergyNbrParticlesPlus2 = 1.0e200;
      double MinEnergyNbrParticles = 1.0e200;
      int OptimizeNbrParticles = Manager.GetInteger("optimize-chargingenergy");
      double AverageNbrParticles = Manager.GetDouble("average-nbrparticles");
      for (int i = 0; i < TotalNbrLevels; ++i)
	{
	  if (TwoLayerNbrParticles[i] == OptimizeNbrParticles)
	    {
	      if (MinEnergyNbrParticles > TwoLayerBareEnergies[i])
		MinEnergyNbrParticles = TwoLayerBareEnergies[i];
	    }
	  else
	    {
	      if (TwoLayerNbrParticles[i] == (OptimizeNbrParticles - 2))
		{
		  if (MinEnergyNbrParticlesMinus2  > TwoLayerBareEnergies[i])
		    MinEnergyNbrParticlesMinus2 = TwoLayerBareEnergies[i];
		}
	      else
		{
		  if (TwoLayerNbrParticles[i] == (OptimizeNbrParticles + 2))
		    { 
		      if (MinEnergyNbrParticlesPlus2  > TwoLayerBareEnergies[i])
			MinEnergyNbrParticlesPlus2 = TwoLayerBareEnergies[i];
		    }
		}
	    }
	} 
      ChargingEnergy = (MinEnergyNbrParticlesMinus2 - MinEnergyNbrParticlesPlus2) / (8.0 * (((double) OptimizeNbrParticles) - AverageNbrParticles));
      cout << "optimal charging energy = " << ChargingEnergy << endl;
      if (Manager.GetBoolean("only-optimize") == true)
	{
	  return 0;
	}
      for (int i = 0; i < TotalNbrLevels; ++i)
	{
	  TwoLayerEnergies[i] = TwoLayerBareEnergies[i] + (ChargingEnergy * (((double) TwoLayerNbrParticles[i]) - AverageNbrParticles) * 
							   (((double) TwoLayerNbrParticles[i]) - AverageNbrParticles));
	}
    }


  SortArrayUpOrdering<int>(TwoLayerSortedEnergies, TwoLayerIndices, TotalNbrLevels);

  if (Manager.GetBoolean("write-fullspectrum"))
    {
      char* OutputFileName = ReplaceString(TmpOutputFileName, ".dat", "_fullspectrum.dat");
      ofstream File;  
      File.open(OutputFileName, ios::binary | ios::out); 
      File.precision(14); 
      File << "# N Sz Lz N_u N_d Lz_u Lz_d E" << endl;
      for (int i = 0; i < TotalNbrLevels; ++i)
	{
	  int UpLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
	  int DownLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
	  int DownLayerLzValue = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
	  File << TwoLayerNbrParticles[TwoLayerIndices[i]] << " " << TwoLayerSzSector << " " << TwoLayerLzSector
	       << " " << UpLayerNbrParticles << " " << DownLayerNbrParticles << " " 
	       << TwoLayerUpLzValues[TwoLayerIndices[i]] << " " << DownLayerLzValue << " " << TwoLayerEnergies[TwoLayerIndices[i]] << endl;
	}     
      File.close();
    }

  int EffectiveSubspaceDimension = Manager.GetInteger("nbr-states");
  double Error = Manager.GetDouble("degeneracy-error");
  if (EffectiveSubspaceDimension >= TotalNbrLevels)
    {
      EffectiveSubspaceDimension = TotalNbrLevels;
    }
  else
    {
      while ((EffectiveSubspaceDimension > 0) && (fabs (TwoLayerEnergies[EffectiveSubspaceDimension - 1] - TwoLayerEnergies[EffectiveSubspaceDimension]) < Error))
	{
	  --EffectiveSubspaceDimension;
	}
    }
  int* TmpTwoLayerNbrParticles = new int [EffectiveSubspaceDimension];
  for (int i = 0; i < EffectiveSubspaceDimension; ++i)
    {
      TmpTwoLayerNbrParticles[i] = TwoLayerNbrParticles[TwoLayerIndices[i]];
    }
  SortArrayDownOrdering<int>(TmpTwoLayerNbrParticles, TwoLayerIndices, EffectiveSubspaceDimension);
  delete[] TmpTwoLayerNbrParticles;
  
  int** LargestUsedEigenstate = new int* [MaxNbrParticlesPerLayer + 1];
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      LargestUsedEigenstate[i] = new int [NbrLzSectors[i]];
      for (int j = 0; j < NbrLzSectors[i]; ++j)
	{
	  LargestUsedEigenstate[i][j] = -1;
	}
    }
  for (int i = 0; i < EffectiveSubspaceDimension; ++i)
    {
      int UpLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
      int DownLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
      int DownLayerLzValue = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
      if (LargestUsedEigenstate[UpLayerNbrParticles][(TwoLayerUpLzValues[TwoLayerIndices[i]] - MinLzValues[UpLayerNbrParticles]) / 2] < TwoLayerUpIndices[TwoLayerIndices[i]])
	{
	  LargestUsedEigenstate[UpLayerNbrParticles][(TwoLayerUpLzValues[TwoLayerIndices[i]] - MinLzValues[UpLayerNbrParticles]) / 2] = TwoLayerUpIndices[TwoLayerIndices[i]];
	}
      if (LargestUsedEigenstate[DownLayerNbrParticles][(DownLayerLzValue - MinLzValues[DownLayerNbrParticles]) / 2] < TwoLayerDownIndices[TwoLayerIndices[i]])
	{
	  LargestUsedEigenstate[DownLayerNbrParticles][(DownLayerLzValue - MinLzValues[DownLayerNbrParticles]) / 2] = TwoLayerDownIndices[TwoLayerIndices[i]];
	}
    }
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      for (int j = 0; j < NbrLzSectors[i]; ++j)
	{
	  if (LargestUsedEigenstate[i][j] >= 0)
	    {
	      cout << "require " << (LargestUsedEigenstate[i][j] + 1) << " eigenstates for N=" << i << " and Ky=" << ((j * 2) + MinLzValues[i]) << endl;
	    }
	}
    }
  int*** UsedEigenstateGlobalIndex = new int** [MaxNbrParticlesPerLayer + 1];
  int NbrGlobalIndices = 0;
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      UsedEigenstateGlobalIndex[i] = new int* [NbrLzSectors[i]];
      for (int j = 0; j < NbrLzSectors[i]; ++j)
	{
	  UsedEigenstateGlobalIndex[i][j] = new int[LargestUsedEigenstate[i][j] + 1];
	  for (int k = 0; k <= LargestUsedEigenstate[i][j]; ++k)
	    {
	      UsedEigenstateGlobalIndex[i][j][k] = NbrGlobalIndices;
	      ++NbrGlobalIndices;
	    }
	}
    }



  char* OutputFileNameExtension = new char[128];
  sprintf (OutputFileNameExtension, "_effective_%d.basis", EffectiveSubspaceDimension);
  char* OutputFileName = ReplaceString(TmpOutputFileName, ".dat", OutputFileNameExtension);
  ofstream File;  
  File.open(OutputFileName, ios::binary | ios::out); 
  File.precision(14); 
  File << "# N Sz Lz N_u N_d Lz_u Lz_d E E-E_c file_index_up file_index_down" << endl;
  for (int i = 0; i < EffectiveSubspaceDimension; ++i)
    {
      int UpLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
      int DownLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
      int DownLayerLzValue = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
      File << TwoLayerNbrParticles[TwoLayerIndices[i]] << " " << TwoLayerSzSector << " " << TwoLayerLzSector
	   << " " << UpLayerNbrParticles << " " << DownLayerNbrParticles << " " 
	   << TwoLayerUpLzValues[TwoLayerIndices[i]] << " " << DownLayerLzValue << " " << TwoLayerEnergies[TwoLayerIndices[i]] << " " << TwoLayerBareEnergies[TwoLayerIndices[i]]
	   << " " << UsedEigenstateGlobalIndex[UpLayerNbrParticles][(TwoLayerUpLzValues[TwoLayerIndices[i]] - MinLzValues[UpLayerNbrParticles]) / 2][TwoLayerUpIndices[TwoLayerIndices[i]]]
	   << " " << UsedEigenstateGlobalIndex[DownLayerNbrParticles][(DownLayerLzValue - MinLzValues[DownLayerNbrParticles]) / 2][TwoLayerDownIndices[TwoLayerIndices[i]]] << endl;
    }
  File.close();

  char* OutputFileNameVectorList = ReplaceExtensionToFileName(OutputFileName, "basis", "vectors");
  File.open(OutputFileNameVectorList, ios::binary | ios::out); 
  File.precision(14); 
  char* TmpExtension = new char[128];
  char** UsedEigenstateFileName = new char* [NbrGlobalIndices];
  NbrGlobalIndices = 0; 
  for (int i = 0; i <= MaxNbrParticlesPerLayer; ++i)
    {
      for (int j = 0; j < NbrLzSectors[i]; ++j)
	{
	  for (int k = 0; k <= LargestUsedEigenstate[i][j]; ++k)
	    {
	      sprintf(TmpExtension, "_%d.%d.vec", ((2 * j) + MinLzValues[i]), k);
	      UsedEigenstateFileName[NbrGlobalIndices] = ReplaceExtensionToFileName(SpectrumFileNames[i], ".dat", TmpExtension);
	      File << NbrGlobalIndices << " " << UsedEigenstateFileName[NbrGlobalIndices] << endl;
	      ++NbrGlobalIndices;
	    }
	}
    }
  File.close();

  QuasiholeOnSphereWithSpinAndPairing* InputSpace = 0;
  char* FilePrefix = new char[512];      
  if (Perimeter > 0.0)	
    {
      if (Statistics == true)
	{
	  sprintf (FilePrefix, "fermions_cylinder_perimeter_%.6f", Perimeter);
	}
      else
	{
	  sprintf (FilePrefix, "bosons_cylinder_perimeter_%.6f", Perimeter);
	}
    }
  else
    {
      if (Statistics == true)
	{
	  sprintf (FilePrefix, "fermions_cylinder_ratio_%.6f", Ratio);
	}
      else
	{
	  sprintf (FilePrefix, "bosons_cylinder_ratio_%.6f", Ratio);
	}
    }
  timeval TotalEndingTime;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
  cout << "Effective basis generated in " << Dt << "s" << endl;

  if (Manager.GetBoolean("build-eigenstates"))
    {
      if (Manager.GetBoolean("build-effectivehamiltonian"))
	{
	  InputSpace = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, TwoLayerLzSector, NbrFluxQuanta, TwoLayerSzSector, 
								Manager.GetString("directory"), FilePrefix);
	}
      else
	{
	  InputSpace = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, TwoLayerLzSector, NbrFluxQuanta, TwoLayerSzSector, 
								Manager.GetString("directory"), FilePrefix, true);
	}
      char* TmpEigenstateBasisFile = 0;
      if (Manager.GetBoolean("write-eigenstatebasis") == true)
	{
	  char* OutputVectorFileNameExtension = new char[128];
	  sprintf (OutputVectorFileNameExtension, "_effective_%d_n_0_", EffectiveSubspaceDimension);
	  char* OutputVectorFileName1 = ReplaceString(TmpOutputFileName, "_n_0_", OutputVectorFileNameExtension);
	  TmpEigenstateBasisFile = ReplaceExtensionToFileName(OutputVectorFileName1, "dat", "basis");
	  ofstream File;  
	  File.open(TmpEigenstateBasisFile, ios::binary | ios::out); 
	  File << "Basis=";
	  File.close();	  
	  delete[] OutputVectorFileNameExtension;
	  delete[] OutputVectorFileName1;
	}
      for (int i = 0; i < EffectiveSubspaceDimension; ++i)
	{
	  int UpLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
	  int DownLayerNbrParticles = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
	  int DownLayerLzValue = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
	  sprintf(TmpExtension, "_%d.%d.vec", TwoLayerUpLzValues[TwoLayerIndices[i]], TwoLayerUpIndices[TwoLayerIndices[i]]);
	  char* TmpUpLayerVectorFileName = ReplaceExtensionToFileName(SpectrumFileNames[UpLayerNbrParticles], ".dat", TmpExtension);
	  RealVector TmpUpLayerVector;
	  if (TmpUpLayerVector.ReadVector(TmpUpLayerVectorFileName) == false)
	    {
	      cout << "error, can't read " << TmpUpLayerVectorFileName << endl;
	  return -1;
	    }
	  sprintf(TmpExtension, "_%d.%d.vec", DownLayerLzValue, TwoLayerDownIndices[TwoLayerIndices[i]]);
	  char* TmpDownLayerVectorFileName = ReplaceExtensionToFileName(SpectrumFileNames[UpLayerNbrParticles], ".dat", TmpExtension);
	  RealVector TmpDownLayerVector;
	  if (TmpDownLayerVector.ReadVector(TmpDownLayerVectorFileName) == false)
	    {
	      cout << "error, can't read " << TmpUpLayerVectorFileName << endl;
	      return -1;
	    }
	  char* OutputVectorFileNameExtension = new char[128];
	  sprintf (OutputVectorFileNameExtension, "_effective_%d_fixedn_%d_n_0_", EffectiveSubspaceDimension, (UpLayerNbrParticles + DownLayerNbrParticles));
	  char* OutputVectorFileName1 = ReplaceString(TmpOutputFileName, "_n_0_", OutputVectorFileNameExtension);
	  sprintf (OutputVectorFileNameExtension, "%d.vec", i);
	  char* OutputVectorFileName2 = ReplaceExtensionToFileName(OutputVectorFileName1, "dat", OutputVectorFileNameExtension);
	  RealVector TmpTwoLayerVector = InputSpace->BuildFromTwoSingleLayerEigenstates(TmpUpLayerVector, UpLayerNbrParticles, TwoLayerUpLzValues[TwoLayerIndices[i]], 
											TmpDownLayerVector, DownLayerNbrParticles, DownLayerLzValue);
	  if (TmpTwoLayerVector.WriteVector(OutputVectorFileName2) == false)
	    {
	      cout << "error, can't write " << OutputVectorFileName2 << endl;
	      return -1;
	    }	  
	  cout << "generating eigenstate " << OutputVectorFileName2 << endl;
	  if (Manager.GetBoolean("write-eigenstatebasis") == true)
	    {
	      ofstream File;  
	      File.open(TmpEigenstateBasisFile, ios::binary | ios::out | ios::app); 
	      File << " " << OutputVectorFileName2;
	      File.close();	  
	    }
	  delete[] OutputVectorFileName1;
	  delete[] OutputVectorFileName2;
	  delete[] TmpUpLayerVectorFileName;
	  delete[] TmpDownLayerVectorFileName;
	}
      if (Manager.GetBoolean("write-eigenstatebasis") == true)
	{
	  ofstream File;  
	  File.open(TmpEigenstateBasisFile, ios::binary | ios::out | ios::app); 
	  File << endl;
	  File.close();	  
	  delete[] TmpEigenstateBasisFile;
	}
    }
      
  if (Manager.GetBoolean("build-effectivehamiltonian"))
    {
      if (InputSpace == 0)
	{
	  InputSpace = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, TwoLayerLzSector, NbrFluxQuanta, TwoLayerSzSector, 
								Manager.GetString("directory"), FilePrefix);
	}
      gettimeofday (&(TotalStartingTime), 0);
      RealVector* SingleLayerVectors = new RealVector[NbrGlobalIndices];
      for (int i = 0; i < NbrGlobalIndices; ++i)
	{
	  if (SingleLayerVectors[i].ReadVector(UsedEigenstateFileName[i]) == false)
	    {
	      cout << "can't open " << UsedEigenstateFileName[i] << endl;
	    }
	}  

      char* OutputMatrixFileNameExtension = new char[128];
      sprintf (OutputMatrixFileNameExtension, "_effective_%d_n_0_", EffectiveSubspaceDimension);
      char* OutputMatrixFileNameExtension1 = ReplaceString(TmpOutputFileName, "_n_0_", OutputMatrixFileNameExtension);
	  
      char* OutputPairingMatrixFileNameList = ReplaceString(OutputMatrixFileNameExtension1, ".dat", ".cddcdu.list");
      char* OutputConfiningUpMatrixFileNameList = ReplaceString(OutputMatrixFileNameExtension1, ".dat", ".cducu.list");
      char* OutputConfiningDownMatrixFileNameList = ReplaceString(OutputMatrixFileNameExtension1, ".dat", ".cddcd.list");
      ofstream FilePairing;
      ofstream FileConfiningUp;
      ofstream FileConfiningDown;
      FilePairing.open(OutputPairingMatrixFileNameList, ios::binary | ios::out); 
      FileConfiningUp.open(OutputConfiningUpMatrixFileNameList, ios::binary | ios::out); 
      FileConfiningDown.open(OutputConfiningDownMatrixFileNameList, ios::binary | ios::out); 
    
      for (int OrbitalIndex = 0; OrbitalIndex <= NbrFluxQuanta; ++OrbitalIndex)
	{
	  RealSymmetricMatrix TmpAdAd(EffectiveSubspaceDimension, true);
	  int TmpMomentumShift = (2 * OrbitalIndex) - NbrFluxQuanta;
	  for (int i = 0; i < EffectiveSubspaceDimension; ++i)
	    {
	      int UpLayerNbrParticlesRight = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
	      int DownLayerNbrParticlesRight  = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
	      int UpLayerLzValueRight  = TwoLayerUpLzValues[TwoLayerIndices[i]];
	      int DownLayerLzValueRight  = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
	      int UpLayerGlobalIndexRight = UsedEigenstateGlobalIndex[UpLayerNbrParticlesRight][(UpLayerLzValueRight- MinLzValues[UpLayerNbrParticlesRight]) / 2][TwoLayerUpIndices[TwoLayerIndices[i]]];
	      int DownLayerGlobalIndexRight = UsedEigenstateGlobalIndex[DownLayerNbrParticlesRight][(DownLayerLzValueRight - MinLzValues[DownLayerNbrParticlesRight]) / 2][TwoLayerDownIndices[TwoLayerIndices[i]]];	      
	      int TmpCount = 0;
	      int TmpUpLayerGlobalIndexLeft;
	      int TmpDownLayerGlobalIndexLeft;
	      for (int j = 0; j < EffectiveSubspaceDimension; ++j)
		{
		  int UpLayerNbrParticlesLeft = (TwoLayerNbrParticles[TwoLayerIndices[j]] + TwoLayerSzSector) / 2;
		  int DownLayerNbrParticlesLeft  = (TwoLayerNbrParticles[TwoLayerIndices[j]] - TwoLayerSzSector) / 2;  
		  int UpLayerLzValueLeft  = TwoLayerUpLzValues[TwoLayerIndices[j]];
		  int DownLayerLzValueLeft  = TwoLayerUpLzValues[TwoLayerIndices[j]] - TwoLayerLzSector;
		  int UpLayerGlobalIndexLeft = UsedEigenstateGlobalIndex[UpLayerNbrParticlesLeft][(UpLayerLzValueLeft- MinLzValues[UpLayerNbrParticlesLeft]) / 2][TwoLayerUpIndices[TwoLayerIndices[j]]];
		  int DownLayerGlobalIndexLeft = UsedEigenstateGlobalIndex[DownLayerNbrParticlesLeft][(DownLayerLzValueLeft - MinLzValues[DownLayerNbrParticlesLeft]) / 2][TwoLayerDownIndices[TwoLayerIndices[j]]];	      
		  if (((UpLayerNbrParticlesRight - 1) == UpLayerNbrParticlesLeft) && ((UpLayerLzValueRight - TmpMomentumShift) == UpLayerLzValueLeft) &&
		      ((DownLayerNbrParticlesRight - 1) == DownLayerNbrParticlesLeft) && ((DownLayerLzValueRight - TmpMomentumShift) == DownLayerLzValueLeft))
		    {
		      ++TmpCount;
		      TmpUpLayerGlobalIndexLeft = UpLayerGlobalIndexLeft;
		      TmpDownLayerGlobalIndexLeft = DownLayerGlobalIndexLeft;
		    }		  
		}
	      if (TmpCount > 0)
		{
		  RealVector TmpVectorUp (SingleLayerVectors[TmpUpLayerGlobalIndexLeft].GetVectorDimension());
		  RealVector TmpVectorDown (SingleLayerVectors[TmpDownLayerGlobalIndexLeft].GetVectorDimension());
		  InputSpace->Au(OrbitalIndex, SingleLayerVectors[UpLayerGlobalIndexRight], TmpVectorUp, UpLayerNbrParticlesRight, UpLayerLzValueRight);
//		  cout << TmpVectorUp << endl;
		  InputSpace->Ad(OrbitalIndex, SingleLayerVectors[DownLayerGlobalIndexRight], TmpVectorDown, DownLayerNbrParticlesRight, DownLayerLzValueRight);		  
		  for (int j = 0; j < EffectiveSubspaceDimension; ++j)
		    {
		      int UpLayerNbrParticlesLeft = (TwoLayerNbrParticles[TwoLayerIndices[j]] + TwoLayerSzSector) / 2;
		      int DownLayerNbrParticlesLeft  = (TwoLayerNbrParticles[TwoLayerIndices[j]] - TwoLayerSzSector) / 2;  
		      int UpLayerLzValueLeft  = TwoLayerUpLzValues[TwoLayerIndices[j]];
		      int DownLayerLzValueLeft  = TwoLayerUpLzValues[TwoLayerIndices[j]] - TwoLayerLzSector;
		      int UpLayerGlobalIndexLeft = UsedEigenstateGlobalIndex[UpLayerNbrParticlesLeft][(UpLayerLzValueLeft - MinLzValues[UpLayerNbrParticlesLeft]) / 2][TwoLayerUpIndices[TwoLayerIndices[j]]];
		      int DownLayerGlobalIndexLeft = UsedEigenstateGlobalIndex[DownLayerNbrParticlesLeft][(DownLayerLzValueLeft - MinLzValues[DownLayerNbrParticlesLeft]) / 2][TwoLayerDownIndices[TwoLayerIndices[j]]];	      
		      if (((UpLayerNbrParticlesRight - 1) == UpLayerNbrParticlesLeft) && ((UpLayerLzValueRight - TmpMomentumShift) == UpLayerLzValueLeft) &&
			  ((DownLayerNbrParticlesRight - 1) == DownLayerNbrParticlesLeft) && ((DownLayerLzValueRight - TmpMomentumShift) == DownLayerLzValueLeft))
			{
			  if (SingleLayerVectors[UpLayerGlobalIndexLeft].GetVectorDimension() != TmpVectorUp.GetVectorDimension())
			    {
			      cout << "error 1" << endl;
			    }
			  if (SingleLayerVectors[DownLayerGlobalIndexLeft].GetVectorDimension() != TmpVectorDown.GetVectorDimension())
			    {
			      cout << "error 2 : " << SingleLayerVectors[DownLayerGlobalIndexLeft].GetVectorDimension() << " " <<  TmpVectorDown.GetVectorDimension() << endl;
			    }
			  double Tmp = (SingleLayerVectors[DownLayerGlobalIndexLeft] * TmpVectorDown) * (SingleLayerVectors[UpLayerGlobalIndexLeft] * TmpVectorUp);		 
			  TmpAdAd.AddToMatrixElement(i, j, Tmp);
			}
		    }
		}
	    }
	  char* TmpExtension = new char [128];
	  sprintf (TmpExtension, "_cddcdu_%d.mat", OrbitalIndex);
	  char* OutputMatrixFileNameExtension2 = ReplaceString(OutputMatrixFileNameExtension1, ".dat", TmpExtension);
	  FilePairing << OutputMatrixFileNameExtension2 << " 0.0" << endl;
	  if (TmpAdAd.WriteMatrix(OutputMatrixFileNameExtension2) == false)
	    {
	      cout << "can't write " << OutputMatrixFileNameExtension2 << endl;
	      return -1;
	    }
	  delete[] TmpExtension;
	  delete[] OutputMatrixFileNameExtension2;
	}
      FilePairing.close();
      for (int OrbitalIndex = 0; OrbitalIndex <= NbrFluxQuanta; ++OrbitalIndex)
	{
	  RealSymmetricMatrix TmpAduAu(EffectiveSubspaceDimension, true);
	  RealSymmetricMatrix TmpAddAd(EffectiveSubspaceDimension, true);
	  for (int i = 0; i < EffectiveSubspaceDimension; ++i)
	    {
	      int UpLayerNbrParticlesRight = (TwoLayerNbrParticles[TwoLayerIndices[i]] + TwoLayerSzSector) / 2;
	      int DownLayerNbrParticlesRight  = (TwoLayerNbrParticles[TwoLayerIndices[i]] - TwoLayerSzSector) / 2;  
	      int UpLayerLzValueRight  = TwoLayerUpLzValues[TwoLayerIndices[i]];
	      int DownLayerLzValueRight  = TwoLayerUpLzValues[TwoLayerIndices[i]] - TwoLayerLzSector;
	      int UpLayerGlobalIndexRight =  UsedEigenstateGlobalIndex[UpLayerNbrParticlesRight][(UpLayerLzValueRight- MinLzValues[UpLayerNbrParticlesRight]) / 2][TwoLayerUpIndices[TwoLayerIndices[i]]];
	      int DownLayerGlobalIndexRight = UsedEigenstateGlobalIndex[DownLayerNbrParticlesRight][(DownLayerLzValueRight - MinLzValues[DownLayerNbrParticlesRight]) / 2][TwoLayerDownIndices[TwoLayerIndices[i]]];	      
	      RealVector TmpVectorUp (SingleLayerVectors[UpLayerGlobalIndexRight].GetVectorDimension());
	      RealVector TmpVectorDown (SingleLayerVectors[DownLayerGlobalIndexRight].GetVectorDimension());
	      InputSpace->AduAu(OrbitalIndex, SingleLayerVectors[UpLayerGlobalIndexRight], TmpVectorUp, UpLayerNbrParticlesRight, UpLayerLzValueRight);
	      InputSpace->AddAd(OrbitalIndex, SingleLayerVectors[DownLayerGlobalIndexRight], TmpVectorDown, DownLayerNbrParticlesRight, DownLayerLzValueRight);
	      for (int j = i; j < EffectiveSubspaceDimension; ++j)
		{
		  int UpLayerNbrParticlesLeft = (TwoLayerNbrParticles[TwoLayerIndices[j]] + TwoLayerSzSector) / 2;
		  int DownLayerNbrParticlesLeft  = (TwoLayerNbrParticles[TwoLayerIndices[j]] - TwoLayerSzSector) / 2;  
		  int UpLayerLzValueLeft  = TwoLayerUpLzValues[TwoLayerIndices[j]];
		  int DownLayerLzValueLeft  = TwoLayerUpLzValues[TwoLayerIndices[j]] - TwoLayerLzSector;
		  int UpLayerGlobalIndexLeft = UsedEigenstateGlobalIndex[UpLayerNbrParticlesLeft][(UpLayerLzValueLeft - MinLzValues[UpLayerNbrParticlesLeft]) / 2][TwoLayerUpIndices[TwoLayerIndices[j]]];
		  int DownLayerGlobalIndexLeft = UsedEigenstateGlobalIndex[DownLayerNbrParticlesLeft][(DownLayerLzValueLeft - MinLzValues[DownLayerNbrParticlesLeft]) / 2][TwoLayerDownIndices[TwoLayerIndices[j]]];	      
		  if ((UpLayerNbrParticlesRight == UpLayerNbrParticlesLeft) && (UpLayerLzValueRight == UpLayerLzValueLeft) &&
		      (DownLayerNbrParticlesRight == DownLayerNbrParticlesLeft) && (DownLayerLzValueRight == DownLayerLzValueLeft) &&
		      (TwoLayerDownIndices[TwoLayerIndices[i]] == TwoLayerDownIndices[TwoLayerIndices[j]]))
		    {			  
		      TmpAduAu.AddToMatrixElement(j, i, SingleLayerVectors[UpLayerGlobalIndexLeft] * TmpVectorUp);
		    }      		  
		  if ((UpLayerNbrParticlesRight == UpLayerNbrParticlesLeft) && (UpLayerLzValueRight == UpLayerLzValueLeft) &&
		      (DownLayerNbrParticlesRight == DownLayerNbrParticlesLeft) && (DownLayerLzValueRight == DownLayerLzValueLeft) &&
		      (TwoLayerUpIndices[TwoLayerIndices[i]] == TwoLayerUpIndices[TwoLayerIndices[j]]))
		    {			  
		      TmpAddAd.AddToMatrixElement(j, i, SingleLayerVectors[DownLayerGlobalIndexLeft] * TmpVectorDown);
		    }      		  
		}
	    }
	  char* TmpExtension = new char [128];
	  sprintf (TmpExtension, "_cducu_%d.mat", OrbitalIndex);
	  char* OutputMatrixFileNameExtension2 = ReplaceString(OutputMatrixFileNameExtension1, ".dat", TmpExtension);
	  FileConfiningUp << OutputMatrixFileNameExtension2 << " 0.0" << endl;
	  if (TmpAduAu.WriteMatrix(OutputMatrixFileNameExtension2) == false)
	    {
	      cout << "can't write " << OutputMatrixFileNameExtension2 << endl;
	      return -1;
	    }
	  delete[] OutputMatrixFileNameExtension2;
	  sprintf (TmpExtension, "_cddcd_%d.mat", OrbitalIndex);
	  OutputMatrixFileNameExtension2 = ReplaceString(OutputMatrixFileNameExtension1, ".dat", TmpExtension);
	  FileConfiningDown << OutputMatrixFileNameExtension2 << " 0.0" << endl;
	  if (TmpAddAd.WriteMatrix(OutputMatrixFileNameExtension2) == false)
	    {
	      cout << "can't write " << OutputMatrixFileNameExtension2 << endl;
	      return -1;
	    }
	  delete[] TmpExtension;
	  delete[] OutputMatrixFileNameExtension2;
	}
      FileConfiningDown.close();
      FileConfiningUp.close();
      delete[]  SingleLayerVectors;   
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));               
      cout << "Effective hamiltonian generated in " << Dt << "s" << endl;
    }
  File.close();
  delete[] TwoLayerEnergies;
  delete[] TwoLayerNbrParticles;
  delete[] TwoLayerIndices;
  delete[] TwoLayerUpIndices;
  delete[] TwoLayerDownIndices;
  return 0;
}

// extract the spectrum from a file
//
// spectrumFile = spectrum file name
// nbrLzSectors = reference on the number of momentum sectors
// minLzSector = reference on the minimal momentum
// nbrEnergies =  reference on the array where the number of levels per momentum sector
// energies =  reference on the array where the emergies will be stored
// return value = true if no error occured

bool FQHECylinderQuasiholesWithSpinTimeReversalSymmetryExtractSpectrum(char* spectrumFile, int& nbrLzSectors, int& minLzSector,
								       int*& nbrEnergies, double**& energies)
{
  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(spectrumFile) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return false;
    }
  if (SpectrumFile.GetNbrColumns() < 2)
    {
      cout << "too few columns in " << spectrumFile << endl;
      return false;
    }
  int TmpNbrValues = SpectrumFile.GetNbrLines();
  int* TmpLzSectors = SpectrumFile.GetAsIntegerArray(0);
  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(1);
  minLzSector = TmpLzSectors[0];
  int CurrentLzValue = minLzSector;
  nbrLzSectors = 1;
  for (int i = 1; i < TmpNbrValues; ++i)
    {
      if (CurrentLzValue != TmpLzSectors[i])
	{
	  CurrentLzValue = TmpLzSectors[i];
	  ++nbrLzSectors;
	  if (minLzSector > CurrentLzValue)
	    {
	      minLzSector = CurrentLzValue;
	    }	  
	}
    }
  nbrEnergies = new int [nbrLzSectors];
  energies = new double* [nbrLzSectors];
  for (int i = 0; i < nbrLzSectors; ++i)
    {
      nbrEnergies[i] = 0;
    }
  for (int i = 0; i < TmpNbrValues; ++i)
    {
      ++nbrEnergies[(TmpLzSectors[i] - minLzSector) / 2];
    }
  for (int i = 0; i < nbrLzSectors; ++i)
    {
      energies[i] = new double[nbrEnergies[i]];
      nbrEnergies[i] = 0;
    }
  for (int i = 0; i < TmpNbrValues; ++i)
    {
      energies[(TmpLzSectors[i] - minLzSector) / 2][nbrEnergies[(TmpLzSectors[i] - minLzSector) / 2]] = TmpEnergies[i];
      ++nbrEnergies[(TmpLzSectors[i] - minLzSector) / 2];
    }
  delete[] TmpEnergies;
  delete[] TmpLzSectors;
  return true;
}
