#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/Endian.h"

#include "MathTools/NumericalAnalysis/Constant1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Linear1DRealFunction.h"

#include "MathTools/NumericalAnalysis/Constant1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Linear1DRealFunction.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates);


// perform the density of states on a parsed spectrum with a cut-off on the entanglement energies
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// minEnergyCutOff = mininum entanglement energy cut-off
// minEnergyCutOff = maxinum entanglement energy cut-off
void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates, double minEnergyCutOff, double maxEnergyCutOff);

// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = reference on the minimum level spacing 
// maxAverageSpacing = reference on the maximum level spacing 
// averageSpacing = reference on the average level spacing 
// nbrSpacings = reference on the number of level spacings
// return value = true if no error occured
bool LevelStatisticsParseSpectrumFile(double* spectrum, int spectrumSize,
				      double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings, Abstract1DRealFunction* densityOfStates);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainMultipleEntanglementSpectrumLevelStatistics" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectra", "name of the file that contains the binary entanglement spectra");
  (*SystemGroup) += new SingleIntegerOption('\n', "window-min", "set the index of the first entanglement spectrum to consider", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "window-max", " set the index of the last entanglement spectrum to consider (negative if up to the last available entanglement spectrum)", -1);
  (*SystemGroup) += new SingleDoubleOption('\n', "min-entenergy", "reject all entanglement energies below a given energy", 0.0);
  (*SystemGroup) += new SingleDoubleOption('\n', "max-entenergy", " reject all entanglement energies above a given energy (negative if no cut-off should be applied)", -1.0);
  (*SystemGroup) += new BooleanOption('\n', "fixed-sza", "focus on a single SzA sector");
  (*SystemGroup) += new SingleIntegerOption('\n', "sza-value", "if the option --fixed-sza is turned on, indicate which SzA sector should be considered", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the output file where the density of states will be stored (if none, try to deduce it from either --spectra replacing the dat extension with dos/levelstat)");
  (*OutputGroup) += new SingleIntegerOption ('\n', "nbr-bins", "number of bins for the density of states", 100);
  (*OutputGroup) += new SingleDoubleOption ('\n', "bin-spacing", "spacing range for each bin for the level statistics", 0.1);
  (*OutputGroup) += new BooleanOption ('\n', "discard-outputfiles", "do not save any results on disk, just display the summary using the standard output");
  (*OutputGroup) += new SingleIntegerOption ('\n', "extract-singlespectrum", "extract a single entanglement spectrum from --spectra instead of computing level statistics", -1);
  (*OutputGroup) += new BooleanOption ('\n', "only-rvalues", "only compute the average ratio of adjacent gaps for each entanglement spectrum");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainMultipleEntanglementSpectrumLevelStatistics -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if (Manager.GetString("spectra") == 0)
    {
      cout << "error, no entanglement spectrum was provided" << endl;
      return 0;
    }
  
  
  double MaxEnergy = 0.0;
  double MinEnergy = 1e300;
  double MinEnergyCutOff = 0.0;
  double MaxEnergyCutOff = 1e300;
  bool CutOffFlag = false;
  if ((Manager.GetDouble("min-entenergy") > 0.0) || (Manager.GetDouble("max-entenergy") > 0.0))
    {
      CutOffFlag  = true;
      MinEnergyCutOff = Manager.GetDouble("min-entenergy");
      if (Manager.GetDouble("max-entenergy") > 0.0)
	{
	  MaxEnergyCutOff = Manager.GetDouble("max-entenergy");
	}
    }
  long NbrLevels = 0l;
  int* SpectrumSize;
  double** Spectrum;
  int* SpectrumWeight;
  int NbrSectors;
  ifstream File;
  File.open(Manager.GetString("spectra"), ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "cannot open the file " << Manager.GetString("spectra") << endl;
      return 0;
    }

  int NbrEntanglementSpectra = 0;
  ReadLittleEndian(File, NbrEntanglementSpectra);  
  int NbrSzASectors = 0;
  ReadLittleEndian(File, NbrSzASectors);  
  int* SzaSectors = new int [NbrSzASectors];
  ReadBlockLittleEndian(File, SzaSectors, NbrSzASectors);

  int MinEntanglementSpectrumIndex = Manager.GetInteger("window-min");
  int MaxEntanglementSpectrumIndex = Manager.GetInteger("window-max");
  if (MaxEntanglementSpectrumIndex < MinEntanglementSpectrumIndex)
    {
      MaxEntanglementSpectrumIndex = NbrEntanglementSpectra - 1;
    }
//  int TotalNbrEntanglementSpectra = (MaxEntanglementSpectrumIndex - MinEntanglementSpectrumIndex + 1);
  int TotalNbrEntanglementSpectra = NbrSzASectors * (MaxEntanglementSpectrumIndex - MinEntanglementSpectrumIndex + 1);
  if ((NbrSzASectors > 1) && (Manager.GetBoolean("fixed-sza") == true))
    {
      TotalNbrEntanglementSpectra = (MaxEntanglementSpectrumIndex - MinEntanglementSpectrumIndex + 1);
    }
  if (Manager.GetInteger("extract-singlespectrum") >= 0)
    {
      MaxEntanglementSpectrumIndex = Manager.GetInteger("extract-singlespectrum");
      MinEntanglementSpectrumIndex = MaxEntanglementSpectrumIndex;
    }
  cout << "number of entanglement spectra = " << NbrEntanglementSpectra << endl;
  cout << "number of SzA sectors = " << NbrSzASectors << endl;
  cout << "will process the entanglement spectra between " << MinEntanglementSpectrumIndex << " and " << MaxEntanglementSpectrumIndex << endl;
  Spectrum = new double*[TotalNbrEntanglementSpectra];
  SpectrumSize = new int[TotalNbrEntanglementSpectra];
  int Index = 0;

  for (int i = 0; i < MinEntanglementSpectrumIndex; ++i)
    {
      for (int j = 0; j < NbrSzASectors; ++j)
	{
	  int TmpSize = 0;
	  ReadLittleEndian(File, TmpSize);
	  File.seekg (TmpSize * sizeof(double), ios::cur);
	}
    }
  if (Manager.GetInteger("extract-singlespectrum") >= 0)
    {
      char* SingleSpectrumExtension = new char[128];
      sprintf (SingleSpectrumExtension, "spec_%d.ent", MaxEntanglementSpectrumIndex);
      char* SingleSpectrumOutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", SingleSpectrumExtension);
      delete[] SingleSpectrumExtension;
      ofstream File2;
      File2.open(SingleSpectrumOutputFileName, ios::binary | ios::out);
      File2.precision(14);
      for (int j = 0; j < NbrSzASectors; ++j)
	{
	  int TmpSize = 0;
	  ReadLittleEndian(File, TmpSize);
	  double* TmpArray = new double[TmpSize];
	  ReadBlockLittleEndian(File, TmpArray, TmpSize); 
	  for (int k = 0; k <  TmpSize; ++k)
	    File2 << SzaSectors[j] << " " << TmpArray[k] << " " << (-log(TmpArray[k])) << endl;
	  delete[] TmpArray;
	}      
      File2.close();
      File.close();
      return 0;
    }
  if ((NbrSzASectors == 1) || (Manager.GetBoolean("fixed-sza") == true))
    {
      int FixedSzSector = Manager.GetInteger("sza-value");
      if (Manager.GetBoolean("fixed-sza") == false)
	FixedSzSector = SzaSectors[0];
       for (int i = MinEntanglementSpectrumIndex; i <= MaxEntanglementSpectrumIndex; ++i)
	{
	  for (int j = 0; j < NbrSzASectors; ++j)
	    {
	      if (SzaSectors[j] == FixedSzSector)
		{
		  if (CutOffFlag == false)
		    DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels);
		  else
		    DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels, MinEnergyCutOff, MaxEnergyCutOff);
		  ++Index;
		}
	      else
		{
		  int TmpSize = 0;
		  ReadLittleEndian(File, TmpSize);
		  File.seekg (TmpSize * sizeof(double), ios::cur);
		}
	    }
	}
   }
  else
    {
      for (int i = MinEntanglementSpectrumIndex; i <= MaxEntanglementSpectrumIndex; ++i)
	{
	  for (int j = 0; j < NbrSzASectors; ++j)
	    {
	      if (CutOffFlag == false)
		DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels);
	      else
		DensityOfStatesParseSpectrumFile(File, Spectrum[Index], SpectrumSize[Index], MinEnergy, MaxEnergy, NbrLevels, MinEnergyCutOff, MaxEnergyCutOff);
	      ++Index;
	    }
	}
    }
  File.close();

  char* DensityOfStateOutputFileName = 0;
  char* LevelStatisticOutputFileName = 0;
  char* ROutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      char* DensityOfStateExtension = new char[128];
      sprintf (DensityOfStateExtension, "range_%d_%d.ent.dos", MinEntanglementSpectrumIndex, MaxEntanglementSpectrumIndex);
      DensityOfStateOutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", DensityOfStateExtension);
      delete[] DensityOfStateExtension;
      char*  LevelStatisticExtension = new char[128];
      sprintf (LevelStatisticExtension, "range_%d_%d.ent.levelstat", MinEntanglementSpectrumIndex, MaxEntanglementSpectrumIndex);
      LevelStatisticOutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", LevelStatisticExtension);
      delete[] LevelStatisticExtension;
      char* ROutputExtension = new char[128];
      sprintf (ROutputExtension, "range_%d_%d.ent.rvalues", MinEntanglementSpectrumIndex, MaxEntanglementSpectrumIndex);
      ROutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectra"), "ent", ROutputExtension);
      delete[] ROutputExtension;
    }
  else
    {
      LevelStatisticOutputFileName = new char [strlen (Manager.GetString("output-file")) + 1];
      strcpy (LevelStatisticOutputFileName, Manager.GetString("output-file"));
    }

  
  int NbrBins = Manager.GetInteger("nbr-bins");
  MaxEnergy +=  0.001 * (MaxEnergy - MinEnergy) / ((double) NbrBins);
  double BinSize = (MaxEnergy - MinEnergy) / ((double) NbrBins);
  long* NbrStatePerBin = new long [NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    NbrStatePerBin[i] = 0l;
  long NbrRejectedStates = 0l;
  long NbrAcceptedStates = 0l;

  for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
    {
      for (int j = 0; j < SpectrumSize[i]; ++j)
	{
	  int TmpIndex = int ((Spectrum[i][j] - MinEnergy) / BinSize);
	  if ((TmpIndex < NbrBins) && (TmpIndex >= 0))
	    {
	      NbrStatePerBin[TmpIndex]++;
	      ++NbrAcceptedStates;
	    }
	  else
	    {
	      ++NbrRejectedStates;
	    }      
	}
    }


  long MaxNbrStatePerBin = 0l;
  long MinNbrStatePerBin = NbrLevels;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      if (NbrStatePerBin[i] > MaxNbrStatePerBin)
	{
	  MaxNbrStatePerBin = NbrStatePerBin[i];
	}
      if (NbrStatePerBin[i] < MinNbrStatePerBin)
	{
	  MinNbrStatePerBin = NbrStatePerBin[i];
	}
    }
  
  double Sum = 0.0;
  double* DensityOfStateEnergies =  new double[NbrBins];
  double* DensityOfStateValues =  new double[NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double DOS = (((double) NbrStatePerBin[i]) / ((double) NbrLevels));
      DensityOfStateEnergies[i] = (MinEnergy + (BinSize * ((double) i)));
      DensityOfStateValues[i] = (DOS / BinSize);
      Sum += DOS;
    }

  if ((DensityOfStateOutputFileName != 0) && (Manager.GetBoolean("discard-outputfiles") == false) 
      && (Manager.GetBoolean("only-rvalues") == false))
    {
      ofstream File;
      File.open(DensityOfStateOutputFileName, ios::binary | ios::out);
      File.precision(14);
      File << "# Min energy " << MinEnergy << endl;
      File << "# Max energy " << MaxEnergy << endl;
      File << "# Nbr of states = " << NbrLevels << endl;
      File << "# Nbr rejected states = " << NbrRejectedStates << endl;
      File << "# int dE rho(E) = " << Sum << endl;
      File << "# Min nbr of states per bin " << MinNbrStatePerBin << endl;
      File << "# Max nbr of states per bin " << MaxNbrStatePerBin << endl;
      File << "# Bin size " << BinSize << endl;
      File << "# E P(E)" << endl;
      for (int i = 0 ; i < NbrBins; ++i)
	{
	  double DOS = (((double) NbrStatePerBin[i]) / ((double) NbrLevels));
	  File << (MinEnergy + (BinSize * ((double) i))) << " " << (DOS / BinSize) << endl;
	}  
      File.close();
    }

  cout << "Min energy " << MinEnergy << endl;
  cout << "Max energy " << MaxEnergy << endl;
  cout << "Nbr of states = " << NbrLevels << endl;
  cout << "Nbr rejected states = " << NbrRejectedStates << endl;
  cout << "int dE rho(E) = " << Sum << endl;
  cout << "Min nbr of states per bin " << MinNbrStatePerBin << endl;
  cout << "Max nbr of states per bin " << MaxNbrStatePerBin << endl;
  cout << "Bin size " << BinSize << endl;

  Abstract1DRealFunction* DensityOfStates = 0;
  double DensityOfStatesThreshold = 1e-5;
  if (Manager.GetBoolean("only-rvalues") == false)
    {
      Tabulated1DRealFunction TmpFunction(DensityOfStateEnergies, DensityOfStateValues, NbrBins);
      DensityOfStates = TmpFunction.GetPrimitive();
    }
  else
    {
      DensityOfStates = new Linear1DRealFunction(1.0);
    }
  double AverageSpacing = 0.0;
  double MinAverageSpacing = 1.0e300;
  double MaxAverageSpacing = 0.0;
  long NbrSpacings = 0l;
  for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
    {
      LevelStatisticsParseSpectrumFile(Spectrum[i], SpectrumSize[i],
				       MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings, DensityOfStates);
    }

  if (Manager.GetBoolean("only-rvalues") == true)
    {
      double* AverageRValues =  new double[TotalNbrEntanglementSpectra];
      double* VarianceAverageRValues =  new double[TotalNbrEntanglementSpectra];
      for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
	{
	  if (SpectrumSize[i] > 1)
	    {
	      double TmpAverageR = 0.0;
	      double TmpNbrRatios = 0.0;
	      double TmpVarianceAverageR = 0.0;
	      int Lim = SpectrumSize[i] - 1;	  
	      double TmpInfDiff;
	      double TmpSupDiff;
	      for (int j = 1; j < Lim; ++j)
		{
		  double TmpDiff = (Spectrum[i][j] - Spectrum[i][j - 1]) / AverageSpacing;
		  double TmpDiff2 = (Spectrum[i][j + 1] - Spectrum[i][j]) / AverageSpacing;
		  if (TmpDiff2 > TmpDiff)
		    {
		      if (TmpDiff2 > MACHINE_PRECISION)
			{
			  double Tmp = TmpDiff / TmpDiff2;
			  TmpAverageR += Tmp;	
			  TmpVarianceAverageR += Tmp * Tmp;	      
			  TmpNbrRatios += 1.0;
			}
		    }
		  else
		    {
		      if (TmpDiff > MACHINE_PRECISION)
			{
			  double Tmp =  TmpDiff2 / TmpDiff;
			  TmpAverageR += Tmp;	
			  TmpVarianceAverageR += Tmp * Tmp;	      
			  TmpNbrRatios += 1.0;
			}
		    }
		}
	      AverageRValues[i] = (TmpAverageR / TmpNbrRatios);
	      VarianceAverageRValues[i] = sqrt (((TmpVarianceAverageR / TmpNbrRatios) - ((TmpAverageR / TmpNbrRatios) * (TmpAverageR / TmpNbrRatios))) / (TmpNbrRatios - 1));
	    }
	}
      ofstream File2;
      File2.open(ROutputFileName, ios::binary | ios::out);
      File2.precision(14);
      File2 << "# state_index <r> D<r>" << endl;
      for (int i = 0 ; i < TotalNbrEntanglementSpectra; ++i)
	{
	  if (SpectrumSize[i] > 1)
	    {
	      File2 << i << " " << AverageRValues[i] << " " <<  VarianceAverageRValues[i] << endl;
	    }  
	}
      delete[] AverageRValues;
      delete[] VarianceAverageRValues;
      File2.close();
      return 0;
    }


  cout << MinAverageSpacing << " " << MaxAverageSpacing << " " << AverageSpacing << endl;
  AverageSpacing /= (double) NbrSpacings;
  MinAverageSpacing /= AverageSpacing;
  MaxAverageSpacing /= AverageSpacing;

  double SpacingBinSize = Manager.GetDouble("bin-spacing");
  int SpacingNbrBins = ((int) (MaxAverageSpacing / SpacingBinSize)) + 1;
  long* NbrSpacingPerBin = new long [SpacingNbrBins];
  for (int i = 0 ; i < SpacingNbrBins; ++i)
    NbrSpacingPerBin[i] = 0l;
  long NbrRejectedSpacings = 0l;
  long NbrAcceptedSpacings = 0l;
  double AverageR = 0.0;
  double VarianceAverageR = 0.0;

  for (int i = 0; i < TotalNbrEntanglementSpectra; ++i)
    {
      if (SpectrumSize[i] > 1)
	{
	  double TmpAverageR = 0.0;
	  double TmpNbrRatios = 0.0;
	  double TmpVarianceAverageR = 0.0;
	  int Lim = SpectrumSize[i] - 1;	  
	  double TmpInfDiff;
	  double TmpSupDiff;
	  for (int j = 1; j < Lim; ++j)
	    {
	      double TmpDiff = (Spectrum[i][j] - Spectrum[i][j - 1]) / AverageSpacing;
	      int TmpIndex = int (TmpDiff / SpacingBinSize);
	      if (TmpIndex < SpacingNbrBins)
		{
		  NbrSpacingPerBin[TmpIndex]++;
		  ++NbrAcceptedSpacings;
		}
	      else
		{
		  ++ NbrRejectedSpacings;
		}
	      double TmpDiff2 = (Spectrum[i][j + 1] - Spectrum[i][j]) / AverageSpacing;
	      if (TmpDiff2 > TmpDiff)
		{
		  if (TmpDiff2 > MACHINE_PRECISION)
		    {
		      double Tmp = TmpDiff / TmpDiff2;
		      TmpAverageR += Tmp;	
		      TmpVarianceAverageR += Tmp * Tmp;	      
		      TmpNbrRatios += 1.0;
		    }
		}
	      else
		{
		  if (TmpDiff > MACHINE_PRECISION)
		    {
		      double Tmp =  TmpDiff2 / TmpDiff;
		      TmpAverageR += Tmp;	
		      TmpVarianceAverageR += Tmp * Tmp;	      
		      TmpNbrRatios += 1.0;
		    }
		}
	    }
	  double TmpDiff = (Spectrum[i][Lim] - Spectrum[i][Lim - 1]) / AverageSpacing;
	  int TmpIndex = int (TmpDiff / SpacingBinSize);
	  if (TmpIndex < SpacingNbrBins)
	    {
	      NbrSpacingPerBin[TmpIndex]++;
	      ++NbrAcceptedSpacings;
	    }
	  else
	    {
	      ++NbrRejectedSpacings;
	    }
	  if (TmpNbrRatios > 0.0)
	    {
	      AverageR += (TmpAverageR / TmpNbrRatios);
	      VarianceAverageR += (TmpAverageR / TmpNbrRatios) * (TmpAverageR / TmpNbrRatios);
	    }
	}
    }

  Sum = 0.0;
  for (int i = 0 ; i < SpacingNbrBins; ++i)
    {
      double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
      Sum += PSpacing;
    }

  long MaxNbrSpacingPerBin = 0l;
  long MinNbrSpacingPerBin = NbrSpacings;
  for (int i = 0 ; i < SpacingNbrBins; ++i)
    {
      if (NbrSpacingPerBin[i] > MaxNbrSpacingPerBin)
	{
	  MaxNbrSpacingPerBin = NbrSpacingPerBin[i];
	}
      if (NbrSpacingPerBin[i] < MinNbrSpacingPerBin)
	{
	  MinNbrSpacingPerBin = NbrSpacingPerBin[i];
	}
    }

  if (Manager.GetBoolean("discard-outputfiles") == false)
    {
      ofstream File2;
      File2.open(LevelStatisticOutputFileName, ios::binary | ios::out);
      File2.precision(14);
      File2 << "# Min spacing " << MinAverageSpacing << endl;
      File2 << "# Max spacing " << MaxAverageSpacing << endl;
      File2 << "# Average spacing = " << AverageSpacing << endl;
      File2 << "# Nbr points = " << NbrSpacings << endl;
      File2 << "# Nbr rejected points = " << NbrRejectedSpacings << endl;
      File2 << "# int ds P(s) = " << Sum << endl;
      File2 << "# Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
      File2 << "# Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
      File2 << "# Max spacing " << MaxAverageSpacing << endl;
      File2 << "# s P(s)" << endl;
      File2 << "# <r>=" << (AverageR / TotalNbrEntanglementSpectra) << " " 
	    << sqrt (((VarianceAverageR / TotalNbrEntanglementSpectra) - ((AverageR / TotalNbrEntanglementSpectra) * (AverageR / TotalNbrEntanglementSpectra))) / (TotalNbrEntanglementSpectra - 1)) << endl;
      for (int i = 0 ; i < SpacingNbrBins; ++i)
	{
	  double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
	  File2 << (SpacingBinSize * (((double) i) + 0.5)) << " " << (PSpacing / SpacingBinSize) << endl;
	}  
      File2.close();
    }

  cout << "Min spacing " << MinAverageSpacing << endl;
  cout << "Max spacing " << MaxAverageSpacing << endl;
  cout << "Average spacing = " << AverageSpacing << endl;
  cout << "Nbr points = " << NbrSpacings << endl;
  cout << "Nbr rejected points = " << NbrRejectedSpacings << endl;
  cout << "int ds P(s) = " << Sum << endl;
  cout << "Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
  cout << "Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
  cout << "<r>=" << (AverageR / TotalNbrEntanglementSpectra) << " " 
       << sqrt (((VarianceAverageR / TotalNbrEntanglementSpectra) - ((AverageR / TotalNbrEntanglementSpectra) * (AverageR / TotalNbrEntanglementSpectra))) / (TotalNbrEntanglementSpectra - 1)) << endl;
  delete[] NbrSpacingPerBin;
  delete[] NbrStatePerBin;

  return 0;
}

// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states

void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates)
{
  ReadLittleEndian(inputFile, spectrumSize);
  spectrum = new double[spectrumSize];
  ReadBlockLittleEndian(inputFile, spectrum, spectrumSize);  
  nbrStates += (long) spectrumSize;
  for (int i = 0; i < spectrumSize; ++i)
    {     
      if (spectrum[i] > 0.0)
	{
	  spectrum[i] = -log(spectrum[i]);
	  if (spectrum[i] < minEnergy)
	    minEnergy = spectrum[i];
	  if (spectrum[i] > maxEnergy)
	    maxEnergy = spectrum[i];
	}
      else
	{
	  spectrum[i] = 600.0;
	}
     }  
}

// perform the density of states on a parsed spectrum with a cut-off on the entanglement energies
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// minEnergyCutOff = mininum entanglement energy cut-off
// minEnergyCutOff = maxinum entanglement energy cut-off

void DensityOfStatesParseSpectrumFile(ifstream& inputFile, double*& spectrum, int& spectrumSize,
				      double& minEnergy, double& maxEnergy, long& nbrStates, double minEnergyCutOff, double maxEnergyCutOff)
{
  int TmpSpectrumSize = 0;
  ReadLittleEndian(inputFile, TmpSpectrumSize);
  double* TmpSpectrum = new double[TmpSpectrumSize];
  ReadBlockLittleEndian(inputFile, TmpSpectrum, TmpSpectrumSize);  
  spectrumSize = 0;
  for (int i = 0; i < TmpSpectrumSize; ++i)
    {     
      if (TmpSpectrum[i] > 0.0)
	{
	  TmpSpectrum[i] = -log(TmpSpectrum[i]);
	}
      else
	{
	  TmpSpectrum[i] = 600.0;
	}
      if ((TmpSpectrum[i] >= minEnergyCutOff) && (TmpSpectrum[i] <= maxEnergyCutOff))
	++spectrumSize;
    }
  spectrum = new double[spectrumSize];
  nbrStates += (long) spectrumSize;
  spectrumSize = 0;
   for (int i = 0; i < TmpSpectrumSize; ++i)
    {     
      if ((TmpSpectrum[i] >= minEnergyCutOff) && (TmpSpectrum[i] <= maxEnergyCutOff))
	{
	  if (TmpSpectrum[i] < minEnergy)
	    minEnergy = TmpSpectrum[i];
	  if (TmpSpectrum[i] > maxEnergy)
	    maxEnergy = TmpSpectrum[i];
	  spectrum[spectrumSize] = TmpSpectrum[i];
	  ++spectrumSize;
	}  
    }
   delete[] TmpSpectrum;
}

// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = reference on the minimum level spacing 
// maxAverageSpacing = reference on the maximum level spacing 
// averageSpacing = reference on the average level spacing 
// nbrSpacings = reference on the number of level spacings
// return value = true if no error occured

bool LevelStatisticsParseSpectrumFile(double* spectrum, int spectrumSize,
				      double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings, Abstract1DRealFunction* densityOfStates)
{
  for (int i = 0; i < spectrumSize; ++i)
    spectrum[i] = (*densityOfStates)(spectrum[i]);
  nbrSpacings += spectrumSize - 1;
  for (int j = 1; j < spectrumSize; ++j)
    {
      double TmpDiff = spectrum[j] - spectrum[j - 1];
      averageSpacing += TmpDiff;
      if (TmpDiff > maxAverageSpacing)
	{
	  maxAverageSpacing = TmpDiff;
	}
      if (TmpDiff < minAverageSpacing)
	{
	  minAverageSpacing = TmpDiff;
	}
    }
  return true;
}

