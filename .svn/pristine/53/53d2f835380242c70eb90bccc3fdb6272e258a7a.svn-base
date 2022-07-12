#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/Options.h"

#include "MathTools/NumericalAnalysis/Constant1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"
#include "MathTools/NumericalAnalysis/Linear1DRealFunction.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


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
bool LevelStatisticsParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings, Abstract1DRealFunction* densityOfStates);


// perform level statistics on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = minimum level spacing 
// maxAverageSpacing = maximum level spacing 
// averageSpacing =average level spacing 
// nbrSpacings = number of level spacings
// nbrBins = number of bins for th level statistics
// binSize = level spacing range for each bin
// nbrSpacingPerBin = array that contains the number of level spacing per bin
// nbrRejectedSpacings = reference on the number of rejected level spacings (i.e. that cannot be stored in any bin)
// nbrAcceptedSpacings = reference on the number of accepted level spacings (i.e. that can be stored in a bin)
// ratio = reference on the ratio of the adjacent gaps
// varianceRatio = reference on the variance of ratio of the adjacent gaps
// nbrRatios = reference on the number of ratios
void LevelStatisticsPerformLevelStatistics(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					   double minAverageSpacing, double maxAverageSpacing, double averageSpacing, long nbrSpacings,
					   int nbrBins, double binSize, long* nbrSpacingPerBin, long& nbrRejectedSpacings, long& nbrAcceptedSpacings, 
					   Abstract1DRealFunction* densityOfStates, double densityOfStatesThreshold,
					   double& ratio, double& varianceRatio, double& nbrRatios);



int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("LevelStatistics" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectrum", "name of the file that contains the spectrum");
  (*SystemGroup) += new SingleStringOption ('\n', "multiple-spectra", "column formated ASCII file that lists all the spectra to analyze");
  (*SystemGroup) += new SingleIntegerOption ('c', "energy-column", "index of the column that contains the energies (0 being the first column)", 1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-column", "index of the first column that gives the quantum numbers (0 being the first column)", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-column", "index of the last column that gives the quantum numbers (negative if it is the one just before the energy column)", -1);
  (*SystemGroup) += new SingleIntegerOption ('d', "degeneracy-column", "index of the optional column that contains the energies degeneracy not appearing explicitly in the spectrum (must be larger than energy-column)", 0);
  (*SystemGroup) += new BooleanOption('\n', "z2-symmetry", "assume that the spectrum is invariant under the transformation Q<->-Q and that only the Q>=0 are available in the spectrum");
  (*SystemGroup) += new SingleStringOption ('\n', "unfold-spectrum", "unfold the spectrum using the tabulated density of states in an ASCII file");
  (*SystemGroup) += new SingleIntegerOption('\n', "filter-column", "if equal or greater than 0, apply a filter on the corresponding quantum number column", -1);
  (*SystemGroup) += new SingleIntegerOption('\n', "filter-value", "value of the quantum number that should be kept if --filter-value is applied", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "window-min", "when having a single quantum number sector, set the index of the first eigenvalue to consider", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "window-max", "when having a single quantum number sector, set the index of the last eigenvalue to consider (negative if up to the last available eigenvalue)", -1);
  (*SystemGroup) += new BooleanOption('\n', "discard-quantumnumbers", "do not perform the level statistics per quantum number sector");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the output file where the level statistics will be stored (if none, try to deduce it from either --spectrum or --multiple-spectra replacing the dat extension with levelstat)");
  (*OutputGroup) += new SingleDoubleOption ('b', "bin-spacing", "spacing range for each bin", 0.1);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type LevelStatistics -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("multiple-spectra") == 0) && (Manager.GetString("spectrum") == 0))
    {
      cout << "error, at least one spectrum should be provided" << endl;
      return 0;
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") == 0)
    {
      if (Manager.GetString("multiple-spectra") == 0)
	{  
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectrum"), "dat", "levelstat");
	  if (OutputFileName == 0)
	    {
	      cout << "error, can't guess output file name from " << Manager.GetString("spectrum") << endl;
	      return 0;
	    }
	}
      else
	{
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-spectra"), "dat", "levelstat");
	  if (OutputFileName == 0)
	    {
	      cout << "error, can't guess output file name from " << Manager.GetString("multiple-spectra") << endl;
	      return 0;
	    }
	}
    }
  else
    {
      OutputFileName = new char [strlen (Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }


  double AverageSpacing = 0.0;
  double MinAverageSpacing = 1.0e300;
  double MaxAverageSpacing = 0.0;
  long NbrSpacings = 0l;
  int* SpectrumSize;
  double** Spectrum;
  int* SpectrumWeight;
  int NbrSectors;

  Abstract1DRealFunction* DensityOfStates = 0;
  double DensityOfStatesThreshold = 1e-5;
  if (Manager.GetString("unfold-spectrum") != 0)
    {
      MultiColumnASCIIFile DensityOfStateFile;
      if (DensityOfStateFile.Parse(Manager.GetString("unfold-spectrum")) == false)
	{
	  DensityOfStateFile.DumpErrors(cout);
	  return false;
	}
      Tabulated1DRealFunction TmpFunction(DensityOfStateFile.GetAsDoubleArray(0), DensityOfStateFile.GetAsDoubleArray(1), 
					  DensityOfStateFile.GetNbrLines());
      DensityOfStates = TmpFunction.GetPrimitive();
    }
  else
    {
      DensityOfStates = new Linear1DRealFunction(1.0);
    }

  if (Manager.GetString("multiple-spectra") == 0)
    {
      if (LevelStatisticsParseSpectrumFile(Manager.GetString("spectrum"), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					   MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings, DensityOfStates) == false)
	{
	  return 0;
	}
      for (int i = 0; i < NbrSectors; ++i)
	{    
	  if (SpectrumSize[i] > 0)
	    delete[] Spectrum[i];      
	}
      delete[] Spectrum;
      delete[] SpectrumSize;
      delete[] SpectrumWeight;
    }
  else
    {
      MultiColumnASCIIFile SpectraFile;
      if (SpectraFile.Parse(Manager.GetString("multiple-spectra")) == false)
	{
	  SpectraFile.DumpErrors(cout);
	  return false;
	}
      for (int i = 0; i < SpectraFile.GetNbrLines(); ++i)
	{
	  cout << "checking file " << SpectraFile(0, i) << endl;
	  if (LevelStatisticsParseSpectrumFile(SpectraFile(0, i), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					       MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings, DensityOfStates) == false)
	    {
	      return 0;
	    }
	  for (int i = 0; i < NbrSectors; ++i)
	    {    
	      if (SpectrumSize[i] > 0)
		delete[] Spectrum[i];      
	    }
	  delete[] Spectrum;
	  delete[] SpectrumSize;
	  delete[] SpectrumWeight;	  
	}
    }
  
  cout << MinAverageSpacing << " " << MaxAverageSpacing << " " << AverageSpacing << endl;
  AverageSpacing /= (double) NbrSpacings;
  MinAverageSpacing /= AverageSpacing;
  MaxAverageSpacing /= AverageSpacing;

  double BinSize = Manager.GetDouble("bin-spacing");
  int NbrBins = ((int) (MaxAverageSpacing / BinSize)) + 1;
  long* NbrSpacingPerBin = new long [NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    NbrSpacingPerBin[i] = 0l;
  long NbrRejectedSpacings = 0l;
  long NbrAcceptedSpacings = 0l;

  if (Manager.GetString("multiple-spectra") == 0)
    {
      double DummyAverageSpacing = 0.0;
      double DummyMinAverageSpacing = 1.0e300;
      double DummyMaxAverageSpacing = 0.0;
      long DummyNbrSpacings = 0l; 
      if (LevelStatisticsParseSpectrumFile(Manager.GetString("spectrum"), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					   DummyMinAverageSpacing, DummyMaxAverageSpacing, DummyAverageSpacing, DummyNbrSpacings, DensityOfStates) == false)
	{
	  return 0;
	}

      double AverageR = 0.0;
      double NbrRatios = 0.0;
      double VarianceAverageR = 0.0;
      LevelStatisticsPerformLevelStatistics(Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					    MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings,
					    NbrBins, BinSize, NbrSpacingPerBin, NbrRejectedSpacings, NbrAcceptedSpacings,
					    DensityOfStates, DensityOfStatesThreshold,
					    AverageR, VarianceAverageR, NbrRatios);
      cout << "<r>=" << (AverageR / NbrRatios) << " 0.0" << endl;
      for (int i = 0; i < NbrSectors; ++i)
	{    
	  if (SpectrumSize[i] > 0)
	    delete[] Spectrum[i];      
	}
      delete[] Spectrum;
      delete[] SpectrumSize;
      delete[] SpectrumWeight;	  
    }
  else
    {
      MultiColumnASCIIFile SpectraFile;
      if (SpectraFile.Parse(Manager.GetString("multiple-spectra")) == false)
	{
	  SpectraFile.DumpErrors(cout);
	  return false;
	}
      double AverageR = 0.0;
      double VarianceAverageR = 0.0;
      for (int i = 0; i < SpectraFile.GetNbrLines(); ++i)
	{
	  double DummyAverageSpacing = 0.0;
	  double DummyMinAverageSpacing = 1.0e300;
	  double DummyMaxAverageSpacing = 0.0;
	  long DummyNbrSpacings = 0l; 
	  double TmpAverageR = 0.0;
	  double TmpNbrRatios = 0.0;
	  double TmpVarianceAverageR = 0.0;
	  cout << "processing file " << SpectraFile(0, i) << endl;
	  if (LevelStatisticsParseSpectrumFile(SpectraFile(0, i), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					       DummyMinAverageSpacing, DummyMaxAverageSpacing, DummyAverageSpacing, DummyNbrSpacings, DensityOfStates) == false)
	    {
	      return 0;
	    }
	  
	  LevelStatisticsPerformLevelStatistics(Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
						MinAverageSpacing, MaxAverageSpacing, AverageSpacing, NbrSpacings,
						NbrBins, BinSize, NbrSpacingPerBin, NbrRejectedSpacings, NbrAcceptedSpacings,
						DensityOfStates, DensityOfStatesThreshold,
						TmpAverageR, TmpVarianceAverageR, TmpNbrRatios);
	  AverageR += (TmpAverageR / TmpNbrRatios);
	  VarianceAverageR +=  (TmpAverageR / TmpNbrRatios) * (TmpAverageR / TmpNbrRatios);
	  for (int i = 0; i < NbrSectors; ++i)
	    {    
	      if (SpectrumSize[i] > 0)
		delete[] Spectrum[i];      
	    }
	  delete[] Spectrum;
	  delete[] SpectrumSize;
	  delete[] SpectrumWeight;	  
	}
      char* OutputRFileName = ReplaceExtensionToFileName(OutputFileName, "levelstat", "rvalue");
      ofstream File;
      File.open(OutputRFileName, ios::binary | ios::out);
      File.precision(14);
      File << "# Average r = " << (AverageR / SpectraFile.GetNbrLines()) << " " 
	   << sqrt (((VarianceAverageR / SpectraFile.GetNbrLines()) - ((AverageR / SpectraFile.GetNbrLines()) * (AverageR / SpectraFile.GetNbrLines()))) / (SpectraFile.GetNbrLines() - 1)) << endl;
      File << "# file name <r>" << endl;
      for (int i = 0; i < SpectraFile.GetNbrLines(); ++i)
	{
	  File << SpectraFile(0, i) << endl;
	}
      File.close();
      cout << "<r>=" << (AverageR / SpectraFile.GetNbrLines()) << " " 
	   << sqrt (((VarianceAverageR / SpectraFile.GetNbrLines()) - ((AverageR / SpectraFile.GetNbrLines()) * (AverageR / SpectraFile.GetNbrLines()))) / (SpectraFile.GetNbrLines() - 1)) << endl;
    }

  double Sum = 0.0;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
      Sum += PSpacing;
    }

  long MaxNbrSpacingPerBin = 0l;
  long MinNbrSpacingPerBin = NbrSpacings;
  for (int i = 0 ; i < NbrBins; ++i)
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

  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);
  File << "# Min spacing " << MinAverageSpacing << endl;
  File << "# Max spacing " << MaxAverageSpacing << endl;
  File << "# Average spacing = " << AverageSpacing << endl;
  File << "# Nbr points = " << NbrSpacings << endl;
  File << "# Nbr rejected points = " << NbrRejectedSpacings << endl;
  File << "# int ds P(s) = " << Sum << endl;
  File << "# Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
  File << "# Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
  File << "# Max spacing " << MaxAverageSpacing << endl;
  File << "# s P(s)" << endl;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double PSpacing = (((double) NbrSpacingPerBin[i]) / ((double) NbrSpacings));
      File << (BinSize * (((double) i) + 0.5)) << " " << (PSpacing / BinSize) << endl;
    }  
  File.close();

  cout << "Min spacing " << MinAverageSpacing << endl;
  cout << "Max spacing " << MaxAverageSpacing << endl;
  cout << "Average spacing = " << AverageSpacing << endl;
  cout << "Nbr points = " << NbrSpacings << endl;
  cout << "Nbr rejected points = " << NbrRejectedSpacings << endl;
  cout << "int ds P(s) = " << Sum << endl;
  cout << "Min nbr spacings per bin " << MinNbrSpacingPerBin << endl;
  cout << "Max nbr spacings per bin " << MaxNbrSpacingPerBin << endl;
  delete[] NbrSpacingPerBin;

  return 0;
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

bool LevelStatisticsParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minAverageSpacing, double& maxAverageSpacing, double& averageSpacing,
				      long& nbrSpacings, Abstract1DRealFunction* densityOfStates)
{
  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(spectrumFileName) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return false;
    }
  nbrSectors = 1;
  if ((manager->GetInteger("energy-column") == 0) || (manager->GetBoolean("discard-quantumnumbers") == true))
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      spectrumSize = new int [1];
      spectrum = new double*[1];
      spectrumWeight = new int[1];
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
      if (manager->GetBoolean("discard-quantumnumbers") == true)
	{
	  SortArrayUpOrdering<double>(TmpSpectrum, TmpSize);
	}
      if ((manager->GetInteger("window-min") > 0) || (manager->GetInteger("window-max") >= 0))
	{
	  int MaxIndex = manager->GetInteger("window-max");
	  if (MaxIndex < 0)
	    MaxIndex = TmpSize - 1;
	  int MinIndex = manager->GetInteger("window-min");
	  spectrumSize[0] = MaxIndex - MinIndex + 1;
	  spectrumWeight[0] = 1;
	  spectrum[0] = new double[spectrumSize[0]];
	  for (int i = MinIndex; i <= MaxIndex; ++i)
	    spectrum[0][i - MinIndex] = (*densityOfStates)(TmpSpectrum[i]);
	}
      else
	{
	  spectrumSize[0] = TmpSize;
	  spectrum[0] = new double [TmpSize];
	  for (int i = 0; i < TmpSize; ++i)
	    {
	      spectrum[0][i] =  (*densityOfStates)(TmpSpectrum[i]);
	    }
	  spectrumWeight[0] = 1;
	}
      delete[] TmpSpectrum;
    }
  else
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      int QuantumNumberShift = manager->GetInteger("min-column");
      int NbrQuantumNumber = manager->GetInteger("energy-column");
      if ((manager->GetInteger("max-column") >= 0) && (manager->GetInteger("max-column") < manager->GetInteger("energy-column")))
	NbrQuantumNumber = manager->GetInteger("max-column") - QuantumNumberShift + 1;
      int** QuantumNumbers =  new int*[NbrQuantumNumber];
      double* TmpSpectrum = 0;
      int* TmpDegeneracy = 0;
      if (manager->GetInteger("filter-column") >= 0)
	{
	  int FilterColumn = manager->GetInteger("filter-column");
	  int FilterValue = manager->GetInteger("filter-value");
	  int TmpActualSize = 0;
	  int* TmpFilterArray = SpectrumFile.GetAsIntegerArray(FilterColumn);
	  for (int i = 0; i < TmpSize; ++i)
	    if (TmpFilterArray[i] == FilterValue)
	      ++TmpActualSize;
	  if (TmpActualSize == 0)
	    {
	      cout << "error, applying filter leads to an empty spectrum" << endl;
	      return false;
	    }
	  for (int i = 0; i < NbrQuantumNumber; ++i)
	    {
	      QuantumNumbers[i] = new int [TmpActualSize];	  
	      int* TmpArray = SpectrumFile.GetAsIntegerArray(i);
	      TmpActualSize = 0;
	      for (int j = 0; j < TmpSize; ++j)
		{
		  if (TmpFilterArray[j] == FilterValue)
		    {
		      QuantumNumbers[i][TmpActualSize] = TmpArray[j];
		      ++TmpActualSize;
		    }
		}
	      delete[] TmpArray;
	    }
	  double* TmpSpectrum2 = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
	  TmpSpectrum = new double[TmpActualSize];
	  TmpActualSize = 0;
	  for (int j = 0; j < TmpSize; ++j)
	    {
	      if (TmpFilterArray[j] == FilterValue)
		{
		  TmpSpectrum[TmpActualSize] = TmpSpectrum2[j];
		  ++TmpActualSize;
		}
	    }	  
	  if (manager->GetInteger("degeneracy-column") >= 0)
	    {
	      TmpDegeneracy = new int[TmpActualSize];
	      int* TmpDegeneracy2 = SpectrumFile.GetAsIntegerArray(manager->GetInteger("degeneracy-column"));
	      TmpActualSize = 0;
	      for (int j = 0; j < TmpSize; ++j)
		{
		  if (TmpFilterArray[j] == FilterValue)
		    {
		      TmpDegeneracy[TmpActualSize] = TmpDegeneracy2[j];
		      ++TmpActualSize;
		    }
		}	  
	      delete[] TmpDegeneracy2;
	    }
	  delete[] TmpSpectrum2;
	  delete[] TmpFilterArray;
	  TmpSize = TmpActualSize;
	}     
      else
	{
	  for (int i = 0; i < NbrQuantumNumber; ++i)
	    QuantumNumbers[i] = SpectrumFile.GetAsIntegerArray(i + QuantumNumberShift);
	  TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
	  if (manager->GetInteger("degeneracy-column") >= 0)
	    {
	      TmpDegeneracy = SpectrumFile.GetAsIntegerArray(manager->GetInteger("degeneracy-column"));
	    }
	}
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  j = NbrQuantumNumber;
		  ++nbrSectors;
		}
	    }
	}
      spectrum = new double*[nbrSectors];
      spectrumSize = new int[nbrSectors];
      spectrumWeight = new int[nbrSectors];
      nbrSectors = 0;
      int CurrentIndex = 0;
      spectrumWeight[0] = 1;
      if (manager->GetBoolean("z2-symmetry"))
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][0] > 0)
		{
		  spectrumWeight[0] *= 2;
		}
	    }
	}
      else
	{
	  if (TmpDegeneracy != 0)
	    {
	      spectrumWeight[0] = TmpDegeneracy[0];
	    }
	}
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  spectrum[nbrSectors] = new double [i - CurrentIndex];
		  spectrumSize[nbrSectors] = i - CurrentIndex;
		  CurrentIndex = i;
		  j = NbrQuantumNumber;
		  ++nbrSectors;
		  spectrumWeight[nbrSectors] = 1;
		  if (manager->GetBoolean("z2-symmetry"))
		    {
		      for (int k = 0; k < NbrQuantumNumber; ++k)
			{
			  if (QuantumNumbers[k][i] > 0)
			    {
			      spectrumWeight[nbrSectors] *= 2;
			    }
			}
		    }
		  else
		    {
		      if (TmpDegeneracy != 0)
			{
			  spectrumWeight[nbrSectors] = TmpDegeneracy[i];
			}
		    }
		}
	    }
	}
      spectrum[nbrSectors] = new double [TmpSize - CurrentIndex];
      spectrumSize[nbrSectors]= TmpSize - CurrentIndex;
      ++nbrSectors;            
      spectrum[0][0] = (*densityOfStates)(TmpSpectrum[0]);
      nbrSectors = 0;
      CurrentIndex = 0;
      for (int i = 1; i < TmpSize; ++i)
	{
	  for (int j = 0; j < NbrQuantumNumber; ++j)
	    {
	      if (QuantumNumbers[j][i - 1] != QuantumNumbers[j][i])
		{
		  CurrentIndex = i;
		  j = NbrQuantumNumber;
		  ++nbrSectors;
		}
	    }
	  spectrum[nbrSectors][i - CurrentIndex] = (*densityOfStates)(TmpSpectrum[i]);
	}
      ++nbrSectors;            
      for (int i = 0; i < NbrQuantumNumber; ++i)
	delete[] QuantumNumbers[i];
      delete[] QuantumNumbers;
    }

  for (int i = 0; i < nbrSectors; ++i)
    {     
      if (spectrumSize[i] > 1)
	{
	  nbrSpacings += spectrumSize[i] - 1;
	  int Lim = spectrumSize[i];	  
	  for (int j = 1; j < Lim; ++j)
	    {
	      double TmpDiff = spectrum[i][j] - spectrum[i][j - 1];
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
	}
    }
  return true;
}

// perform level statistics on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minAverageSpacing = minimum level spacing 
// maxAverageSpacing = maximum level spacing 
// averageSpacing =average level spacing 
// nbrSpacings = number of level spacings
// nbrBins = number of bins for th level statistics
// binSize = level spacing range for each bin
// nbrSpacingPerBin = array that contains the number of level spacing per bin
// nbrRejectedSpacings = reference on the number of rejected level spacings (i.e. that cannot be stored in any bin)
// nbrAcceptedSpacings = reference on the number of accepted level spacings (i.e. that can be stored in a bin)
// ratio = reference on the ratio of the adjacent gaps
// varianceRatio = reference on the variance of ratio of the adjacent gaps
// nbrRatios = reference on the number of ratios

void LevelStatisticsPerformLevelStatistics(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					     double minAverageSpacing, double maxAverageSpacing, double averageSpacing, long nbrSpacings,
					     int nbrBins, double binSize, long* nbrSpacingPerBin, long& nbrRejectedSpacings, long& nbrAcceptedSpacings, 
					     Abstract1DRealFunction* densityOfStates, double densityOfStatesThreshold,
					     double& ratio, double& varianceRatio, double& nbrRatios)
{
  for (int i = 0; i < nbrSectors; ++i)
    {     
      if (spectrumSize[i] > 1)
	{
	  int Lim = spectrumSize[i] - 1;	  
	  double TmpInfDiff;
	  double TmpSupDiff;
	  for (int j = 1; j < Lim; ++j)
	    {
	      double TmpDiff = (spectrum[i][j] - spectrum[i][j - 1]) / averageSpacing;
	      int TmpIndex = int (TmpDiff / binSize);
	      if (TmpIndex < nbrBins)
		{
		  nbrSpacingPerBin[TmpIndex]++;
		  ++nbrAcceptedSpacings;
		}
	      else
		{
		  ++ nbrRejectedSpacings;
		}
	      double TmpDiff2 = (spectrum[i][j + 1] - spectrum[i][j]) / averageSpacing;
	      if (TmpDiff2 > TmpDiff)
		{
		  if (TmpDiff2 > MACHINE_PRECISION)
		    {
		      double Tmp = TmpDiff / TmpDiff2;
		      ratio += Tmp;	
		      varianceRatio += Tmp * Tmp;	      
		      nbrRatios += 1.0;
		    }
		}
	      else
		{
		  if (TmpDiff > MACHINE_PRECISION)
		    {
		      double Tmp =  TmpDiff2 / TmpDiff;
		      ratio += Tmp;	
		      varianceRatio += Tmp * Tmp;	      
		      nbrRatios += 1.0;
		    }
		}
	    }
	  double TmpDiff = (spectrum[i][Lim] - spectrum[i][Lim - 1]) / averageSpacing;
	  int TmpIndex = int (TmpDiff / binSize);
	  if (TmpIndex < nbrBins)
	    {
	      nbrSpacingPerBin[TmpIndex]++;
	      ++nbrAcceptedSpacings;
	    }
	  else
	    {
	      ++nbrRejectedSpacings;
	    }
	}
    }
}
