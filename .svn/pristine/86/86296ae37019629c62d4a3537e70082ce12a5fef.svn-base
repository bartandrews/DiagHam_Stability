#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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


// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of sattes per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = reference on the minimum energy
// maxEnergy = reference on the maximum energy
// nbrStates = reference on the number of states
// return value = true if no error occured
bool DensityOfStatesParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minEnergy, double& maxEnergy, long& nbrStates, Abstract1DRealFunction* densityOfStates);


// perform the density of states on a parsed spectrum
//
// spectrum = two dimensional array where the spectrum is stored
// spectrumSize = number of levels per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = minimum level spacing 
// maxEnergy = maximum level spacing 
// nbrStates = number of states
// nbrBins = number of bins for the density of states
// binSize = density of states range for each bin
// nbrStatePerBin = array that contains the number of level spacing per bin
// nbrRejectedStates = reference on the number of states (i.e. that cannot be stored in any bin)
// nbrAcceptedStates = reference on the number of states (i.e. that can be stored in a bin)
void DensityOfStatesPerformDensityOfStates(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					   double minEnergy, double maxEnergy, long nbrStates,
					   int nbrBins, double binSize, long* nbrStatePerBin, long& nbrRejectedStates, long& nbrAcceptedStates);



int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("DensityOfStates" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectrum", "name of the file that contains the spectrum");
  (*SystemGroup) += new SingleStringOption ('\n', "multiple-spectra", "column formated ASCII file that lists all the spectra to analyze");
  (*SystemGroup) += new SingleIntegerOption ('c', "energy-column", "index of the column that contains the energies (0 being the first column)", 1);
  (*SystemGroup) += new SingleIntegerOption ('d', "degeneracy-column", "index of the optional column that contains the energies degeneracy not appearing explicitly in the spectrum (must be larger than energy-column)", 0);
  (*SystemGroup) += new BooleanOption('\n', "z2-symmetry", "assume that the spectrum is invariant under the transformation Q<->-Q and that only the Q>=0 are available in the spectrum");
  (*SystemGroup) += new SingleStringOption ('\n', "unfold-spectrum", "unfold the spectrum using the tabulated density of states in an ASCII file");
  (*SystemGroup) += new SingleIntegerOption('\n', "filter-column", "if equal or greater than 0, apply a filter on the corresponding quantum number column", -1);
  (*SystemGroup) += new SingleIntegerOption('\n', "filter-value", "value of the quantum number that should be kept if --filter-value is applied", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "window-min", "when having a single quantum number sector, set the index of the first eigenvalue to consider", 0);
  (*SystemGroup) += new SingleIntegerOption('\n', "window-max", "when having a single quantum number sector, set the index of the last eigenvalue to consider (negative if up to the last available eigenvalue)", -1);
  (*SystemGroup) += new BooleanOption('\n', "discard-quantumnumbers", "do not use any information (such as the degeneragy) from the quantum numbers");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the output file where the density of states will be stored (if none, try to deduce it from either --spectrum or --multiple-spectra replacing the dat extension with dos)");
  (*OutputGroup) += new SingleIntegerOption ('b', "nbr-bins", "number of bins", 100);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type DensityOfStates -h" << endl;
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
	  if (Manager.GetString("unfold-spectrum") != 0)
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectrum"), "dat", "unfoldeddos");
	    }
	  else
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("spectrum"), "dat", "dos");
	    }
	  if (OutputFileName == 0)
	    {
	      cout << "error, can't guess output file name from " << Manager.GetString("spectrum") << endl;
	      return 0;
	    }
	}
      else
	{
	  if (Manager.GetString("unfold-spectrum") != 0)
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-spectra"), "dat", "unfoldeddos");
	    }
	  else
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("multiple-spectra"), "dat", "dos");
	    }
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


  double MinEnergy = 1.0e300;
  double MaxEnergy = -1.0e300;
  long NbrStates = 0l;
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
      if (DensityOfStatesParseSpectrumFile(Manager.GetString("spectrum"), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					   MinEnergy, MaxEnergy, NbrStates, DensityOfStates) == false)
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
	  if (DensityOfStatesParseSpectrumFile(SpectraFile(0, i), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					       MinEnergy, MaxEnergy, NbrStates, DensityOfStates) == false)
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

  int NbrBins = Manager.GetInteger("nbr-bins");
  MaxEnergy +=  0.001 * (MaxEnergy - MinEnergy) / ((double) NbrBins);
  double BinSize = (MaxEnergy - MinEnergy) / ((double) NbrBins);
  long* NbrStatePerBin = new long [NbrBins];
  for (int i = 0 ; i < NbrBins; ++i)
    NbrStatePerBin[i] = 0l;
  long NbrRejectedStates = 0l;
  long NbrAcceptedStates = 0l;


  if (Manager.GetString("multiple-spectra") == 0)
    {
      double DummyAverageSpacing = 0.0;
      double DummyMinEnergy = 1.0e300;
      double DummyMaxEnergy = -1.0e300;
      long DummyNbrStates = 0l; 
      if (DensityOfStatesParseSpectrumFile(Manager.GetString("spectrum"), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					   DummyMinEnergy, DummyMaxEnergy, DummyNbrStates, DensityOfStates) == false)
	{
	  return 0;
	}

      DensityOfStatesPerformDensityOfStates(Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					    MinEnergy, MaxEnergy, NbrStates,
					    NbrBins, BinSize, NbrStatePerBin, NbrRejectedStates, NbrAcceptedStates);
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
	  double DummyAverageSpacing = 0.0;
	  double DummyMinEnergy = 1.0e300;
	  double DummyMaxEnergy = 0.0;
	  long DummyNbrStates = 0l; 
	  cout << "processing file " << SpectraFile(0, i) << endl;
	  if (DensityOfStatesParseSpectrumFile(SpectraFile(0, i), &Manager, Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
					       DummyMinEnergy, DummyMaxEnergy, DummyNbrStates, DensityOfStates) == false)
	    {
	      return 0;
	    }
	  
	  DensityOfStatesPerformDensityOfStates(Spectrum, SpectrumSize, SpectrumWeight, NbrSectors,
						MinEnergy, MaxEnergy, NbrStates,
						NbrBins, BinSize, NbrStatePerBin, NbrRejectedStates, NbrAcceptedStates);
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

  double Sum = 0.0;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double DOS = (((double) NbrStatePerBin[i]) / ((double) NbrStates));
      Sum += DOS;
    }

  long MaxNbrStatePerBin = 0l;
  long MinNbrStatePerBin = NbrStates;
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
  
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);
  File << "# Min energy " << MinEnergy << endl;
  File << "# Max energy " << MaxEnergy << endl;
  File << "# Nbr of states = " << NbrStates << endl;
  File << "# Nbr rejected states = " << NbrRejectedStates << endl;
  File << "# int dE rho(E) = " << Sum << endl;
  File << "# Min nbr of states per bin " << MinNbrStatePerBin << endl;
  File << "# Max nbr of states per bin " << MaxNbrStatePerBin << endl;
  File << "# Bin size " << BinSize << endl;
  File << "# E P(E)" << endl;
  for (int i = 0 ; i < NbrBins; ++i)
    {
      double DOS = (((double) NbrStatePerBin[i]) / ((double) NbrStates));
      File << (MinEnergy + (BinSize * ((double) i))) << " " << (DOS / BinSize) << endl;
    }  
  File.close();

  cout << "Min energy " << MinEnergy << endl;
  cout << "Max energy " << MaxEnergy << endl;
  cout << "Nbr of states = " << NbrStates << endl;
  cout << "Nbr rejected states = " << NbrRejectedStates << endl;
  cout << "int dE rho(E) = " << Sum << endl;
  cout << "Min nbr of states per bin " << MinNbrStatePerBin << endl;
  cout << "Max nbr of states per bin " << MaxNbrStatePerBin << endl;
  cout << "Bin size " << BinSize << endl;
  delete[] NbrStatePerBin;

  return 0;
}

// parse the information contained in a single spectrum file
//
// spectrumFileName = specrtum file name
// manager = pointer to the option manager
// spectrum = reference on the two dimensional array where the spectrum will be stored
// spectrumSize = number of sattes per quantum number sector
// spectrumWeight = degeneracy associated to each quantum number sector
// nbrSectors = number of quantum number sectors
// minEnergy = reference on the minimum energy
// maxEnergy = reference on the maximum energy
// nbrStates = reference on the number of states
// return value = true if no error occured

bool DensityOfStatesParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
				      int*& spectrumWeight, int& nbrSectors, double& minEnergy, double& maxEnergy, long& nbrStates, Abstract1DRealFunction* densityOfStates)
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
      spectrumSize = new int [1];
      spectrum = new double*[1];
      spectrumWeight = new int[1]; 
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
      int TmpSize = SpectrumFile.GetNbrLines();
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
	    spectrum[0][i - MinIndex] = TmpSpectrum[i];
	  delete[] TmpSpectrum;
	}
      else
	{
	  spectrum[0] = TmpSpectrum;
	  spectrumSize[0] = TmpSize;
	  spectrumWeight[0] = 1;
	}
   }
  else
    {
      int TmpSize = SpectrumFile.GetNbrLines();
      int NbrQuantumNumber = manager->GetInteger("energy-column");
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
	    QuantumNumbers[i] = SpectrumFile.GetAsIntegerArray(i);
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
      if (spectrumSize[i] > 0)
	{
	  nbrStates += spectrumSize[i];
	  int Lim = spectrumSize[i];	  
	  for (int j = 0; j < Lim; ++j)
	    {
	      double TmpEnergy =  spectrum[i][j];
	      if (TmpEnergy > maxEnergy)
		{
		  maxEnergy = TmpEnergy;
		}
	      if (TmpEnergy < minEnergy)
		{
		  minEnergy = TmpEnergy;
		}
	    }
	}
    }
  return true;
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
// nbrBins = number of bins for the density of states
// binSize = density of states range for each bin
// nbrStatePerBin = array that contains the number of level spacing per bin
// nbrRejectedStates = reference on the number of states (i.e. that cannot be stored in any bin)
// nbrAcceptedStates = reference on the number of states (i.e. that can be stored in a bin)

void DensityOfStatesPerformDensityOfStates(double** spectrum, int* spectrumSize, int* spectrumWeight, int nbrSectors, 
					   double minEnergy, double maxEnergy, long nbrStates,
					   int nbrBins, double binSize, long* nbrStatePerBin, long& nbrRejectedStates, long& nbrAcceptedStates)
{
  for (int i = 0; i < nbrSectors; ++i)
    {     
      if (spectrumSize[i] > 0)
	{
	  int Lim = spectrumSize[i];	  
	  for (int j = 0; j < Lim; ++j)
	    {
	      double TmpDiff = (spectrum[i][j] - minEnergy);
	      if (TmpDiff < 0.0)
		TmpDiff = 0.0;
	      int TmpIndex = int (TmpDiff / binSize);
	      if (TmpIndex < nbrBins)
		{
		  nbrStatePerBin[TmpIndex]++;
		  ++nbrAcceptedStates;
		}
	      else
		{
		  ++ nbrRejectedStates;
		}
	    }
	}
    }
}
