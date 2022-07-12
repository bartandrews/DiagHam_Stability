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
// spectrumLowestEnergy = reference on the lowest energy in the spectrum
// return value = true if no error occured
bool SpectrumSimilarityParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
					 int*& spectrumWeight, int& nbrSectors, double& spectrumLowestEnergy);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("ComputeSpectrumSimilarity" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('1', "spectrum-1", "first spectrum to compare");
  (*SystemGroup) += new SingleStringOption ('2', "spectrum-2", "second spectrum to compare");
  (*SystemGroup) += new SingleIntegerOption ('c', "energy-column", "index of the column that contains the energies (0 being the first column)", 1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-column", "index of the first column that gives the quantum numbers (0 being the first column)", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-column", "index of the last column that gives the quantum numbers (negative if it is the one just before the energy column)", -1);
  (*SystemGroup) += new SingleIntegerOption ('d', "degeneracy-column", "index of the optional column that contains the energies degeneracy not appearing explicitly in the spectrum (must be larger than energy-column)", 0);
  (*SystemGroup) += new BooleanOption('\n', "discard-quantumnumbers", "do not perform the level statistics per quantum number sector");
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

  if ((Manager.GetString("spectrum-1") == 0) || (Manager.GetString("spectrum-2") == 0))
    {
      cout << "error, two spectra should be provided" << endl;
      return 0;
    }


  int* FirstSpectrumSize;
  double** FirstSpectrum;
  int* FirstSpectrumWeight;
  int FirstSpectrumNbrSectors;
  double FirstSpectrumLowestEnegry;
  if (SpectrumSimilarityParseSpectrumFile(Manager.GetString("spectrum-1"), &Manager, FirstSpectrum, FirstSpectrumSize, FirstSpectrumWeight, 
					  FirstSpectrumNbrSectors, FirstSpectrumLowestEnegry) == false)
    {
      return 0;
    }
  int* SecondSpectrumSize;
  double** SecondSpectrum;
  int* SecondSpectrumWeight;
  int SecondSpectrumNbrSectors;
  double SecondSpectrumLowestEnegry;
  if (SpectrumSimilarityParseSpectrumFile(Manager.GetString("spectrum-1"), &Manager, SecondSpectrum, SecondSpectrumSize, SecondSpectrumWeight, 
					  SecondSpectrumNbrSectors, SecondSpectrumLowestEnegry) == false)
    {
      return 0;
    }

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
// spectrumLowestEnergy = reference on the lowest energy in the spectrum
// return value = true if no error occured

bool SpectrumSimilarityParseSpectrumFile(char* spectrumFileName, OptionManager* manager, double**& spectrum, int*& spectrumSize, 
					 int*& spectrumWeight, int& nbrSectors, double& spectrumLowestEnergy)
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
      spectrum[0] = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
      if (manager->GetBoolean("discard-quantumnumbers") == true)
	{
	  SortArrayUpOrdering<double>(spectrum[0], TmpSize);
	}
      spectrumSize[0] = TmpSize;
      spectrumWeight[0] = 1;
      spectrumLowestEnergy = spectrum[0][0];
      for (int i = 0; i < TmpSize; ++i)
	spectrum[0][i] -= spectrumLowestEnergy;
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
      for (int i = 0; i < NbrQuantumNumber; ++i)
	QuantumNumbers[i] = SpectrumFile.GetAsIntegerArray(i + QuantumNumberShift);
      TmpSpectrum = SpectrumFile.GetAsDoubleArray(manager->GetInteger("energy-column"));
      if (manager->GetInteger("degeneracy-column") >= 0)
	{
	  TmpDegeneracy = SpectrumFile.GetAsIntegerArray(manager->GetInteger("degeneracy-column"));
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
      if (TmpDegeneracy != 0)
	{
	  spectrumWeight[0] = TmpDegeneracy[0];
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
		  if (TmpDegeneracy != 0)
		    {
		      spectrumWeight[nbrSectors] = TmpDegeneracy[i];
		    }
		}
	    }
	}
      spectrum[nbrSectors] = new double [TmpSize - CurrentIndex];
      spectrumSize[nbrSectors]= TmpSize - CurrentIndex;
      ++nbrSectors;            
      spectrum[0][0] = TmpSpectrum[0];
      spectrumLowestEnergy = TmpSpectrum[0];
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
	  spectrum[nbrSectors][i - CurrentIndex] = TmpSpectrum[i];
	  if (spectrumLowestEnergy > TmpSpectrum[i])
	    spectrumLowestEnergy = TmpSpectrum[i];
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
	  int Lim = spectrumSize[i];	  
	  for (int j = 1; j < Lim; ++j)
	    {
	      spectrum[i][j] -= spectrumLowestEnergy;
	    }
	}
    }
  return true;
}

