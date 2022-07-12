#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "Options/Options.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinSystemSzToS" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectrum", "name of the file that contains the spectrum");
  (*SystemGroup) += new SingleIntegerOption ('z', "sz-column", "index of the column that contains the Sz values (0 being the first column)", 0);
  (*SystemGroup) += new SingleIntegerOption ('c', "energy-column", "index of the column that contains the energies (0 being the first column)", 1);
  (*SystemGroup) += new SingleIntegerOption ('k', "momentum-column", "if positive, index of the column that contains the momentum (0 being the first column)", -1);
  (*OutputGroup) += new  SingleStringOption ('\n', "output-file", "use a specific for the S-sorted spectrum instead of the default name deduced from the input file name (replacing .dat with _s.dat)");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinSystemSzToS -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  MultiColumnASCIIFile SpectrumFile;
  if (SpectrumFile.Parse(Manager.GetString("spectrum")) == false)
    {
      SpectrumFile.DumpErrors(cout);
      return -1;
    }
  double** Spectrum;
  int* SpectrumSize;
  int* SzValues;
  int** Momenta;
  int NbrSzSectors = 1;

  int TmpSize = SpectrumFile.GetNbrLines();
  int* TmpSzValues = SpectrumFile.GetAsIntegerArray(Manager.GetInteger("sz-column"));
  for (int i = 1; i < TmpSize; ++i)
    {
      if (TmpSzValues[i - 1] != TmpSzValues[i])
	{
	  ++NbrSzSectors;
	}
    }

  if (Manager.GetInteger("momentum-column") < 0)
    {
      Spectrum = new double*[NbrSzSectors];
      SpectrumSize = new int[NbrSzSectors];
      SzValues = new int[NbrSzSectors];
      NbrSzSectors = 0;
      int CurrentIndex = 0;  
      SzValues[0] = TmpSzValues[0];
      for (int i = 1; i < TmpSize; ++i)
	{
	  if (TmpSzValues[i - 1] != TmpSzValues[i])
	    {
	      Spectrum[NbrSzSectors] = new double [i - CurrentIndex];
	      SpectrumSize[NbrSzSectors] = i - CurrentIndex;
	      CurrentIndex = i;
	      ++NbrSzSectors;
	      SzValues[NbrSzSectors] = TmpSzValues[i];
	    }
	}
      Spectrum[NbrSzSectors] = new double [TmpSize - CurrentIndex];
      SpectrumSize[NbrSzSectors]= TmpSize - CurrentIndex;
      ++NbrSzSectors; 
           
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(Manager.GetInteger("energy-column"));
      Spectrum[0][0] = TmpSpectrum[0];
      NbrSzSectors = 0;
      CurrentIndex = 0;
      for (int i = 1; i < TmpSize; ++i)
	{
	  if (TmpSzValues[i - 1] != TmpSzValues[i])
	    {
	      CurrentIndex = i;
	      ++NbrSzSectors;
	    }
	  Spectrum[NbrSzSectors][i - CurrentIndex] = TmpSpectrum[i];
	}
      ++NbrSzSectors;            
    }
  else
    {
      Spectrum = new double*[NbrSzSectors];
      SpectrumSize = new int[NbrSzSectors];
      SzValues = new int[NbrSzSectors];
      Momenta = new int*[NbrSzSectors];
      NbrSzSectors = 0;
      int CurrentIndex = 0;  
      SzValues[0] = TmpSzValues[0];
      for (int i = 1; i < TmpSize; ++i)
	{
	  if (TmpSzValues[i - 1] != TmpSzValues[i])
	    {
	      Spectrum[NbrSzSectors] = new double [i - CurrentIndex];
	      Momenta[NbrSzSectors] = new int [i - CurrentIndex];
	      SpectrumSize[NbrSzSectors] = i - CurrentIndex;
	      CurrentIndex = i;
	      ++NbrSzSectors;
	      SzValues[NbrSzSectors] = TmpSzValues[i];
	    }
	}
      Spectrum[NbrSzSectors] = new double [TmpSize - CurrentIndex];
      Momenta[NbrSzSectors] = new int [TmpSize - CurrentIndex];
      SpectrumSize[NbrSzSectors]= TmpSize - CurrentIndex;
      ++NbrSzSectors; 
           
      double* TmpSpectrum = SpectrumFile.GetAsDoubleArray(Manager.GetInteger("energy-column"));
      int* TmpMomenta = SpectrumFile.GetAsIntegerArray(Manager.GetInteger("momentum-column"));
      Spectrum[0][0] = TmpSpectrum[0];
      Momenta[0][0] = TmpMomenta[0];
      NbrSzSectors = 0;
      CurrentIndex = 0;
      for (int i = 1; i < TmpSize; ++i)
	{
	  if (TmpSzValues[i - 1] != TmpSzValues[i])
	    {
	      CurrentIndex = i;
	      ++NbrSzSectors;
	    }
	  Spectrum[NbrSzSectors][i - CurrentIndex] = TmpSpectrum[i];
	  Momenta[NbrSzSectors][i - CurrentIndex] = TmpMomenta[i];
	}
      ++NbrSzSectors;            
   }
 
  char* OutputFileName;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen (Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = ReplaceString(Manager.GetString("spectrum"), ".dat", "_s.dat");
    }  
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out); 
  File.precision(14); 

  double MaxError = 0.0;
  if (Manager.GetInteger("momentum-column") < 0)
    {
      File << "# 2S Energy deg." << endl;
      for (int i = 0; i < (NbrSzSectors - 1); ++i)
	{
	  bool* TmpFlags = new bool[SpectrumSize[i]];
	  for (int j = 0; j <  SpectrumSize[i]; ++j)
	    {
	      TmpFlags[j] = true;	  
	    }
	  for (int j = 0; j <  SpectrumSize[i + 1]; ++j)
	    {
	      double MinError = 1e14;
	      int CurrentIndex = 0;
	      for (int k = 0; k < SpectrumSize[i]; ++k)
		{
		  if (TmpFlags[k] == true)
		    {
		      double Tmp = fabs(Spectrum[i + 1][j] - Spectrum[i][k]);
		      if (Tmp < MinError)
			{
			  CurrentIndex = k;
			  MinError = Tmp;
			}
		    }
		}	  
	      if (MinError > MaxError)
		{
		  MaxError = MinError;
		}
	      TmpFlags[CurrentIndex] = false;
	    }
	  for (int k = 0; k < SpectrumSize[i]; ++k)
	    {
	      if (TmpFlags[k] == true)
		{
		  File << SzValues[i] << " " << Spectrum[i][k] << " " << (SzValues[i] + 1) << endl;
		}
	    }
	  delete[] TmpFlags;
	}
      for (int k = 0; k < SpectrumSize[NbrSzSectors - 1]; ++k)
	{
	  File << SzValues[NbrSzSectors - 1] << " " << Spectrum[NbrSzSectors - 1][k] << " " << (SzValues[NbrSzSectors - 1] + 1) << endl;
	}
    }
  else
    {
      File << "# 2S K Energy deg." << endl;
      for (int i = 0; i < (NbrSzSectors - 1); ++i)
	{
	  bool* TmpFlags = new bool[SpectrumSize[i]];
	  for (int j = 0; j <  SpectrumSize[i]; ++j)
	    {
	      TmpFlags[j] = true;	  
	    }
	  for (int j = 0; j <  SpectrumSize[i + 1]; ++j)
	    {
	      double MinError = 1e14;
	      int CurrentIndex = 0;
	      for (int k = 0; k < SpectrumSize[i]; ++k)
		{
		  if (TmpFlags[k] == true)
		    {
		      double Tmp = fabs(Spectrum[i + 1][j] - Spectrum[i][k]);
		      if ((Tmp < MinError) && (Momenta[i + 1][j] == Momenta[i][k]))
			{
			  CurrentIndex = k;
			  MinError = Tmp;
			}
		    }
		}	  
	      if (MinError > MaxError)
		{
		  MaxError = MinError;
		}
	      TmpFlags[CurrentIndex] = false;
	    }
	  for (int k = 0; k < SpectrumSize[i]; ++k)
	    {
	      if (TmpFlags[k] == true)
		{
		  File << SzValues[i] << " " << Momenta[i][k] << " " << Spectrum[i][k] << " " << (SzValues[i] + 1) << endl;
		}
	    }
	  delete[] TmpFlags;
	}
      for (int k = 0; k < SpectrumSize[NbrSzSectors - 1]; ++k)
	{
	  File << SzValues[NbrSzSectors - 1] << " " << Momenta[NbrSzSectors - 1][k] << " " << Spectrum[NbrSzSectors - 1][k] << " " << (SzValues[NbrSzSectors - 1] + 1) << endl;
	}
    }
  File.close();
  cout << "Maximum error = " << MaxError << endl;
  delete[] TmpSzValues;

  return 0;
}
