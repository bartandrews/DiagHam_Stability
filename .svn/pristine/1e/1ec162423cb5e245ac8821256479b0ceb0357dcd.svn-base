#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
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

  OptionManager Manager ("FQHES2xS2LzKzToLK" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption ('s', "spectrum", "name of the file that contains the spectrum");
  (*OutputGroup) += new  SingleStringOption ('\n', "output-file", "use a specific for the (L,K)-sorted spectrum instead of the default name deduced from the input file name (replacing .dat with _lk.dat)");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHES2xS2LzKzToLK -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  MultiColumnASCIIFile Spectrum;
  if (Spectrum.Parse(Manager.GetString("spectrum")) == false)
    {
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }

  if (Spectrum.GetNbrLines() < 1)
    {
      cout << "file should contain at least one line" << endl;
      return -1;
    }
  if (Spectrum.GetNbrColumns() < 3)
    {
      cout << "file should contain at least three columns" << endl;
      return -1;
    }
  bool FermionFlag = false;
  if (strstr(Manager.GetString("spectrum"), "fermions") != 0)
    {
      FermionFlag = true;
    }
  char* OutputFileName;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen (Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = ReplaceString(Manager.GetString("spectrum"), ".dat", "_lk.dat");
    }  
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out); 
  File.precision(14); 
  File << "# 2L 2K E" << endl;
  int TotalSize = Spectrum.GetNbrLines();
  int* FullLzValues = Spectrum.GetAsIntegerArray(0);
  if (FullLzValues == 0)
    {
      cout << "error while retrieving Lz values" << endl;
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }
  int* FullKzValues = Spectrum.GetAsIntegerArray(1);
  if (FullKzValues == 0)
    {
      cout << "error while retrieving Kz values" << endl;
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }
  double* FullEigenvalues = Spectrum.GetAsDoubleArray(2);
  if (FullEigenvalues == 0)
    {
      cout << "error while retrieving energy values" << endl;
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }

  int NbrValues = 1;
  int CurrentLzValue = FullLzValues[0];
  int CurrentKzValue = FullKzValues[0];
  int TotalPos = 1;
  int MaxNbrLzKzValues = 1;
  int MaxLzValue = CurrentLzValue;
  int MaxKzValue = CurrentKzValue;
  int CurrentNbrLzKzValues = 1;
  while (TotalPos < TotalSize)
    {
      if ((FullLzValues[TotalPos] != CurrentLzValue) || (FullKzValues[TotalPos] != CurrentKzValue))
	{
	  CurrentLzValue = FullLzValues[TotalPos];
	  CurrentKzValue = FullKzValues[TotalPos];
	  ++NbrValues;
	  if (CurrentLzValue > MaxLzValue)
	    MaxLzValue = CurrentLzValue;
	  if (CurrentKzValue > MaxKzValue)
	    MaxKzValue = CurrentKzValue;
	}
      ++TotalPos;
    }

  int* Dimensions = new int [NbrValues];
  double** Eigenvalues = new double* [NbrValues];
  bool** EigenvalueFlags = new bool* [NbrValues];
  int** SectorPositions = new int*[MaxKzValue + 1];
  for (CurrentKzValue = MaxKzValue; CurrentKzValue >= 0; CurrentKzValue -= 2)
    {
      SectorPositions[CurrentKzValue] = new int [MaxLzValue + 1];
      for (CurrentLzValue = MaxLzValue; CurrentLzValue >= 0; CurrentLzValue -= 2)
	{
	  SectorPositions[CurrentKzValue][CurrentLzValue] = -1;
	}
    }
  int Pos = 0;
  TotalPos = 1;
  CurrentLzValue = FullLzValues[0];
  CurrentKzValue = FullKzValues[0];
  SectorPositions[CurrentKzValue][CurrentLzValue] = Pos;
  Dimensions[Pos] = 1;  
  while (TotalPos < TotalSize)
    {
      if ((FullLzValues[TotalPos] != CurrentLzValue) || (FullKzValues[TotalPos] != CurrentKzValue))
	{
	  CurrentLzValue = FullLzValues[TotalPos];
	  CurrentKzValue = FullKzValues[TotalPos];
	  Eigenvalues[Pos] = new double [Dimensions[Pos]];
	  EigenvalueFlags[Pos] = new bool [Dimensions[Pos]];
	  ++Pos;
	  SectorPositions[CurrentKzValue][CurrentLzValue] = Pos;
	  Dimensions[Pos] = 0;
	}
      ++Dimensions[Pos];
      ++TotalPos;
    }

  Eigenvalues[Pos] = new double [Dimensions[Pos]];      
  EigenvalueFlags[Pos] = new bool [Dimensions[Pos]];
  Pos = 0;
  TotalPos = 1;
  int Pos2 = 0;
  Eigenvalues[Pos][Pos2] = FullEigenvalues[0];
  EigenvalueFlags[Pos][Pos2] = true;
  ++Pos2;
  CurrentLzValue = FullLzValues[0];
  CurrentKzValue = FullKzValues[0];
  while (TotalPos < TotalSize)
    {
      if ((FullLzValues[TotalPos] != CurrentLzValue) || (FullKzValues[TotalPos] != CurrentKzValue))
	{
	  ++Pos;
	  CurrentLzValue = FullLzValues[TotalPos];
	  CurrentKzValue = FullKzValues[TotalPos];
	  Pos2 = 0;
	}
      Eigenvalues[Pos][Pos2] = FullEigenvalues[TotalPos];
      EigenvalueFlags[Pos][Pos2] = true;
      ++Pos2;
      ++TotalPos;
    }
  
  double MaxError = 0.0;

  for (CurrentKzValue = MaxKzValue; CurrentKzValue >= 0; CurrentKzValue -= 2)
    {
      for (CurrentLzValue = MaxLzValue; CurrentLzValue >= 0; CurrentLzValue -= 2)
	{
	  int TmpPos1 = SectorPositions[CurrentKzValue][CurrentLzValue];
	  if (TmpPos1 < 0)
	    {
	      if (FermionFlag == false)
		{
		  cout << "error, missing sector (Lz,Kz) = (" << CurrentLzValue << ", " << CurrentKzValue << ")" << endl;
		  return 0;
		}
	    }
	  else
	    {
	      for (int i = 0; i < Dimensions[TmpPos1]; ++i)
		{
		  if (EigenvalueFlags[TmpPos1][i] == true)
		    {
		      for (int CurrentKzValue2 = CurrentKzValue; CurrentKzValue2 >= 0; CurrentKzValue2 -= 2)
			{
			  for (int CurrentLzValue2 = CurrentLzValue; CurrentLzValue2 >= 0; CurrentLzValue2 -= 2)
			    {
			      int TmpPos2 = SectorPositions[CurrentKzValue2][CurrentLzValue2];
			      if (TmpPos2 < 0)
				{
				  cout << "error, missing sector (Lz,Kz) = (" << CurrentLzValue2 << ", " << CurrentKzValue2 << ")" << endl;
				  return 0;
				}		  
			      if ((CurrentKzValue2 != CurrentKzValue) || (CurrentLzValue2 != CurrentLzValue))
				{
				  //			      cout << CurrentKzValue << " " << CurrentLzValue << " " <<  CurrentKzValue2 << " " << CurrentLzValue2 << endl;
				  double MinError = 1.0e300;
				  int TmpPos3 = -1;
				  for (int j = 0; j < Dimensions[TmpPos2]; ++j)
				    {
				      if (EigenvalueFlags[TmpPos2][j] == true)
					{
					  if (fabs(Eigenvalues[TmpPos1][i] - Eigenvalues[TmpPos2][j]) < MinError)
					    {
					      TmpPos3 = j;
					      MinError = fabs(Eigenvalues[TmpPos1][i] - Eigenvalues[TmpPos2][j]);
					    }
					}
				    }			      
				  if (TmpPos3 >= 0)
				    {
				      EigenvalueFlags[TmpPos2][TmpPos3] = false;
				      if (MinError > MaxError)
					{
					  MaxError = MinError;
					}
				    }
				}
			    }
			}
		      File << CurrentLzValue << " " << CurrentKzValue << " " << Eigenvalues[TmpPos1][i] << endl;
		      EigenvalueFlags[TmpPos1][i] = false;
		    }
		}
	    }
	}      
    }
  File.close();

  for (CurrentKzValue = MaxKzValue; CurrentKzValue >= 0; CurrentKzValue -= 2)
    {
      for (CurrentLzValue = MaxLzValue; CurrentLzValue >= 0; CurrentLzValue -= 2)
	{
	  int TmpPos = SectorPositions[CurrentKzValue][CurrentLzValue];
	  for (int i = 0; i < Dimensions[TmpPos]; ++i)
	    {
	      if (EigenvalueFlags[TmpPos][i] == true)
		{
		  cout << "orphan eigenvalue at (Lz,Kz) = (" << CurrentLzValue << ", " << CurrentKzValue << ") with energy " << Eigenvalues[TmpPos][i] << endl;
		}
	    }
	}
    }
  cout << "Maximum error = " << MaxError << endl;


  return 0;
}

