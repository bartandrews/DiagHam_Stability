#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleIntegerOption.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFileOption ('\0', "lz-file", "name of the file containing the lz sorted spectrum", 0);
  SingleDoubleOption PrecisionOption ('\n', "precision", "precision used to compare two energy values", 1e-7);
  BooleanOption FullOption ('f', "full", "indicates that all Lz value are contained in the file (if not L output will have one less value that the Lz counterpart)");
  SingleIntegerOption LColumnOption ('\n', "meanl-column", "index of the column that contains mean L value. If negative, use only the energy to guess the L value.", -1);
  BooleanOption LValidityOption ('\n', "check-meanl", "check if the L mean value is valid (i.e. is an integer up to a given error bar and compatible with the Lz value).");
  SingleDoubleOption LErrorOption ('\n', "meanl-error", "allowed error on the L mean value", 1e-10);
  BooleanOption LSortOption ('\n', "sort-l", "sort l values from the smallest to the largest one");
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFileOption;
  OptionList += &PrecisionOption;
  OptionList += &FullOption;
  OptionList += &LColumnOption;
  OptionList += &LValidityOption;
  OptionList += &LErrorOption;
  OptionList += &LSortOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type LzToL -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  char* InputName = InputFileOption.GetString();
  double Precision = PrecisionOption.GetDouble();

  MultiColumnASCIIFile Spectrum;
  if (Spectrum.Parse(InputName) == false)
    {
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }
  if (Spectrum.GetNbrColumns() < 2)
    {
      cout << "file should contain at least two columns" << endl;
      return -1;
    }
  if (Spectrum.GetNbrLines() < 1)
    {
      cout << "file should contain at least one line" << endl;
      return -1;
    }
  int TotalSize = Spectrum.GetNbrLines();
  int* FullLzValues = Spectrum.GetAsIntegerArray(0);
  if (FullLzValues == 0)
    {
      cout << "error while retrieving Lz values" << endl;
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }
  double* FullEigenvalues = Spectrum.GetAsDoubleArray(1);
  if (FullEigenvalues == 0)
    {
      cout << "error while retrieving energy values" << endl;
      Spectrum.DumpErrors(cout) << endl;
      return -1;
    }

  int NbrValue = 1;
  int LColumn = LColumnOption.GetInteger();
  bool LCheck = LValidityOption.GetBoolean();
  double LError = LErrorOption.GetDouble();
  int MaxNbrLzValues = 1;
  int CurrentNbLzValues = 1;
  int CurrentLzValue = FullLzValues[0];
  int TotalPos = 1;
  while (TotalPos < TotalSize)
    {
      if (FullLzValues[TotalPos] != CurrentLzValue)
	{
	  CurrentLzValue = FullLzValues[TotalPos];
	  ++NbrValue;
	  if (CurrentNbLzValues > MaxNbrLzValues)
	    MaxNbrLzValues = CurrentNbLzValues;
	  CurrentNbLzValues = 0;
	}
      ++TotalPos;
      ++CurrentNbLzValues;
    }
  int* Dimensions = new int [NbrValue];
  double** Eigenvalues = new double* [NbrValue];
  int Pos = 0;
  if (LColumn < 0)
    {
      TotalPos = 1;
      CurrentLzValue = FullLzValues[0];
      Dimensions[Pos] = 1;  
      while (TotalPos < TotalSize)
	{
	  if (FullLzValues[TotalPos] != CurrentLzValue)
	    {
	      CurrentLzValue = FullLzValues[TotalPos];
	      Eigenvalues[Pos] = new double [Dimensions[Pos]];
	      ++Pos;
	      Dimensions[Pos] = 0;
	    }
	  ++Dimensions[Pos];
	  ++TotalPos;
	}
      Eigenvalues[Pos] = new double [Dimensions[Pos]];      
      Pos = 0;
      TotalPos = 1;
      int Pos2 = 0;
      Eigenvalues[Pos][Pos2] = FullEigenvalues[0];
      ++Pos2;
      CurrentLzValue = FullLzValues[0];
      while (TotalPos < TotalSize)
	{
	  if (FullLzValues[TotalPos] != CurrentLzValue)
	    {
	      ++Pos;
	      CurrentLzValue = FullLzValues[TotalPos];
	      Pos2 = 0;
	    }
	  Eigenvalues[Pos][Pos2] = FullEigenvalues[TotalPos];
	  ++Pos2;
	  ++TotalPos;
	}
      
      int SpectrumSize = 0;
      double* Spectrum = new double [TotalSize];
      double* TmpSpectrum = new double [TotalSize];
      bool* Degeneracy = new bool [TotalSize];
      bool Flag;
      double TmpEigenvalue;
      int L = NbrValue;
      if (FullOption.GetBoolean() == false)
	{
	  for (int i = 0; i < Dimensions[NbrValue - 1]; ++i)
	    {
	      Spectrum[i] = Eigenvalues[NbrValue - 1][i];
	      ++SpectrumSize;
	    }
	  --L;
	}
      else
	{
	  for (int i = 0; i < Dimensions[NbrValue - 1]; ++i)
	    {
	      Spectrum[i] = Eigenvalues[NbrValue - 1][i];
	      cout << (NbrValue - 1) << " " << Eigenvalues[NbrValue - 1][i] << endl;
	      ++SpectrumSize;
	    }
	}
      for (int L = NbrValue - 2; L >= 0; --L)
	{
	  for (int j = 0; (j < SpectrumSize); ++j)
	    Degeneracy[j] = false;
	  int KeptValues = 0;
	  for (int i = 0; i < Dimensions[L]; ++i)
	    {
	      Flag = false;
	      TmpEigenvalue = Eigenvalues[L][i];
	      for (int j = 0; ((j < SpectrumSize) && (Flag == false)); ++j)
		if ((Degeneracy[j] == false) && (fabs((TmpEigenvalue - Spectrum[j]) / TmpEigenvalue) < Precision))
		  {
		    Flag = true;
		    Degeneracy[j] = true;
		  }
	      if (Flag == false)
		{
		  TmpSpectrum[KeptValues] = TmpEigenvalue;
		  ++KeptValues;
		  cout << L << " " << TmpEigenvalue << endl;
		}
	    }
	  for (int i = 0; i < KeptValues; ++i)
	    {
	      Spectrum[SpectrumSize] = TmpSpectrum[i];
	      Degeneracy[SpectrumSize]  = false;
	      ++SpectrumSize;
	    }
	  delete[] Eigenvalues[L];
	}
      delete[] Degeneracy;
      delete[] Spectrum;
    }
  else
    {
      if (Spectrum.GetNbrColumns() < LColumn)
	{
	  cout << "file should contain at least " << LColumn << " columns" << endl;
	  return -1;
	}
      double* FullLValues = Spectrum.GetAsDoubleArray(LColumn - 1);
      if (FullLValues == 0)
	{
	  cout << "error while retrieving mean L values" << endl;
	  Spectrum.DumpErrors(cout) << endl;
	  return -1;
	}
      for (int i = 0; i < NbrValue; ++i)
	{
	  Dimensions[i] = 0;
	  Eigenvalues[i] = new double [MaxNbrLzValues];
	}
      int LineIndex = 1;
      int RoundedLValue;
      TotalPos = 0;
      if (LCheck == true)
	{
	  while (TotalPos < TotalSize)
	    {
 	      RoundedLValue = (int) round(FullLValues[TotalPos]);
	      if ((fabs(FullLValues[TotalPos] - ((double) RoundedLValue)) > (LError * fabs((double) RoundedLValue))) || (RoundedLValue < FullLzValues[TotalPos]))
		{
		  cout << "invalid mean L value at line " << LineIndex << endl;
		  cout << "found " << FullLValues[TotalPos] << ", should be " << RoundedLValue << endl;
		  exit(1);
		}
	      if (RoundedLValue == FullLzValues[TotalPos])
		{
		  Eigenvalues[RoundedLValue][Dimensions[RoundedLValue]] = FullEigenvalues[TotalPos];
		  ++Dimensions[RoundedLValue];  
		}
	      ++LineIndex;
	    }
	  ++TotalPos;
	}
      else
	{
	  while (TotalPos < TotalSize)
	    {
	      RoundedLValue = (int) round(FullLValues[TotalPos]);
	      if (RoundedLValue == FullLzValues[TotalPos])
		{
		  Eigenvalues[RoundedLValue][Dimensions[RoundedLValue]] = FullEigenvalues[TotalPos];
		  ++Dimensions[RoundedLValue];  
		}
	      ++TotalPos;
	    }
	}      
      if (LSortOption.GetBoolean() == false)
	{
	  for (int i = NbrValue - 1; i >= 0; --i)
	    {
	      if (Dimensions[i] > 0)
		for (int j = 0; j < Dimensions[i]; ++j)
	      cout << i << " " << Eigenvalues[i][j] << endl;
	      delete[] Eigenvalues[i];
	    }     
	}
      else
	{
	  for (int i = 0; i < NbrValue; --i)
	    {
	      if (Dimensions[i] > 0)
		for (int j = 0; j < Dimensions[i]; ++j)
	      cout << i << " " << Eigenvalues[i][j] << endl;
	      delete[] Eigenvalues[i];
	    }     
	}
      delete[] FullLValues;
    }
  delete[] FullEigenvalues;
  delete[] FullLzValues;
  delete[] Eigenvalues;
  delete[] Dimensions;

  return 0;
}

