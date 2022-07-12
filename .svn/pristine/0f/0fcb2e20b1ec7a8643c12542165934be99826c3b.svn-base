#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Vector/RealVector.h"

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using std::cout;
using std::endl;
using std::ifstream;

void QuickSort(double* data, int* indice, int N);

int main(int argc, char** argv)
{
  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", "");
  SingleIntegerOption XCell('X', "Xcells", "number of the cells in X direction", 50);
  SingleIntegerOption YCell('Y', "Ycells", "number of the cells in Y direction", 50);
  SingleIntegerOption ZCell('Z', "Zcells", "number of the cells in Z direction", 30);
  SingleIntegerOption Division('D', "division", "number of the sub-divisions in each direction", 2);
  SingleStringOption Output('r', "output", "name of the output file", "default_output.txt");

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &XCell;
  OptionList += &YCell;
  OptionList += &ZCell;
  OptionList += &Division;
  OptionList += &Output;

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "bad options" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }
  char* FileName = InputFile.GetString();
  int M = XCell.GetInteger();
  int N = YCell.GetInteger();
  int P = ZCell.GetInteger();
  int Number = Division.GetInteger();
  char * FileOutput = Output.GetString();

  int NbrStateX = M; int NbrStateY = N; int NbrStateZ = P;  
  int Dim = NbrStateX * NbrStateY * NbrStateZ;
  int NbrX = M * Number; int NbrY = N * Number; int NbrZ = P * Number; 

  RealVector* SinX = new RealVector[NbrX];
  RealVector* SinY = new RealVector[NbrY];
  RealVector* SinZ = new RealVector[NbrZ];

  for (int i = 0; i < NbrX; ++i)
    {
      SinX[i] = RealVector(NbrStateX);
      for (int j = 0; j < NbrStateX; ++j)
	SinX[i][j] = sin((M_PI * (j + 1.0) * (i + 0.5)) / NbrX);
    }
  for (int i = 0; i < NbrY; ++i)
    {
      SinY[i] = RealVector(NbrStateY);
      for (int j = 0; j < NbrStateY; ++j)
	SinY[i][j] = sin((M_PI * (j + 1.0) * (i + 0.5)) / NbrY);
    }
  for (int i = 0; i < NbrZ; ++i)
    {
      SinZ[i] = RealVector(NbrStateZ);
      for (int j = 0; j < NbrStateZ; ++j)
	SinZ[i][j] = sin((M_PI * (j + 1.0) * (i + 0.5)) / NbrZ);
    }

  ifstream input(FileName);
  if (! input.is_open())
    {
      cout << "Error when open the eigenvector file: " << FileName << " . Exit now" << endl;
      exit(1);
    }

  double*** Coefficient = new double** [NbrStateX];
  for (int m = 0; m < NbrStateX; ++m)
    {
      Coefficient[m] = new double* [NbrStateY];
      for (int n = 0; n < NbrStateY; ++n)
	{
	  Coefficient[m][n] = new double [NbrStateZ];
	  for (int p = 0; p < NbrStateZ; ++p)
	    input >> Coefficient[m][n][p];
	}
    }
  input.close();

  double*** Inter1 = new double** [NbrX]; 
  double Tmp = 0.0;
  for (int i = 0; i < NbrX; ++i)
    {
      Inter1[i] = new double* [NbrStateZ];
      for (int p = 0; p < NbrStateZ; ++p)
	{
	  Inter1[i][p] = new double [NbrStateY];
	  for (int n = 0; n < NbrStateY; ++n)
	    {
	      Tmp = 0.0;	    
	      for (int m = 0; m < NbrStateX; ++m)		
		Tmp += SinX[i][m] * Coefficient[m][n][p];		  
	      Inter1[i][p][n] = Tmp;	      
	    }
	}
    }

  double*** Inter2 = new double** [NbrY];
  for (int j = 0; j < NbrY; ++j)
    {
      Inter2[j] = new double* [NbrX];
      for (int i = 0; i < NbrX; ++i)
	{
	  Inter2[j][i] = new double [NbrStateZ];
	  for (int p = 0; p < NbrStateZ; ++p)
	    {
	      Tmp = 0.0;	  
	      for (int n = 0; n < NbrStateY; ++n)
		Tmp += SinY[j][n] * Inter1[i][p][n];	    
	      Inter2[j][i][p] = Tmp;
	    }
	}
    }

  ofstream output(FileOutput);
  if (! output.is_open())
    {
      cout << "Error when open the output file: " << FileOutput << " . Exit now" << endl;
      exit(1);
    }

  int NbrTotal = NbrX * NbrY * NbrZ;
  double* Function = new double[NbrTotal];
  double Factor = 8.0 / (NbrTotal);

  int* Index = new int [NbrTotal];
  for (int i = 0; i < NbrTotal; ++i)
    Index[i] = i;

  int tmpIndice = 0;
  for (int k = 0; k < NbrZ; ++k)    
    for (int j = 0; j < NbrY; ++j)
      for (int i = 0; i < NbrX; ++i)
	{
	  Tmp = 0.0;
	  for (int p = 0; p < NbrStateZ; ++p)
	    Tmp += SinZ[k][p] * Inter2[j][i][p];	
	  Function[tmpIndice] = Tmp * Tmp * Factor;
	  ++tmpIndice;
	}

  QuickSort(Function, Index, NbrTotal);
  
  double Isovalue = 0.7;
  double Down_Isovalue = Isovalue - 0.05;
  double Up_Isovalue = Isovalue + 0.05;
  double TmpTotal = 0.0;

  int LowIndex = 0; int HighIndex = 0; int MainIndex = 0;
  int n = 0;
  for (; TmpTotal < Down_Isovalue; ++n)
    TmpTotal += Function[n];
  LowIndex = n;
  for (; TmpTotal < Isovalue; ++n)
    TmpTotal += Function[n];
  MainIndex = n;
  for (; TmpTotal < Up_Isovalue; ++n)
    TmpTotal += Function[n];
  HighIndex = n;

  cout << Function[MainIndex] << endl;
  int tmpX = 0; int tmpY = 0; int tmpZ = 0; 
  output << "#Isovalue corresponds to 0.7: " << Function[MainIndex] << '\n';
  int tmpNbr = NbrY * NbrZ; int Tmp1 = 0;
  for (n = LowIndex; n < HighIndex; ++n)
    {
      tmpX = Index[n] / tmpNbr;
      Tmp1 = Index[n] - tmpNbr * tmpX;
      tmpY = Tmp1 / NbrZ;
      tmpZ = Tmp1 - tmpY * NbrZ;

      output << (tmpX + 0.5) / Number << '\t' << (tmpY + 0.5) / Number << '\t' << (tmpZ + 0.5) / Number << '\t' << Function[n] << '\n';
    }

  output.close();
  return 0;
}

void QuickSort(double* data, int* indice, int N)
{
  int i, j, tmpIndice;
  double v, t;

  if (N <= 1) 
    return;

  // Partition elements
  v = data[0];
  i = 0;
  j = N;
  for( ; ; )
    {
      while(data[++i] > v && i < N) 
	;
      while(data[--j] < v) 
	;
      if(i >= j) 
	break;
      t = data[i]; data[i] = data[j]; data[j] = t;
      tmpIndice = indice[i]; indice[i] = indice[j]; indice[j] = tmpIndice;
    }
  t = data[i - 1]; data[i - 1] = data[0]; data[0] = t;
  tmpIndice = indice[i - 1]; indice[i - 1] = indice[0]; indice[0] = tmpIndice;
  QuickSort(data, indice, i - 1);
  QuickSort(data + i, indice + i, N - i);
}
