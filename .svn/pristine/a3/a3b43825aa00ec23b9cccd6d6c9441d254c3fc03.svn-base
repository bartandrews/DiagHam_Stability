#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

int M = 10;//x
int N = 10;//y
int P = 14;//z

char * FileName = "eigenvector.2";

int Number = 5;

using namespace std;

int main(int argc, char** argv)
{
  // some running options and help 
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption XCell('X', "Xcells", "number of the cells in X direction", 10);
  SingleIntegerOption YCell('Y', "Ycells", "number of the cells in Y direction", 10);
  SingleIntegerOption ZCell('Z', "Zcells", "number of the cells in Z direction", 10);
  SingleIntegerOption Division('D', "division", "number of the sub-divisions in each direction", 5);
  SingleIntegerOption Orientation('o', "orientation", "(1) X | (2) Y | (3) Z", 1);
  SingleStringOption Output('r', "output", "name of the output file", "default_output.txt");

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &XCell;
  OptionList += &YCell;
  OptionList += &ZCell;
  OptionList += &Division;
  OptionList += &Orientation;
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
  int M1, N1, P1;
  FileName = InputFile.GetString();
  M = XCell.GetInteger(); M1 = M;
  N = YCell.GetInteger(); N1 = N;
  P = ZCell.GetInteger(); P1 = P;
  Number = Division.GetInteger();
  int choice = Orientation.GetInteger();
  char * out = Output.GetString();

  int exchange;
  switch(choice)
    {
    case 1: cout << "The output format is: X for lines and Y for columns\n"; break;
    case 2: cout << "The output format is: Z for lines and Y for columns\n";
      exchange = M; M = N; N = exchange; break;
    case 3: cout << "The output format is: X for lines and Z for columns\n";
      exchange = M; M = P; P = exchange; break;
    default: cout << "Bad choice\n"; exit(1);
    }  
  
  int NX = M * Number;
  
  double ** C; double*** f;
  f = new double** [M];
  for (int i = 0; i < M; ++i)
    {
      f[i] = new double * [M];
      for (int j = 0; j < M ; ++j)
	f[i][j] = new double [NX + 1];
    }

  //read data from the input file
  double *** Coef;
  Coef = new double ** [M];
  for (int i = 0; i < M; ++i)
    {
      Coef[i] = new double * [N];
      for (int j = 0; j < N ; ++j)
	{
	  Coef[i][j] = new double [P];
	}
    }
  
  ifstream Coefficients(FileName,ios::in|ios::binary);
  for (int i = 0; i < M1; ++i)
    for (int j = 0; j < N1; ++j)
      for (int k = 0; k < P1; ++k)
	{
	  switch(choice)
	    {
	    case 1: Coefficients >> Coef[i][j][k]; break;
	    case 2: Coefficients >> Coef[j][i][k]; break;
	    default: Coefficients >> Coef[k][j][i]; break;  
	    }
	}
  Coefficients.close();
  
  C = new double * [M];
  for (int i = 0; i < M; ++i)
    {
      C[i] = new double [M];
      for (int j = 0; j < i; ++j)
	{
	  C[i][j] = 0.0;
	  for (int k = 0; k < N; ++k)
	    for (int l = 0; l < P; ++l)
	      {
		C[i][j] += Coef[i][k][l]*Coef[j][k][l];
	      }
	  C[j][i] = C[i][j];
	}
      C[i][i] = 0.0;
      for (int k = 0; k < N; ++k)
	for (int l = 0; l < P; ++l)
	  {
	    C[i][i] += Coef[i][k][l]*Coef[i][k][l];
	  }     
    }
  

  double ** Sin = new double * [M];
  for (int i = 0 ; i < M; ++i)
    Sin[i] = new double [NX + 1];

  for (int i = 1 ; i <= M; ++i)
    for (int j = 0; j <= NX; ++j)
      Sin[i - 1][j] = sin((M_PI * i * j) / NX);

  for (int i = 0;  i < M; ++i)
    {
      for (int j = 0; j < i; ++j)
	for (int k = 0; k <= NX; ++k)
	  {
	    f[i][j][k] = Sin[i][k] * Sin[j][k];
	    f[j][i][k] = f[i][j][k];
	    //cout << f[i][j][k] << endl;
	  }
      for (int k = 0; k <= NX; ++k)
	f[i][i][k] = Sin[i][k] * Sin[i][k];
    }

  ofstream Visual(out,ios::out);
  double tmp = 0.0;
  Visual << "#x\tFunction" << endl;
  for (int i = 0; i <= NX; ++i)
    {
      Visual << double(i) / Number << "\t";
      tmp = 0.0;
      for (int m1 = 0; m1 < M; ++m1)
	{
	  for (int m2 = 0; m2 < m1; ++m2)
	  {
	    tmp += (2 * f[m1][m2][i] * C[m1][m2]);
	  }
	  tmp += f[m1][m1][i] * C[m1][m1];
	}
      Visual << tmp;
      Visual << '\n';
    }
  
  
  Visual.close();
  return -1;
}
