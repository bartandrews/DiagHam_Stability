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
//using namespace ;

int main(int argc, char** argv)
{
  // some running options and help 
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption XCell('X', "Xcells", "number of the cells in X direction", 10);
  SingleIntegerOption YCell('Y', "Ycells", "number of the cells in Y direction", 10);
  SingleIntegerOption ZCell('Z', "Zcells", "number of the cells in Z direction", 10);
  SingleIntegerOption Division('D', "division", "number of the sub-divisions in each direction", 5);
  SingleIntegerOption Orientation('o', "orientation", "(1) X-Y | (2) Z-Y | (3) X-Z", 1);
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
      exchange = M; M = P; P = exchange; break;
    case 3: cout << "The output format is: X for lines and Z for columns\n";
      exchange = N; N = P; P = exchange; break;
    default: cout << "Bad choice\n"; exit(1);
    }  
  
  int NX = M * Number;
  int NY = N * Number;
  
  double **** C; double*** f; double*** g;
  f = new double** [M];
  for (int i = 0; i < M; ++i)
    {
      f[i] = new double * [M];
      for (int j = 0; j < M ; ++j)
	f[i][j] = new double [NX - 1];
    }
  g = new double** [N];
  for (int i = 0 ; i < N; ++i)
    {
      g[i] = new double * [N];
      for (int j = 0; j < N; ++j)
	g[i][j] = new double [NY - 1];
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
	    case 2: Coefficients >> Coef[k][j][i]; break;
	    default: Coefficients >> Coef[i][k][j]; break;  
	    }
	}
  
  Coefficients.close();
  C = new double *** [M];
  for (int i = 0; i < M; ++i)
    {
      C[i] = new double ** [M];
      for (int j = 0; j < M; ++j)
	{
	  C[i][j] = new double * [N];
	  for (int k = 0; k < N; ++k)
	    {
	      C[i][j][k] = new double [N];
	      for (int l = 0; l < N; ++l)
		{
		  C[i][j][k][l] = 0.0;
		  for (int p = 0 ; p < P; ++p){
		    // cout << "Before " << i << " " << j << " " << k << " " << l << " " << p << endl;
		    C[i][j][k][l] += Coef[i][k][p] * Coef[j][l][p];
		    // cout << "After " << i << " " << j << " " << k << " " << l << " " << p << endl;
		  }
		}
	    }   
	}
    }
  
  double ** SinX = new double * [M];
  double ** SinY = new double * [N];
  for (int i = 0 ; i < M; ++i)
    SinX[i] = new double [NX - 1];
  for (int i = 0 ; i < N; ++i)
    SinY[i] = new double [NY - 1];  

  for (int i = 1 ; i <= M; ++i)
    for (int j = 1; j < NX; ++j)
      SinX[i - 1][j - 1] = sin((M_PI * i * j) / NX);
  
  for (int i = 1 ; i <= N; ++i)
    for (int j = 1; j < NY; ++j)
      SinY[i - 1][j - 1] = sin((M_PI * i * j) / NY);

  for (int i = 0;  i < M; ++i)
    for (int j = 0; j < M; ++j)
      for (int k = 0; k < NX - 1; ++k)
	{
	  f[i][j][k] = SinX[i][k] * SinX[j][k];
	}
  
  for (int i = 0;  i < N; ++i)
    for (int j = 0; j < N; ++j)
      for (int k = 0; k < NY - 1; ++k)
	g[i][j][k] = SinY[i][k] * SinY[j][k];

  ofstream Visual(out,ios::out);

  Visual << "#x\ty\tFunction" << endl;
  
  for (int i = 0; i < NX + 1; ++i)
    Visual << double(i) / Number << '\t' <<  0.0 << '\t' << 0.0 << '\n';
    //Visual << 0.0 << " ";
  //Visual << ' \n';

  double tmp = 0.0;
  double** Temp = new double* [M];
  for (int i = 0; i < M; ++i)
    Temp[i] = new double [M];
  
  for (int j = 0; j < NY - 1; ++j)
    {
      Visual << 0.0 << '\t' << double (j) / Number << '\t' << 0.0 << '\n';
      //Visual << 0.0 << " ";
      
      for (int m1 = 0; m1 < M; ++m1)
	for (int m2 = 0; m2 < M; ++m2)
	  {
	    Temp[m1][m2] = 0.0;
	    for (int n1 = 0; n1 < N; ++n1)
	      for (int n2 = 0; n2 < N; ++n2)
		Temp[m1][m2] += g[n1][n2][j] * C[m1][m2][n1][n2];
	  }
      for (int i = 0; i < NX - 1; ++i)
	{
	  tmp = 0.0;
	  for (int m1 = 0; m1 < M; ++m1)
	    for (int m2 = 0; m2 < M; ++m2)
	      tmp += f[m1][m2][i] * Temp[m1][m2];
	  Visual << double(i) / Number << '\t' << double(j) / Number << '\t' << tmp << '\n';
	  //Visual << tmp << " ";
	}
      Visual << M << '\t' << double(j) / Number << '\t' << 0.0 << '\n';
      //Visual << 0.0 << '\n';
    }
  
  for (int i = 0; i < M; ++i)
    delete[] Temp[i];
  delete[] Temp;

  for (int i = 0; i < NX + 1; ++i)
    Visual << double(i)/Number  << '\t' <<  N << '\t' << 0.0 << '\n';
    //Visual << 0.0 << " ";
  Visual << '\n';
  
  Visual.close();
  char x = 'X', y = 'Y';
  if (choice == 2)
    {
      x = 'Z';
      y = 'Y';
    }
  if (choice == 3)
    y = 'Z';
  
  ofstream Maple("Maple.mws",ios::out|ios::app);
  
  Maple << ">with(plots);\n";
  Maple << ">with(linalg);\n";
  Maple << ">flux:=readdata(\"" << out << "\",float," << NX + 1 << ");\n";
  Maple << ">Function:=Matrix(" << NY + 1 << "," << NX + 1 << ",flux);\n";
  Maple << ">matrixplot(Function,axes=frame,labels=[\"" << y << "\",\"" << x <<  "\",\"\"]);\n";
  Maple.close();
  
  return -1;
}
