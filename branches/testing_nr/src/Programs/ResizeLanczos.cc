#include "Matrix/RealTriDiagonalSymmetricMatrix.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "new number of eigenvalues", 0);
  SingleStringOption InputFileOption ('\n', "input", "input file name", "lanczos.dat");
  SingleStringOption OutputFileOption ('o', "output", "output file name (same as input file if not defined)",0);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &InputFileOption;
  OptionList += &OutputFileOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  int NewNbrEigenvalue = NbrEigenvaluesOption.GetInteger();  
  char* InputFile = InputFileOption.GetString();
  char* OutputFile = OutputFileOption.GetString();
  if (OutputFile == 0)
    {
      cout << "The input file will be overwritten" << endl;
      OutputFile = InputFile;
    }

  ifstream File;
  File.open(InputFile, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the log file: " << InputFile << endl;
      return false;
    }

  int Index;
  double PreviousLastWantedEigenvalue;
  double EigenvaluePrecision;
  int NbrEigenvalue;

  File.read((char*) (&Index), sizeof(int));
  File.read((char*) (&PreviousLastWantedEigenvalue), sizeof(double));
  File.read((char*) (&EigenvaluePrecision), sizeof(double));
  File.read((char*) (&NbrEigenvalue), sizeof(int));
  int TmpDimension;
  File.read((char*) (&TmpDimension), sizeof(int));
  cout << "Number of done iterations: " << TmpDimension << endl;
  if (TmpDimension < NewNbrEigenvalue)
    {
      cout << "The number of wanted eigenvalues must be smaller than the number of done iterations. Please reduce the number of wanted eigenvalues." << endl;
      File.close();
      return 1;
    }  
  cout << "Previous last wanted eigenvalue is: " << PreviousLastWantedEigenvalue << endl;
  cout << "Old number of wanted eigenvalues is: " << NbrEigenvalue << endl;
  cout << "New number of wanted eigenvalues is: " << NewNbrEigenvalue << endl;

  RealTriDiagonalSymmetricMatrix TridiagonalizedMatrix (TmpDimension, true);
  for (int i = 0; i <= (Index + 1); ++i)
    {
      File.read((char*) (&TridiagonalizedMatrix.DiagonalElement(i)), sizeof(double));
    }
  for (int i = 0; i <= Index; ++i)
    {
      File.read((char*) (&TridiagonalizedMatrix.UpperDiagonalElement(i)), sizeof(double));
    }
  File.close();  

  double* PreviousWantedEigenvalues = new double [NewNbrEigenvalue];

  ofstream File2;
  File2.open(OutputFile, ios::binary | ios::out);
  File2.write((char*) (&Index), sizeof(int));
  File2.write((char*) (&PreviousLastWantedEigenvalue), sizeof(double));
  File2.write((char*) (&EigenvaluePrecision), sizeof(double));
  File2.write((char*) (&NewNbrEigenvalue), sizeof(int));
  TmpDimension = TridiagonalizedMatrix.GetNbrRow();
  File2.write((char*) (&TmpDimension), sizeof(int));
  for (int i = 0; i <= (Index + 1); ++i)    
    {    
      File2.write((char*) (&TridiagonalizedMatrix.DiagonalElement(i)), sizeof(double));
    }
  for (int i = 0; i <= Index; ++i)
    {
      File2.write((char*) (&TridiagonalizedMatrix.UpperDiagonalElement(i)), sizeof(double));
    }
  for (int i = 0; i < NewNbrEigenvalue; ++i)
    {
      File2.write((char*) (&PreviousWantedEigenvalues[i]), sizeof(double));
    }
  File2.close();
  return 0;  
}
