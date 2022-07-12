#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Options/OptionManager.h"
#include "Options/Options.h"

#include <iostream>
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

  // some running options and help
  OptionManager Manager ("MatrixExtractColumns" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-matrix", "name of the file containing the binary matrix");
  (*SystemGroup) += new SingleStringOption  ('o', "output-prefix", "prefix to use of each extracted vector", "vector");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "first-column", "index of the first column to extract", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-columns", "number of consecutive columns to extract", 1);
  (*SystemGroup) += new BooleanOption  ('c', "complex", "indicate that the input matrix has complex coefficients");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type MatrixExtractColumns -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-matrix") == 0)
    {
      cout << "MatrixExtractColumns requires an input file" << endl << "see man page for option syntax or type MatrixExtractColumns -h" << endl;
      return -1;
    }


  if (Manager.GetBoolean("complex") == false)
    {
      RealMatrix InputMatrix;
      if (InputMatrix.ReadMatrix (Manager.GetString("input-matrix")) == false)
	{
	  cout << "can't open matrix file " << Manager.GetString("input-matrix") << endl;
	  return -1;      
	}

      int FirstIndex = Manager.GetInteger("first-column");
      int LastIndex = FirstIndex + Manager.GetInteger("nbr-columns") - 1;
      for (; FirstIndex <= LastIndex; ++FirstIndex)
	{
	  char* TmpFileName = new char [strlen(Manager.GetString("output-prefix")) + 128];
	  sprintf (TmpFileName, "%s.%d.vec", Manager.GetString("output-prefix"), FirstIndex);
	  if (InputMatrix[FirstIndex].WriteVector(TmpFileName) == false)
	    {
	      cout << "error, can't write " << TmpFileName << endl;
	      return -1;
	    }
	}

      return 0;
    }

  ComplexMatrix InputMatrix;
  if (InputMatrix.ReadMatrix (Manager.GetString("input-matrix")) == false)
    {
      cout << "can't open matrix file " << Manager.GetString("input-matrix") << endl;
      return -1;      
    }
  
  int FirstIndex = Manager.GetInteger("first-column");
  int LastIndex = FirstIndex + Manager.GetInteger("nbr-columns") - 1;
  for (; FirstIndex <= LastIndex; ++FirstIndex)
    {
      char* TmpFileName = new char [strlen(Manager.GetString("output-prefix")) + 128];
      sprintf (TmpFileName, "%s.%d.vec", Manager.GetString("output-prefix"), FirstIndex);
      if (InputMatrix[FirstIndex].WriteVector(TmpFileName) == false)
	{
	  cout << "error, can't write " << TmpFileName << endl;
	  return -1;
	}
    }
  return 0;
}
