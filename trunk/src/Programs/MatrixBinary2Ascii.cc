#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

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
  OptionManager Manager ("MatrixBinary2Ascii" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-matrix", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-matrix", "name of the file where the vector will be stored in ASCII mode");
  (*SystemGroup) += new BooleanOption  ('s', "std-output", "use standard output instead of an output file");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "indicate that the input vector has complex coefficients");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type MatrixBinary2Ascii -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-matrix") == 0)
    {
      cout << "MatrixBinary2Ascii requires an input file" << endl << "see man page for option syntax or type MatrixBinary2Ascii -h" << endl;
      return -1;
    }

  if ((Manager.GetString("output-matrix") == 0) && (Manager.GetBoolean("std-output") == false))
    {
      cout << "MatrixBinary2Ascii requires an output file" << endl << "see man page for option syntax or type MatrixBinary2Ascii -h" << endl;
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

      if (Manager.GetBoolean("std-output") == false)
	{
	  ofstream File;
	  File.open(Manager.GetString("output-matrix"), ios::out);
	  File.precision(14);
	  for (int i = 0; i < InputMatrix.GetNbrRow(); ++i)
	    for (int j = 0; j < InputMatrix.GetNbrColumn(); ++j)
	      File << i << " " << j << " " << InputMatrix[j][i] << endl;
	  File.close();
	}
      else
	{
	  cout.precision(14);
	  for (int i = 0; i < InputMatrix.GetNbrRow(); ++i)
	    for (int j = 0; j < InputMatrix.GetNbrColumn(); ++j)
	      cout << i << " " << j << " " << InputMatrix[j][i] << endl;
	}
      return 0;
    }

  ComplexMatrix InputMatrix;
  if (InputMatrix.ReadMatrix (Manager.GetString("input-matrix")) == false)
    {
      cout << "can't open matrix file " << Manager.GetString("input-matrix") << endl;
      return -1;      
    }
  
  if (Manager.GetBoolean("std-output") == false)
    {
      ofstream File;
      File.open(Manager.GetString("output-matrix"), ios::out);
      File.precision(14);
      for (int i = 0; i < InputMatrix.GetNbrRow(); ++i)
	for (int j = 0; j < InputMatrix.GetNbrColumn(); ++j)
	  File << i << " " << j << " " << InputMatrix[j][i] << endl;
      File.close();
    }
  else
    {
      cout.precision(14);
      for (int i = 0; i < InputMatrix.GetNbrRow(); ++i)
	for (int j = 0; j < InputMatrix.GetNbrColumn(); ++j)
	  cout << i << " " << j << " " << InputMatrix[j][i] << endl;
    }
  return 0;
}
