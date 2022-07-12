#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

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
  OptionManager Manager ("VectorRealImaginary2Complex" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the input vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the output has to be stored");
  (*SystemGroup) += new BooleanOption  ('\n', "set-imag", "set input vector to be the imaginary part, otherwise it will be real", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorRealImaginary2Complex -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "no input file" << endl;
    }
  if (Manager.GetString("output-vector") == 0)
    {
      cout << "no output file" << endl;
    }
  bool SetImag = Manager.GetBoolean("set-imag"); 

  RealVector InputVector;
  if (InputVector.ReadVector(Manager.GetString("input-vector")) == false)
        {
	  return -1;
	}

  ComplexVector OutputVector (InputVector.GetVectorDimension(), true);
  for (int i = 0; i < InputVector.GetVectorDimension(); ++i)
    {
      if (SetImag == false)
	{
	  OutputVector.Re(i) = InputVector[i];
	  OutputVector.Im(i) = 0;
	}
      else
	{
	  OutputVector.Re(i) = 0;
	  OutputVector.Im(i) = InputVector[i];
	}
    }

 OutputVector.WriteVector(Manager.GetString("output-vector"));
 return 0;
}
