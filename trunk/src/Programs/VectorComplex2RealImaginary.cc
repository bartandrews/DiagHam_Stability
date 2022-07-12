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
  OptionManager Manager ("VectorComple2RealImaginary" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the input vector (assumed to be complex)");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-real", "name of the file where the real part has to be stored");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-imag", "name of the file where the imaginary part has to be stored");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorComplex2RealImaginary -h" << endl;
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
  if (Manager.GetString("output-real") == 0)
    {
      cout << "no output file" << endl;
    }
  if (Manager.GetString("output-imag") == 0)
    {
      cout << "no output file" << endl;
    }

  ComplexVector InputVector;
  if (InputVector.ReadVectorTest(Manager.GetString("input-vector")))  
    {
      if (InputVector.ReadVector(Manager.GetString("input-vector")) == false)
	{
	  return -1;
	}
      RealVector OutputVectorRe (InputVector.GetVectorDimension(), true);
      RealVector OutputVectorIm (InputVector.GetVectorDimension(), true); 
      Complex Comp;
      for (int i = 0; i < InputVector.GetVectorDimension(); ++i)
        {
         Comp = InputVector[i];
         OutputVectorRe[i] = Comp.Re;
         OutputVectorIm[i] = Comp.Im;
        }

      OutputVectorRe.WriteVector(Manager.GetString("output-real"));
      OutputVectorIm.WriteVector(Manager.GetString("output-imag"));
    }
 return 0;
}
