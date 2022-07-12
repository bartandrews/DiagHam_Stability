#include "Vector/RealVector.h"

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
  OptionManager Manager ("VectorBinary2Ascii" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the vector will be stored in ASCII mode");
  (*SystemGroup) += new BooleanOption  ('\n', "add-index", "write ascii vector in a two-column formatted output (first column for the component index, second for the component value)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-range", "display vector starting from a given component (is negative, start couting from the last component) from ", 0l);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-range", "display vector up to a given component (0 if up to the end)", 0l);
  (*SystemGroup) += new BooleanOption  ('s', "std-output", "use standard output instead of an output file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorBinary2Ascii -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-vector"])->GetString() == 0)
    {
      cout << "VectorBinary2Ascii requires an input file" << endl << "see man page for option syntax or type VectorBinary2Ascii -h" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["input-vector"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-vector"])->GetString() << endl;
      return -1;      
    }
  if ((((SingleStringOption*) Manager["output-vector"])->GetString() == 0) && (Manager.GetBoolean("std-output") == false))
    {
      cout << "VectorBinary2Ascii requires an output file" << endl << "see man page for option syntax or type VectorBinary2Ascii -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("std-output") == false)
    {
      ofstream File;
      File.open(((SingleStringOption*) Manager["output-vector"])->GetString(), ios::out);
      File.precision(14);
      long MinValue = 0l;
      if ((Manager.GetInteger("min-range") > 0l) && (Manager.GetInteger("min-range") < State.GetLargeVectorDimension()))
	MinValue = Manager.GetInteger("min-range");
      long MaxValue = State.GetLargeVectorDimension(); 
      if ((Manager.GetInteger("max-range") < State.GetLargeVectorDimension()) && (Manager.GetInteger("max-range") > MinValue))
	MaxValue = Manager.GetInteger("max-range");
      
      if (((BooleanOption*) Manager["add-index"])->GetBoolean() == true)
	for (long i = MinValue; i < MaxValue; ++i)
	  File << i << " " << State[i] << endl;
      else
	for (long i = MinValue; i < MaxValue; ++i)
	  File << State[i] << endl;
      File.close();
    }
  else
    {
      cout.precision(14);
      long MinValue = 0l;
      if (Manager.GetInteger("min-range") > 0l)
	{
	  if (Manager.GetInteger("min-range") < State.GetLargeVectorDimension())
	    MinValue = Manager.GetInteger("min-range");      
	}
      else
	{
	  long Tmp = State.GetLargeVectorDimension() + Manager.GetInteger("min-range");
	  if ((Tmp >= 0) && (Tmp < State.GetLargeVectorDimension()))
	    MinValue = Tmp;
	}
      long MaxValue = State.GetLargeVectorDimension();
      if ((Manager.GetInteger("max-range") < State.GetLargeVectorDimension()) && (Manager.GetInteger("max-range") > MinValue))
        MaxValue = Manager.GetInteger("max-range");
      if (Manager.GetBoolean("add-index") == true)
        for (long i = MinValue; i < MaxValue; ++i)
          cout << i << " " << State[i] << endl;
      else
        for (long i = MinValue; i < MaxValue; ++i)
          cout << State[i] << endl;
    }
  return 0;
}
