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
  OptionManager Manager ("ZeroingVector" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the vector will be stored in ASCII mode");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-range", "zero vector starting from a given component (is negative, start couting from the last component) from ", 0l);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-range", "zero vector up to a given component (0 if up to the end)", 0l);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ZeroingVector -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "ZeroingVector requires an input file" << endl << "see man page for option syntax or type ZeroingVector -h" << endl;
      return -1;
    }
  RealVector State;
  if (Manager.GetString("output-vector") == 0)
    {
      cout << "ZeroingVector requires an output file" << endl << "see man page for option syntax or type ZeroingVector -h" << endl;
      return -1;
    }
  if (State.ReadVector (Manager.GetString("input-vector")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
      return -1;      
    }

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
  for (long i = MinValue; i < MaxValue; ++i)
    State[i] = 0.0;

  State.WriteVector(Manager.GetString("output-vector"));

  return 0;
}
