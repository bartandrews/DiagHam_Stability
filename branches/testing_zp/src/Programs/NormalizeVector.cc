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
  OptionManager Manager ("NormalizeVector" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "vector", "vector to normalize");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type NormalizeVector -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  RealVector State;
  if (State.ReadVector (Manager.GetString("vector")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("vector") << endl;
      return -1;      
    }
  State /= State.Norm();
  if (State.WriteVector (Manager.GetString("vector")) == false)
    {
      cout << "can't overwrite vector file " << Manager.GetString("vector") << endl;
      return -1;
    }
  return 0;
}
