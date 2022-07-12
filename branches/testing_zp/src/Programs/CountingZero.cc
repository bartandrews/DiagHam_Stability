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
  OptionManager Manager ("CountingZero" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "input-vector", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleDoubleOption  ('e', "error", "rounding error", 1e-14);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type CountingZero -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  double Error = Manager.GetDouble("error");

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "CountingZero requires an input file" << endl << "see man page for option syntax or type CountingZero -h" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (Manager.GetString("input-vector")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
      return -1;      
    }

  long Count = 0l;

  for (long i = 0; i < State.GetLargeVectorDimension(); ++i)
    if (fabs(State[i]) < Error)
      ++Count;

  cout << Count << " / " << State.GetLargeVectorDimension() << " (" << ((((double) Count) * 100.0) / ((double) State.GetLargeVectorDimension())) << "%)" << endl; 

  return 0;
}
