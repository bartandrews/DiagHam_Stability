#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

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
  OptionManager Manager ("VectorRational2Double" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the vector will be stored floating point mode");
#ifdef __GMP__
  (*SystemGroup) += new BooleanOption  ('\n', "use-gmp", "use arbitrary precision integers instead of fixed precision integers in rational mode");
#else
  (*SystemGroup) += new BooleanOption  ('\n', "use-longlong", "use 128bit(64bits) integers instead of 64bits(32bits) integers in rational mode");
#endif

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorRational2Double -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "VectorRational2Double requires an input file" << endl << "see man page for option syntax or type VectorRational2Double -h" << endl;
      return -1;
    }

#ifdef __GMP__
  if (Manager.GetBoolean("use-gmp") == false)
#else
    if (Manager.GetBoolean("use-longlong") == false)	
#endif
      {
	RationalVector State;
	if (State.ReadVector (Manager.GetString("input-vector")) == false)
	  {
	    cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
	    return -1;      
	  }
	if (Manager.GetString("output-vector") == 0)
	  {
	    cout << "VectorRational2Double requires an output file" << endl << "see man page for option syntax or type VectorRational2Double -h" << endl;
	    return -1;
	  }
	RealVector OutputState (State.GetLargeVectorDimension());	
	for (long i = 0l; i < State.GetLargeVectorDimension(); ++i)
	  OutputState[i] = State[i].GetNumericalValue();
	OutputState.WriteVector(Manager.GetString("output-vector"));
      }
    else
      {
	LongRationalVector State;
	if (State.ReadVector (Manager.GetString("input-vector")) == false)
	  {
	    cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
	    return -1;      
	  }
	if (Manager.GetString("output-vector") == 0)
	  {
	    cout << "VectorRational2Double requires an output file" << endl << "see man page for option syntax or type VectorRational2Double -h" << endl;
	    return -1;
	  }
	RealVector OutputState (State.GetLargeVectorDimension());	
	for (long i = 0l; i < State.GetLargeVectorDimension(); ++i)
	  OutputState[i] = State[i].GetNumericalValue();
	OutputState.WriteVector(Manager.GetString("output-vector"));
      }
  return 0;
}
