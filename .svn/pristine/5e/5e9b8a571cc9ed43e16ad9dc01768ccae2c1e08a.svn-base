#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
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
  OptionManager Manager ("CountingZero" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleDoubleOption  ('e', "error", "rounding error (for floattig point vectors)", 1e-14);
  (*SystemGroup) +=  new BooleanOption  ('r', "rational", "input vectors are rational vectors");
  (*SystemGroup) +=  new BooleanOption  ('c', "complex", "input vectors are complex vectors");
  (*SystemGroup) +=  new BooleanOption  ('\n', "histogram", "compute the number of components whithin each order of magnitude");

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
  int MaxMagnitude = -500;
  int MinMagnitude = 499;
  int MagnitudeShift = 500;
  long* NbrComponentPerMagnitude = new long [1000];
  double* WeightPerMagnitude = new double [1000];
  for (int i = 0; i < 1000; ++i)
    {
      NbrComponentPerMagnitude[i] = 0;
      WeightPerMagnitude[i] = 0.0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "CountingZero requires an input file" << endl << "see man page for option syntax or type CountingZero -h" << endl;
      return -1;
    }

  long Count = 0l;
  long Dimension = 0l;
  if (Manager.GetBoolean("rational"))
    {
      LongRationalVector State;
      if (State.ReadVector (Manager.GetString("input-vector")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
	  return -1;      
	}
      Dimension = State.GetLargeVectorDimension();
      for (long i = 0; i < Dimension; ++i)
	if (State[i] == 0l)
	  ++Count;
     }
  else
    {
      if (Manager.GetBoolean("complex"))
	{
	  ComplexVector State;
	  if (State.ReadVector (Manager.GetString("input-vector")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
	      return -1;      
	    }
	  Dimension = State.GetLargeVectorDimension();
	  for (long i = 0; i < Dimension; ++i)
	    {
	      if (Norm(State[i]) < Error)
		{
		  ++Count;
		}
	      int TmpMagnitude = (int) (log10(Norm(State[i])));
	      if (Norm(State[i]) < 1.0)
		{
		  TmpMagnitude--;
		}
	      if (MaxMagnitude < TmpMagnitude)
		{
		  MaxMagnitude = TmpMagnitude;
		}
	      if (MinMagnitude > TmpMagnitude)
		{
		  MinMagnitude = TmpMagnitude;
		}
	      NbrComponentPerMagnitude[TmpMagnitude + MagnitudeShift]++;
	      WeightPerMagnitude[TmpMagnitude + MagnitudeShift] += SqrNorm(State[i]);	  
	    }
	}
      else
	{
	  RealVector State;
	  if (State.ReadVector (Manager.GetString("input-vector")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-vector") << endl;
	      return -1;      
	    }
	  Dimension = State.GetLargeVectorDimension();
	  for (long i = 0; i < Dimension; ++i)
	    {
	      if (fabs(State[i]) < Error)
		{
		  ++Count;
		}
	      int TmpMagnitude = (int) (log10(fabs(State[i])));
	      if (fabs(State[i]) < 1.0)
		{
		  TmpMagnitude--;
		}
	      if (MaxMagnitude < TmpMagnitude)
		{
		  MaxMagnitude = TmpMagnitude;
		}
	      if (MinMagnitude > TmpMagnitude)
		{
		  MinMagnitude = TmpMagnitude;
		}
	      NbrComponentPerMagnitude[TmpMagnitude + MagnitudeShift]++;
	      WeightPerMagnitude[TmpMagnitude + MagnitudeShift] += State[i] * State[i];	  
	    }
	}
    } 
  if (Manager.GetBoolean("histogram"))
    {
      cout << "# magnitude nbr_components sum_nbr_components weight sum_weight" << endl;
      double TmpSum = 0.0;
      long TmpSumNbr = 0l;
      for (int i = MaxMagnitude; i >= MinMagnitude; --i)
	{
	  TmpSum += WeightPerMagnitude[MagnitudeShift + i];
	  TmpSumNbr += NbrComponentPerMagnitude[MagnitudeShift + i];
	  cout << i << " " << NbrComponentPerMagnitude[MagnitudeShift + i] << " " << TmpSumNbr << " " 
	       << WeightPerMagnitude[MagnitudeShift + i] << " " << TmpSum << endl;
	}
    }
  else
    {
      cout << Count << " / " << Dimension << " (" << ((((double) Count) * 100.0) / ((double) Dimension)) << "%)" << endl; 
    }
  delete[] NbrComponentPerMagnitude;
  delete[] WeightPerMagnitude;

  return 0;
}
