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
  OptionManager Manager ("GenericOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of the vector files obtained using exact diagonalization");
  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "conjugate", "Conjugate the second (complex) number");
  (*SystemGroup) += new BooleanOption  ('\n', "discard-sign", "compute sum_i |v1_i * v2_i| instead of sum_i v1_i * v2_i");
  (*SystemGroup) += new BooleanOption  ('x', "no-cross", "calculate only overlap of 1st vector with all others");
  (*SystemGroup) += new BooleanOption  ('\n', "quiet", "discard any output except the overlaps");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool QuietFlag = Manager.GetBoolean("quiet");

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);

  if (NbrVectors<2)
    {
      cout << "At least two vector files are required!"<<endl;
      exit(1);
    }

  if (QuietFlag == false)
    for (int i=0; i<NbrVectors; ++i)    
      cout << "File "<<i<<"  "<<VectorFiles[i]<<endl;

  Complex sp=0.0;

  int MaxVectors=(Manager.GetBoolean("no-cross")?1:NbrVectors);

  bool HaveComplex=Manager.GetBoolean("complex");
  
  if (HaveComplex)
    {      
      ComplexVector State1, State2;
      for (int i=0; i<MaxVectors; ++i)
	{
	  if (State1.ReadVector (VectorFiles[i]) == false)
	    {
	      cout << "can't open vector file " << VectorFiles[i] << endl;
	      return -1;      
	    }
	  for (int j=i+1; j<NbrVectors; ++j)
	    {	      
	      if (State2.ReadVector (VectorFiles[j]) == false)
		{
		  cout << "can't open vector file " << VectorFiles[j] << endl;
		  return -1;      
		}
	      if (State1.GetVectorDimension() != State2.GetVectorDimension() )
		{
		  cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
		  return -2;
		}
	      sp=0.0;
	      if (Manager.GetBoolean("discard-sign"))
		for (int i=0; i<State1.GetVectorDimension(); ++i)
		  sp+=Norm(State1[i]*State2[i]);
	      else
		if (Manager.GetBoolean("conjugate"))
		  for (int i=0; i<State1.GetVectorDimension(); ++i)
		    sp+= State1[i]*State2[i];
		else
		  for (int i=0; i<State1.GetVectorDimension(); ++i)
		    sp+= Conj(State1[i])*State2[i];
	      if (QuietFlag == false)
		cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
	      else
		cout << SqrNorm(sp) << endl;
	    }
	}
    }
  else // real vectors
    {
      RealVector State1, State2;
      for (int i=0; i<MaxVectors; ++i)
	{
	  if (State1.ReadVector (VectorFiles[i]) == false)
	    {
	      cout << "can't open vector file " << VectorFiles[i] << endl;
	      return -1;      
	    }
	  for (int j=i+1; j<NbrVectors; ++j)
	    {	      
	      if (State2.ReadVector (VectorFiles[j]) == false)
		{
		  cout << "can't open vector file " << VectorFiles[j] << endl;
		  return -1;      
		}
	      if (State1.GetVectorDimension() != State2.GetVectorDimension() )
		{
		  cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
		  return -2;
		}
	      
	      sp=0.0;
	      if (Manager.GetBoolean("discard-sign"))
		for (int i=0; i<State1.GetVectorDimension(); ++i)
		  sp+=fabs(State1[i]*State2[i]);
	      else
		for (int i=0; i<State1.GetVectorDimension(); ++i)
		  sp+= State1[i]*State2[i];
	      
	      if (QuietFlag == false)
		cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
	      else
		cout << SqrNorm(sp) << endl;
	    }
	}
    }
}
