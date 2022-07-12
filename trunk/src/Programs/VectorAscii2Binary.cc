#include "MathTools/Complex.h"

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

#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("VectorAscii2Binary" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the ASCII vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the vector will be stored in binary");
  (*SystemGroup) += new BooleanOption  ('r', "rational", "indicate that the input vector has rational coefficients");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "indicate that the input vector has complexl coefficients");
#ifdef __GMP__
  (*SystemGroup) += new BooleanOption  ('\n', "use-gmp", "use arbitrary precision integers instead of fixed precision integers in rational mode");
#else
  (*SystemGroup) += new BooleanOption  ('\n', "use-longlong", "use 128bit(64bits) integers instead of 64bits(32bits) integers in rational mode");
#endif

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorAscii2Binary -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "VectorAscii2Binary requires an input file" << endl << "see man page for option syntax or type VectorAscii2Binary -h" << endl;
      return -1;
    }

  if ((Manager.GetString("output-vector") == 0) && (Manager.GetBoolean("std-output") == false))
    {
      cout << "VectorAscii2Binary requires an output file" << endl << "see man page for option syntax or type VectorAscii2Binary -h" << endl;
      return -1;
    }

  MultiColumnASCIIFile AsciiVector;
  if (AsciiVector.Parse(Manager.GetString("input-vector")) == false)
    {
      AsciiVector.DumpErrors(cout) << endl;
      return -1;
    }

  if (Manager.GetBoolean("rational") == false && Manager.GetBoolean("complex")  == false )
    {
      double* TmpData = AsciiVector.GetAsDoubleArray(0);
      if (TmpData == 0)
	{
	  AsciiVector.DumpErrors(cout) << endl;
	  return -1;     
	}
      RealVector BinaryVector(TmpData, AsciiVector.GetNbrLines());
      BinaryVector.WriteVector(Manager.GetString("output-vector"));
    }
  else if (Manager.GetBoolean("rational") == false && Manager.GetBoolean("complex")  == true )
    {
      Complex* TmpData = AsciiVector.GetAsComplexArray(0);
      if (TmpData == 0)
	{
	  AsciiVector.DumpErrors(cout) << endl;
	  return -1;     
	}
      ComplexVector BinaryVector(TmpData, AsciiVector.GetNbrLines());
      BinaryVector.WriteVector(Manager.GetString("output-vector"));      
    }
  else
    {
#ifdef __GMP__
      if (Manager.GetBoolean("use-gmp") == false)
#else
      if (Manager.GetBoolean("use-longlong") == false)	
#endif
	{
	  long* TmpNumerators = AsciiVector.GetAsLongArray(0);
	  long* TmpDenominators = AsciiVector.GetAsLongArray(1);
	  if ((TmpNumerators == 0) || (TmpDenominators == 0))
	    {
	      AsciiVector.DumpErrors(cout) << endl;
	      return -1;     
	    }
	  RationalVector BinaryVector(TmpNumerators, TmpDenominators, AsciiVector.GetNbrLines());
	  BinaryVector.WriteVector(Manager.GetString("output-vector"));
	}
      else 
	{
	  LongRational* TmpCoefficients = AsciiVector.GetAsLongRationalArray(0);
	  if (TmpCoefficients == 0)
	    {
	      AsciiVector.DumpErrors(cout) << endl;
	      return -1;     
	    }
	  LongRationalVector BinaryVector(TmpCoefficients, AsciiVector.GetNbrLines());
	  BinaryVector.WriteVector(Manager.GetString("output-vector"));
	}
    }

 
  return 0;
}
