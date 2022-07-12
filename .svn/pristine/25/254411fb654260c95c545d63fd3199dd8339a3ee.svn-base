
#include "MathTools/IntegerAlgebraTools.h"
#include "MathTools/LongRational.h"

#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// evalute an element of the S matrix
//
// pValue = p value of the M(p,q) minimal model
// qValue = q value of the M(p,q) minimal model
//
// return value = S matrix element
double SMatrixElement (long pValue, long qValue, long rTransformedField, long sTransformedField,  long rField, long sField);

int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHETorusMinimalModelQuantumDimensions" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "p-value", "p value of the M(p,q) minimal model", 4);
  (*SystemGroup) += new SingleIntegerOption  ('q', "q-value", "q value of the M(p,q) minimal model", 3);
  (*OutputGroup) += new BooleanOption  ('\n', "latex", "output the result in a latex ready format");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMinimalModelQuantumDimensions -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  long PValue = Manager.GetInteger("p-value"); 
  long QValue = Manager.GetInteger("q-value"); 

  int SMatrixDimension = (QValue - 1) * (PValue -1);
  RealMatrix SMatrix (SMatrixDimension, SMatrixDimension, true);
  for (long r1 = 1l; r1 < QValue; ++r1)
    {
      for (long s1 = 1l; s1 < PValue; ++s1)
	{
	  for (long r2 = 1l; r2 < QValue; ++r2)
	    {
	      for (long s2 = 1l; s2 < PValue; ++s2)
		{
		  SMatrix.SetMatrixElement((r1 - 1) * (PValue - 1) + (s1 - 1), (r2 - 1) * (PValue - 1) + (s2 - 1), SMatrixElement(PValue, QValue, r1, s1, r2, s2));
		}
	    }
	}
    }

  int PerronFroebeniusIndex = -1;
  for (int i = 0; (i < SMatrixDimension) && (PerronFroebeniusIndex < 0); ++i)
    {
      int j = 0;
      while ((j < SMatrixDimension) && (SMatrix[j][i] > 0.0))
	++j;
      if (j == SMatrixDimension)
	PerronFroebeniusIndex = i;
    }
  int PerronFroebeniusRIndex = 1 + (PerronFroebeniusIndex / (PValue - 1));
  int PerronFroebeniusSIndex = 1 + (PerronFroebeniusIndex % (PValue - 1));
  cout << "Highest weight in r=" << PerronFroebeniusRIndex << ",s=" << PerronFroebeniusSIndex << endl;
  double TotalQuantumDimension  = 1.0 / SMatrixElement(PValue, QValue, PerronFroebeniusRIndex, PerronFroebeniusSIndex, PerronFroebeniusRIndex, PerronFroebeniusSIndex);
  LongRational Shift ((PValue - QValue) * (PValue - QValue), 4l * PValue * QValue);
  for (long r = 1l; r < QValue; ++r)
    {
      for (long s = 1l; s < PValue; ++s)
	{
	  LongRational Tmp (((PValue * r) - (QValue * s)) * ((PValue * r) - (QValue * s)), 4l * PValue * QValue);
	  Tmp -= Shift;	
	  cout << r << "," << s << " : h=" << Tmp << " d=" << (TotalQuantumDimension * SMatrixElement(PValue, QValue, PerronFroebeniusRIndex, PerronFroebeniusSIndex, r, s)) << endl;
	}
    }
  cout << "total quantum dimension  = " << fabs(TotalQuantumDimension) << endl;
  cout << "S matrix : " << endl;
  for (long r1 = 1l; r1 < QValue; ++r1)
    {
      for (long s1 = 1l; s1 < PValue; ++s1)
	{
	  for (long r2 = 1l; r2 < QValue; ++r2)
	    {
	      for (long s2 = 1l; s2 < PValue; ++s2)
		{
		  cout << r1 << "," << s1 << " ; " << r2 << "," << s2 << " : " 
		       << (SMatrixElement(PValue, QValue, r1, s1, r2, s2) * TotalQuantumDimension) << endl;
		}
	    }
	}
    }
  return 0;

}

// evalute an element of the S matrix
//
// pValue = p value of the M(p,q) minimal model
// qValue = q value of the M(p,q) minimal model
//
// return value = S matrix element

double SMatrixElement (long pValue, long qValue, long rTransformedField, long sTransformedField,  long rField, long sField)
{
  long Tmp = rTransformedField * sField + sTransformedField * rField + 1;
  double SMatrixElement = 1.0;
  if ((Tmp & 1l) != 0l)
    SMatrixElement = 1.0;
  SMatrixElement *= sqrt (8.0 / (pValue * qValue)) * sin ((M_PI * (pValue * rTransformedField * rField)) / (qValue)) * sin ((M_PI * (qValue * sTransformedField * sField)) / (pValue));
  return SMatrixElement;
}
