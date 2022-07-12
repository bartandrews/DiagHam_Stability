#include "Vector/RealVector.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "GeneralTools/ConfigurationParser.h"

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
  OptionManager Manager ("QHEGenericCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");


  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleStringOption  ('e', "expansion-file", "file containing expansion coefficients", 0);
  (*SystemGroup) += new SingleStringOption  ('i', "interaction-name", "name of the interaction (used for output file name)", "generic");
  (*SystemGroup) += new SingleStringOption ('a', "add-filename", "add a string with additional informations to the output file name(just before the .dat extension)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 200);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEGenericCorrelation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();

  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  double* ExpansionCoeff = 0;
  if (((SingleStringOption*) Manager["expansion-file"])->GetString() == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(((SingleStringOption*) Manager["expansion-file"])->GetString()) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      int TmpNbrCoefficients;
      if (InteractionDefinition.GetAsDoubleArray("ExpansionCoeff", ' ', ExpansionCoeff, TmpNbrCoefficients) == false)
	{
	  cout << "Expansion coefficients are not defined or as a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;
	}
      if (TmpNbrCoefficients != (LzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials" << endl;
	  return -1;	  
	}
    }

  char* OutputNameCorr = new char [256 + strlen (((SingleStringOption*) Manager["interaction-name"])->GetString())];
  if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
    {
      sprintf (OutputNameCorr, "corr_%s_2s_%d.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), LzMax);
    }
  else
    {
      sprintf (OutputNameCorr, "fermions_%s_2s_%d_%s.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), LzMax,
	       ((SingleStringOption*) Manager["add-filename"])->GetString());
    }

  ParticleOnSphereFunctionBasis Basis(LzMax);

  double Sum;
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);

  ofstream File;
  File.precision(14);
  File.open(OutputNameCorr, ios::binary | ios::out);
  double Factor2;
  if (((BooleanOption*) Manager["radians"])->GetBoolean() == true)
    Factor2 = 1.0;
  else
    Factor2 = sqrt (0.5 * LzMax );
  for (int x = 0; x < NbrPoints; ++x)
    {
      Value[0] = X;
      Sum = 0.0;
      for (int i = 0; i <= LzMax; ++i)
	{
	  Basis.GetFunctionValue(Value, TmpValue, i);
	  Sum += ExpansionCoeff[i] * Norm(TmpValue* TmpValue);
	}
      File << (X * Factor2) << " " << Norm(Sum)  << endl;
      X += XInc;
    }
  File.close();

  delete[] OutputNameCorr;	  
  delete[] ExpansionCoeff;

  return 0;
}
