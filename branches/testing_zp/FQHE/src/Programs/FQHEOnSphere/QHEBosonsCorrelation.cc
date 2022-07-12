#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphere.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

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
  OptionManager Manager ("QHEBosonsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector", 0, true, 0);
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleStringOption  ('i', "interaction-name", "name of the interaction (used for output file name)", "delta");
  (*SystemGroup) += new SingleStringOption ('a', "add-filename", "add a string with additional informations to the output file name(just before the .dat extension)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsCorrelation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEBosonsCorrelation requires a state" << endl;
      return -1;
    }
  RealVector State;
//   if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
//     {
//       cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
//       return -1;      
//     }
   char* OutputNameCorr = new char [256 + strlen (((SingleStringOption*) Manager["interaction-name"])->GetString())];
//   if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
//     {
//       sprintf (OutputNameCorr, "bosons_%s_n_%d_2s_%d.rho_rho.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrBosons, LzMax);
//     }
//   else
//     {
//       sprintf (OutputNameCorr, "bosons_%s_n_%d_2s_%d_%s.rho_rho.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrBosons, LzMax,
// 	       ((SingleStringOption*) Manager["add-filename"])->GetString());
//     }

  BosonOnSphere Space (NbrBosons, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);

  for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
    {
      Space.PrintState(cout, i) << endl;
    }
  return 0;

  Complex Sum (0.0, 0.0);
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);
  Complex* PrecalculatedValues = new Complex [LzMax + 1];
	  
  for (int i = 0; i <= LzMax; ++i)
    {
      Basis.GetFunctionValue(Value, TmpValue, LzMax);
      ParticleOnSphereDensityDensityOperator Operator (&Space, i, LzMax, i, LzMax);
      PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }
  ofstream File;
  File.precision(14);
  File.open(OutputNameCorr, ios::binary | ios::out);
  double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrBosons * NbrBosons));
  double Factor2;
    if (((BooleanOption*) Manager["radians"])->GetBoolean() == true)
    Factor2 = 1.0;
  else
    Factor2 = sqrt (0.5 * LzMax );
  for (int x = 0; x < NbrPoints; ++x)
    {
      Value[0] = X;
      int Pos = 0;
      Sum = 0.0;
      for (int i = 0; i <= LzMax; ++i)
	{
	  Basis.GetFunctionValue(Value, TmpValue, i);
	  Sum += PrecalculatedValues[Pos] * (Conj(TmpValue) * TmpValue);
	  ++Pos;
	}
      File << (X * Factor2) << " " << Norm(Sum)  * Factor1 << endl;
      X += XInc;
    }
  File.close();

  delete[] OutputNameCorr;	  
  delete[] PrecalculatedValues;

  return 0;
}


