#include "config.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

#include "Operator/ParticleOnSphereWithSU4DensityDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;

using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEFermionsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin (i.e valley SU(2) degeneracy) of the system", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-entanglement", "use a define value for the spin-isopsin entanglement of the system");
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsSU4Correlation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  int IsoSzTotal = ((SingleIntegerOption*) Manager["total-isosz"])->GetInteger();
  int TotalEntanglement = ((SingleIntegerOption*) Manager["total-entanglement"])->GetInteger();
  int LandauLevel = ((SingleIntegerOption*) Manager["landau-level"])->GetInteger();
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  bool ChordFlag = ((BooleanOption*) Manager["chord"])->GetBoolean();
  unsigned long MemorySpace = 9ul << 20;
  bool Statistics = true;

  if (NbrParticles==0)
    if (FQHEOnSphereWithSU4SpinFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["state"])->GetString(), NbrParticles, LzMax, TotalLz, 
								SzTotal, IsoSzTotal, TotalEntanglement, Statistics) == false)
      {
	cout << "error while retrieving system informations from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
	return -1;
      }
  cout << NbrParticles << " " << TotalLz << " " << SzTotal << " " << IsoSzTotal << " " << TotalEntanglement << endl;

  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEFermionsCorrelation requires a state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
      return -1;      
    }

  int NbrUp = (NbrParticles + SzTotal) >> 1;
  int NbrDown = (NbrParticles - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int NbrPlus = (NbrParticles + IsoSzTotal) >> 1;
  int NbrMinus = (NbrParticles - IsoSzTotal) >> 1;
  if ((NbrPlus < 0) || (NbrMinus < 0))
    {
      cout << "This value of the isospin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int NbrEntanglementPlus = (NbrParticles + TotalEntanglement) >> 1;
  int NbrEntanglementMinus = (NbrParticles - TotalEntanglement) >> 1;
  if ((((BooleanOption*) Manager["use-entanglement"])->GetBoolean()) && 
      ((NbrEntanglementPlus < 0) || (NbrEntanglementMinus < 0)))
    {
      cout << "This value of the entanglement projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

  ParticleOnSphereWithSU4Spin* Space;
  if (((BooleanOption*) Manager["use-entanglement"])->GetBoolean())
    Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, MemorySpace);
  else
    Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, LzMax, SzTotal, IsoSzTotal, MemorySpace);

  AbstractFunctionBasis* Basis;
  if (LandauLevel == 0)
    Basis = new ParticleOnSphereFunctionBasis(LzMax);
  else
    Basis = new ParticleOnSphereGenericLLFunctionBasis(LzMax - (2 * LandauLevel), LandauLevel);


  Complex* Sum = new Complex [10];
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue = 0.0;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);

  Complex** PrecalculatedValues = new Complex* [10];
  for (int i = 0; i < 10; ++i)
    PrecalculatedValues[i] = new Complex [LzMax + 1];
  Basis->GetFunctionValue(Value, TmpValue, LzMax);
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator (Space, i, 0, LzMax, 0, i, 0, LzMax, 0);
      PrecalculatedValues[0][i] =   Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator (Space, i, 1, LzMax, 1, i, 1, LzMax, 1);
      PrecalculatedValues[4][i] =  Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator (Space, i, 2, LzMax, 2, i, 2, LzMax, 2);
      PrecalculatedValues[7][i] =  Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator (Space, i, 3, LzMax, 3, i, 3, LzMax, 3);
      PrecalculatedValues[9][i] =  Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator1 (Space, i, 0, LzMax, 1, i, 0, LzMax, 1);
      PrecalculatedValues[1][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
//       ParticleOnSphereWithSU4DensityDensityOperator Operator2 (Space, i, 1, LzMax, 0, i, 1, LzMax, 0);
//       PrecalculatedValues[1][i] =  (Operator1.MatrixElement(State, State) + Operator2.MatrixElement(State, State))* TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator1 (Space, i, 0, LzMax, 2, i, 0, LzMax, 2);
      PrecalculatedValues[2][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
//       ParticleOnSphereWithSU4DensityDensityOperator Operator2 (Space, i, 2, LzMax, 0, i, 2, LzMax, 0);
//       PrecalculatedValues[2][i] =  (Operator1.MatrixElement(State, State) + Operator2.MatrixElement(State, State))* TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator1 (Space, i, 0, LzMax, 3, i, 0, LzMax, 3);
      PrecalculatedValues[3][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
//       ParticleOnSphereWithSU4DensityDensityOperator Operator2 (Space, i, 3, LzMax, 0, i, 3, LzMax, 0);
//       PrecalculatedValues[3][i] =  (Operator1.MatrixElement(State, State) + Operator2.MatrixElement(State, State))* TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator1 (Space, i, 1, LzMax, 2, i, 1, LzMax, 2);
      PrecalculatedValues[5][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
//       ParticleOnSphereWithSU4DensityDensityOperator Operator2 (Space, i, 2, LzMax, 1, i, 2, LzMax, 1);
//       PrecalculatedValues[5][i] =  (Operator1.MatrixElement(State, State) + Operator2.MatrixElement(State, State))* TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator1 (Space, i, 1, LzMax, 3, i, 1, LzMax, 3);
      PrecalculatedValues[6][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
//       ParticleOnSphereWithSU4DensityDensityOperator Operator2 (Space, i, 3, LzMax, 1, i, 3, LzMax, 1);
//       PrecalculatedValues[6][i] =  (Operator1.MatrixElement(State, State) + Operator2.MatrixElement(State, State))* TmpValue * Conj(TmpValue);
    }
  for (int i = 0; i <= LzMax; ++i)
    {
      ParticleOnSphereWithSU4DensityDensityOperator Operator1 (Space, i, 2, LzMax, 3, i, 2, LzMax, 3);
      PrecalculatedValues[8][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
//       ParticleOnSphereWithSU4DensityDensityOperator Operator2 (Space, i, 3, LzMax, 2, i, 3, LzMax, 2);
//       PrecalculatedValues[8][i] =  (Operator1.MatrixElement(State, State) + Operator2.MatrixElement(State, State))* TmpValue * Conj(TmpValue);
    }
  double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
  double Factor2;
  ofstream File;
  File.precision(14);
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(((SingleStringOption*) Manager["state"])->GetString(), "vec", "rhorho");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << ((SingleStringOption*) Manager["state"])->GetString() << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  if (((BooleanOption*) Manager["radians"])->GetBoolean() == true)
    {
      Factor2 = 1.0;
      File << "# dist (rad)";
    }
  else
    {
      Factor2 = sqrt (0.5 * LzMax);
      if (ChordFlag == false)
	File << "# dist (sqrt{S})";
      else
	File << "# dist (chord)";
    }
//  File << " rho_{up,up} rho_{up,um}+rho_{um,up} rho_{up,dp}+rho_{dp,up} rho_{up,dm}+rho_{dm,up} rho_{um,um} rho_{um,dp}+rho_{dp,um} rho_{um,dm}+rho_{dm,um} rho_{dp,dp} rho_{dp,dm}+rho_{dm,dp} rho_{dm,dm}" << endl;
  File << " rho_{up,up} rho_{up,um} rho_{up,dp} rho_{up,dm} rho_{um,um} rho_{um,dp} rho_{um,dm} rho_{dp,dp} rho_{dp,dm} rho_{dm,dm}" << endl;
  for (int x = 0; x < NbrPoints; ++x)
    {
      Value[0] = X;
      for (int j = 0; j < 10; ++j)
	Sum[j] = 0.0;
      for (int i = 0; i <= LzMax; ++i)
 	{
 	  Basis->GetFunctionValue(Value, TmpValue, i);
	  for (int j = 0; j < 10; ++j)	    
	    Sum[j] += PrecalculatedValues[j][i] * (Conj(TmpValue) * TmpValue);
 	}
      if (ChordFlag == false)
	File << (X * Factor2);
      else
	File << (2.0 * Factor2 * sin (X * 0.5));
      for (int j = 0; j < 10; ++j)
	File << " " << (Norm(Sum[j])  * Factor1);
      File << endl;
      X += XInc;
    }
  for (int i = 0; i < 10; ++i)
    delete[] PrecalculatedValues[i];
  delete[] PrecalculatedValues;
  delete[] Sum;

}

