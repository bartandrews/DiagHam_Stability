#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"


#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
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
  OptionManager Manager ("FQHESphereBosonsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;
  ParticleOnSphereManager ParticleManager(false, true, 1);
  ParticleManager.AddOptionGroup(&Manager);

  Architecture.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int LandauLevel = Manager.GetInteger("landau-level");
  int NbrPoints = Manager.GetInteger("nbr-points");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  bool DensityFlag = Manager.GetBoolean("density");
  bool ChordFlag = Manager.GetBoolean("chord");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("eigenstate"),
						   NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("eigenstate") << endl;
      return -1;
    }

  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHESphereFermionsCorrelation requires a state" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("eigenstate")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }

  ParticleOnSphere* Space = ParticleManager.GetHilbertSpace(TotalLz);
  cout << "N=" << NbrParticles << ", LzMax=" << LzMax << ", TotalLz=" << TotalLz << endl;
  cout << "dim = " << Space->GetHilbertSpaceDimension() << endl;


  AbstractFunctionBasis* Basis;
  if (LandauLevel == 0)
    Basis = new ParticleOnSphereFunctionBasis(LzMax);
  else
    Basis = new ParticleOnSphereGenericLLFunctionBasis(LzMax - (2 * LandauLevel), LandauLevel);

  Complex Sum (0.0, 0.0);
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);

  Complex* PrecalculatedValues = new Complex [LzMax + 1];	  
  RealVector State;
  if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
    {
      if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	  return -1;      
	}
      if (DensityFlag == false)
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, LzMax);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax, i, LzMax);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	  }
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    ParticleOnSphereDensityOperator Operator (Space, i);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State);
	  }
    }
  else
    {
      ComplexVector ComplexState;
      if (ComplexState.ReadVector (Manager.GetString("eigenstate")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	  return -1;      
	}
      if (DensityFlag == false)
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, LzMax);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax, i, LzMax);
	    PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState) * TmpValue * Conj(TmpValue);
	  }
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    ParticleOnSphereDensityOperator Operator (Space, i);
	    PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState);
	  }
    }

  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  if (Manager.GetBoolean("coefficients-only"))
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho-c");
	  else
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho");
	}
      else
	{
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("eigenstate") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  if (DensityFlag == true)      
    File << "# density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  else
    File << "# density-density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  File << "#" << endl << "# (l+S)    n_l" << endl;
  if (CoefficientOnlyFlag == false)
    {
      for (int i = 0; i <= LzMax; ++i)
	File << "# " << i << " " << PrecalculatedValues[i]<< endl;
    }
  else
    {
      for (int i = 0; i <= LzMax; ++i)
	File << i << " " << PrecalculatedValues[i]<< endl;
    }
  if (CoefficientOnlyFlag == false)
    {
      double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      if (DensityFlag == true)
	Factor1 = 1.0;//4.0 * M_PI;
      double Factor2;
      if (Manager.GetBoolean("radians") == true)
	Factor2 = 1.0;
      else
	Factor2 = sqrt (0.5 * LzMax);
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  Sum = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis->GetFunctionValue(Value, TmpValue, i);
	      Sum += PrecalculatedValues[i] * (Conj(TmpValue) * TmpValue);
	    }
	  if (ChordFlag == false)
	    File << (X * Factor2) << " " << (Norm(Sum)  * Factor1) << endl;
	  else
	    File << (2.0 * Factor2 * sin (X * 0.5)) << " " << Norm(Sum)  * Factor1 << endl;
	  X += XInc;
	}
    }
  File.close();
 
  delete[] PrecalculatedValues;

  return 0;
}


