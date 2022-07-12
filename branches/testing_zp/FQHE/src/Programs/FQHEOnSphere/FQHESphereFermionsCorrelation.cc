#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

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
  OptionManager Manager ("FQHESphereFermionsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsCorrelation -h" << endl;
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
  int LandauLevel = ((SingleIntegerOption*) Manager["landau-level"])->GetInteger();
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  bool DensityFlag = ((BooleanOption*) Manager["density"])->GetBoolean();
  bool ChordFlag = ((BooleanOption*) Manager["chord"])->GetBoolean();
  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  bool CoefficientOnlyFlag = ((BooleanOption*) Manager["coefficients-only"])->GetBoolean();
  bool Statistics = true;
 if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["eigenstate"])->GetString(),
						  NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["eigenstate"])->GetString() << endl;
      return -1;
    }

  if (((SingleStringOption*) Manager["eigenstate"])->GetString() == 0)
    {
      cout << "FQHESphereFermionsCorrelation requires a state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["eigenstate"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["eigenstate"])->GetString() << endl;
      return -1;      
    }

  ParticleOnSphere* Space = 0;
  if (HaldaneBasisFlag == false)
    {
#ifdef __64_BITS__
      if (LzMax <= 62)
#else
	if (LzMax <= 30)
#endif
	  if ((SymmetrizedBasis == false) || (TotalLz != 0))
	    {
	      Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
	    }
	  else
	    {
	      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		Space = new FermionOnSphereSymmetricBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
	      else
		Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
	      if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		{
		  ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
		  return 0;
		}
	    }
	else
#ifdef __128_BIT_LONGLONG__
	  if (LzMax <= 126)
#else
	    if (LzMax <= 62)
#endif
	      {
		if ((SymmetrizedBasis == false) || (TotalLz != 0))
		  Space = new FermionOnSphereLong(NbrParticles, TotalLz, LzMax, MemorySpace);
		else
		  {
		    if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		      Space = new FermionOnSphereSymmetricBasisLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		    else
		      Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		    if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		      {
			((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			return 0;
		      }
		  }
	      }
	    else
	      Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax, MemorySpace);
    }
  else
    {
      int* ReferenceState = 0;
      if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
	{
	  ReferenceState = new int[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    ReferenceState[i] = 0;
	  if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "laughlin") == 0)
	    for (int i = 0; i <= LzMax; i += 3)
	      ReferenceState[i] = 1;
	  else
	    if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "pfaffian") == 0)
	      for (int i = 0; i <= LzMax; i += 4)
		{
		  ReferenceState[i] = 1;
		  ReferenceState[i + 1] = 1;
		}
	    else
	      if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "readrezayi3") == 0)
		for (int i = 0; i <= LzMax; i += 5)
		  {
		    ReferenceState[i] = 1;
		    ReferenceState[i + 1] = 1;
		    ReferenceState[i + 2] = 1;
		  }
	      else
		{
		  cout << "unknown reference state " << ((SingleStringOption*) Manager["reference-state"])->GetString() << endl;
		  return -1;
		}
	}
      else
	{
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
	      return -1;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return -1;     
	    }
	}	
      if (SymmetrizedBasis == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		  Space = new FermionOnSphereHaldaneBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		else
		  Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
		if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		  {
		    ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
		    return 0;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		      Space = new FermionOnSphereHaldaneBasisLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		    else
		      Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
		    if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		      {
			((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			return 0;
		      }
		  }	       
	}
      else
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		  Space = new FermionOnSphereHaldaneSymmetricBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		else
		  Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		  {
		    ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			 return 0;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		      Space = new FermionOnSphereHaldaneSymmetricBasisLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		    else
		      Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
		    if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
		      {
			((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
			return 0;
		      }
		  }
	}
    }

  cout << Space->GetHilbertSpaceDimension() << endl;

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
  ofstream File;
  File.precision(14);
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
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
  if (CoefficientOnlyFlag == true)
    {
      if (DensityFlag == true)      
	File << "# density correlation coefficients for " << ((SingleStringOption*) Manager["eigenstate"])->GetString() << endl;
      else
	File << "# density-density correlation coefficients for " << ((SingleStringOption*) Manager["eigenstate"])->GetString() << endl;
      File << "#" << endl << "# (l+S)    n_l" << endl;
      for (int i = 0; i <= LzMax; ++i)
	File << i << " " << PrecalculatedValues[i]<< endl;
    }
  else
    {
      double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      if (DensityFlag == true)
	Factor1 = 1.0;//4.0 * M_PI;
      double Factor2;
      if (((BooleanOption*) Manager["radians"])->GetBoolean() == true)
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


