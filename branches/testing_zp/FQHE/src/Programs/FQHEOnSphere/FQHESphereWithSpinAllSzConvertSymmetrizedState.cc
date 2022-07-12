#include "Vector/RealVector.h"


#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereWithSpinConvertSymmetrizedState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of _*sym_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-file"])->GetString() == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereWithSpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger(); 
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  bool SymmetrizeFlag = ((BooleanOption*) Manager["symmetrize"])->GetBoolean();
  bool LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
  bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();
  long MemorySpace = 9l << 20;

/*
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(), NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetrizedBasis, SzMinusParity, 
							   LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
 */     
  if (((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean() == true)
    {
      LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
      LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();      
    }

  if (((NbrParticles * LzMax) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
    }

  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;      
    }


      RealVector OutputState;
#ifdef __64_BITS__
      if (LzMax <= 31)
#else
	if (LzMax <= 15)
#endif
	  {
	    FermionOnSphereWithSpinAllSzLzSymmetry* InitialSpace = 0;
	      if (LzSymmetrizedBasis == true)
	       InitialSpace = new FermionOnSphereWithSpinAllSzLzSymmetry(NbrParticles, LzMax, LzMinusParity, MemorySpace);

	    FermionOnSphereWithSpinAllSz TargetSpace(NbrParticles, TotalLz, LzMax);
	    if (SymmetrizeFlag)
	      {
		if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		OutputState = InitialSpace->ConvertToSymmetricNbodyBasis(State, TargetSpace);
	      }
	    else
	      {
		if (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		OutputState = InitialSpace->ConvertToNbodyBasis(State, TargetSpace);
		delete InitialSpace;
	      }
	  }
	else
	  {
		cout<<"Not implemented"<<endl;
	    /*FermionOnSphereWithSpinLzSzSymmetryLong* InitialSpace = 0;
	    if (SzSymmetrizedBasis == true) 
	      if (LzSymmetrizedBasis == false)
		InitialSpace = new FermionOnSphereWithSpinSzSymmetryLong(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
	      else
		InitialSpace = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrParticles, LzMax, SzMinusParity, LzMinusParity, MemorySpace);
	    else
	      InitialSpace = new FermionOnSphereWithSpinLzSymmetryLong(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
	    FermionOnSphereWithSpinLong TargetSpace(NbrParticles, TotalLz, LzMax, TotalSz);
	    if (SymmetrizeFlag)
	      {
		if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		OutputState = InitialSpace->ConvertToSymmetricNbodyBasis(State, TargetSpace);
	      }
	    else
	      {
		if (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between Hilbert space and input state" << endl;
		    return -1;
		  }
		OutputState = InitialSpace->ConvertToNbodyBasis(State, TargetSpace);
		delete InitialSpace;
	      }
		*/
	  }
      if (OutputState.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString()) == false)
	{
	  cout << "error while writing output state " << ((SingleStringOption*) Manager["output-file"])->GetString() << endl;
	  return -1;
	}

}

