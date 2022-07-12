#include "Vector/RealVector.h"


#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
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
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('f', "fermion", "use fermionic statistic (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*SystemGroup) += new BooleanOption  ('\n', "conjugate-down", "particle-hole conjugate down spins only");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of _*sym_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinConvertHaldaneSymmetrizedState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-file"])->GetString() == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereWithSpinHaldaneConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger(); 
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  int TotalSz = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  bool SymmetrizeFlag = ((BooleanOption*) Manager["symmetrize"])->GetBoolean();
  bool SzSymmetrizedBasis = ((BooleanOption*) Manager["szsymmetrized-basis"])->GetBoolean();
  bool SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();
  bool LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
  bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();
  bool Statistics = true;
  long MemorySpace = 9l << 20;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(), NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetrizedBasis, SzMinusParity, 
							   LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
      
  if (((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean() == true)
    {
      LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
      LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();      
    }
  if (((BooleanOption*) Manager["szsymmetrized-basis"])->GetBoolean() == true)
    {
      SzSymmetrizedBasis = ((BooleanOption*) Manager["szsymmetrized-basis"])->GetBoolean();
      SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();
    }
  if ((((BooleanOption*) Manager["boson"])->GetBoolean() == true) || (((BooleanOption*) Manager["fermion"])->GetBoolean() == true))
    {
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	Statistics = false;
      else
	Statistics = true;
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

  if (Manager.GetBoolean("conjugate-down"))
    {
      #ifdef __64_BITS__
      if (LzMax <= 31)
#else
	if (LzMax <= 15)
#endif
	  {
	    FermionOnSphereWithSpin MainSpace(NbrParticles, TotalLz, LzMax, TotalSz);
	    RealVector Destination;
	    MainSpace.ParticleHoleConjugateDownSpins(State,Destination);
	    char *OutputName=Manager.GetString("output-file");
	    if (OutputName==NULL)
	      {
		OutputName= new char[strlen(Manager.GetString("input-file"))+10];
		sprintf(OutputName,"%s.qhd",Manager.GetString("input-file"));
	      }
	    Destination.WriteVector(OutputName);
	    cout << "Down particle conjugated state written to "<<OutputName<<endl;
	    exit(0);
	  }
	else
	  {
	    cout << "Cannot represent state in a single word"<<endl;
	    exit(1);
	  }
    } 

  if (Statistics == true)
    {
      RealVector OutputState;
#ifdef __64_BITS__
      if (LzMax <= 31)
#else
	if (LzMax <= 15)
#endif
	  {
	    FermionOnSphereWithSpinLzSzSymmetry* InitialSpace = 0;
	    if (SzSymmetrizedBasis == true) 
	      if (LzSymmetrizedBasis == false)
		InitialSpace = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
	      else
		{		  
		  if (Manager.GetString("reference-file") != 0)
		    {		      
		      int** ReferenceStates = 0;
		      int NbrReferenceStates;		    
		      bool TexturelessFlag;
		      if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates, TexturelessFlag) == false)
			{
			  cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
			  return 0;
			}      		  		  
		      if (TexturelessFlag == false ) 
			{	
			   InitialSpace = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, ReferenceStates, NbrReferenceStates, (unsigned long)MemorySpace);
			}
		      else
			{			  			    
			    int **TexturelessReferenceState = new int*[NbrReferenceStates];
			    for ( int j = 0 ; j < NbrReferenceStates ; j++ ) 
			      {
				TexturelessReferenceState[j] = new int[LzMax+1];
				for ( int i = 0 ; i < (LzMax + 1) ; i++ )
				  {
				    if ( ReferenceStates[j][i] == 3 ) 
				      {
					  TexturelessReferenceState[j][i] = 2;
				      }
				    else if ( (ReferenceStates[j][i] == 1) || (ReferenceStates[j][i] == 2) ) 
				      {
					  TexturelessReferenceState[j][i] = 1;
				      }
				    else
				      {
					  TexturelessReferenceState[j][i] = 0;
				      }
				  }
			      }	
			  InitialSpace = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, TexturelessReferenceState, NbrReferenceStates, true, MemorySpace);
			}
		    }
		  else
		    {
		      InitialSpace = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, MemorySpace);
		    }
		}
	    else
	      InitialSpace = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
	    FermionOnSphereWithSpin TargetSpace(NbrParticles, TotalLz, LzMax, TotalSz);
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
	    FermionOnSphereWithSpinLzSzSymmetryLong* InitialSpace = 0;
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
	  }
      if (OutputState.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString()) == false)
	{
	  cout << "error while writing output state " << ((SingleStringOption*) Manager["output-file"])->GetString() << endl;
	  return -1;
	}
    }
}

