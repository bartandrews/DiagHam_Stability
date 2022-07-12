#include "config.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSpinLValue" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "no-spin", "do not compute the S^2 value of the state");
//  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
//  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "show-extracted", "show values extracted from file name");
//  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinLValue -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereLValue -h" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  int TotalSz = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  bool SzSymmetrizedBasis = ((BooleanOption*) Manager["szsymmetrized-basis"])->GetBoolean();
  bool SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();
  bool LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
  bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (NbrParticles==0)
    if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["state"])->GetString(), NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetrizedBasis, SzMinusParity, 
							     LzSymmetrizedBasis, LzMinusParity, FermionFlag) == false)
      {
	return -1;
      }
    else
      {
	if (((BooleanOption*) Manager["show-extracted"])->GetBoolean() == true)
	  {
	    cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << "  TotalSz=" << TotalSz;
	    if (LzSymmetrizedBasis == true)
	      {
		cout << "  Lz symmetrized basis ";
		if (LzMinusParity == true)
		  cout << "(minus parity) ";
		else
		  cout << "(plus parity) ";
	      }
	    if (SzSymmetrizedBasis == true)
	      {
		cout << "  Sz symmetrized basis ";
		if (SzMinusParity == true)
		  cout << "(minus parity) ";
		else
		  cout << "(plus parity) ";
	      }
	    cout << endl;
	  }
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
  if (((SingleStringOption*) Manager["statistics"])->GetString() != 0)
    if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
      {
	FermionFlag = true;
      }
    else
      if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
	{
	  FermionFlag = false;
	}
      else
	{
	  cout << ((SingleStringOption*) Manager["statistics"])->GetString() << " is an undefined statistics" << endl;
	}  
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* StateFileName = ((SingleStringOption*) Manager["state"])->GetString();
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  RealVector State;
  if (State.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }


  long MemorySpace = 9l << 20;
  ParticleOnSphereWithSpin* Space;
  if (FermionFlag == true)
    {
      if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (LzMax <= 63)
#else
		  if (LzMax <= 31)
#endif
		    {
		      Space = new FermionOnSphereWithSpinLong(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return -1;
		    }	
	      }
	}
      else
	{
#ifdef __128_BIT_LONGLONG__
	  if (LzMax >= 61)
#else
	    if (LzMax >= 29)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  if (SzSymmetrizedBasis == true) 
	    if (LzSymmetrizedBasis == false)
	      {
#ifdef __64_BITS__
		if (LzMax <= 28)
#else
		  if (LzMax <= 13)
#endif
		    {
		      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetry(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		    }
		  else
		    {
		      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinSzSymmetryLong(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetryLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		    }
		  }
	    else
#ifdef __64_BITS__
	      if (LzMax <= 28)
#else
		if (LzMax <= 13)
#endif
		  {
		    if (((SingleStringOption*) Manager["load-hilbert"])->GetString() == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, LzMax, SzMinusParity,
									LzMinusParity, MemorySpace);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetry(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		  }
		else
		  {
		    if (((SingleStringOption*) Manager["load-hilbert"])->GetString() == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrParticles, LzMax, SzMinusParity,
									    LzMinusParity, MemorySpace);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetryLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
		    
		  }
	      else
#ifdef __64_BITS__
		if (LzMax <= 28)
#else
		  if (LzMax <= 13)
#endif
		    {
		      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetry(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);	      
		    }
		  else
		    {
		      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinLzSymmetryLong(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetryLong(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);	      
		    }
	}      
    }
  else
    {
      Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
    }
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  ParticleOnSphereWithSpinL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0);
  RealVector TmpState(Space->GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &State, &TmpState);
  Operation.ApplyOperation(Architecture.GetArchitecture());
  double L2Value = TmpState * State;
  double RawTmpAngularMomentum = 0.5 * (sqrt ((4.0 * L2Value) + 1.0) - 1.0);
  cout << "<L^2> = " << L2Value << endl
       << "<L> = " << RawTmpAngularMomentum << endl;
  if (((BooleanOption*) Manager["no-spin"])->GetBoolean() == false)
    {
      ParticleOnSphereWithSpinS2Hamiltonian Hamiltonian2 (Space, NbrParticles, LzMax, TotalLz, TotalSz, Architecture.GetArchitecture(), 1.0, 0);
      VectorHamiltonianMultiplyOperation Operation2 (&Hamiltonian2, &State, &TmpState);
      Operation2.ApplyOperation(Architecture.GetArchitecture());
      L2Value = TmpState * State;
      RawTmpAngularMomentum = 0.5 * (sqrt ((4.0 * L2Value) + 1.0) - 1.0);
      cout << "<S^2> = " << L2Value << endl
	   << "<S> = " << RawTmpAngularMomentum << endl;
    }
  delete Space;
  return 0;
}

