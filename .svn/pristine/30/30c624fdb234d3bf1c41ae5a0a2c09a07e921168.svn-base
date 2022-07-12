#include "Vector/RealVector.h"


#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "input state file name");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
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
      cout << "see man page for option syntax or type FQHESphereWithSpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereWithSpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  bool SymmetrizeFlag = Manager.GetBoolean("symmetrize");
  bool SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
  bool SzMinusParity = Manager.GetBoolean("minus-szparity");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool LzMinusParity = Manager.GetBoolean("minus-lzparity");
  bool Statistics = true;
  long MemorySpace = 9l << 20;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetrizedBasis, SzMinusParity, 
							   LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-file") << endl;
      return -1;
    }
      
  if (Manager.GetBoolean("lzsymmetrized-basis") == true)
    {
      LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
      LzMinusParity = Manager.GetBoolean("minus-lzparity");      
    }
  if (Manager.GetBoolean("szsymmetrized-basis") == true)
    {
      SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
      SzMinusParity = Manager.GetBoolean("minus-szparity");
    }
  if ((Manager.GetBoolean("boson") == true) || (Manager.GetBoolean("fermion") == true))
    {
      if (Manager.GetBoolean("boson") == true)
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
  if (State.ReadVector (Manager.GetString("input-file")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-file") << endl;
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
      char* OutputFileName = 0;
      if (Manager.GetString("output-file") != 0)
	{
	  OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFileName, Manager.GetString("output-file"));
	}
      else
	{
	  char* TmpOldString = new char [8];
	  sprintf (TmpOldString, "_su2_");
	  char* TmpSymString = new char [128];
	  if (LzSymmetrizedBasis == true)
	    {
	      if (SzSymmetrizedBasis == true)
		{
		  if (LzMinusParity == false)
		    {
		      if (SzMinusParity == false)
			sprintf (TmpSymString, "_su2_lzsym_1_szsym_1_");
		      else
			sprintf (TmpSymString, "_su2_lzsym_1_szsym_-1_");
		    }
		  else
		    {
		      if (SzMinusParity == false)
			sprintf (TmpSymString, "_su2_lzsym_-1_szsym_1_");
		      else
			sprintf (TmpSymString, "_su2_lzsym_-1_szsym_-1_");
		    }
		}
	      else
		{		 
		  if (LzMinusParity == false)
		    sprintf (TmpSymString, "_su2_lzsym_1_");
		  else
		    sprintf (TmpSymString, "_su2_lzsym_-1_");
		}
	    }
	  else
	    {
	      if (SzMinusParity == false)
		sprintf (TmpSymString, "_su2_szsym_1_");
	      else
		sprintf (TmpSymString, "_su2_szsym_-1_");
	    }
	  if (SymmetrizeFlag)
	    {
	      OutputFileName = ReplaceString(Manager.GetString("input-file"), TmpOldString, TmpSymString);
	    }
	  else
	    {
	      OutputFileName = ReplaceString(Manager.GetString("input-file"), TmpSymString, TmpOldString);
	    }
	  if (OutputFileName == 0)
	    {
	      cout << "can't guess output file name from " << Manager.GetString("input-file") << endl;
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
		InitialSpace = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, MemorySpace);
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
      if (OutputState.WriteVector(OutputFileName) == false)
	{
	  cout << "error while writing output state " << OutputFileName << endl;
	  return -1;
	}
    }
  else
    {
      RealVector OutputState;
      BosonOnSphereWithSU2Spin* InitialSpace = 0;
      if (LzSymmetrizedBasis == true)
	{
	  if (SzSymmetrizedBasis == true)
	    {
	      InitialSpace = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrParticles, LzMax, TotalSz, SzMinusParity, LzMinusParity);
	    }
	  else
	    {
	      InitialSpace = new BosonOnSphereWithSU2SpinLzSymmetry(NbrParticles, LzMax, TotalSz, LzMinusParity);
	    }
	}
      else
	{
	  if (SzSymmetrizedBasis == true)
	    {
	      InitialSpace = new BosonOnSphereWithSU2SpinSzSymmetry(NbrParticles, TotalLz, LzMax, TotalSz, SzMinusParity);
	    }
	  else
	    {
	      InitialSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, LzMax, TotalSz);
	    }
	}
      BosonOnSphereWithSU2Spin* TargetSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, LzMax, TotalSz);
      if (SymmetrizeFlag)
	{
	  if (TargetSpace->GetHilbertSpaceDimension() != State.GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	  OutputState = InitialSpace->ConvertToNbodyBasis(State, TargetSpace);
	}
      else
	{
	  if (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	  OutputState = InitialSpace->ConvertFromNbodyBasis(State, TargetSpace);
	}
      if (OutputState.WriteVector(OutputFileName) == false)
	{
	  cout << "error while writing output state " << OutputFileName << endl;
	  return -1;
	}
      delete[] OutputFileName;
      delete InitialSpace;
      delete TargetSpace;
    }
}

