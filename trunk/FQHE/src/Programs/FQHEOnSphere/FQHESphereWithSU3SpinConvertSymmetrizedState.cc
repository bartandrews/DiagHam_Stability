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
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzZ3Symmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinZ3Symmetry.h"

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
  OptionManager Manager ("FQHESphereWithSU3SpinConvertSymmetrizedState" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('t', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
  (*SystemGroup) += new BooleanOption  ('f', "fermion", "use fermionic statistic (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*SystemGroup) += new BooleanOption  ('\n', "tzsymmetrized-basis", "use Tz <-> -Tz symmetrized version of the basis (only valid if total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "z3symmetrized-basis", "use Z3 symmetrized version of the basis (only valid if total-y=0 and total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of _*sym_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU3SpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereWithSU3SpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  bool SymmetrizeFlag = Manager.GetBoolean("symmetrize");
  bool TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
  bool Z3SymmetrizedBasis = Manager.GetBoolean("z3symmetrized-basis");
  bool TzMinusParity = Manager.GetBoolean("minus-tzparity");
  bool Statistics = true;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  long MemorySpace = 9l << 20;
  if (FQHEOnSphereWithSU3SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), 
							      NbrParticles, LzMax, TotalLz, TotalTz, TotalY,
							      TzSymmetrizedBasis, TzMinusParity, Z3SymmetrizedBasis,
							      LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-file") << endl;
      return -1;
    }
      
  if (Manager.GetBoolean("z3symmetrized-basis") == true)
    {
      Z3SymmetrizedBasis = Manager.GetBoolean("z3symmetrized-basis");
    }
  if (Manager.GetBoolean("tzsymmetrized-basis") == true)
    {
      TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
      TzMinusParity = Manager.GetBoolean("minustszparity");
    }
  if ((Manager.GetBoolean("boson") == true) || (Manager.GetBoolean("fermion") == true))
    {
      if (Manager.GetBoolean("boson") == true)
	Statistics = false;
      else
	Statistics = true;
    }
  if (((NbrParticles * LzMax) & 1) != (TotalLz != 0))
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


  if (Statistics == true)
    {
      RealVector OutputState;
      FermionOnSphereWithSU3SpinTzSymmetry* InitialSpace = 0;
      if (TzSymmetrizedBasis == true) 
	if (Z3SymmetrizedBasis == false)
	  InitialSpace = new FermionOnSphereWithSU3SpinTzSymmetry(NbrParticles, TotalLz, LzMax, TotalY, TzMinusParity, MemorySpace);
	else
	  InitialSpace = new FermionOnSphereWithSU3SpinTzZ3Symmetry(NbrParticles, TotalLz, LzMax, TzMinusParity, MemorySpace);
      else
	InitialSpace = new FermionOnSphereWithSU3SpinZ3Symmetry(NbrParticles, TotalLz, LzMax, MemorySpace);
      FermionOnSphereWithSU3Spin TargetSpace(NbrParticles, TotalLz, LzMax, TotalTz, TotalY, MemorySpace);
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
      if (OutputState.WriteVector(Manager.GetString("output-file")) == false)
	{
	  cout << "error while writing output state " << Manager.GetString("output-file") << endl;
	  return -1;
	}
    }
}

