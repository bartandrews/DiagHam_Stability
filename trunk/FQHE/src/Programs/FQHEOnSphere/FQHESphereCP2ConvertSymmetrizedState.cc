#include "Vector/RealVector.h"


#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/BosonOnCP2.h"
#include "HilbertSpace/BosonOnCP2TzSymmetry.h"
#include "HilbertSpace/BosonOnCP2TzZ3Symmetry.h"

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
  OptionManager Manager ("FQHESphereCP2ConvertSymmetrizedState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('t', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "tzsymmetrized-basis", "use Tz <-> -Tz symmetrized version of the basis (only valid if total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "tzZ3symmetrized-basis", "use Tz <-> -Tz and Z3 permutations symmetrized version of the basis (only valid if total-tz=0 and total-y = 0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
//   (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of _*sym_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereCP2ConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-file"])->GetString() == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereWithSU3SpinConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
//   bool SymmetrizeFlag = Manager.GetBoolean("symmetrize");
  bool TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
  bool TzZ3SymmetrizedBasis = Manager.GetBoolean("tzZ3symmetrized-basis");
  bool TzMinusParity = Manager.GetBoolean("minus-tzparity");
  bool Statistics = true;
  long MemorySpace = 9l << 20;
  
  if (FQHEOnCP2FindSystemInfoFromVectorFileName(Manager.GetString("input-file"), 
							      NbrParticles, NbrFluxQuanta, TotalTz, TotalY,
							      TzSymmetrizedBasis, TzMinusParity, TzZ3SymmetrizedBasis,
							      Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-file") << endl;
      return -1;
    }
      
  if (Manager.GetBoolean("tzZ3symmetrized-basis") == true)
    {
      TzZ3SymmetrizedBasis = Manager.GetBoolean("TzZ3symmetrized-basis");
    }
  if (Manager.GetBoolean("tzsymmetrized-basis") == true)
    {
      TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
      TzMinusParity = Manager.GetBoolean("minus-tzparity");
    }
  
  if ((((TotalY + 3*TotalTz + 2*NbrParticles*NbrFluxQuanta) % 6) != 0) || (((TotalY - 3*TotalTz + 2*NbrParticles*NbrFluxQuanta) % 6) != 0))
    {
      cout << "incompatible values for nbr-particles, nbr-flux, total-tz and total-y" << endl;
      return -1;
    }

  RealVector State;
  if (State.ReadVector (Manager.GetString("input-file")) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;      
    }


  if (Statistics == false)
    {
      RealVector OutputState;
      BosonOnCP2TzSymmetry* InitialSpace = 0;
      if (TzSymmetrizedBasis == true) 
	{
	  if (TzZ3SymmetrizedBasis == false)
	    InitialSpace = new BosonOnCP2TzSymmetry(NbrParticles, NbrFluxQuanta, TotalTz, TotalY, TzMinusParity, MemorySpace);
	  else
	    InitialSpace = new BosonOnCP2TzZ3Symmetry(NbrParticles, NbrFluxQuanta, TotalTz, TotalY, TzMinusParity, MemorySpace);
	}
      BosonOnCP2 TargetSpace(NbrParticles, NbrFluxQuanta, TotalTz, TotalY, MemorySpace);
//       if (SymmetrizeFlag)
// 	{
// 	  if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
// 	    {
// 	      cout << "dimension mismatch between Hilbert space and input state" << endl;
// 	      return -1;
// 	    }
// 	  OutputState = InitialSpace->ConvertToSymmetricNbodyBasis(State, TargetSpace);
// 	}
//       else
// 	{
	  if (InitialSpace->GetHilbertSpaceDimension() != State.GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	  OutputState = InitialSpace->ConvertToNbodyBasis(State, TargetSpace);
	  delete InitialSpace;
// 	}
      if (OutputState.WriteVector(Manager.GetString("output-file")) == false)
	{
	  cout << "error while writing output state " << Manager.GetString("output-file") << endl;
	  return -1;
	}
    }
}

