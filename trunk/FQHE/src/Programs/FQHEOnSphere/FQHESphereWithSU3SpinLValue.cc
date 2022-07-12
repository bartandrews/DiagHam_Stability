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

#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinZ3Symmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzZ3Symmetry.h"

#include "Hamiltonian/ParticleOnSphereWithSU3SpinL2Hamiltonian.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSU3SpinLValue" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('t', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "tzsymmetrized-basis", "use Tz <-> -Tz symmetrized version of the basis (only valid if total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "z3symmetrized-basis", "use Z3 symmetrized version of the basis (only valid if total-y=0 and total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "show-extracted", "show values extracted from file name");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU3SpinLValue -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereWithSU3SpinLValue -h" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
  bool Z3SymmetrizedBasis = Manager.GetBoolean("z3symmetrized-basis");
  bool TzMinusParity = Manager.GetBoolean("minus-tzparity");
  bool LzMinusParity = Manager.GetBoolean("minus-lzparity");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (NbrParticles==0)
    {
      if (FQHEOnSphereWithSU3SpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, TotalTz, TotalY,
								  TzSymmetrizedBasis, TzMinusParity, Z3SymmetrizedBasis,
								  LzSymmetrizedBasis, LzMinusParity, FermionFlag) == false)
	{
	  return -1;
	}
      else
	{
	  if (Manager.GetBoolean("show-extracted") == true)
	    {
	      cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << "  TotalTz=" << TotalTz;
	      if (LzSymmetrizedBasis == true)
		{
		  cout << "  Lz symmetrized basis ";
		  if (LzMinusParity == true)
		    cout << "(minus parity) ";
		  else
		    cout << "(plus parity) ";
		}
	      if (TzSymmetrizedBasis == true)
		{
		  cout << "  Tz symmetrized basis ";
		  if (TzMinusParity == true)
		    cout << "(minus parity) ";
		  else
		    cout << "(plus parity) ";
		}
	      if (Z3SymmetrizedBasis == true)
		{
		  cout << "  Z3 symmetrized basis ";
		}
	      cout << endl;
	    }
	}
    }
  if (Manager.GetBoolean("lzsymmetrized-basis") == true)
    {
      LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
      LzMinusParity = Manager.GetBoolean("minus-lzparity");      
    }
  if (Manager.GetBoolean("tzsymmetrized-basis") == true)
    {
      TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
      TzMinusParity = Manager.GetBoolean("minus-tzparity");
    }
  if (Manager.GetBoolean("z3symmetrized-basis") == true)
    {
      Z3SymmetrizedBasis = Manager.GetBoolean("z3symmetrized-basis");
    }
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	  {
	    FermionFlag = false;
	}
	else
	  {
	    cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	  }  
    }
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* StateFileName = Manager.GetString("state");
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
  ParticleOnSphereWithSU3Spin* Space=0;
  if (FermionFlag == true)
    {
      if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == false))
	{
#ifdef __64_BITS__
	  if (LzMax <= 20)
#else
	    if (LzMax <= 9)
#endif
	      {
		Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, LzMax, TotalTz, TotalY, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	}
      else
	{
#ifdef __64_BITS__
	  if (LzMax > 20)
#else
	    if (LzMax > 9)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == false))
	    {
	      if (Manager.GetString("load-hilbert") == 0)
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(NbrParticles, TotalLz, LzMax, TotalY, Manager.GetBoolean("minus-tzparity"), MemorySpace);
		  if (Manager.GetString("save-hilbert") != 0)
		    {
		      ((FermionOnSphereWithSU3SpinTzSymmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		}
	      else
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
		}
	    }
	  else
	    if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == true))
	      {
		if (Manager.GetString("load-hilbert") == 0)
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(NbrParticles, TotalLz, LzMax, TotalTz, MemorySpace);
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereWithSU3SpinZ3Symmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
		else
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(Manager.GetString("load-hilbert"), MemorySpace);
		  }
	      }
	    else
	      if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == true))
		{
		  if (Manager.GetString("load-hilbert") == 0)
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(NbrParticles, TotalLz, LzMax, Manager.GetBoolean("minus-tzparity"), MemorySpace);
		      if (Manager.GetString("save-hilbert") != 0)
			{
			  ((FermionOnSphereWithSU3SpinTzZ3Symmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			  return 0;
			}
		    }
		  else
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(Manager.GetString("load-hilbert"), MemorySpace);
		    }
		}
	}
    }
  else
    {
      Space = 0;
    }
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  ParticleOnSphereWithSU3SpinL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0);
  RealVector TmpState(Space->GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &State, &TmpState);
  Operation.ApplyOperation(Architecture.GetArchitecture());
  double L2Value = TmpState * State;
  double RawTmpAngularMomentum = 0.5 * (sqrt ((4.0 * L2Value) + 1.0) - 1.0);
  cout << "<L^2> = " << L2Value << endl
       << "<L> = " << RawTmpAngularMomentum << endl;
  delete Space;
  return 0;
}

