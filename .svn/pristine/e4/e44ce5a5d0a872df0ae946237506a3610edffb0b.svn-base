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

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

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
  OptionManager Manager ("FQHESphereSUKToU1" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file that contains the SU(K) state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin  of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin (override symmetry found from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(3) spin (override symmetry found from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin (override symmetry found from file name)");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSUKToU1 -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state " << endl << "see man page for option syntax or type FQHESphereSUKToU1 -h" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  int TotalIsoSz = Manager.GetInteger("total-isosz");
  int TotalEntanglement = Manager.GetInteger("total-entanglement");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool TzSymmetrizedBasis = false;
  bool TzMinusParity = false;
  bool YSymmetrizedBasis = false;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  bool FermionFlag = false;
  bool SU2SymmetryFlag = false;
  bool SU3SymmetryFlag = false;
  bool SU4SymmetryFlag = false;

  char* StateFileName = Manager.GetString("state");
  FQHEOnSphereFindInternalSymmetryGroupFromFileName(StateFileName, SU2SymmetryFlag, SU3SymmetryFlag, SU4SymmetryFlag);

  if (Manager.GetBoolean("su2-spin") == true)
    SU2SymmetryFlag = true;
  else
    if (Manager.GetBoolean("su3-spin") == true)
      SU3SymmetryFlag = true;
    else
      if (Manager.GetBoolean("su4-spin") == true)
	SU4SymmetryFlag = true;

  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (NbrParticles==0)
    {
      if (SU2SymmetryFlag == true)
	if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetrizedBasis, SzMinusParity, 
								 LzSymmetrizedBasis, LzMinusParity, FermionFlag) == false)
	  {
	    cout << "error while retrieving system informations from file name " << Manager.GetString("state") << endl;
	    return -1;
	  }
      if (SU3SymmetryFlag == true)
	if (FQHEOnSphereWithSU3SpinFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, LzMax, TotalLz, TotalTz, TotalY, TzSymmetrizedBasis, TzMinusParity, 
								    YSymmetrizedBasis, LzSymmetrizedBasis, LzMinusParity, FermionFlag) == false)
	  {
	    cout << "error while retrieving system informations from file name " << Manager.GetString("state") << endl;
	    return -1;
	  }
      if (SU4SymmetryFlag == true)
	if (FQHEOnSphereWithSU4SpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, 
								    TotalLz, TotalIsoSz, TotalEntanglement, FermionFlag) == false)
	  {
	    cout << "error while retrieving system informations from file name " << Manager.GetString("state") << endl;
	    return -1;
	  }
    }
  cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << endl;
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	{
	  if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	    {
	      FermionFlag = false;
	    }
	  else
	    {
	      cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	    }
	}
    }
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

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
  char* OutputName = new char [512 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputName, "fermions_sphere_%s_n_%d_2s_%d_lz_%d.0.vec", Manager.GetString("interaction-name"), 
	   NbrParticles, LzMax, TotalLz);
  if (FermionFlag == true)
    {
      FermionOnSphere* U1Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
      if (SU2SymmetryFlag == true)
	{
	  FermionOnSphereWithSpin* Space;
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
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  RealVector OutputState = Space->ForgeU1FromSU2(State, *U1Space);
	  OutputState.WriteVector(OutputName);
	  delete Space;
	}
      if (SU3SymmetryFlag == true)
	{
	  FermionOnSphereWithSU3Spin* Space;
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
	  RealVector OutputState = Space->ForgeU1FromSU3(State, *U1Space);
	  OutputState.WriteVector(OutputName);
	  delete Space;
	}
      if (SU4SymmetryFlag == true)
	{
	  FermionOnSphereWithSU4Spin* Space;
#ifdef __64_BITS__
	  if (LzMax <= 15)
#else
	    if (LzMax <= 7)
#endif
	      {
		Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, LzMax, TotalSz, TotalIsoSz, TotalEntanglement, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  RealVector OutputState = Space->ForgeU1FromSU4(State, *U1Space);
	  OutputState.WriteVector(OutputName);
	  delete Space;
	}
      delete U1Space;
    }
  return 0;
}

