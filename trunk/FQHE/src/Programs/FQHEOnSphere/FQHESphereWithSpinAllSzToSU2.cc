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
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"


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
  OptionManager Manager ("FQHESphereFermionsAllSzToSU2" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file that contains state from the tunneling space");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-sz", "what is the desired Sz value of the output state?", 0);
  (*SystemGroup) += new SingleDoubleOption  ('t', "tunneling-amp", "tunneling amplitude", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsAllSzToSU2 -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state " << endl << "see man page for option syntax or type FQHESphereFermionsAllSzToSU2 -h" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  int PairParity = Manager.GetInteger("pair-parity");
  bool SzSymmetrizedBasis;
  bool SzMinusParity;  
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool LzMinusParity = Manager.GetBoolean("minus-lzparity");
  double tunneling = Manager.GetDouble("tunneling-amp");
  bool FermionFlag = false;

  char* StateFileName = Manager.GetString("state");
 
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;

  int TmpTotalSz=-1;
  if (NbrParticles==0)
    if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, TmpTotalSz, SzSymmetrizedBasis, SzMinusParity, 
							     LzSymmetrizedBasis, LzMinusParity, FermionFlag) == false)
      {
	return -1;
      }
  cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << "  FermionFlag="<< FermionFlag <<endl;
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	{
	  if ((strcmp ("bosons", Manager.GetString("statistics")) == 0))
	    {
	      FermionFlag = false;
	    }
	  else
	    {
	      cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	    }
	}
    }
  if (NbrParticles==0)
    {
      cout<<"Please provide the number of particles!"<<endl;
      exit(0);
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

  unsigned long MemorySpace = 9l << 20;
  char* OutputName;
  if (FermionFlag==true)
    {
      OutputName = new char [512 + strlen(Manager.GetString("interaction-name"))];
      sprintf (OutputName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_t_%f_lz_%d.0.vec", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalSz, tunneling, TotalLz);
    }
  else
    {
      char *Insertion = new char[10];
      sprintf(Insertion,"_sz_%d",TotalSz);
      OutputName=AddSegmentInFileName(Manager.GetString("state"), Insertion, "_2s_", true);
      delete [] Insertion;
    }
      
  if (FermionFlag == true)
    {
      if (LzSymmetrizedBasis == false)
	{
	  FermionOnSphereWithSpin* SU2Space = 0;
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		SU2Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }

	  FermionOnSphereWithSpinAllSz* Space;
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
    
	  RealVector OutputState = Space->ForgeSU2FromTunneling(State, *SU2Space, TotalSz);
	  OutputState.WriteVector(OutputName);	
	  delete Space;
	  delete SU2Space;
	}
      else
	{
	  FermionOnSphereWithSpinLzSymmetry* SU2Space = 0;
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		SU2Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	

	  FermionOnSphereWithSpinAllSzLzSymmetry* Space;
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpinAllSzLzSymmetry(NbrParticles, LzMax, LzMinusParity, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
    
	  RealVector OutputState = Space->ForgeSU2FromTunneling(State, *SU2Space, TotalSz);
	  OutputState.WriteVector(OutputName);	
	  delete Space;
	  delete SU2Space;
	}
    }
  else
    {
      if (LzSymmetrizedBasis == false)
	{
	  BosonOnSphereWithSpin* SU2Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);

	  BosonOnSphereWithSpinAllSz* Space;
	  
	  if ( PairParity >=0 ) 
	    Space = new BosonOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, PairParity, MemorySpace);
	  else
	    Space = new BosonOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);
    
	  RealVector OutputState = Space->ForgeSU2FromTunneling(State, *SU2Space, TotalSz);
	  OutputState.WriteVector(OutputName);	
	  delete Space;
	  delete SU2Space;
	}
      else
	{
	  cout << "Lz-symmetrized states not available for Bosons with Spin."<<endl;
	  return -1;
	}
    }


  return 0;
}

