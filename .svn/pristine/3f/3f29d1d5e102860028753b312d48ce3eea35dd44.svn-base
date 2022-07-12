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
  OptionManager Manager ("FQHESphereFermionsWithSpinAllSzToSz" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state");
  (*SystemGroup) += new SingleIntegerOption  ('s', "sz", "Sz sector of the full Hilbert space on which the state has to be projected", 0);
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinAllSzToSz -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state " << endl << "see man page for option syntax or type FQHESphereFermionsWithSpinAllSzToSz -h" << endl;
      return -1;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int Sz = Manager.GetInteger("sz");
  int TotalLz = 0;

  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  bool FermionFlag = true;

  char* StateFileName = Manager.GetString("state");

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, LzMax, TotalLz, FermionFlag) == false)
    {
      cout << "error while retrieving system informations from file name " << Manager.GetString("state") << endl;
      return -1;
    }

  cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << endl;

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
  char* OldExtension = new char [512];
  sprintf (OldExtension, "lz_%d.0.vec", TotalLz);  
  char* NewExtension = new char [512];
  sprintf (OldExtension, "sz_%d_lz_%d.0.vec", Sz, LzMax);  
  char* OutputName = ReplaceExtensionToFileName(Manager.GetString("interaction-name"), OldExtension, NewExtension);

  if (FermionFlag == true)
    {
      FermionOnSphereWithSpin* SzSpace = 0;
#ifdef __64_BITS__
      if (LzMax <= 63)
#else
	if (LzMax <= 31)
#endif
	  {
	    SzSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, Sz, MemorySpace);
	  }
	else
	  {
	    cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	    return -1;
	  }	
      
      FermionOnSphereWithSpinAllSz* Space = 0;
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
      
      RealVector OutputState = SzSpace->ConvertFromNbodyBasis(State, *Space);
      OutputState.WriteVector(OutputName);
      
      delete Space;
      delete SzSpace;
    }
  return 0;
}

