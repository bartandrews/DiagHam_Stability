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
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsWithSpinAllSzSxValue" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file that contains the SU2WithTunneling state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('t', "tunneling-amp", "tunneling amplitude", 0.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinAllSzSxValue -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "no input state " << endl << "see man page for option syntax or type FQHESphereFermionsWithSpinAllSzSxValue -h" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  double tunneling = Manager.GetDouble("tunneling-amp");
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  bool FermionFlag = false;

  char* StateFileName = ((SingleStringOption*) Manager["state"])->GetString();
 
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;


  if (NbrParticles==0)
   {
    cout<<"Please provide the number of particles!"<<endl;
    exit(0);
   }

  //cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << endl;
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
  //char* OutputName = new char [512 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  //sprintf (OutputName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_t_%f_lz_%d.0.vec", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, TotalSz, tunneling, TotalLz);
  if (FermionFlag == true)
    {
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
    
	  double meanSzvalue = Space->MeanSzValue(State);
	  cout<< "Mean value <Sz> in the state: "<<endl;
	  cout<< meanSzvalue<<endl;

	  double meanSxvalue = Space->MeanSxValue(State);
	  cout<< "Mean value <Sx> in the state: "<<endl;
	  cout<< meanSxvalue<<endl;
	
	delete Space;
    }
  return 0;
}

