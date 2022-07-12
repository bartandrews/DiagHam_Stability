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
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"

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
  OptionManager Manager ("FQHESphereU1ToTwoLandauLevels" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('u', "up-state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleStringOption  ('d', "down-state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "up-nbrparticles", "number of particles with spin up (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "down-nbrparticles", "number of particleswith spin down (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lzmaxup", "twice the maximum momentum for a single particle with spin up (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lzmaxdown", "twice the maximum momentum for a single particle with spin down (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "up-totallz", "twice the total lz value for spin up particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "down-totallz", "twice the total lz value for spin down particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");


  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereU1ToTwoLandauLevels -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("up-state") == 0)
    {
      cout << "no input state for the spin up part" << endl << "see man page for option syntax or type FQHESphereU1ToTwoLandauLevels -h" << endl;
      return -1;
    }
  if(Manager.GetString("down-state") == 0)
    {
      cout << "no input state for the spin down part" << endl << "see man page for option syntax or type FQHESphereU1ToTwoLandauLevels -h" << endl;
      return -1;
    }

  int UpNbrParticles = Manager.GetInteger("up-nbrparticles");
  int DownNbrParticles = Manager.GetInteger("down-nbrparticles");
  int LzMaxUp = Manager.GetInteger("lzmaxup");
  int LzMaxDown = Manager.GetInteger("lzmaxdown");
  int UpTotalLz = Manager.GetInteger("up-totallz");
  int DownTotalLz = Manager.GetInteger("down-totallz");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (UpNbrParticles == 0)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("up-state"), UpNbrParticles, LzMaxUp, UpTotalLz, FermionFlag) == false)
	{
	  return -1;
	}      
    }
  if (DownNbrParticles == 0)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("down-state"), DownNbrParticles, LzMaxDown, DownTotalLz, FermionFlag) == false)
	{
	  return -1;
	}
      
    }
  if (Manager.GetBoolean("statistics") != 0)
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

  int UpParity = UpTotalLz & 1;
  if (UpParity != ((UpNbrParticles * LzMaxUp) & 1))
    {
      cout << "Lz and (UpNbrParticles * LzMaxUp) must have the parity" << endl;
      return -1;           
    }
  int DownParity = DownTotalLz & 1;
  if (DownParity != ((DownNbrParticles * LzMaxDown) & 1))
    {
      cout << "Lz and (DownNbrParticles * LzMaxDown) must have the parity" << endl;
      return -1;           
    }

  char* UpStateFileName = Manager.GetString("up-state");
  if (IsFile(UpStateFileName) == false)
    {
      cout << "state " << UpStateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  RealVector UpState;
  if (UpState.ReadVector(UpStateFileName) == false)
    {
      cout << "error while reading " << UpStateFileName << endl;
      return -1;
    }
  char* DownStateFileName = Manager.GetString("down-state");
  if (IsFile(DownStateFileName) == false)
    {
      cout << "state " << DownStateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  RealVector DownState;
  if (DownState.ReadVector(DownStateFileName) == false)
    {
      cout << "error while reading " << DownStateFileName << endl;
      return -1;
    }


  long MemorySpace = 9l << 20;
  ParticleOnSphere* UpSpace = 0;
  ParticleOnSphere* DownSpace = 0;
  ParticleOnSphereWithSpin* SU2Space = 0;
  if (FermionFlag == true)
    {
      UpSpace = new FermionOnSphere(UpNbrParticles, UpTotalLz, LzMaxUp);
      DownSpace = new FermionOnSphere(DownNbrParticles, DownTotalLz, LzMaxDown);
      SU2Space = new FermionOnSphereTwoLandauLevels(UpNbrParticles + DownNbrParticles, UpTotalLz + DownTotalLz, LzMaxUp, LzMaxDown);
    }
  else
    {
      cout << "bosonic statistics is not yet supported" << endl;
    }
  
  if (UpSpace->GetHilbertSpaceDimension() != UpState.GetVectorDimension())
    {
      cout << "dimension mismatch between the up state (" << UpState.GetVectorDimension() << ") and the Hilbert space (" << UpSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  if (DownSpace->GetHilbertSpaceDimension() != DownState.GetVectorDimension())
    {
      cout << "dimension mismatch between the down state (" << DownState.GetVectorDimension() << ") and the Hilbert space (" << DownSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  RealVector OutputState = ((FermionOnSphereWithSpin*) SU2Space)->ForgeSU2FromU1(UpState, *(FermionOnSphere*)UpSpace, DownState, *(FermionOnSphere*)DownSpace);
  char* OutputName = new char [512 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputName, "fermions_sphere_ll0_ll_%d_%s_n_%d_2s_%d_lz_%d.0.vec", ((LzMaxUp - LzMaxDown) >> 1),
	   Manager.GetString("interaction-name"), 
	   (UpNbrParticles + DownNbrParticles), LzMaxDown, (UpTotalLz + DownTotalLz));

  OutputState.WriteVector(OutputName);
  delete UpSpace;
  delete DownSpace;
  delete SU2Space;
  return 0;
}

