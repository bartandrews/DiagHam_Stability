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
  OptionManager Manager ("FQHESphereU1ToSU2" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "up-totallz", "twice the total lz value for spin up particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "down-totallz", "twice the total lz value for spin down particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereU1ToSU2 -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(((SingleStringOption*) Manager["up-state"])->GetString() == 0)
    {
      cout << "no input state for the spin up part" << endl << "see man page for option syntax or type FQHESphereU1ToSU2 -h" << endl;
      return -1;
    }
  if(((SingleStringOption*) Manager["down-state"])->GetString() == 0)
    {
      cout << "no input state for the spin down part" << endl << "see man page for option syntax or type FQHESphereU1ToSU2 -h" << endl;
      return -1;
    }

  int UpNbrParticles = ((SingleIntegerOption*) Manager["up-nbrparticles"])->GetInteger();
  int DownNbrParticles = ((SingleIntegerOption*) Manager["down-nbrparticles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int UpTotalLz = ((SingleIntegerOption*) Manager["up-totallz"])->GetInteger();
  int DownTotalLz = ((SingleIntegerOption*) Manager["down-totallz"])->GetInteger();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (UpNbrParticles==0)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["up-state"])->GetString(), UpNbrParticles, LzMax, UpTotalLz, FermionFlag) == false)
	{
	  return -1;
	}      
    }
  if (DownNbrParticles==0)
    {
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["down-state"])->GetString(), DownNbrParticles, LzMax, DownTotalLz, FermionFlag) == false)
	{
	  return -1;
	}
      
    }
  int SzTotal = UpNbrParticles - DownNbrParticles;
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
  int UpParity = UpTotalLz & 1;
  if (UpParity != ((UpNbrParticles * LzMax) & 1))
    {
      cout << "Lz and (UpNbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }
  int DownParity = DownTotalLz & 1;
  if (DownParity != ((DownNbrParticles * LzMax) & 1))
    {
      cout << "Lz and (DownNbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* UpStateFileName = ((SingleStringOption*) Manager["up-state"])->GetString();
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
  char* DownStateFileName = ((SingleStringOption*) Manager["down-state"])->GetString();
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
#ifdef  __64_BITS__
      if (((LzMax + UpNbrParticles - 1) < 31) && ((LzMax + DownNbrParticles - 1) < 31))
#else
	if (((LzMax + DownNbrParticles - 1) < 15)	&& ((LzMax + DownNbrParticles - 1) < 15))
#endif
	{
	  UpSpace = new FermionOnSphere(UpNbrParticles, UpTotalLz, LzMax);
	  DownSpace = new FermionOnSphere(DownNbrParticles, DownTotalLz, LzMax);
	  SU2Space = new FermionOnSphereWithSpin(UpNbrParticles + DownNbrParticles, UpTotalLz + DownTotalLz, LzMax, SzTotal);
	}
      else
	{
	  cout << "LzMax value exceeds integer precision" << endl;
	  return -1;
	}
    }
  else
    {
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
  char* OutputName = new char [512 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf (OutputName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d.0.vec", ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
	   (UpNbrParticles + DownNbrParticles), LzMax, SzTotal, (UpTotalLz + DownTotalLz));

  OutputState.WriteVector(OutputName);
  delete UpSpace;
  delete DownSpace;
  delete SU2Space;
  return 0;
}

