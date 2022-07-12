#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereFermionsParticleHoleSymmetrize.cc" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
//  (*SystemGroup) += new BooleanOption  ('\n', "v2", "use a different scheme to compute particle/hole conjugason");

  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (changing N and 2S)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsParticleHoleSymmetrize.cc -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-file"])->GetString() == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereFermionsParticleHoleSymmetrize -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int LzMax = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  bool Statistics = true;
  int TotalLz = 0;
  bool BergholtzFlag = false;//Manager.GetBoolean("v2");

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(),
						  NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
  int NbrHoles = LzMax + 1 - NbrParticles;

  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;      
    }

  RealVector HoleState;
  if (((BooleanOption*) Manager["haldane"])->GetBoolean() == false)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
#else
      if (LzMax <= 31)
#endif
	{
	  if ((SymmetrizedBasis == false) || (TotalLz != 0))
	    {
	      FermionOnSphere InputSpace (NbrParticles, TotalLz, LzMax, MemorySpace);
	      FermionOnSphere OutputSpace (NbrHoles, -TotalLz, LzMax, MemorySpace);
	      HoleState = InputSpace.ParticleHoleSymmetrize(State, OutputSpace);
	    }
	  else
	    {
	      FermionOnSphereSymmetricBasis InputSpace(NbrParticles, LzMax, MemorySpace);
	      FermionOnSphereSymmetricBasis OutputSpace (NbrHoles, LzMax, MemorySpace);
	      HoleState = InputSpace.ParticleHoleSymmetrize(State, OutputSpace);
	    }
	}
    }
  else
    {
      int* ReferenceState = 0;
      if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
	{
	  cout << "a --reference-file is needed when using the --haldane option" << endl;
	  return -1;
	}
      else
	{
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
	      return -1;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return -1;     
	    }
	}
      if (SymmetrizedBasis == false)
	{
	  if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
	    {
	      FermionOnSphereHaldaneBasis InputSpace (((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
	      FermionOnSphere OutputSpace (NbrHoles, TotalLz, LzMax, MemorySpace);
	      HoleState = InputSpace.ParticleHoleSymmetrize(State, OutputSpace);
	    }
	  else
	    {
	      FermionOnSphereHaldaneBasis InputSpace (NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
	      FermionOnSphere OutputSpace (NbrHoles, -TotalLz, LzMax, MemorySpace);
	      HoleState = InputSpace.ParticleHoleSymmetrize(State, OutputSpace);
	    }
	}
      else
	{
	  if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
	    {
	      FermionOnSphereHaldaneSymmetricBasis InputSpace(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
	      FermionOnSphereSymmetricBasis OutputSpace (NbrHoles, LzMax, MemorySpace);
	      HoleState = InputSpace.ParticleHoleSymmetrize(State, OutputSpace);
	    }
	  else
	    {
	      FermionOnSphereHaldaneSymmetricBasis InputSpace (NbrParticles, LzMax, ReferenceState, MemorySpace);
	      FermionOnSphereSymmetricBasis OutputSpace (NbrHoles, LzMax, MemorySpace);
	      HoleState = InputSpace.ParticleHoleSymmetrize(State, OutputSpace);
	    }
	}
    }


  if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
    {
      char* InputFileName = ((SingleStringOption*) Manager["input-file"])->GetString(); 
      char* TagPosition = strcasestr(InputFileName, "fermions_");
      if (TagPosition != InputFileName)
	{
	  cout << "no default output name can be built from " << InputFileName << endl;
	  return -1;
	}
      char* OutputFileName = new char [strlen(InputFileName) + 8];
      strcpy (OutputFileName, "fermions_holes_");
      TagPosition = strcasestr(InputFileName, "_n_");
      if (TagPosition == 0)
	{
	  cout << "no default output name can be built from " << InputFileName << endl;
	  return -1;
	} 
      strncpy (OutputFileName + 15, InputFileName + 9, (TagPosition - InputFileName - 9));
      sprintf (OutputFileName + 6 + (TagPosition - InputFileName), "_n_%d_2s_%d_lz_%d.0.vec", NbrHoles, LzMax, -TotalLz);
      HoleState.WriteVector(OutputFileName);
    }
  else
    {
      HoleState.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString());
    }
  return 0;
}
