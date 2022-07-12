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
  OptionManager Manager ("FQHESphereLzMinusLz.cc" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
  

  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (changing Lz)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLzToMinusLz -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereFermionsParticleHoleSymmetrize -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("nbr-flux"); 
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  bool Statistics = true;
  int TotalLz = 0;

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-file"),
						  NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-file") << endl;
      return -1;
    }
    
  RealVector State;
  if (State.ReadVector (Manager.GetString("input-file")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-file") << endl;
      return -1;      
    }
  
  RealVector MinusLzState;
  if (Manager.GetBoolean("haldane") == false)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
#else
	if (LzMax <= 31)
#endif
	  {
	    FermionOnSphere InputSpace (NbrParticles, TotalLz, LzMax, MemorySpace);
	    FermionOnSphere OutputSpace (NbrParticles, -TotalLz, LzMax, MemorySpace);
	    MinusLzState = InputSpace.GetLzSymmetricVector(&OutputSpace,State);
	  }
    }
  else
    {
      int* ReferenceState = 0;
      if (Manager.GetString("reference-file") == 0)
	{
	  ReferenceState = new int[LzMax + 1];
	  for (int i = 0; i <= LzMax; ++i)
	    ReferenceState[i] = 0;
	  if (strcasecmp(Manager.GetString("reference-state"), "laughlin") == 0)
	    for (int i = 0; i <= LzMax; i += 3)
	      ReferenceState[i] = 1;
	  else
	    if (strcasecmp(Manager.GetString("reference-state"), "pfaffian") == 0)
	      for (int i = 0; i <= LzMax; i += 4)
		{
		  ReferenceState[i] = 1;
		  ReferenceState[i + 1] = 1;
		}
	    else
	      if (strcasecmp(Manager.GetString("reference-state"), "readrezayi3") == 0)
		for (int i = 0; i <= LzMax; i += 5)
		  {
		    ReferenceState[i] = 1;
		    ReferenceState[i + 1] = 1;
		    ReferenceState[i + 2] = 1;
		  }
	      else
		{
		  cout << "unknown reference state " << Manager.GetString("reference-state") << endl;
		  return -1;
		}
	}
      else
	{
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
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
	      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	      return -1;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return -1;     
	    }
	}
      if (Manager.GetString("load-hilbert") != 0)
	{
	  FermionOnSphereHaldaneBasis InputSpace (Manager.GetString("load-hilbert"), MemorySpace);
	  FermionOnSphere OutputSpace (NbrParticles, -TotalLz, LzMax, MemorySpace);
	  MinusLzState = InputSpace.GetLzSymmetricVector(&OutputSpace, State);
	}
      else
	{
	  FermionOnSphereHaldaneBasis InputSpace (NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
	  FermionOnSphere OutputSpace (NbrParticles, -TotalLz, LzMax, MemorySpace);
	  MinusLzState = InputSpace.GetLzSymmetricVector(&OutputSpace, State);
	}
    }
  
  if (Manager.GetString("output-file") == 0)
    {
      char* InputFileName = Manager.GetString("input-file"); 
      char* TagPosition = strcasestr(InputFileName, "fermions_");
      if (TagPosition != InputFileName)
	{
	  cout << "no default output name can be built from " << InputFileName << endl;
	  return -1;
	}
      char* OutputFileName = new char [strlen(InputFileName) + 1];
      strcpy (OutputFileName, "fermions_");
      TagPosition = strcasestr(InputFileName, "_lz_");
      strncpy (OutputFileName + 9, InputFileName + 9, (TagPosition - InputFileName - 9));
      sprintf (OutputFileName + (TagPosition - InputFileName), "_lz_%d", -TotalLz);
      long TmpPos = strlen (OutputFileName);
      TagPosition = strcasestr(InputFileName, ".");
      if (TagPosition == 0)
        {
          cout << "no default output name can be built from " << InputFileName << endl;
          return -1;
        }      
      strcpy (OutputFileName + TmpPos, TagPosition);
      MinusLzState.WriteVector(OutputFileName);
    }
  else
    {
      MinusLzState.WriteVector(Manager.GetString("output-file"));
    }
  return 0;
}
