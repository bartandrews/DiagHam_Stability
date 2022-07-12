#include "config.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnCP2.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"

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
  OptionManager Manager ("FQHESphereMultipleCP2ToCP2" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
	
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('1', "state-1", "vector file that corresponds to the first component");
  (*SystemGroup) += new SingleStringOption  ('2', "state-2", "vector file that corresponds to the second component");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file names");
  (*OutputGroup) += new BooleanOption  ('u', "unnormalized-basis", "indicates that calculations and data are in the unnormalized basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMultipleCP2ToCP2 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHESphereMultipleCP2ToCP2 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if (Manager.GetString("state-2") == 0)
    {
      cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHESphereMultipleCP2ToCP2 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-2")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }

  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalTz1 = 0;
  int TotalY1 = 0;
  bool TzSymmetry1 = false;
  bool MinusTzSymmetry1 = false;
  bool TzZ3Symmetry1 = false;
  
  bool Statistics = true;
  if (FQHEOnCP2FindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						   NbrParticles1, NbrFluxQuanta1, TotalTz1, TotalY1, TzSymmetry1, MinusTzSymmetry1, TzZ3Symmetry1, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalTz2 = 0;
  int TotalY2 = 0;
  bool TzSymmetry2 = false;
  bool MinusTzSymmetry2 = false;
  bool TzZ3Symmetry2 = false;
  Statistics = true;
  if (FQHEOnCP2FindSystemInfoFromVectorFileName(Manager.GetString("state-2"),
						   NbrParticles2, NbrFluxQuanta2, TotalTz2, TotalY2, TzSymmetry2, MinusTzSymmetry2, TzZ3Symmetry2, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
      return -1;
    }
  /*  if (NbrParticles1 != NbrParticles2)
      {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of particles" << endl;
      return -1;
      }*/
  if (NbrFluxQuanta2 != NbrFluxQuanta1)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
      return -1;
    }
  

  RealVector State1;
  if (State1.ReadVector (Manager.GetString("state-1")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
      return -1;      
    }
  RealVector State2;
  if (State2.ReadVector (Manager.GetString("state-2")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
      return -1;      
    }

  BosonOnCP2* Space1 = new BosonOnCP2(NbrParticles1, NbrFluxQuanta1, TotalTz1, TotalY1);
  BosonOnCP2* Space2 = new BosonOnCP2(NbrParticles2, NbrFluxQuanta2, TotalTz2, TotalY2);
  
  
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [512];
      if (Manager.GetBoolean("unnormalized-basis") == false)
	{
	  sprintf (OutputFileName, "bosons_cp2_delta_symmetrized_n_%d_2s_%d_tz_%d_y_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalTz1 + TotalTz2), (TotalY1 + TotalY2));
	}
      else
	{
	  sprintf (OutputFileName, "bosons_cp2_delta_unnormalized_symmetrized_n_%d_2s_%d_tz_%d_y_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalTz1 + TotalTz2), (TotalY1 + TotalY2));
	}
    }


  BosonOnCP2* TargetSpace = new BosonOnCP2(NbrParticles1 + NbrParticles2, NbrFluxQuanta1, TotalTz1 + TotalTz2, TotalY1 + TotalY2);
  
    
  RealVector OutputState = TargetSpace->SymmetrizeU1U1State (State1 , State2, Space1 , Space2 , Manager.GetBoolean("unnormalized-basis") , Architecture.GetArchitecture());
    
	
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }
  delete Space1;
  delete Space2;
  delete TargetSpace;
}

