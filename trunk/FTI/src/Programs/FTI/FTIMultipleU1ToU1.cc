#include "config.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESquareLatticeSymmetrizeU1U1StateOperation.h"

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
  OptionManager Manager ("FQHETopInsulatorMultipleU1ToU1" , "0.01");
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
      cout << "see man page for option syntax or type FQHETopInsulatorMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHETopInsulatorMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if (Manager.GetString("state-2") == 0)
    {
      cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHETopInsulatorMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-2")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }

  int NbrParticles1 = 0; 
  int NbrSiteX1 = 0;
  int NbrSiteY1 = 0;
  int MomentumX1 = 0;
  int MomentumY1 = 0;
  bool Statistics = true;
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state-1"), NbrParticles1, NbrSiteX1, NbrSiteY1, MomentumX1, MomentumY1, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  int NbrParticles2 = 0; 
  int NbrSiteX2 = 0;
  int NbrSiteY2 = 0;
  int MomentumX2 = 0;
  int MomentumY2 = 0;
  Statistics = true;
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state-2"), NbrParticles2, NbrSiteX2, NbrSiteY2, MomentumX2, MomentumY2, Statistics) == false)

    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
      return -1;
    }
  if (NbrSiteX1 != NbrSiteX2)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same lattice geometry in the x-direction" << endl;
      return -1;
    }
  if (NbrSiteY1 != NbrSiteY2)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same lattice geometry in the y-direction" << endl;
      return -1;
    }


  ComplexVector State1;
  if (State1.ReadVector (Manager.GetString("state-1")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
      return -1;      
    }
  ComplexVector State2;
  if (State2.ReadVector (Manager.GetString("state-2")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
      return -1;      
    }

  BosonOnSquareLatticeMomentumSpace* Space1 = new BosonOnSquareLatticeMomentumSpace(NbrParticles1, NbrSiteX1, NbrSiteY1, MomentumX1, MomentumY1);	       
  
  BosonOnSquareLatticeMomentumSpace* Space2 = new BosonOnSquareLatticeMomentumSpace(NbrParticles2, NbrSiteX2, NbrSiteY2, MomentumX2, MomentumY2);	       
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
	  sprintf (OutputFileName, "bosons_symmetrized_n_%d_x_%d_y_%d_kx_%d_ky_%d.0.vec", NbrParticles1 + NbrParticles2, NbrSiteX1, NbrSiteY1, (MomentumX1+MomentumX2)%NbrSiteX1, (MomentumY1+MomentumY2)%NbrSiteY1);
	}
      else
	{
	  sprintf (OutputFileName, "bosons_unnormalized_symmetrized_n_%d_x_%d_y_%d_kx_%d_ky_%d.0.vec", NbrParticles1 + NbrParticles2, NbrSiteX1, NbrSiteY1, (MomentumX1+MomentumX2)%NbrSiteX1, (MomentumY1+MomentumY2)%NbrSiteY1);
	}
    }

  BosonOnSquareLatticeMomentumSpace* TargetSpace = new BosonOnSquareLatticeMomentumSpace(NbrParticles1 + NbrParticles2, NbrSiteX1, NbrSiteY1, (MomentumX1+MomentumX2)%NbrSiteX1, (MomentumY1+MomentumY2)%NbrSiteY1);	       

    
  ComplexVector OutputState = TargetSpace->SymmetrizeU1U1State (State1 , State2, Space1 , Space2 , Manager.GetBoolean("unnormalized-basis") , Architecture.GetArchitecture());
    
	
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }
  delete Space1;
  delete Space2;
  delete TargetSpace;
}

