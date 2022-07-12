#include "config.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"


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
  OptionManager Manager ("FQHELatticeMultipleU1ToU1" , "0.01");
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
      cout << "see man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if (Manager.GetString("state-2") == 0)
    {
      cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-2")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }
  bool HardCore = false;
  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalLz1 = 0;
  int Lx=0;
  int Ly=0;
  bool Statistics = true;
  double TmpI=1.0;
  if (FQHEOnLatticeFindSystemInfoFromFileName(Manager.GetString("state-1"), NbrParticles1,Lx,Ly,TmpI, NbrFluxQuanta1, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalLz2 = 0;
  Statistics = true;
  
  if (FQHEOnLatticeFindSystemInfoFromFileName(Manager.GetString("state-2"), NbrParticles2,Lx,Ly,TmpI,  NbrFluxQuanta2, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
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

  BosonOnLattice * Space1 = 0;

   Space1 = new BosonOnLattice(NbrParticles1, Lx, Ly, NbrFluxQuanta1);	       
   
  
  BosonOnLattice* Space2 = 0;
  
      Space2 = new BosonOnLattice(NbrParticles2, Lx, Ly, NbrFluxQuanta2); 
  
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
	  if (Manager.GetBoolean("haldane-output") == true)
	    sprintf (OutputFileName, "bosons_haldane_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	  else
	    sprintf (OutputFileName, "bosons_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	}
      else
	{
	  if (Manager.GetBoolean("haldane-output") == true)
	    sprintf (OutputFileName, "bosons_haldane_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	  else
	    sprintf (OutputFileName, "bosons_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
	}
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

  BosonOnLattice* TargetSpace = 0;
  
  TargetSpace = new BosonOnLattice(NbrParticles1 + NbrParticles2, Lx, Ly, NbrFluxQuanta2);    

  ComplexVector OutputState(TargetSpace->GetHilbertSpaceDimension(),true);
  
   TargetSpace->SymmetrizeU1U1State (OutputState,State1 , State2, Space1 , Space2 , Manager.GetBoolean("unnormalized-basis") ,0ul,(unsigned long)Space1->GetHilbertSpaceDimension());
    
	
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }
  delete Space1;
  delete Space2;
  delete TargetSpace;
}

