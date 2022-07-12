#include "config.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

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
  OptionManager Manager ("FQHESphereWithSU2SpinMultipleSU2ToSU2" , "0.01");
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
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinMultipleSU2ToSU2 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("state-1") == 0) || (Manager.GetString("state-2") == 0))
    {
      cout << "error, two input files should be provided. See man page for option syntax or type FQHESphereWithSU2SpinMultipleSU2ToSU2 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if (IsFile(Manager.GetString("state-2")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }

  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalLz1 = 0;
  int TotalSz1 = 0;
  int LzSymmetry1 = 0;
  int SzSymmetry1 = 0;
  bool Statistics = true;
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalLz2 = 0;
  int TotalSz2 = 0;
  int LzSymmetry2 = 0;
  int SzSymmetry2 = 0;
  
    
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state-1"), NbrParticles1, NbrFluxQuanta1, TotalLz1, 
							   TotalSz1, LzSymmetry1, SzSymmetry1, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state-2"), NbrParticles2, NbrFluxQuanta2, TotalLz2, 
							   TotalSz2, LzSymmetry2, SzSymmetry2, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
      return -1;
    }

  if (NbrFluxQuanta2 != NbrFluxQuanta1)
    {
      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
      return -1;
    }
  
  RealVector State1;
  RealVector State2;
  if (State1.ReadVector (Manager.GetString("state-1")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
      return -1;      
    }
   
  if (State2.ReadVector (Manager.GetString("state-2")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
      return -1;      
    }

  ParticleOnSphereWithSpin* Space1 = 0;
  ParticleOnSphereWithSpin* Space2 = 0;
  if (Statistics == false)
    {
      Space1 = new BosonOnSphereWithSU2Spin(NbrParticles1, TotalLz1, NbrFluxQuanta1, TotalSz1);	       
    }
  else
    {
      Space1 = new FermionOnSphereWithSpin(NbrParticles1, TotalLz1, NbrFluxQuanta1, TotalSz1);	       
    }
            
  
  if (Statistics == false)
    {
      Space2 = new BosonOnSphereWithSU2Spin(NbrParticles2, TotalLz2, NbrFluxQuanta2, TotalSz2);	       
    }
  else
    {
      Space2 = new FermionOnSphereWithSpin(NbrParticles2, TotalLz2, NbrFluxQuanta2, TotalSz2);
    }
  
  char* OutputFileName = 0;
  int NbrFluxQuanta = NbrFluxQuanta1;
  
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [512];
      if (Statistics == false)
	{
	  sprintf (OutputFileName, "bosons_sphere_su2_symmetrized_n_%d_2s_%d_sz_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta, 
		   (TotalSz1 + TotalSz2), (TotalLz1 + TotalLz2));
	}
      else
	{
	  sprintf (OutputFileName, "fermions_sphere_su2_symmetrized_n_%d_2s_%d_sz_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta, 
		   (TotalSz1 + TotalSz2), (TotalLz1 + TotalLz2));
	}
    }
  
  
  ParticleOnSphereWithSpin* TargetSpace = 0;
  
  if (Statistics == false)
    {
      TargetSpace = new BosonOnSphereWithSU2Spin(NbrParticles1 + NbrParticles2, TotalLz1 + TotalLz2, NbrFluxQuanta2, TotalSz1 + TotalSz2);	       
    }
  else
    {
      TargetSpace = new FermionOnSphereWithSpin(NbrParticles1 + NbrParticles2, TotalLz1 + TotalLz2, NbrFluxQuanta2, TotalSz1 + TotalSz2);	  
    }

  RealVector OutputState = TargetSpace->SymmetrizeSU2SU2State (State1 , State2, Space1, Space2, false, Architecture.GetArchitecture());  
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }
  
  if (Space1 != 0)
    delete Space1;
  if (Space2 != 0)
    delete Space2;
  if (TargetSpace != 0)
    delete TargetSpace;
}

