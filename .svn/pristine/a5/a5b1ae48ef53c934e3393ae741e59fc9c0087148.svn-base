#include "Operator/FermionParityOperator.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("Z2InteractingChainComputeInvariant" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('p', "periodic-state", "name of the file that contains the ground state with periodic boundary conditions");
  (*SystemGroup) += new  SingleStringOption ('a', "antiperiodic-state", "name of the file that contains the ground state with antiperiodic boundary conditions");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Z2InteractingChainComputeInvariant -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrSpinsPeriodic = 0;
  if (SpinFindSystemInfoFromFileName(Manager.GetString("periodic-state"), NbrSpinsPeriodic) == false)
    {
      return 0;
    }
  int NbrSpinsAntiPeriodic = 0;
  if (SpinFindSystemInfoFromFileName(Manager.GetString("antiperiodic-state"), NbrSpinsAntiPeriodic) == false)
    {
      return 0;
    }
  if (NbrSpinsPeriodic != NbrSpinsAntiPeriodic)
    {
      cout << "the number of spins is different for " << Manager.GetString("periodic-state") 
	   << " and " << Manager.GetString("antiperiodic-state") << endl;
      return 0;
    }

  RealVector PeriodicGroundState;
  if (PeriodicGroundState.ReadVector(Manager.GetString("periodic-state")) == false)
    {
      cout << "can't read " << Manager.GetString("periodic-state") << endl;
    }
  RealVector AntiPeriodicGroundState;
  if (AntiPeriodicGroundState.ReadVector(Manager.GetString("antiperiodic-state")) == false)
    {
      cout << "can't read " << Manager.GetString("antiperiodic-state") << endl;
    }

  Spin1_2Chain* Chain = new Spin1_2ChainFull (NbrSpinsPeriodic);

  FermionParityOperator Operator(Chain, NbrSpinsPeriodic);
  
  Complex Tmp1 = Operator.MatrixElement(PeriodicGroundState, PeriodicGroundState);
  cout << "<GS+|FP|GS+> = " << Tmp1.Re << endl;
  Complex Tmp2 = Operator.MatrixElement(AntiPeriodicGroundState, AntiPeriodicGroundState);
  cout << "<GS-|FP|GS-> = " << Tmp2.Re << endl;
  cout << "Z2 = " << (Tmp1.Re * Tmp2.Re) << endl;
  delete Chain;
  return 0;
}
