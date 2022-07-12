#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceFixedParity.h"

#include "Hamiltonian/ParticleOnLatticeRealSpacePairingHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelKitaevChain.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  cout.precision(14);
  OptionManager Manager ("MultipleKitaevChainsToSpin1" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('1', "state-1", "name of the first state");
  (*SystemGroup) += new SingleStringOption  ('2', "state-2", "name of the second state");
  (*SystemGroup) += new SingleStringOption  ('3', "state-3", "name of the third state");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type MultipleKitaevChainsToSpin1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  
  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, --state-1 is not set" << endl;
    }
  if (Manager.GetString("state-2") == 0)
    {
      cout << "error, --state-1 is not set" << endl;
    }
  if (Manager.GetString("state-3") == 0)
    {
      cout << "error, --state-1 is not set" << endl;
    }
  int NbrSites = 0; 
  int FermionParity1 = 0;
  int FermionParity2 = 0;
  int FermionParity3 = 0;
  bool Statistics = true;

  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("state-1"), NbrSites, Statistics, FermionParity1) == false)
    {
      cout << "error, can't guess the system information from " << Manager.GetString("state-1") << endl;
    }
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("state-2"), NbrSites, Statistics, FermionParity2) == false)
    {
      cout << "error, can't guess the system information from " << Manager.GetString("state-2") << endl;
    }
  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("state-3"), NbrSites, Statistics, FermionParity3) == false)
    {
      cout << "error, can't guess the system information from " << Manager.GetString("state-3") << endl;
    }


  ComplexVector State1;
  if (State1.ReadVector(Manager.GetString("state-1")) == false)
    {
      cout << "can't read " << Manager.GetString("state-1") <<  endl;
    }
  ComplexVector State2;
  if (State2.ReadVector(Manager.GetString("state-2")) == false)
    {
      cout << "can't read " << Manager.GetString("state-2") <<  endl;
    }
  ComplexVector State3;
  if (State3.ReadVector(Manager.GetString("state-3")) == false)
    {
      cout << "can't read " << Manager.GetString("state-3") <<  endl;
    }

  
  FermionOnLatticeRealSpaceFixedParity* Space1 = new FermionOnLatticeRealSpaceFixedParity(NbrSites, FermionParity1);
  FermionOnLatticeRealSpaceFixedParity* Space2 = new FermionOnLatticeRealSpaceFixedParity(NbrSites, FermionParity2);
  FermionOnLatticeRealSpaceFixedParity* Space3 = new FermionOnLatticeRealSpaceFixedParity(NbrSites, FermionParity3);
  Space1->TripleTensorProductAndGutzwillerProjection(State1, Space1, State2, Space2, State3, Space3);

  return 0;
}
