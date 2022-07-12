#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/Potts3ChainHamiltonian.h"
#include "Hamiltonian/Potts3ChainHamiltonianWithTranslations.h"
#include "Hamiltonian/Potts3ChainNaturalBoundaryTermHamiltonian.h"

#include "HilbertSpace/Potts3Chain.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

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
  OptionManager Manager ("Potts3ChainModel" , "0.01");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-q", "initial q sector that has to computed (can be either 0, 1 or 2)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-q", "number of q value to evaluate (0 for all q sectors)", 0);
  (*SystemGroup) += new  SingleDoubleOption ('j', "coupling", "magnitude of the nearest nieighbor coupling", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('f', "flip", "magnitude of the on-site flip term", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "phi-coupling", "phase (in 2 \\pi units) of the nearest nieighbor coupling", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "phi-flip", "phase (in 2 \\pi units) of the on-site flip term", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  SingleIntegerOption ('b', "boundary-conditions", "type of boundary conditions (0 for 1, 1 for exp(2i \\pi / 3) and -1 for exp(-2i \\pi / 3)", 0);
  (*SystemGroup) += new  BooleanOption  ('\n', "natural-boundaryterms", "use a natural boundary term instead of the translation invariant boundary term");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "boundaryterm-order", "perturbation order for the edge mode development involved in the natural boundary term", 0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "filter-0", "first factor coming from the filter function when using the first order correction", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "filter-1", "second factor coming from the filter function when using the first order correction", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "filter-2", "third factor coming from the filter function when using the first order correction", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "use-momentum", "use the momentum quantum number");
  (*SystemGroup) += new  SingleIntegerOption  ('k', "k-sector", "look at a given momentum sector (-1 if all momentum sectors have to be computed)", -1);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Potts3ChainModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  double JValue = Manager.GetDouble("coupling");
  double FValue = Manager.GetDouble("flip");
  double PhiJ = Manager.GetDouble("phi-coupling");
  double PhiF = Manager.GetDouble("phi-flip");
  int BoundaryCondition = Manager.GetInteger("boundary-conditions");
  bool UseMomentumFlag = Manager.GetBoolean("use-momentum");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  if (Manager.GetBoolean("periodic") == false)
    {
      UseMomentumFlag = false;
      sprintf (OutputFileName, "potts3_openchain_j_%.6f_phij_%.6f_f_%.6f_phif_%.6f_n_%d", JValue, PhiJ, FValue, PhiF, NbrSpins);
      sprintf (CommentLine, " open potts 3 chain with %d sites and J=%.6f, PhiJ=%.6f, F=%.6f, PhiF=%.6f \n# Q", NbrSpins, JValue, PhiJ, FValue, PhiF);
    }
  else
    {
      if (Manager.GetBoolean("natural-boundaryterms") == false)
	{
	  sprintf (OutputFileName, "potts3_closedchain_j_%.6f_phij_%.6f_f_%.6f_phif_%.6f_b_%d_n_%d", JValue, PhiJ, FValue, PhiF, BoundaryCondition, NbrSpins);
	  if (UseMomentumFlag == false)
	    {
	      sprintf (CommentLine, " close potts 3 chain with %d sites  and J=%.6f, PhiJ=%.6f, F=%.6f, PhiF=%.6f, B=%d \n#\n# Q", NbrSpins, JValue, PhiJ, FValue, PhiF, BoundaryCondition);
	    }
	  else
	    {
	      sprintf (CommentLine, " close potts 3 chain with %d sites  and J=%.6f, PhiJ=%.6f, F=%.6f, PhiF=%.6f, B=%d \n#\n# Q K", NbrSpins, JValue, PhiJ, FValue, PhiF, BoundaryCondition);
	    }
	}
      else
	{
	  sprintf (OutputFileName, "potts3_closedchain_naturalboundaryterms_%d_j_%.6f_phij_%.6f_f_%.6f_phif_%.6f_b_%d_n_%d", (int) Manager.GetInteger("boundaryterm-order"), JValue, PhiJ, FValue, PhiF, BoundaryCondition, NbrSpins);
	  sprintf (CommentLine, " close potts 3 chain an natural boundray conditions at the order %d with %d sites  and J=%.6f, PhiJ=%.6f, F=%.6f, PhiF=%.6f, B=%d \n#\n# Q", (int) Manager.GetInteger("boundaryterm-order"), NbrSpins, JValue, PhiJ, FValue, PhiF, BoundaryCondition);
	}
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  int MaxQValue = 2;
  int InitialQValue = 0;
  if (Manager.GetInteger("initial-q") >= 0)
    {
      InitialQValue = Manager.GetInteger("initial-q") % 3;
    }
  if (Manager.GetInteger("nbr-q") > 0)
    {
      MaxQValue = InitialQValue + Manager.GetInteger("nbr-q") - 1;
      if (MaxQValue >= 3)
	MaxQValue = 2;
    }
  bool FirstRun = true;

  if (UseMomentumFlag == true)
    {
      for (; InitialQValue <= MaxQValue; ++InitialQValue)
	{
	  int NbrMomentumSector = NbrSpins;
	  int MomentumSector = 0;
	  if ((Manager.GetInteger("k-sector") >= 0) && (Manager.GetInteger("k-sector") < NbrMomentumSector))
	    NbrMomentumSector = Manager.GetInteger("k-sector");
	  for (; MomentumSector < NbrMomentumSector; ++MomentumSector)
	    {
	      Potts3ChainWithTranslations* Chain = new Potts3ChainWithTranslations (NbrSpins, InitialQValue, MomentumSector, 1000000);      
	      Potts3ChainHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, MomentumSector, JValue, PhiJ, FValue, PhiF, BoundaryCondition);
	      char* TmpQString = new char[64];
	      sprintf (TmpQString, "%d %d", InitialQValue, MomentumSector);
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_q_%d_k_%d", OutputFileName, InitialQValue, MomentumSector);
	      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
					  FirstRun, TmpEigenstateString);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      FirstRun = false;
	      delete Chain;
	      delete[] TmpQString;
	    }
	}
    }
  else
    {
      for (; InitialQValue <= MaxQValue; ++InitialQValue)
	{
	  Potts3Chain* Chain = new Potts3Chain (NbrSpins, InitialQValue, 1000000);      
	  Potts3ChainHamiltonian* Hamiltonian = 0;
	  if (Manager.GetBoolean("natural-boundaryterms") == false)
	    {
	      Hamiltonian = new Potts3ChainHamiltonian (Chain, NbrSpins, JValue, PhiJ, FValue, PhiF, Manager.GetBoolean("periodic"), 
							BoundaryCondition, ((long) Manager.GetInteger("memory")) << 20);
	    }
	  else
	    {
	      Hamiltonian = new Potts3ChainNaturalBoundaryTermHamiltonian(Chain, NbrSpins, JValue, PhiJ, FValue, PhiF, BoundaryCondition, 
									  (int) Manager.GetInteger("boundaryterm-order"), 
									  Manager.GetDouble("filter-0"), Manager.GetDouble("filter-1"), Manager.GetDouble("filter-2"),
									  ((long) Manager.GetInteger("memory")) << 20);
	    }
	  char* TmpQString = new char[64];
	  sprintf (TmpQString, "%d", InitialQValue);
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpEigenstateString, "%s_q_%d", OutputFileName, InitialQValue);
	  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
				      FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Chain;
	  delete[] TmpQString;
	  delete Hamiltonian;
	}
    }
  return 0;
}
