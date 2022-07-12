#include "Hamiltonian/SpinChainXYZHamiltonian.h"
#include "Hamiltonian/SpinChainXYZNaturalBoundaryTermHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFixedParity.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

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
  OptionManager Manager ("SpinChainXYZ" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 8);
  (*SystemGroup) += new  SingleDoubleOption('x', "jx-value", "coupling constant along the x axis", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('y', "jy-value", "coupling constant along the y axis", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('z', "jz-value", "coupling constant along the z axis", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('f', "h-value", "Zeeman strength along the z axis", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-q", "initial parity sector that has to computed (can be either 0, 1)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-q", "number of parity sectors to evaluate (0 for all parity sectors)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('b', "boundary-conditions", "boundary conditions (0 for open, 1 for periodic, -1 for antiperiodic)", 0);
  (*SystemGroup) += new  BooleanOption  ('\n', "no-parity", "do not take into account the parity when computing the spectrum");
  (*SystemGroup) += new  BooleanOption  ('\n', "natural-boundaryterms", "use a natural boundary term instead of the translation invariant boundary term");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "boundaryterm-order", "perturbation order for the edge mode development involved in the natural boundary term", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "jy-dominated", "the 0-th order of the natural boundary term is Jy dominated instead of being Jx dominated");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*OutputGroup) += new  SingleStringOption ('\n', "output-suffix", "apprend and extra suffix to the string describing the system in the output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainXYZ -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  double JxValue = Manager.GetDouble("jx-value");
  double JyValue = Manager.GetDouble("jy-value");
  double JzValue = Manager.GetDouble("jz-value");
  double HValue = Manager.GetDouble("h-value");
  int BValue = Manager.GetInteger("boundary-conditions");
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* FileNamePrefix = new char [256];
  if ((Manager.GetBoolean("natural-boundaryterms") == false) || (BValue == 0))
    {
      sprintf (FileNamePrefix, "spin_1_2");
    }
  else
    {
      if (Manager.GetBoolean("jy-dominated") == false)
	{
	  sprintf (FileNamePrefix, "spin_1_2_naturalboundaryterms_%d", (int) Manager.GetInteger("boundaryterm-order"));
	}
      else
	{
	  sprintf (FileNamePrefix, "spin_1_2_naturalboundaryterms_jydominated_%d", (int) Manager.GetInteger("boundaryterm-order"));
	}
    }
  if (Manager.GetString("output-suffix") != 0)
    {
      char* TmpString = new char [strlen(FileNamePrefix) + strlen(Manager.GetString("output-suffix")) + 2];
      sprintf (TmpString, "%s_%s", FileNamePrefix, Manager.GetString("output-suffix"));
      delete[] FileNamePrefix;
      FileNamePrefix = TmpString;
    }
  if (HValue == 0.0)
    {
      sprintf (OutputFileName, "%s_x_%.6f_y_%.6f_z_%.6f_b_%d_n_%d", FileNamePrefix, JxValue, JyValue, JzValue, BValue, NbrSpins);
      sprintf (CommentLine, " XYZ chain with %d sites, Jx=%.6f, Jy= %.6f, Jz=%.6f and boundary conditions B=%d\n# ", NbrSpins, JxValue, JyValue, JzValue, BValue);
    }
  else
    {
      sprintf (OutputFileName, "%s_x_%.6f_y_%.6f_z_%.6f_h_%.6f_b_%d_n_%d", FileNamePrefix, JxValue, JyValue, JzValue, HValue, BValue, NbrSpins);
      sprintf (CommentLine, " XYZ chain with %d sites, Jx=%.6f, Jy= %.6f, Jz=%.6f, H=%.6f and boundary conditions B=%d\n# ", NbrSpins, JxValue, JyValue, JzValue, HValue, BValue);
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  bool FirstRun = true;
  if (Manager.GetBoolean("no-parity") == true)
    { 
      Spin1_2Chain* Chain = new Spin1_2ChainFull (NbrSpins);
      
      if (Chain->GetHilbertSpaceDimension() > 0)
	{
	  SpinChainXYZHamiltonian* Hamiltonian = 0;
	  if (Manager.GetBoolean("natural-boundaryterms") == false)
	    Hamiltonian = new SpinChainXYZHamiltonian (Chain, NbrSpins, JxValue, JyValue, JzValue, HValue, (double) BValue);
	  else
	    Hamiltonian = new SpinChainXYZNaturalBoundaryTermHamiltonian (Chain, NbrSpins, JxValue, JyValue, JzValue, HValue, (double) BValue, 
									  (int) Manager.GetInteger("boundaryterm-order"), false, 0, Manager.GetBoolean("jy-dominated"));
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpEigenstateString, "%s", OutputFileName);
	  char TmpEntry = '\0';
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, &TmpEntry, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete[] TmpEigenstateString;
	  delete Hamiltonian;
	}
      delete Chain;
    }
  else
    {
      int MaxQValue = 1;
      int InitialQValue = 0;
      if (Manager.GetInteger("initial-q") >= 0)
	{
	  InitialQValue = Manager.GetInteger("initial-q") % 2;
	}
      if (Manager.GetInteger("nbr-q") > 0)
	{
	  MaxQValue = InitialQValue + Manager.GetInteger("nbr-q") - 1;
	  if (MaxQValue >= 1)
	    MaxQValue = 1;
	}
      for (; InitialQValue <= MaxQValue; ++InitialQValue)
	{
	  Spin1_2Chain* Chain = new Spin1_2ChainFixedParity (NbrSpins, InitialQValue);
	  if (Chain->GetHilbertSpaceDimension() > 0)
	    {	     
	      SpinChainXYZHamiltonian* Hamiltonian = 0;
	      if (Manager.GetBoolean("natural-boundaryterms") == false)
		Hamiltonian = new SpinChainXYZHamiltonian (Chain, NbrSpins, JxValue, JyValue, JzValue, HValue, (double) BValue);
	      else
		Hamiltonian = new SpinChainXYZNaturalBoundaryTermHamiltonian (Chain, NbrSpins, JxValue, JyValue, JzValue, HValue, (double) BValue, 
									      (int) Manager.GetInteger("boundaryterm-order"), true, InitialQValue, 
									      Manager.GetBoolean("jy-dominated"));
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_q_%d", OutputFileName, InitialQValue);
	      char* TmpQString = new char[64];
	      sprintf (TmpQString, "%d", InitialQValue);
	      GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
				       FirstRun, TmpEigenstateString);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      FirstRun = false;
	      delete[] TmpEigenstateString;
	      delete Hamiltonian;
	    }
	  delete Chain;
	}
    }
  return 0;
}
