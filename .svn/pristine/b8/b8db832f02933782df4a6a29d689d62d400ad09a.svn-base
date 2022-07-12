#include "Hamiltonian/HaldaneShastryHamiltonianWithTranslations.h"

#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"

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
  OptionManager Manager ("HaldaneShastrySpinChain" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HaldaneShastrySpinChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "spin_%d_haldaneshastrychain_n_%d", (SpinValue / 2), NbrSpins);
      sprintf (CommentLine, " Haldane Shastry spin %d chain with %d sites \n# 2Sz K", (SpinValue / 2), NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_haldaneshastrychain_n_%d", SpinValue, NbrSpins);
      sprintf (CommentLine, " Haldane Shastry spin %d/2 chain with %d sites \n# 2Sz K ", SpinValue, NbrSpins);
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if (Manager.GetInteger("initial-sz") > 1)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      for (int Momentum = 0; Momentum < NbrSpins; ++Momentum)
	{
	  AbstractSpinChainWithTranslations* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      Chain = new Spin1_2ChainWithTranslations (NbrSpins, Momentum, 1, InitalSzValue, 1000000, 1000000);
	      break;
	    case 2 :
	      Chain = new Spin1ChainWithTranslations (NbrSpins, Momentum, InitalSzValue);
	      break;
	    default :
	      {
		if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		else 
		  cout << "spin " << SpinValue << "/2 are not available" << endl;
		return -1;
	      }
	    }

	  if (Chain->GetHilbertSpaceDimension() > 0)
	    {
	      cout << "Sz = " << InitalSzValue << ", K = " << Momentum << endl; 
	      HaldaneShastryHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins);
	      char* TmpSzString = new char[64];
	      sprintf (TmpSzString, "%d %d", InitalSzValue, Momentum);
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_sz_%d_k_%d", OutputFileName, InitalSzValue, Momentum);
	      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					  FirstRun, TmpEigenstateString);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      FirstRun = false;
	      delete[] TmpSzString;
	    }
	  delete Chain;
	}
    }
  return 0;
}
