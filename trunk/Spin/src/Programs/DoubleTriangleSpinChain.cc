#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/DoubleTrianglePeriodicSpinChainHamiltonian.h"
#include "Hamiltonian/DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
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
  OptionManager Manager ("DoubleTriangleSpinChain" , "0.01");
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
  (*SystemGroup) += new  SingleIntegerOption ('\n', "fixed-k", "compute only a given single k sector (all if negative)", -1);
  (*SystemGroup) += new  SingleDoubleOption ('j', "j1-value", "nearest neighbourg coupling constant value", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('g', "j2-value", "second nearest neighbourg coupling constant value", 0.5);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "djz1-value", "delta added to the nearest neighbourg coupling constant value along z", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "djz2-value", "delta added to the second nearest neighbourg coupling constant value along z", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-translations", "do not use translation symmetry");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type DoubleTriangleSpinChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");
  if ((NbrSpins & 1) != 0)
    {
      cout << "the number of spin has to be even" << endl;
      return -1;
    }
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  double J1Value = Manager.GetDouble("j1-value");
  double J2Value = Manager.GetDouble("j2-value");
  double DJz1Value = Manager.GetDouble("djz1-value");
  double DJz2Value = Manager.GetDouble("djz2-value");
  if (Manager.GetBoolean("no-translations") == true)
    {
      Lanczos.SetRealAlgorithms();
      if ((SpinValue & 1) == 0)
	{
	  sprintf (OutputFileName, "spin_%d_doubletrianglechain_j1_%f_j2_%f_djz1_%f_djz2_%f_n_%d", (SpinValue / 2), J1Value, J2Value, DJz1Value, DJz2Value, NbrSpins);
	  sprintf (CommentLine, " double triangle spin %d / 2 chain with %d sites j1=%f j2=%f djz1=%f djz2=%f \n# 2Sz ", (SpinValue / 2), NbrSpins, J1Value, J2Value, DJz1Value, DJz2Value);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_2_doubletrianglechain_j1_%f_j2_%f_djz1_%f_djz2_%f_n_%d", SpinValue, J1Value, J2Value, DJz1Value, DJz2Value, NbrSpins);
	  sprintf (CommentLine, " double triangle spin %d / 2 chain with %d sites j1=%f j2=%f djz1=%f djz2=%f \n# 2Sz ", SpinValue, NbrSpins, J1Value, J2Value, DJz1Value, DJz2Value);
	}
    }
  else
    {
      Lanczos.SetComplexAlgorithms();
      if ((SpinValue & 1) == 0)
	{
	  sprintf (OutputFileName, "spin_%d_translations_doubletrianglechain_j1_%f_j2_%f_djz1_%f_djz2_%f_n_%d", (SpinValue / 2), J1Value, J2Value, DJz1Value, DJz2Value, NbrSpins);
	  sprintf (CommentLine, " double triangle spin %d chain with tanslations and %d sites j1=%f j2=%f djz1=%f djz2=%f \n# 2Sz K ", (SpinValue / 2), NbrSpins, J1Value, J2Value, DJz1Value, DJz2Value);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_2_translations_doubletrianglechain_j1_%f_j2_%f_djz1_%f_djz2_%f_n_%d", SpinValue, J1Value, J2Value, DJz1Value, DJz2Value, NbrSpins);
	  sprintf (CommentLine, " double triangle spin %d / 2 chain with tanslations and %d sites j1=%f j2=%f djz1=%f djz2=%f \n# 2Sz K ", SpinValue, NbrSpins, J1Value, J2Value, DJz1Value, DJz2Value);
	}
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
  if (Manager.GetBoolean("no-translations") == true)
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  AbstractSpinChain* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      Chain = new Spin1_2Chain (NbrSpins, InitalSzValue, 1000000);
	      break;
	    case 2 :
	      Chain = new Spin1Chain (NbrSpins, InitalSzValue, 1000000);
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
	  
	  DoubleTrianglePeriodicSpinChainHamiltonian Hamiltonian (Chain, NbrSpins, J1Value, J2Value, DJz1Value, DJz2Value);
	  char* TmpSzString = new char[64];
	  sprintf (TmpSzString, "%d", InitalSzValue);
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Chain;
	  delete[] TmpSzString;
	}
    }
  else
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  int MaxKValue = NbrSpins;
	  int InitialKValue = 0;
	  if (Manager.GetInteger("fixed-k") >= 0)
	    {
	      InitialKValue = Manager.GetInteger("fixed-k") % NbrSpins;
	      MaxKValue = InitialKValue;
	    }
	  for (; InitialKValue < MaxKValue; ++InitialKValue)
	    {
	      AbstractSpinChainWithTranslations* Chain = 0;
	      switch (SpinValue)
		{
		case 1 :
		  Chain = new Spin1_2ChainWithTranslations (NbrSpins, InitialKValue, 1, InitalSzValue, 1000000, 1000000);
		  break;
		case 2 :
		  Chain = new Spin1ChainWithTranslations (NbrSpins, InitialKValue, InitalSzValue);
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
		  cout << "2Sz=" <<  InitalSzValue << " K=" << InitialKValue << " dim=" << Chain->GetHilbertSpaceDimension() << endl;
// 		  for (int i = 0; i < Chain->GetHilbertSpaceDimension(); ++i)
// 		    Chain->PrintState(cout, i) << endl;
		  DoubleTrianglePeriodicSpinChainWithTranslationHamiltonian Hamiltonian (Chain, NbrSpins, J1Value, J2Value, DJz1Value, DJz2Value);
		  char* TmpSzString = new char[64];
		  sprintf (TmpSzString, "%d %d", InitalSzValue, InitialKValue);
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_sz_%d_k_%d", OutputFileName, InitalSzValue, InitialKValue);
		  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					      FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  delete[] TmpSzString;
		}
	      FirstRun = false;
	      delete Chain;
	    }
	}
    }

  return 0;
}
