#include "Hamiltonian/SpinChainZ2InteractingHamiltonian.h"
#include "Hamiltonian/SpinChainZ2InteractingHamiltonianWithTranslations.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFixedParity.h"
#include "HilbertSpace/Spin1_2ChainFixedParityWithTranslations.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

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
  OptionManager Manager ("Z2InteractingChain" , "0.01");
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
  (*SystemGroup) += new  SingleDoubleOption('j', "j-value", "", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('f', "f-value", "", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('v', "v-value", "", 1.0);
  (*SystemGroup) += new  SingleStringOption('\n', "full-fvalues", "provide the f value for each site from a text file");
  (*SystemGroup) += new  BooleanOption('\n', "random-fvalues", "use a randome f value for each site, maximum amplitude is set by f-value");
  (*SystemGroup) += new  SingleIntegerOption ('b', "boundary-conditions", "boundary conditions (0 for open, 1 for periodic, -1 for antiperiodic)", 0);
  (*SystemGroup) += new  BooleanOption  ('\n', "no-parity", "do not take into account the parity when computing the spectrum");
  (*SystemGroup) += new  BooleanOption('\n', "use-momentum", "compute the spectrum using the momentum as a good quantum number");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Z2InteractingChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  double JValue = Manager.GetDouble("j-value");
  double FValue = Manager.GetDouble("f-value");
  double VValue = Manager.GetDouble("v-value");
  double* JValues = new double[NbrSpins];
  double* FValues = new double[NbrSpins];
  double* VValues = new double[NbrSpins];
  bool UseMomentumFlag = Manager.GetBoolean("use-momentum");
  if (Manager.GetBoolean("random-fvalues") == false)
    {
      for (int i = 0; i < NbrSpins; ++i)
	{
	  FValues[i] = FValue;
	}
      if (UseMomentumFlag == true)
	{
	  Lanczos.SetComplexAlgorithms();
	}
    }
  else
    {
      if (UseMomentumFlag == true)
	{
	  cout << "error, can't use momentum as a good quantum number when having disorder" << endl;
	  return -1;
	}
      AbstractRandomNumberGenerator* RandomNumber = 0;
      timeval CurrentTime;
      gettimeofday (&(CurrentTime), 0);
       RandomNumber = new StdlibRandomNumberGenerator (CurrentTime.tv_usec);
      for (int i = 0; i < NbrSpins; ++i)
	{
	  FValues[i] = FValue * 2.0 * (RandomNumber->GetRealRandomNumber() - 0.5);
	}
    }

  int BValue = Manager.GetInteger("boundary-conditions");
  char* OutputFileName = new char [512];
  char* CommentLine = new char [1024 + (NbrSpins * 32)];

  if (Manager.GetBoolean("random-fvalues") == false)
    { 
      sprintf (OutputFileName, "z2interactingchain_n_%d_j_%.6f_f_%.6f_v_%.6f_b_%d", NbrSpins, JValue, FValue, VValue, BValue);
    }
  else
    {
      sprintf (OutputFileName, "z2interactingchain_n_%d_j_%.6f_f_random_%.6f_v_%.6f_b_%d", NbrSpins, JValue, FValue, VValue, BValue);
    }
  if (Manager.GetBoolean("no-parity") == true)
    { 
      if (Manager.GetBoolean("random-fvalues") == false)
	{
	  sprintf (CommentLine, " Z2 interacting chain with %d sites, J=%.6f, F= %.6f, V=%.6f and boundary conditions B=%d\n# ", NbrSpins, JValue, FValue, VValue, BValue);
	}
      else
	{
	  sprintf (CommentLine, " Z2 interacting chain with %d sites, J=%.6f, V=%.6f and boundary conditions B=%d\n# F Values = ", NbrSpins, JValue, VValue, BValue);
	  char* TmpCommentLine;
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      TmpCommentLine = CommentLine + strlen(CommentLine);
	      sprintf (TmpCommentLine, " %.14f", FValues[i]);
	    }
	  TmpCommentLine = CommentLine + strlen(CommentLine);
	  sprintf (TmpCommentLine, "\n#");	  
	}
    }
  else
    {
      if (Manager.GetBoolean("random-fvalues") == false)
	{
	  if (UseMomentumFlag == false)
	    {
	      sprintf (CommentLine, " Z2 interacting chain with %d sites, J=%.6f, F= %.6f, V=%.6f and boundary conditions B=%d\n# Q ", NbrSpins, JValue, FValue, VValue, BValue);
	    }
	  else
	    {
	      sprintf (CommentLine, " Z2 interacting chain with %d sites, J=%.6f, F= %.6f, V=%.6f and boundary conditions B=%d\n# Q k ", NbrSpins, JValue, FValue, VValue, BValue);
	    }
	}
      else
	{
	  sprintf (CommentLine, " Z2 interacting chain with %d sites, J=%.6f, V=%.6f and boundary conditions B=%d\n# F Values = ", NbrSpins, JValue, VValue, BValue);
	  char* TmpCommentLine;
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      TmpCommentLine = CommentLine + strlen(CommentLine);
	      sprintf (TmpCommentLine, " %.14f", FValues[i]);
	    }
	  TmpCommentLine = CommentLine + strlen(CommentLine);
	  sprintf (TmpCommentLine, "\n#");	  
	}
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  bool FirstRun = true;

  if (Manager.GetBoolean("no-parity") == true)
    { 
      Spin1_2Chain* Chain = new Spin1_2ChainFull (NbrSpins);
      if (Chain->GetHilbertSpaceDimension() > 0)
	{
	  SpinChainZ2InteractingHamiltonian Hamiltonian (Chain, NbrSpins, JValue, FValues, VValue, (double) BValue);
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpEigenstateString, "%s", OutputFileName);
	  char TmpEntry = '\0';
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, &TmpEntry, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete[] TmpEigenstateString;
	}
      delete Chain;
    }
  else
    {
      if (UseMomentumFlag == false)
	{
	  int InitalQValue = 0;
	  int MaxQValue = 1;
	  for (; InitalQValue <= MaxQValue; ++InitalQValue)
	    {
	      Spin1_2Chain* Chain = new Spin1_2ChainFixedParity (NbrSpins, InitalQValue);
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{	     
		  for (int i = 0 ; i < Chain->GetHilbertSpaceDimension(); ++i)
		    Chain->PrintState(cout, i) << endl; 
		  SpinChainZ2InteractingHamiltonian	Hamiltonian (Chain, NbrSpins, JValue, FValues, VValue, (double) BValue);
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_q_%d", OutputFileName, InitalQValue);
		  char* TmpQString = new char[64];
		  sprintf (TmpQString, "%d", InitalQValue);
		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete[] TmpEigenstateString;
		}
	      delete Chain;
	    }
	}
      else
	{
	  int InitalQValue = 0;
	  int MaxQValue = 1;
	  for (; InitalQValue <= MaxQValue; ++InitalQValue)
	    {
	      for (int Momentum = 0; Momentum < NbrSpins; ++Momentum)
		{
		  Spin1_2ChainWithTranslations* Chain = new Spin1_2ChainFixedParityWithTranslations (NbrSpins, Momentum, 1, InitalQValue);
		  if (Chain->GetHilbertSpaceDimension() > 0)
		    {	     
		      for (int i = 0 ; i < Chain->GetHilbertSpaceDimension(); ++i)
			Chain->PrintState(cout, i) << endl; 
		      SpinChainZ2InteractingHamiltonianWithTranslations	Hamiltonian (Chain, NbrSpins, JValue, FValue, VValue, (double) BValue);
		      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		      sprintf (TmpEigenstateString, "%s_q_%d_k_%d", OutputFileName, InitalQValue, Momentum);
		      char* TmpQString = new char[64];
		      sprintf (TmpQString, "%d %d", InitalQValue, Momentum);
		      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
						  FirstRun, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      FirstRun = false;
		      delete[] TmpEigenstateString;
		    }
		  delete Chain;
		}
	    }
	}
    }
  return 0;
}
