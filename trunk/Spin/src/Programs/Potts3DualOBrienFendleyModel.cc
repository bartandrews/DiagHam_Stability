#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/Potts3ChainDualOBrienFendleyHamiltonian.h"
#include "Hamiltonian/Potts3ChainDualOBrienFendleyHamiltonianWithTranslations.h"
#include "Hamiltonian/Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations.h"

#include "HilbertSpace/Potts3Chain.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"
#include "HilbertSpace/Potts3ChainWithTranslationsAndInversion.h"

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
  OptionManager Manager ("Potts3DualOBrienFendleyModel" , "0.01");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

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
  (*SystemGroup) += new  BooleanOption  ('\n', "periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-inversionsymmetry", "if non zero, set the inversion symmetry sector", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
  (*SystemGroup) += new  BooleanOption  ('\n', "disable-momentum", "disable translation symmetry when using periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-realhamiltonian", "do not use a real Hamiltonian at the inversion symmetric points");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Potts3DualOBrienFendleyModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  bool UseMomentumFlag = true;
  
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  if (Manager.GetBoolean("periodic") == false)
    {
      UseMomentumFlag = false;
      sprintf (OutputFileName, "potts3_opendualobrienfendley_n_%d", NbrSpins);
      sprintf (CommentLine, " Potts 3 dual O'Brien-Fendley model with %d sites and open boundary conditions \n#\n# Q", NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "potts3_periodicdualobrienfendley_n_%d", NbrSpins);
      if (Manager.GetBoolean("disable-momentum") == true)
	{
	  UseMomentumFlag = false;
	  sprintf (CommentLine, " Potts 3 dual O'Brien-Fendley model with %d sites and periodic boundary conditions\n#\n# Q", NbrSpins);
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " Potts 3 dual O'Brien-Fendley model with %d sites and periodic boundary conditions\n#\n# Q k InvSym", NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " Potts 3 dual O'Brien-Fendley model with %d sites and periodic boundary conditions\n#\n# Q k", NbrSpins);
	    }
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
  int MinInversionSymmetrySector = -1;
  int MaxInversionSymmetrySector = 1;
  if (Manager.GetInteger("set-inversionsymmetry") != 0)
    {
      MinInversionSymmetrySector = Manager.GetInteger("set-inversionsymmetry");
      MaxInversionSymmetrySector = MinInversionSymmetrySector;
    }
  bool FirstRun = true;

  if (UseMomentumFlag == true)
    {
       for (; InitialQValue <= MaxQValue; ++InitialQValue)
	{
	  int MaxMomentum = NbrSpins;
	  int MomentumSector = 0;
	  if ((Manager.GetInteger("momentum") >= 0) && (Manager.GetInteger("momentum") < MaxMomentum))
	    {
	      MomentumSector =  Manager.GetInteger("momentum");
	      MaxMomentum = MomentumSector + 1;
	    }
	  for (; MomentumSector < MaxMomentum; ++MomentumSector)
	    {
	      if ((Manager.GetBoolean("disable-inversionsymmetry") == false) &&
		  ((MomentumSector == 0) || (((NbrSpins & 1) == 0) && (MomentumSector == (NbrSpins >> 1)))))
		{
		  for (int InversionSymmetrySector = MinInversionSymmetrySector; InversionSymmetrySector <= MaxInversionSymmetrySector; InversionSymmetrySector += 2)
		    {
		      cout << "Q=" << InitialQValue << " k=" << MomentumSector << " inversion sector=" << InversionSymmetrySector << endl;
		      Potts3ChainWithTranslations* Chain = new Potts3ChainWithTranslationsAndInversion (NbrSpins, InitialQValue, MomentumSector, InversionSymmetrySector, 1000000);      
		      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		      char* TmpQString = new char[64];
		      sprintf (TmpQString, "%d %d %d", InitialQValue, MomentumSector, InversionSymmetrySector);
		      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		      sprintf (TmpEigenstateString, "%s_q_%d_invsym_%d_k_%d", OutputFileName, InitialQValue, InversionSymmetrySector, MomentumSector);
		      if (Manager.GetBoolean("disable-realhamiltonian") == false)
			{
			  Lanczos.SetRealAlgorithms();
			  Potts3ChainDualOBrienFendleyRealHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins);
			  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
						   FirstRun, TmpEigenstateString);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			  Lanczos.SetComplexAlgorithms();
			}
		      else
			{
			  Lanczos.SetComplexAlgorithms();
			  Potts3ChainDualOBrienFendleyHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins);
			  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
						      FirstRun, TmpEigenstateString);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			}
		      FirstRun = false;
		      delete Chain;
		      delete[] TmpQString;
		    }
		}
	      else
		{
		  cout << "Q=" << InitialQValue << " k=" << MomentumSector << endl;
		  Potts3ChainWithTranslations* Chain = new Potts3ChainWithTranslations (NbrSpins, InitialQValue, MomentumSector, 1000000);      
		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		  Potts3ChainDualOBrienFendleyHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins);
		  char* TmpQString = new char[64];
		   if (Manager.GetBoolean("disable-inversionsymmetry") == false)
		     {
		       sprintf (TmpQString, "%d %d 0", InitialQValue, MomentumSector);
		     }
		   else
		     {
		       sprintf (TmpQString, "%d %d", InitialQValue, MomentumSector);
		     }
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_q_%d_k_%d", OutputFileName, InitialQValue, MomentumSector);
		  Lanczos.SetComplexAlgorithms();
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
    }
   else
     {
       for (; InitialQValue <= MaxQValue; ++InitialQValue)
	 {
	   cout << "Q=" << InitialQValue << endl;
	   Potts3Chain* Chain = new Potts3Chain (NbrSpins, InitialQValue, 1000000);      
	   Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
	   Potts3ChainDualOBrienFendleyHamiltonian* Hamiltonian = 0;
	   Hamiltonian = new Potts3ChainDualOBrienFendleyHamiltonian (Chain, NbrSpins, Manager.GetBoolean("periodic"), 
								      ((long) Manager.GetInteger("memory")) << 20);
	   char* TmpQString = new char[64];
	   sprintf (TmpQString, "%d", InitialQValue);
	   char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	   sprintf (TmpEigenstateString, "%s_q_%d", OutputFileName, InitialQValue);
	   
	   GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
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
