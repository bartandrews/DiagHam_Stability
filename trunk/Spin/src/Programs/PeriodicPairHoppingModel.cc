#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"

#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"
#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"

#include "Hamiltonian/PairHoppingHamiltonianWithTranslations.h"
#include "Hamiltonian/PairHoppingRealHamiltonianWithTranslations.h"

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
  OptionManager Manager ("PeriodicPairHoppingModel" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
 OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new  SingleIntegerOption ('p', "p-value", "value that defines the filling factor p/(2p+1)", 1);
  (*SystemGroup) += new  SingleIntegerOption ('n', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-inversionsymmetry", "if non zero, set the inversion symmetry sector", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-realhamiltonian", "do not use a real Hamiltonian at the inversion symmetric points");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "v1-0", "strength of the first electrostatic coupling", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "v2-0", "strength of the second electrostatic coupling", 0.0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PeriodicPairHoppingModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int PValue = Manager.GetInteger("p-value");
  int NbrSpins = Manager.GetInteger("nbr-spin");
  if ((NbrSpins % PValue) != 0)
    {
      cout << "--nbr-spin (here " << NbrSpins << ") should be a multiple of --p-value (here " << PValue << ")" << endl;
      return -1;
    }
  double V10 = Manager.GetDouble("v1-0");
  double V20 = Manager.GetDouble("v2-0");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  sprintf (OutputFileName, "spin_1_periodicpairhopping_p_%d_n_%d", PValue, NbrSpins);
  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
    {
      sprintf (CommentLine, " periodic pair hopping model with p=%d in spin 1 language with %d sites \n# K InvSym ", PValue, NbrSpins);
    }
  else
    {
      sprintf (CommentLine, " periodic pair hopping model with p=%d in spin 1 language with %d sites \n# K ", PValue, NbrSpins);
    }



  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  bool FirstRun = true;
  int InitialMomentum = 0;
  int MaxMomentum = NbrSpins / PValue;
  if ((Manager.GetInteger("momentum") >= 0) && (Manager.GetInteger("momentum") < MaxMomentum))
    {
      InitialMomentum =  Manager.GetInteger("momentum");
      MaxMomentum = InitialMomentum + 1;
    }

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  int MinInversionSymmetrySector = -1;
  int MaxInversionSymmetrySector = 1;
  if (Manager.GetInteger("set-inversionsymmetry") != 0)
    {
      MinInversionSymmetrySector = Manager.GetInteger("set-inversionsymmetry");
      MaxInversionSymmetrySector = MinInversionSymmetrySector;
    }
  if (Manager.GetBoolean("disable-inversionsymmetry") == true)
    {
      MaxInversionSymmetrySector = 0;
      MinInversionSymmetrySector = 0;
    }
    
  for (int Momentum = InitialMomentum; Momentum < MaxMomentum; ++Momentum)
    {
      if  ((Momentum == 0) || ((((NbrSpins / PValue) & 1) == 0) && (Momentum == ((NbrSpins / PValue) >> 1))))
	{
	  for (int InversionSymmetrySector = MinInversionSymmetrySector; InversionSymmetrySector <= MaxInversionSymmetrySector; InversionSymmetrySector += 2)
	    {
	      AbstractSpinChainWithTranslations* Chain = 0;
	      if (Manager.GetBoolean("disable-inversionsymmetry") == false)
		{
		  cout << "K=" << Momentum <<" IP=" << InversionSymmetrySector << endl;
		  switch (PValue)
		    {
		    case 1 :
		      {
			if (NbrSpins <= 32)
			  {
			    Chain = new PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, Momentum, InversionSymmetrySector);
			  }
			else
			  {
			    Chain = new PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (NbrSpins, Momentum, InversionSymmetrySector);
			  }
		      }
		      break;
		    case 2 :
		      {
			if (NbrSpins <= 32)
			  {
			    Chain = new PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, Momentum, InversionSymmetrySector);
			  }
			else
			  {
			    Chain = new PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (NbrSpins, Momentum, InversionSymmetrySector);
			  }
		      }
		      break;
		    default :
		      {
			if (NbrSpins <= 32)
			  {
			    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry (NbrSpins, Momentum, InversionSymmetrySector, PValue);
			  }
			else
			  {
			    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (NbrSpins, Momentum, InversionSymmetrySector, PValue);
			  }
		      }
		    }
		}
	      else
		{
		  cout << "K=" << Momentum << endl;
		  switch (PValue)
		    {
		    case 1 :
		      {
			if (NbrSpins <= 32)
			  {
			    Chain = new PairHoppingP1AsSpin1ChainWithTranslations (NbrSpins, Momentum, 1000000);
			  }
			else
			  {
			    Chain = new PairHoppingP1AsSpin1ChainWithTranslationsLong (NbrSpins, Momentum, 1000000);
			  }
		      }
		      break;
		    case 2 :
		      {
			if (NbrSpins <= 32)
			  {
			    Chain = new PairHoppingP2AsSpin1ChainWithTranslations (NbrSpins, Momentum, 1000000);
			  }
			else
			  {
			    Chain = new PairHoppingP2AsSpin1ChainWithTranslationsLong (NbrSpins, Momentum, 1000000);
			  }
		      }
		      break;
		    default :
		      {
			if (NbrSpins <= 32)
			  {
			    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslations (NbrSpins, Momentum, PValue, 1000000);
			  }
			else
			  {
			    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslationsLong (NbrSpins, Momentum, PValue, 1000000);
			  }
		      }
		    }
		}
	      
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{
		  cout << "Hilbert space dimension = " << Chain->GetLargeHilbertSpaceDimension() << endl;	  
		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());
		  char* TmpString = new char[64];
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
		    {
		      sprintf (TmpString, "%d %d", Momentum, InversionSymmetrySector);
		      sprintf (TmpEigenstateString, "%s_invsym_%d_k_%d", OutputFileName, InversionSymmetrySector, Momentum);
		    }
		  else
		    {
		      sprintf (TmpString, "%d", Momentum);
		      sprintf (TmpEigenstateString, "%s_k_%d", OutputFileName, Momentum);
		    }
		  if (Manager.GetBoolean("disable-realhamiltonian") == false)
		    {
		      Lanczos.SetRealAlgorithms();
		      PairHoppingRealHamiltonianWithTranslations* Hamiltonian = 0;
		      Hamiltonian = new PairHoppingRealHamiltonianWithTranslations(Chain, NbrSpins, PValue);
		      GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
					       FirstRun, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      Lanczos.SetComplexAlgorithms();
		      FirstRun = false;
		      delete Hamiltonian;
		    }
		  else
		    {
		      PairHoppingHamiltonianWithTranslations* Hamiltonian = 0;
		      Hamiltonian = new PairHoppingHamiltonianWithTranslations(Chain, NbrSpins, PValue, V10, V20, Architecture.GetArchitecture(), Memory);
		      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
						  FirstRun, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      FirstRun = false;
		      delete Hamiltonian;
		    }
		}
	    }
	}
      else
	{
	  cout << "K=" << Momentum << endl;
	  AbstractSpinChainWithTranslations* Chain = 0;
	  switch (PValue)
	    {
	    case 1 :
	      {
		if (NbrSpins <= 32)
		  {
		    Chain = new PairHoppingP1AsSpin1ChainWithTranslations (NbrSpins, Momentum, 1000000);
		  }
		else
		  {
		    Chain = new PairHoppingP1AsSpin1ChainWithTranslationsLong (NbrSpins, Momentum, 1000000);
		  }
	      }
	      break;
	    case 2 :
	      {
		if (NbrSpins <= 32)
		  {
		    Chain = new PairHoppingP2AsSpin1ChainWithTranslations (NbrSpins, Momentum, 1000000);
		  }
		else
		  {
		    Chain = new PairHoppingP2AsSpin1ChainWithTranslationsLong (NbrSpins, Momentum, 1000000);
		  }
	      }
	      break;
	    default :
	      {
		if (NbrSpins <= 32)
		  {
		    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslations (NbrSpins, Momentum, PValue, 1000000);
		  }
		else
		  {
		    Chain = new PairHoppingGenericPAsSpin1ChainWithTranslationsLong (NbrSpins, Momentum, PValue, 1000000);
		  }
	      }
	    }
	  
	  if (Chain->GetHilbertSpaceDimension() > 0)
	    {
	      cout << "Hilbert space dimension = " << Chain->GetLargeHilbertSpaceDimension() << endl;	  
	      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
	      PairHoppingHamiltonianWithTranslations* Hamiltonian = 0;
	      Hamiltonian = new PairHoppingHamiltonianWithTranslations(Chain, NbrSpins, PValue, V10, V20, Architecture.GetArchitecture(), Memory);
	      char* TmpString = new char[64];
	      if (Manager.GetBoolean("disable-inversionsymmetry") == false)
		{
		  sprintf (TmpString, "%d 0", Momentum);
		}
	      else
		{
		  sprintf (TmpString, "%d", Momentum);
		}
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_k_%d", OutputFileName, Momentum);
	      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
					  FirstRun, TmpEigenstateString);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      FirstRun = false;
	      delete[] TmpString;
	      delete Hamiltonian;
	    }
	  delete Chain;
	}
    }
  return 0;
}
