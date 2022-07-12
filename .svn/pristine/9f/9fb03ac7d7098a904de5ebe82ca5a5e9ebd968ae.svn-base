#include "Hamiltonian/SpinChainAKLTP3P4HamiltonianWithTranslations.h"
#include "Hamiltonian/SpinChainAKLTP3P4RealHamiltonianWithTranslations.h"

#include "HilbertSpace/Spin2ChainWithTranslations.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzInversionSymmetries.h"

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
  OptionManager Manager ("PeriodicSpinChainGeneralizedAKLT" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 4);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "p3-factor", "scaling factor for the spin 3 projector", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "p4-factor", "scaling factor for the spin 4 projector", 1.0);
  (*SystemGroup) += new  BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-realhamiltonian", "do not use a real Hamiltonian at the inversion symmetric points");
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
      cout << "see man page for option syntax or type PeriodicSpinChainGeneralizedAKLT -h" << endl;
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
  double Spin3ProjectorFactor = Manager.GetDouble("p3-factor");
  double Spin4ProjectorFactor = Manager.GetDouble("p4-factor");
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "spin_%d_periodicaklt_p3_%.6f_p4_%.6f_n_%d", (SpinValue / 2), Spin3ProjectorFactor, Spin4ProjectorFactor, NbrSpins);
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K SzSym InvSym ", (SpinValue / 2), NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K SzSym ", (SpinValue / 2), NbrSpins);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K InvSym ", (SpinValue / 2), NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K ", (SpinValue / 2), NbrSpins);
	    }
	}
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_periodicaklt_p3_%.6f_p4_%.6f_n_%d", SpinValue, Spin3ProjectorFactor, Spin4ProjectorFactor, NbrSpins);
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K SzSym InvSym ", SpinValue, NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K SzSym ", SpinValue, NbrSpins);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K InvSym ", SpinValue, NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K ", SpinValue, NbrSpins);
	    }
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
  int InitialMomentum = 0;
  int MaxMomentum = NbrSpins;
  if ((Manager.GetInteger("momentum") >= 0) && (Manager.GetInteger("momentum") < MaxMomentum))
    {
      InitialMomentum =  Manager.GetInteger("momentum");
      MaxMomentum = InitialMomentum + 1;
    }

  if ((InitalSzValue == 0) && (Manager.GetBoolean("disable-szsymmetry") == false) && (SpinValue == 4))
    {
      for (int Momentum = InitialMomentum; Momentum < MaxMomentum; ++Momentum)
	{
	  for (int SzSymmetrySector = -1; SzSymmetrySector <= 1; SzSymmetrySector += 2)
	    {
	      if ((Manager.GetBoolean("disable-inversionsymmetry") == false)  && (SpinValue == 4) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins >> 1)))))
		{
		  for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
		    {
		      AbstractSpinChainWithTranslations* Chain = 0;
		      switch (SpinValue)
			{
 			case 4 :
 			  Chain = new Spin2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momentum, InversionSymmetrySector, SzSymmetrySector, InitalSzValue);
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
			  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
			  cout << "2Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << ",   inversion sector=" << InversionSymmetrySector << ",  K = " << Momentum << endl; 
			  char* TmpSzString = new char[64];
			  sprintf (TmpSzString, "%d %d %d %d ", InitalSzValue, Momentum, SzSymmetrySector, InversionSymmetrySector);
			  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
			  sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d_szsym_%d_k_%d", OutputFileName, InitalSzValue, InversionSymmetrySector, SzSymmetrySector, Momentum);
			  if (Manager.GetBoolean("disable-realhamiltonian") == false)
			    {
			      Lanczos.SetRealAlgorithms();
			      SpinChainAKLTP3P4RealHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, Spin3ProjectorFactor, Spin4ProjectorFactor);
			      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
							  FirstRun, TmpEigenstateString);
			      MainTaskOperation TaskOperation (&Task);
			      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			      Lanczos.SetComplexAlgorithms();
			    }
			  else
			    {
			      SpinChainAKLTP3P4HamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, Spin3ProjectorFactor, Spin4ProjectorFactor);
			      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
							  FirstRun, TmpEigenstateString);
			      MainTaskOperation TaskOperation (&Task);
			      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			    }
			  FirstRun = false;
			  delete[] TmpSzString;
			}
		      delete Chain;
		    }
		}
	      else
		{
		  AbstractSpinChainWithTranslations* Chain = 0;
		  switch (SpinValue)
		    {
 		    case 4 :
 		      Chain = new Spin2ChainWithTranslationsAndSzSymmetry (NbrSpins, Momentum, SzSymmetrySector, InitalSzValue);
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
		      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		      cout << "2Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << ",  K = " << Momentum << endl; 
		      char* TmpSzString = new char[64];
		      if (Manager.GetBoolean("disable-inversionsymmetry") == false)
			{
			  sprintf (TmpSzString, "%d %d %d 0 ", InitalSzValue, Momentum, SzSymmetrySector);
			}
		      else
			{
			  sprintf (TmpSzString, "%d %d %d ", InitalSzValue, Momentum, SzSymmetrySector);
			}
		      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		      sprintf (TmpEigenstateString, "%s_sz_%d_szsym_%d_k_%d", OutputFileName, InitalSzValue, SzSymmetrySector, Momentum);
		      SpinChainAKLTP3P4HamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, Spin3ProjectorFactor, Spin4ProjectorFactor);
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
	}      
      InitalSzValue +=2;
    }
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      for (int Momentum = InitialMomentum; Momentum < MaxMomentum; ++Momentum)
	{
	  if ((Manager.GetBoolean("disable-inversionsymmetry") == false)  && (SpinValue == 4) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins >> 1)))))
	    {
	      for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
		{
		  AbstractSpinChainWithTranslations* Chain = 0;
		  switch (SpinValue)
		    {
 		    case 4 :
 		      Chain = new Spin2ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momentum, InversionSymmetrySector, InitalSzValue);
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
		      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		      cout << "2Sz = " << InitalSzValue << ", inversion sector=" << InversionSymmetrySector << ",  K = " << Momentum << endl; 
		      char* TmpSzString = new char[64];
		      if (Manager.GetBoolean("disable-szsymmetry") == false)
			{
			  sprintf (TmpSzString, "%d %d 0 %d", InitalSzValue, Momentum, InversionSymmetrySector);
			}
		      else
			{
			  sprintf (TmpSzString, "%d %d %d", InitalSzValue, Momentum, InversionSymmetrySector);
			}
		      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		      sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d_k_%d", OutputFileName, InitalSzValue, InversionSymmetrySector, Momentum);
		      if (Manager.GetBoolean("disable-realhamiltonian") == false)
			{
			  Lanczos.SetRealAlgorithms();
			  SpinChainAKLTP3P4RealHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, Spin3ProjectorFactor, Spin4ProjectorFactor);
			  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
						   FirstRun, TmpEigenstateString);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			  Lanczos.SetComplexAlgorithms();
			}
		      else
			{
			  SpinChainAKLTP3P4HamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, Spin3ProjectorFactor, Spin4ProjectorFactor);
			  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
						      FirstRun, TmpEigenstateString);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			}
		      FirstRun = false;
		      delete[] TmpSzString;
		    }
		  delete Chain;
		}
	    }
	  else
	    {
	      AbstractSpinChainWithTranslations* Chain = 0;
	      switch (SpinValue)
		{
		case 4 :
		  Chain = new Spin2ChainWithTranslations (NbrSpins, Momentum, InitalSzValue);
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
		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		  cout << "2Sz = " << InitalSzValue << ", K = " << Momentum << endl; 
		  SpinChainAKLTP3P4HamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, Spin3ProjectorFactor, Spin4ProjectorFactor);
		  char* TmpSzString = new char[64];
		  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
		    {
		      if (Manager.GetBoolean("disable-szsymmetry") == false)
			{
			  sprintf (TmpSzString, "%d %d 0 0", InitalSzValue, Momentum);
			}
		      else
			{
			  sprintf (TmpSzString, "%d %d 0", InitalSzValue, Momentum);
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean("disable-szsymmetry") == false)
			{
			  sprintf (TmpSzString, "%d %d 0", InitalSzValue, Momentum);
			}
		      else
			{
			  sprintf (TmpSzString, "%d %d", InitalSzValue, Momentum);
			}
		    }
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
    }
  return 0;
}
