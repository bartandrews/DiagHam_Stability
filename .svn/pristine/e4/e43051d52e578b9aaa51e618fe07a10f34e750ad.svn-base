#include "Hamiltonian/SpinChainOBrienFendleyHamiltonian.h"
#include "Hamiltonian/SpinChainOBrienFendleyHamiltonianWithTranslations.h"
#include "Hamiltonian/SpinChainOBrienFendleyRealHamiltonianWithTranslations.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin3_2Chain.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "Matrix/RealMatrix.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"

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


// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
void SpinChainOBrienFendleyComputeCharacteristicPolynomial(SpinChainAKLTRealHamiltonianWithTranslations* hamiltonian, AbstractSpinChainWithTranslations* chain, char* outputFileName, AbstractArchitecture* architecture);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("PeriodicSpinChainOBrienFendley" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 2);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "additional-quadratic", "coefficient in front of the additional quadratic term (0 being the pure OBrienFendley hamiltonian)", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "linear-factor", "if different from 1.0, set the coefficient in front of the Heisenberg term", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "quadratic-factor", "if --linear-factor is different from 1.0, set the coefficient in front of the quadratic term", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-szsymmetry", "if non zero, set the Sz<->-Sz symmetry sector", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-inversionsymmetry", "if non zero, set the inversion symmetry sector", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "disable-momentum", "disable the translation symmetry");
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
  (*ToolsGroup) += new BooleanOption  ('\n', "export-charpolynomial", "export the hamiltonian characteristic polynomial");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PeriodicSpinChainOBrienFendley -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");
  bool UseMomentum = !(Manager.GetBoolean("disable-momentum"));
  
  char* OutputFileName = new char [512];
  char* MomentumFlagString = new char [256];
  if (UseMomentum == true)
    {
      sprintf(MomentumFlagString, "periodicobrienfendley");
    }
  else
    {
      sprintf(MomentumFlagString, "periodicobrienfendley_nomomentum");
    }
  char* CommentLine = new char [512];
  if ((SpinValue & 1) == 0)
    {
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_%s_quadratic_%.6f_n_%d", (SpinValue / 2), MomentumFlagString, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  if (Manager.GetDouble("linear-factor") != 1.0)
	    {
	      sprintf (OutputFileName, "spin_%d_%s_linear_%.6f_quadratic_%.6f_n_%d", (SpinValue / 2), MomentumFlagString, 
		       Manager.GetDouble("linear-factor"),  Manager.GetDouble("quadratic-factor"), NbrSpins);
	    }
	  else
	    {
	      sprintf (OutputFileName, "spin_%d_%s_n_%d", (SpinValue / 2), MomentumFlagString, NbrSpins);
	    }
	}
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K SzSym InvSym ", (SpinValue / 2), NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz SzSym InvSym ", (SpinValue / 2), NbrSpins);
		}
	    }
	  else
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K SzSym ", (SpinValue / 2), NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz SzSym ", (SpinValue / 2), NbrSpins);
		}		  
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K InvSym ", (SpinValue / 2), NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz InvSym ", (SpinValue / 2), NbrSpins);
		}
	    }
	  else
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz K ", (SpinValue / 2), NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d chain with %d sites \n# 2Sz ", (SpinValue / 2), NbrSpins);
		}
	    }
	}
    }
  else
    {
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_2_%s_quadratic_%.6f_n_%d", SpinValue, MomentumFlagString, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_2_%s_n_%d", SpinValue, MomentumFlagString, NbrSpins);
	}
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K SzSym InvSym ", SpinValue, NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz SzSym InvSym ", SpinValue, NbrSpins);
		}
	    }
	  else
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K SzSym ", SpinValue, NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz SzSym ", SpinValue, NbrSpins);
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K InvSym ", SpinValue, NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz InvSym ", SpinValue, NbrSpins);
		}
	    }
	  else
	    {
	      if (UseMomentum == true)
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz K ", SpinValue, NbrSpins);
		}
	      else
		{
		  sprintf (CommentLine, " periodic spin %d/2 chain with %d sites \n# 2Sz ", SpinValue, NbrSpins);
		}
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

  int MinSzSymmetrySector = -1;
  int MaxSzSymmetrySector = 1;
  if (Manager.GetInteger("set-szsymmetry") != 0)
    {
      MinSzSymmetrySector = Manager.GetInteger("set-szsymmetry");
      MaxSzSymmetrySector = MinSzSymmetrySector;
    }
  int MinInversionSymmetrySector = -1;
  int MaxInversionSymmetrySector = 1;
  if (Manager.GetInteger("set-inversionsymmetry") != 0)
    {
      MinInversionSymmetrySector = Manager.GetInteger("set-inversionsymmetry");
      MaxInversionSymmetrySector = MinInversionSymmetrySector;
    }

  if (UseMomentum == true)
    {
      if ((InitalSzValue == 0) && (Manager.GetBoolean("disable-szsymmetry") == false) && (SpinValue == 2))
	{
	  for (int Momentum = InitialMomentum; Momentum < MaxMomentum; ++Momentum)
	    {
	      for (int SzSymmetrySector = MinSzSymmetrySector; SzSymmetrySector <= MaxSzSymmetrySector; SzSymmetrySector += 2)
		{
		  if ((Manager.GetBoolean("disable-inversionsymmetry") == false)  && (SpinValue == 2) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins >> 1)))))
		    {
		      for (int InversionSymmetrySector = MinInversionSymmetrySector; InversionSymmetrySector <= MaxInversionSymmetrySector; InversionSymmetrySector += 2)
			{
			  AbstractSpinChainWithTranslations* Chain = 0;
			  switch (SpinValue)
			    {
			    case 2 :
			      Chain = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momentum, InversionSymmetrySector, SzSymmetrySector, InitalSzValue);
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
				  SpinChainOBrienFendleyRealHamiltonianWithTranslations* Hamiltonian = 0;
				  if (Manager.GetDouble("linear-factor") == 1.0)
				    {
				      Hamiltonian = new SpinChainOBrienFendleyRealHamiltonianWithTranslations(Chain, NbrSpins);
				    }
				  else
				    {
				      Hamiltonian = new SpinChainOBrienFendleyRealHamiltonianWithTranslations(Chain, NbrSpins);
				    }
				  if (Manager.GetBoolean("export-charpolynomial"))
				    {
				      SpinChainOBrienFendleyComputeCharacteristicPolynomial(Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture());
				    }			      				  
				  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
							   FirstRun, TmpEigenstateString);
				  MainTaskOperation TaskOperation (&Task);
				  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
				  Lanczos.SetComplexAlgorithms();
				  delete Hamiltonian;
				}
			      else
				{
				  SpinChainOBrienFendleyHamiltonianWithTranslations* Hamiltonian;
				  if(Manager.GetDouble("linear-factor") == 1.0)
				    { 
				      Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations  (Chain, NbrSpins);
				    }
				  else
				    {
				      Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations (Chain, NbrSpins);
				    }
				  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
							      FirstRun, TmpEigenstateString);
				  MainTaskOperation TaskOperation (&Task);
				  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
				  delete Hamiltonian;
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
			case 2 :
			  Chain = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSpins, Momentum, SzSymmetrySector, InitalSzValue);
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
			  SpinChainOBrienFendleyHamiltonianWithTranslations* Hamiltonian = 0;
			  if(Manager.GetDouble("linear-factor") == 1.0)
			    {			  
			      Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations(Chain, NbrSpins);
			    }
			  else
			    {
			      Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations(Chain, NbrSpins);
			    }
			  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
						      FirstRun, TmpEigenstateString);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			  FirstRun = false;
			  delete[] TmpSzString;
			  delete Hamiltonian;
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
	      if ((Manager.GetBoolean("disable-inversionsymmetry") == false)  && (SpinValue == 2) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins >> 1)))))
		{
		  for (int InversionSymmetrySector = MinInversionSymmetrySector; InversionSymmetrySector <= MaxInversionSymmetrySector; InversionSymmetrySector += 2)
		    {
		      AbstractSpinChainWithTranslations* Chain = 0;
		      switch (SpinValue)
			{
			case 2 :
			  Chain = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momentum, InversionSymmetrySector, InitalSzValue);
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
			      SpinChainOBrienFendleyRealHamiltonianWithTranslations* Hamiltonian = 0;
			      if (Manager.GetDouble("linear-factor") == 1.0)
				{
				  Hamiltonian = new SpinChainOBrienFendleyRealHamiltonianWithTranslations(Chain, NbrSpins);
				}
			      else
				{
				  Hamiltonian = new SpinChainOBrienFendleyRealHamiltonianWithTranslations(Chain, NbrSpins);
				}
			      if (Manager.GetBoolean("export-charpolynomial"))
				{
				  SpinChainOBrienFendleyComputeCharacteristicPolynomial(Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture());
				}			      				  
			      GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
						       FirstRun, TmpEigenstateString);
			      MainTaskOperation TaskOperation (&Task);
			      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			      Lanczos.SetComplexAlgorithms();
			      delete Hamiltonian;
			    }
			  else
			    {
			      SpinChainOBrienFendleyHamiltonianWithTranslations* Hamiltonian ;
			      if (Manager.GetDouble("linear-factor") == 1.0)
				{
				  Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations(Chain, NbrSpins);
				}
			      else
				{
				  Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations(Chain, NbrSpins);
				}
			      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
							  FirstRun, TmpEigenstateString);
			      MainTaskOperation TaskOperation (&Task);
			      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			      delete Hamiltonian;
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
		      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		      cout << "2Sz = " << InitalSzValue << ", K = " << Momentum << endl; 
		      SpinChainOBrienFendleyHamiltonianWithTranslations* Hamiltonian = 0;
		      if (Manager.GetDouble("linear-factor") == 1.0)
			{
			  Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations(Chain, NbrSpins);
			}
		      else
			{
			  Hamiltonian = new SpinChainOBrienFendleyHamiltonianWithTranslations(Chain, NbrSpins);
			}
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
		      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
						  FirstRun, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      FirstRun = false;
		      delete[] TmpSzString;
		      delete Hamiltonian;
		    }
		  delete Chain;
		}
	    }
	}
    }
  else
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  AbstractSpinChain* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      {
		Chain = new Spin1_2Chain (NbrSpins, InitalSzValue, 1000000);
	      }
	      break;
	    case 2 :
	      {
		Chain = new Spin1Chain (NbrSpins, InitalSzValue, 1000000);
	      }
	      break;
	    case 3 :
	      {
		Chain = new Spin3_2Chain (NbrSpins, InitalSzValue, 1000000);
	      }
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
	  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());		  
	  SpinChainOBrienFendleyHamiltonian* Hamiltonian = 0;

	  Hamiltonian = new SpinChainOBrienFendleyHamiltonian(Chain, NbrSpins, true);
	  Lanczos.SetRealAlgorithms();
	  
	  char* TmpSzString = new char[64];
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpSzString, "%d", InitalSzValue);
	  sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);


	  if (Manager.GetBoolean("export-charpolynomial"))
	    {
#ifdef __GMP__
	      LongIntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	      Hamiltonian->GetHamiltonian(TmpMatrix);
	      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();
#else
	      IntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	      Hamiltonian->GetHamiltonian(TmpMatrix);
	      long* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();
#endif	       
	      char* PolynomialOutputFileName = new char[strlen(OutputFileName) + 256];
	      sprintf (PolynomialOutputFileName, "%s_sz_%d.charpol", OutputFileName, InitalSzValue);
	      ofstream OutputFile;
	      OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
	      OutputFile << CharacteristicPolynomial[0];
	      for (int i = 1; i <= Chain->GetHilbertSpaceDimension(); ++i)
		{
		  OutputFile << "," << CharacteristicPolynomial[i];
		}
	      OutputFile << endl;
	      OutputFile.close();
	    }
	  	    
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName, FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Hamiltonian;
	  delete Chain;
	  delete[] TmpSzString;
	  delete[] TmpEigenstateString;
	}
    }
  return 0;
}

// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture

void SpinChainOBrienFendleyComputeCharacteristicPolynomial(SpinChainAKLTRealHamiltonianWithTranslations* hamiltonian, AbstractSpinChainWithTranslations* chain, char* outputFileName, AbstractArchitecture* architecture)
{
#ifdef __GMP__
  RealMatrix TmpRawMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
  hamiltonian->GetHamiltonian(TmpRawMatrix);
  double* TmpNormalizationFactors = chain->GetBasisNormalization();
  for (int i = 0; i < chain->GetHilbertSpaceDimension(); ++i)
    {
      for (int j = 0; j < chain->GetHilbertSpaceDimension(); ++j)
	{
	  double Tmp;
	  TmpRawMatrix.GetMatrixElement(i, j, Tmp);
	  Tmp *= TmpNormalizationFactors[i];
	  Tmp /= TmpNormalizationFactors[j];
	  TmpRawMatrix.SetMatrixElement(i, j, Tmp);
	}
    }
  LongIntegerMatrix TmpMatrix(TmpRawMatrix);
  cout << "Start computing characteristic polynomial (degree " << chain->GetHilbertSpaceDimension() << ")" << endl;
  mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(architecture);
  char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
  sprintf (PolynomialOutputFileName, "%s.charpol", outputFileName);
  ofstream OutputFile;
  OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
  OutputFile << CharacteristicPolynomial[0];
  for (int i = 1; i <= chain->GetHilbertSpaceDimension(); ++i)
    {
      OutputFile << "," << CharacteristicPolynomial[i];
    }
  OutputFile << endl;
  OutputFile.close();
#else
  cout << "GMP library is required for characteristic polynomials" << endl;
#endif	       
}
