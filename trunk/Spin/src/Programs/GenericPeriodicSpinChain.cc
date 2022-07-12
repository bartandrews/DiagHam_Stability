#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"
#include "Hamiltonian/SpinChainRealHamiltonianWithTranslations.h"

#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetryLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetryLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetriesLong.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin2ChainWithTranslations.h"

#include "Operator/SpinWith1DTranslationS2Operator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "Matrix/RealMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
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
// sResolvedFlag = true if the characteristic polynomial should be computed per S sector
// minSValue = twice the minimum S value to consider (when sResolvedFlag is true)
// minSValue = twice the maximum S value to consider (when sResolvedFlag is true)
// sU2Degeneracy = array providing the Hilbert space dimension per S value (when sResolvedFlag is true)
// discardFourFactor = discard a global four factor used to ensure integer numbers
void SpinChainComputeCharacteristicPolynomial(SpinChainRealHamiltonianWithTranslations* hamiltonian, AbstractSpinChainWithTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, bool sResolvedFlag, int minSValue, int maxSValue, int* sU2Degeneracy, bool discardFourFactor = false);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("GenericPeriodicSpinChain" , "0.01");
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
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-szsymmetry", "if non zero, set the Sz<->-Sz symmetry sector", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-inversionsymmetry", "if non zero, set the inversion symmetry sector", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "disable-realhamiltonian", "do not use a real Hamiltonian at the inversion symmetric points");
  (*SystemGroup) += new  SingleDoubleOption ('j', "j-value", "coupling constant value", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "nn-coupling", "add the term to ZZ nearest-neighbour interaction [will be j+nn-coupling]", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "nnn-coupling", "next-nearest-neighbour interaction", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "shift-energies", "shift energies by a constant value", 0.0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "export-charpolynomial", "export the hamiltonian characteristic polynomial");  
  (*ToolsGroup) += new BooleanOption  ('\n', "export-sresolvedcharpolynomial", "export the hamiltonian characteristic polynomial resolved in S quantum number");  
  (*ToolsGroup) += new BooleanOption  ('\n', "charpolynomial-nofour", "when exporting the hamiltonian characteristic polynomial, do not include a global 4 factor");  
   (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericPeriodicSpinChain -h" << endl;
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
      sprintf (OutputFileName, "spin_%d_periodicchain_n_%d", (SpinValue / 2), NbrSpins);
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
      sprintf (OutputFileName, "spin_%d_2_periodicchain_n_%d", SpinValue, NbrSpins);
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
  double JValue =  Manager.GetDouble("j-value");
  double NNCoupling =  Manager.GetDouble("nn-coupling");
  double NNNCoupling =  Manager.GetDouble("nnn-coupling");

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

  bool FirstRun = true;
  if ((InitalSzValue == 0) && (Manager.GetBoolean("disable-szsymmetry") == false) && (SpinValue <= 2))
    {
      for (int Momentum = InitialMomentum; Momentum < MaxMomentum; ++Momentum)
	{
	  for (int SzSymmetrySector = MinSzSymmetrySector; SzSymmetrySector <= MaxSzSymmetrySector; SzSymmetrySector += 2)
	    {
	      if ((Manager.GetBoolean("disable-inversionsymmetry") == false)  && (SpinValue <= 2) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins >> 1)))))
		{
		  for (int InversionSymmetrySector = MinInversionSymmetrySector; InversionSymmetrySector <= MaxInversionSymmetrySector; InversionSymmetrySector += 2)
		    {
		      AbstractSpinChainWithTranslations* Chain = 0;
		      switch (SpinValue)
			{
 			case 1 :
			  {
			    if (NbrSpins < 40)
			      {
				Chain = new Spin1_2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momentum, 1, InversionSymmetrySector, SzSymmetrySector, InitalSzValue, 1000000, 1000000);
			      }
			    else
			      {
				Chain = new Spin1_2ChainWithTranslationsAndSzInversionSymmetriesLong (NbrSpins, Momentum, 1, InversionSymmetrySector, SzSymmetrySector, InitalSzValue, 1000000, 1000000);
			      }
			  }
 			  break;
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
			      SpinChainRealHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, JValue, NNCoupling, NNNCoupling);

			      Hamiltonian.ShiftHamiltonian(Manager.GetDouble("shift-energies"));
			      if (Manager.GetBoolean("export-charpolynomial"))
				{
				  int* SU2Degeneracy = 0;
				  if (Manager.GetBoolean("export-sresolvedcharpolynomial"))
				    {			   
				      SU2Degeneracy = new int [(((SpinValue * NbrSpins) - InitalSzValue) >> 1) + 1];
				      for (int k = InitalSzValue; k <= (SpinValue * NbrSpins); k += 2)
					{
					  switch (SpinValue)
					    {
					    case 1 :
					      {
						if (NbrSpins < 40)
						  {
						    Spin1_2ChainWithTranslationsAndInversionSymmetry TmpHilbert(NbrSpins, Momentum, 1, InversionSymmetrySector, k, 1000000, 1000000);
						    SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
						  }
						else
						  {
						    Spin1_2ChainWithTranslationsAndInversionSymmetryLong TmpHilbert(NbrSpins, Momentum, 1, InversionSymmetrySector, k, 1000000, 1000000);
						    SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
						  }
					      }
					      break;
					    case 2 :
					      {
						Spin1ChainWithTranslationsAndInversionSymmetry TmpHilbert (NbrSpins, Momentum, InversionSymmetrySector, k);
						SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
					      }
					      break;
					    }
					}
				      for (int k = InitalSzValue; k < (SpinValue * NbrSpins); k += 2)
					{
					  SU2Degeneracy[(k - InitalSzValue) >> 1] -= SU2Degeneracy[1 + ((k - InitalSzValue) >> 1)];
					}
				      if (SzSymmetrySector == 1)
					{
					  for (int k = (SpinValue * NbrSpins) - 2; k >= InitalSzValue; k -= 4)
					    {
					      SU2Degeneracy[(k - InitalSzValue) >> 1] = 0;
					    }
					}
				      else
					{
					  for (int k = (SpinValue * NbrSpins); k >= InitalSzValue; k -= 4)
					    {
					      SU2Degeneracy[(k - InitalSzValue) >> 1] = 0;
					    }
					}
				    }
				  SpinChainComputeCharacteristicPolynomial(&Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture(), Manager.GetBoolean("export-sresolvedcharpolynomial"), InitalSzValue, SpinValue * NbrSpins, SU2Degeneracy, Manager.GetBoolean("charpolynomial-nofour"));
				}			      


			      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
							  FirstRun, TmpEigenstateString);
			      MainTaskOperation TaskOperation (&Task);
			      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			      Lanczos.SetComplexAlgorithms();
			    }
			  else
			    {
			      SpinChainHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, JValue, NNCoupling, NNNCoupling);
			      Hamiltonian.ShiftHamiltonian(Manager.GetDouble("shift-energies"));
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
		    case 1 :
		      {
			if (NbrSpins < 40)
			  {
			    Chain = new Spin1_2ChainWithTranslationsAndSzSymmetry (NbrSpins, Momentum, 1, SzSymmetrySector, InitalSzValue, 1000000, 1000000);
			  }
			else
			  {
			    Chain = new Spin1_2ChainWithTranslationsAndSzSymmetryLong (NbrSpins, Momentum, 1, SzSymmetrySector, InitalSzValue, 1000000, 1000000);
			  }
		      }
		      break;
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
		      SpinChainHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, JValue, NNCoupling, NNNCoupling);
		      Hamiltonian.ShiftHamiltonian(Manager.GetDouble("shift-energies"));
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
	  if ((Manager.GetBoolean("disable-inversionsymmetry") == false)  && (SpinValue <= 2) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins >> 1)))))
	    {
	      for (int InversionSymmetrySector = MinInversionSymmetrySector; InversionSymmetrySector <= MaxInversionSymmetrySector; InversionSymmetrySector += 2)
		{
		  AbstractSpinChainWithTranslations* Chain = 0;
		  switch (SpinValue)
		    {
 		    case 1 :
		      {
			if (NbrSpins < 40)
			  {
			    Chain = new Spin1_2ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momentum, 1, InversionSymmetrySector, InitalSzValue, 1000000, 1000000);
			  }
			else
			  {
			    Chain = new Spin1_2ChainWithTranslationsAndInversionSymmetryLong (NbrSpins, Momentum, 1, InversionSymmetrySector, InitalSzValue, 1000000, 1000000);
			  }
		      }
 		      break;
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
			  SpinChainRealHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, JValue, NNCoupling, NNNCoupling);
			  Hamiltonian.ShiftHamiltonian(Manager.GetDouble("shift-energies"));
			  if (Manager.GetBoolean("export-charpolynomial"))
			    {
			      int* SU2Degeneracy = 0;
			      if (Manager.GetBoolean("export-sresolvedcharpolynomial"))
				{			   
				  SU2Degeneracy = new int [(((SpinValue * NbrSpins) - InitalSzValue) >> 1) + 1];
				  for (int k = InitalSzValue; k <= (SpinValue * NbrSpins); k += 2)
				    {
				      switch (SpinValue)
					{
					case 1 :
					  {
					    if (NbrSpins < 40)
					      {
						Spin1_2ChainWithTranslationsAndInversionSymmetry TmpHilbert(NbrSpins, Momentum, 1, InversionSymmetrySector, k, 1000000, 1000000);
						SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
					      }
					    else
					      {
						Spin1_2ChainWithTranslationsAndInversionSymmetryLong TmpHilbert(NbrSpins, Momentum, 1, InversionSymmetrySector, k, 1000000, 1000000);
						SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
					      }
					  }
					  break;
					case 2 :
					  {
					    Spin1ChainWithTranslationsAndInversionSymmetry TmpHilbert (NbrSpins, Momentum, InversionSymmetrySector, k);
					    SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
					  }
					  break;
					}
				    }
				  for (int k = InitalSzValue; k < (SpinValue * NbrSpins); k += 2)
				    {
				      SU2Degeneracy[(k - InitalSzValue) >> 1] -= SU2Degeneracy[1 + ((k - InitalSzValue) >> 1)];
				    }
				}
			      SpinChainComputeCharacteristicPolynomial(&Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture(), Manager.GetBoolean("export-sresolvedcharpolynomial"), InitalSzValue,  (SpinValue * NbrSpins), SU2Degeneracy, Manager.GetBoolean("charpolynomial-nofour"));
			    }			      

			  // RealMatrix TmpTransformation (Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
			  // for (int k = 0; k < (Chain->GetHilbertSpaceDimension() / 2); ++k)
			  //   {
			  //     TmpTransformation.SetMatrixElement(k, k, 1.0);
			  //     TmpTransformation.SetMatrixElement(Chain->GetHilbertSpaceDimension() - 1 - k , k, 1.0);
			  //     TmpTransformation.SetMatrixElement(k, k + (Chain->GetHilbertSpaceDimension() / 2), 1.0);
			  //     TmpTransformation.SetMatrixElement(Chain->GetHilbertSpaceDimension() - 1 - k , k  + (Chain->GetHilbertSpaceDimension() / 2), -1.0);
			  //   }
			  // TmpTransformation /= M_SQRT2;
			  //			  cout << TmpTransformation << endl;
			  // RealMatrix TmpTransformation2 = TmpTransformation.DuplicateAndTranspose();
			  // RealMatrix TmpTransformation3;
			  // TmpTransformation3 = TmpTransformation * TmpTransformation2;
			  // cout << TmpTransformation3 << endl;

			  // RealSymmetricMatrix TmpHamiltonian(Chain->GetHilbertSpaceDimension(), true);
			  // RealSymmetricMatrix TmpConjHamiltonian(Chain->GetHilbertSpaceDimension(), true);
			  // Hamiltonian.GetHamiltonian(TmpHamiltonian);
			  // cout << TmpHamiltonian << endl;			  
			  // TmpHamiltonian.Conjugate(TmpTransformation, 0, 0, TmpConjHamiltonian);
			  // cout << TmpConjHamiltonian << endl;
			  
			  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
						   FirstRun, TmpEigenstateString);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			  Lanczos.SetComplexAlgorithms();
			}
		      else
			{
			  SpinChainHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, JValue, NNCoupling, NNNCoupling);
			  Hamiltonian.ShiftHamiltonian(Manager.GetDouble("shift-energies"));
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
		case 1 :
		  {
		    if (NbrSpins < 40)
		      {
			Chain = new Spin1_2ChainWithTranslations (NbrSpins, Momentum, 1, InitalSzValue, 1000000, 1000000);
		      }
		    else
		      {
			Chain = new Spin1_2ChainWithTranslationsLong (NbrSpins, Momentum, 1, InitalSzValue, 1000000, 1000000);
		      }
		  }
		  break;
		case 2 :
		  Chain = new Spin1ChainWithTranslations (NbrSpins, Momentum, InitalSzValue);
		  break;
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
		  SpinChainHamiltonianWithTranslations Hamiltonian (Chain, NbrSpins, JValue, NNCoupling, NNNCoupling);
		  Hamiltonian.ShiftHamiltonian(Manager.GetDouble("shift-energies"));
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


// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// sResolvedFlag = true if the characteristic polynomial should be computed per S sector
// minSValue = twice the minimum S value to consider (when sResolvedFlag is true)
// minSValue = twice the maximum S value to consider (when sResolvedFlag is true)
// sU2Degeneracy = array providing the Hilbert space dimension per S value (when sResolvedFlag is true)
// discardFourFactor = discard a global four factor used to ensure integer numbers

void SpinChainComputeCharacteristicPolynomial(SpinChainRealHamiltonianWithTranslations* hamiltonian, AbstractSpinChainWithTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, bool sResolvedFlag, int minSValue, int maxSValue, int* sU2Degeneracy, bool discardFourFactor)
{
#ifdef __GMP__
  if (sResolvedFlag == false)
    {
      cout << "Computing the hamiltonian" << endl;
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
      cout << "Converting to integer matrix" << endl;
      LongIntegerMatrix TmpMatrix;
      if (discardFourFactor == false)
	{
	  TmpMatrix = LongIntegerMatrix(TmpRawMatrix, 4.0);
	}
      else
	{
	  TmpMatrix = LongIntegerMatrix(TmpRawMatrix);
	}
      cout << "Start computing characteristic polynomial (degree " << chain->GetHilbertSpaceDimension() << ")" << endl;
      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(architecture);
      char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
      if (discardFourFactor == false)
	{
	  sprintf (PolynomialOutputFileName, "%s.charpol", outputFileName);
	}
      else
	{
	  sprintf (PolynomialOutputFileName, "%s.no4.charpol", outputFileName);
	}
      ofstream OutputFile;
      OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
      OutputFile << CharacteristicPolynomial[0];
      for (int i = 1; i <= chain->GetHilbertSpaceDimension(); ++i)
	{
	  OutputFile << "," << CharacteristicPolynomial[i];
	}
      OutputFile << endl;
      OutputFile.close();
    }
  else
    {
      cout << "Computing the hamiltonian and the S^2 operator" << endl;
      SpinWith1DTranslationS2Operator S2Operator(chain, chain->GetSpinChainLength());
      RealMatrix TmpRawMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      hamiltonian->GetHamiltonian(TmpRawMatrix);
      RealMatrix TmpRawMatrixS2(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      S2Operator.GetOperator(TmpRawMatrixS2);
      LongIntegerMatrix TmpMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      LongIntegerMatrix TmpMatrix2(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      LongIntegerMatrix TmpProjector(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      LongIntegerMatrix TmpMatrixS2(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
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
	      TmpRawMatrixS2.GetMatrixElement(i, j, Tmp);
	      Tmp *= TmpNormalizationFactors[i];
	      Tmp /= TmpNormalizationFactors[j];
	      TmpRawMatrixS2.SetMatrixElement(i, j, Tmp);
	    }
	}

      for (int s = minSValue; s <= maxSValue; s += 2)
	{
	  if (sU2Degeneracy[(s - minSValue) >> 1] > 0)
	    {
	      cout << "computing 2S=" << s << " sector" << endl;
	      cout << "  building projected hamiltonian" << endl;
	      TmpProjector.SetToIdentity();
	      long TmpNormalisation = 1l;
	      for (int j = minSValue; j <= maxSValue; j += 2)
		{
		  if (j != s)
		    {
		      long TmpShift;
		      if ((minSValue & 1) == 0)
			{
			  TmpMatrixS2 = LongIntegerMatrix(TmpRawMatrixS2);
			  TmpShift = ((long) (-j * (j + 2))) >> 2;
			}
		      else
			{
			  TmpMatrixS2 = LongIntegerMatrix(TmpRawMatrixS2, 4.0);
			  TmpShift = ((long) (-j * (j + 2)));
			}
		      for (int k = 0; k < chain->GetHilbertSpaceDimension(); ++k)
			{
			  TmpMatrixS2.AddToMatrixElement(k, k, TmpShift);
			}
		      TmpProjector.Multiply(TmpMatrixS2);
		      if ((minSValue & 1) == 0)
			{
			  TmpNormalisation *= ((long) ((s * (s + 2)) - (j * (j + 2)))) >> 2;
			}
		      else
			{
			  TmpNormalisation *= ((long) ((s * (s + 2)) - (j * (j + 2))));
			}
		    }
		}
	      if (discardFourFactor == false)
		{
		  TmpMatrix = LongIntegerMatrix(TmpRawMatrix, 4.0);
		}
	      else
		{
		  TmpMatrix = LongIntegerMatrix(TmpRawMatrix);
		}
	      TmpMatrix2 = TmpMatrix * TmpProjector;
	      TmpMatrix = TmpProjector * TmpMatrix2;
	      TmpMatrix /= (TmpNormalisation);

	      cout << "Start computing characteristic polynomial (degree " << chain->GetHilbertSpaceDimension() << ")" << endl;
	      mpz_t TmpNormalisation2;
	      mpz_init_set_si(TmpNormalisation2, TmpNormalisation);
	      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(architecture);
	      for (int i = 1; i <= sU2Degeneracy[(s - minSValue) >> 1]; ++i)
		{		  
		  for (int j = 0; j < i; ++j)
		    {
		      mpz_divexact(CharacteristicPolynomial[chain->GetHilbertSpaceDimension() - i], CharacteristicPolynomial[chain->GetHilbertSpaceDimension() - i], TmpNormalisation2);
		    }
		}
	      cout << "  testing polynomial (checking if compatible with 2S=" << s << " with dim=" << sU2Degeneracy[(s - minSValue) >> 1] << ")" << endl;
	      int PolynomialTestingFlag = 0;
	      for (int i = 0; i < (chain->GetHilbertSpaceDimension() - sU2Degeneracy[(s - minSValue) >> 1]); ++i)
		{
		  if (mpz_sgn( CharacteristicPolynomial[i]) != 0)
		    {
		      PolynomialTestingFlag++;
		    }
		}
	      if (PolynomialTestingFlag == 0)
		{
		  cout << "  all clear" << endl;
		}
	      else
		{
		  cout << "  " << PolynomialTestingFlag << " error(s) detected" << endl;
		}
	      
	      char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
	      if (discardFourFactor == false)
		{
		  sprintf (PolynomialOutputFileName, "%s_s_%d.charpol", outputFileName, s);
		}
	      else
		{
		  sprintf (PolynomialOutputFileName, "%s_s_%d.no4.charpol", outputFileName, s);
		}		
	      ofstream OutputFile;
	      OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
	      OutputFile << CharacteristicPolynomial[chain->GetHilbertSpaceDimension() - sU2Degeneracy[(s - minSValue) >> 1]];
	      for (int i = chain->GetHilbertSpaceDimension() - sU2Degeneracy[(s - minSValue) >> 1] + 1; i <= chain->GetHilbertSpaceDimension(); ++i)
		{
		  OutputFile << "," << CharacteristicPolynomial[i];
		}
	      OutputFile << endl;
	      OutputFile.close();
	    }
	}
    }
#else
  cout << "GMP library is required for characteristic polynomials" << endl;
#endif	       
}
