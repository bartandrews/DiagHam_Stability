#include "Hamiltonian/SpinChainAKLTHamiltonianWithTranslations.h"
#include "Hamiltonian/SpinChainAKLTRealHamiltonianWithTranslations.h"

#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"

#include "Operator/SpinWith1DTranslationS2Operator.h"

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
// sResolvedFlag = true if the characteristic polynomial should be computed per S sector
// minSValue = twice the minimum S value to consider (when sResolvedFlag is true)
// minSValue = twice the maximum S value to consider (when sResolvedFlag is true)
// globalMinSValue = twice the minimum S value that can be reached by the system
// globalMaxSValue = twice the maximum S value that can be reached by the system
// sU2Degeneracy = array providing the Hilbert space dimension per S value (when sResolvedFlag is true)
void SpinChainAKLTComputeCharacteristicPolynomial(SpinChainAKLTRealHamiltonianWithTranslations* hamiltonian, AbstractSpinChainWithTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, bool sResolvedFlag, int minSValue = 0, int maxSValue = 0, int globalMinSValue = 0, int globalMaxSValue = 0, int* sU2Degeneracy = 0);



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

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 2);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "additional-quadratic", "coefficient in front of the additional quadratic term (0 being the pure AKLT hamiltonian)", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "linear-factor", "if different from 1.0, set the coefficient in front of the Heisenberg term", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "quadratic-factor", "if --linear-factor is different from 1.0, set the coefficient in front of the quadratic term", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-szsymmetry", "if non zero, set the Sz<->-Sz symmetry sector", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "set-inversionsymmetry", "if non zero, set the inversion symmetry sector", 0);
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
  (*ToolsGroup) += new BooleanOption  ('\n', "export-sresolvedcharpolynomial", "export the hamiltonian characteristic polynomial resolved in S quantum number");  (*ToolsGroup) += new  SingleIntegerOption ('\n', "sresolvedcharpolynomial-svalue", "when exporting the hamiltonian characteristic polynomial resolved in S quantum number, computing only one S sector (compute all of them if negative, should be set to twice S to be valid)", -1);  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PeriodicSpinChainAKLT -h" << endl;
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
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_periodicaklt_quadratic_%.6f_n_%d", (SpinValue / 2), 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  if (Manager.GetDouble("linear-factor") != 1.0)
	    {
	      sprintf (OutputFileName, "spin_%d_periodicaklt_linear_%.6f_quadratic_%.6f_n_%d", (SpinValue / 2), 
		       Manager.GetDouble("linear-factor"),  Manager.GetDouble("quadratic-factor"), NbrSpins);
	    }
	  else
	    {
	      sprintf (OutputFileName, "spin_%d_periodicaklt_n_%d", (SpinValue / 2), NbrSpins);
	    }
	}
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
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_2_periodicaklt_quadratic_%.6f_n_%d", SpinValue, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_2_periodicaklt_n_%d", SpinValue, NbrSpins);
	}
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
			  int* SU2Degeneracy = 0;
			  if (Manager.GetBoolean("export-sresolvedcharpolynomial"))
			    {			   
			      SU2Degeneracy = new int [NbrSpins + 1];
			      for (int k = 0; k <= (2 * NbrSpins); k += 2)
				{
				  Spin1ChainWithTranslationsAndInversionSymmetry TmpHilbert (NbrSpins, Momentum, InversionSymmetrySector, k);
				  SU2Degeneracy[k >> 1] = TmpHilbert.GetHilbertSpaceDimension();
				}
			      for (int k = 0; k < NbrSpins; k++)
				{
				  SU2Degeneracy[k] -= SU2Degeneracy[k + 1];
				}
			      if (SzSymmetrySector == 1)
				{
				  for (int k = NbrSpins - 1; k >= 0; k -= 2)
				    {
				      SU2Degeneracy[k] = 0;
				    }
				}
			      else
				{
				  for (int k = NbrSpins; k >= 0; k -= 2)
				    {
				      SU2Degeneracy[k] = 0;
				    }
				}
			    }
			  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
			  cout << "2Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << ",   inversion sector=" << InversionSymmetrySector << ",  K = " << Momentum << endl; 
			  char* TmpSzString = new char[64];
			  sprintf (TmpSzString, "%d %d %d %d ", InitalSzValue, Momentum, SzSymmetrySector, InversionSymmetrySector);
			  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
			  sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d_szsym_%d_k_%d", OutputFileName, InitalSzValue, InversionSymmetrySector, SzSymmetrySector, Momentum);
			  if (Manager.GetBoolean("disable-realhamiltonian") == false)
			    {
			      Lanczos.SetRealAlgorithms();
			      SpinChainAKLTRealHamiltonianWithTranslations* Hamiltonian = 0;
			      if (Manager.GetDouble("linear-factor") == 1.0)
				{
				  Hamiltonian = new SpinChainAKLTRealHamiltonianWithTranslations(Chain, NbrSpins, 
												 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"));
				}
			      else
				{
				  Hamiltonian = new SpinChainAKLTRealHamiltonianWithTranslations(Chain, NbrSpins, Manager.GetDouble("linear-factor"),
												 Manager.GetDouble("quadratic-factor"));
				}

			      if (Manager.GetBoolean("export-charpolynomial"))
				{
				  int MinSValue = InitalSzValue;
				  int MaxSValue = 2 * NbrSpins;
				  if (Manager.GetInteger("sresolvedcharpolynomial-svalue") >= 0)
				    {
				      MinSValue = Manager.GetInteger("sresolvedcharpolynomial-svalue");
				      if ((MinSValue < InitalSzValue) || (MinSValue > (2 * NbrSpins)))
					{
					  cout << "warning, invalid S value provided by --sresolvedcharpolynomial-svalue" << endl;
					  MinSValue = InitalSzValue;
					}
				      else
					{
					  MaxSValue = MinSValue;
					}
				    }
				  SpinChainAKLTComputeCharacteristicPolynomial(Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture(), Manager.GetBoolean("export-sresolvedcharpolynomial"), MinSValue, MaxSValue, 0, 2 * NbrSpins, SU2Degeneracy);
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
			      SpinChainAKLTHamiltonianWithTranslations* Hamiltonian;
			      if(Manager.GetDouble("linear-factor") == 1.0)
				{ 
				  Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations  (Chain, NbrSpins, 
											       1.0 + 3.0 * Manager.GetDouble("additional-quadratic"));
				}
			      else
				{
				  Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations (Chain, NbrSpins, Manager.GetDouble("linear-factor"),
											      Manager.GetDouble("quadratic-factor"));
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
		      SpinChainAKLTHamiltonianWithTranslations* Hamiltonian = 0;
		      if(Manager.GetDouble("linear-factor") == 1.0)
			{			  
			  Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations(Chain, NbrSpins, 
										     1.0 + 3.0 * Manager.GetDouble("additional-quadratic"));
			}
		      else
			{
			  Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations(Chain, NbrSpins, Manager.GetDouble("linear-factor"),
										     Manager.GetDouble("quadratic-factor"));
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
		      int* SU2Degeneracy = 0;
		      if (Manager.GetBoolean("export-sresolvedcharpolynomial"))
			{			   
			  SU2Degeneracy = new int [((2 * NbrSpins) - InitalSzValue) >> 1 + 1];
			  for (int k = InitalSzValue; k <= (2 * NbrSpins); k += 2)
			    {
			      Spin1ChainWithTranslationsAndInversionSymmetry TmpHilbert (NbrSpins, Momentum, InversionSymmetrySector, k);
			      SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
			    }
			  for (int k = InitalSzValue; k < (2 * NbrSpins); k += 2)
			    {
			      SU2Degeneracy[(k - InitalSzValue) >> 1] -= SU2Degeneracy[1 + ((k - InitalSzValue) >> 1)];
			    }
			}
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
			  SpinChainAKLTRealHamiltonianWithTranslations* Hamiltonian = 0;
			  if (Manager.GetDouble("linear-factor") == 1.0)
			    {
			      Hamiltonian = new SpinChainAKLTRealHamiltonianWithTranslations(Chain, NbrSpins, 
											     1.0 + 3.0 * Manager.GetDouble("additional-quadratic"));
			    }
			  else
			    {
			      Hamiltonian = new SpinChainAKLTRealHamiltonianWithTranslations(Chain, NbrSpins, Manager.GetDouble("linear-factor"),
											     Manager.GetDouble("quadratic-factor"));
			    }
			  
			  if (Manager.GetBoolean("export-charpolynomial"))
			    {
			      int MinSValue = InitalSzValue;
			      int MaxSValue = 2 * NbrSpins;
			      if (Manager.GetInteger("sresolvedcharpolynomial-svalue") >= 0)
				{
				  MinSValue = Manager.GetInteger("sresolvedcharpolynomial-svalue");
				  if ((MinSValue < InitalSzValue) || (MinSValue > (2 * NbrSpins)))
				    {
				      cout << "warning, invalid S value provided by --sresolvedcharpolynomial-svalue" << endl;
				      MinSValue = InitalSzValue;
				    }
				  else
				    {
				      MaxSValue = MinSValue;
				    }
				}
			      SpinChainAKLTComputeCharacteristicPolynomial(Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture(), Manager.GetBoolean("export-sresolvedcharpolynomial"), MinSValue, MaxSValue, 0, 2 * NbrSpins, SU2Degeneracy);
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
			  SpinChainAKLTHamiltonianWithTranslations* Hamiltonian ;
			  if (Manager.GetDouble("linear-factor") == 1.0)
			    {
			      Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations(Chain, NbrSpins, 
											 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"));
			    }
			  else
			    {
			      Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations(Chain, NbrSpins, Manager.GetDouble("linear-factor"),
											 Manager.GetDouble("quadratic-factor"));
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
		  {
		    Chain = new Spin1ChainWithTranslations (NbrSpins, Momentum, InitalSzValue);
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
	      
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{
		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		  cout << "2Sz = " << InitalSzValue << ", K = " << Momentum << endl; 
		  SpinChainAKLTHamiltonianWithTranslations* Hamiltonian = 0;
		  if (Manager.GetDouble("linear-factor") == 1.0)
		    {
		      Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations(Chain, NbrSpins, 
										 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"));
		    }
		  else
		    {
		      Hamiltonian = new SpinChainAKLTHamiltonianWithTranslations(Chain, NbrSpins, Manager.GetDouble("linear-factor"),
										 Manager.GetDouble("quadratic-factor"));
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
// globalMinSValue = twice the minimum S value that can be reached by the system
// globalMaxSValue = twice the maximum S value that can be reached by the system
// sU2Degeneracy = array providing the Hilbert space dimension per S value (when sResolvedFlag is true)

void SpinChainAKLTComputeCharacteristicPolynomial(SpinChainAKLTRealHamiltonianWithTranslations* hamiltonian, AbstractSpinChainWithTranslations* chain, char* outputFileName, AbstractArchitecture* architecture, bool sResolvedFlag, int minSValue, int maxSValue, int globalMinSValue, int globalMaxSValue, int* sU2Degeneracy)
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
      LongIntegerMatrix TmpMatrix(TmpRawMatrix, 3.0);
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
	      for (int j = globalMinSValue; j <= globalMaxSValue; j += 2)
		{
		  if (j != s)
		    {
		      long TmpShift;
		      TmpMatrixS2 = LongIntegerMatrix(TmpRawMatrixS2);
		      TmpShift = ((long) (-j * (j + 2))) >> 2;
		      for (int k = 0; k < chain->GetHilbertSpaceDimension(); ++k)
			{
			  TmpMatrixS2.AddToMatrixElement(k, k, TmpShift);
			}
		      TmpProjector.Multiply(TmpMatrixS2);
		      TmpNormalisation *= ((long) ((s * (s + 2)) - (j * (j + 2)))) >> 2;
		    }
		}
	      TmpMatrix = LongIntegerMatrix (TmpRawMatrix, 3.0);
	      TmpMatrix2 = TmpMatrix * TmpProjector;
	      TmpMatrix = TmpProjector * TmpMatrix2;
	      TmpMatrix /= (TmpNormalisation);

	      // cout << "  simplifying projected hamiltonian" << endl;
	      // cout << "    nbr null columns=" << TmpMatrix.NbrNullColumns() << endl;
	      // int* TmpProportionalColums = new int [TmpMatrix.GetNbrColumn()];
	      // for (int i = 0; i < TmpMatrix.GetNbrColumn(); ++i)
	      // 	{
	      // 	  TmpProportionalColums[i] = -1;
	      // 	}
	      // int TmpNbrProportionalColums = 0;
	      // for (int i = 0; i < (TmpMatrix.GetNbrColumn() - 1); ++i)
	      // 	{
	      // 	  if (TmpProportionalColums[i] == -1)
	      // 	    {
	      // 	      for (int j = i + 1; j < TmpMatrix.GetNbrColumn(); ++j)
	      // 		{
	      // 		  if (TmpProportionalColums[j] == -1)
	      // 		    {
	      // 		      if (TmpMatrix[i].IsProportional(TmpMatrix[j]) == true)
	      // 			{
	      // 			  TmpProportionalColums[j] = i;
	      // 			  ++TmpNbrProportionalColums;
	      // 			}
	      // 		    }
	      // 		}
	      // 	    }
	      // 	}
	      // cout << "    nbr proportional columns=" << TmpNbrProportionalColums << endl;

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
	      sprintf (PolynomialOutputFileName, "%s_s_%d.charpol", outputFileName, s);
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
