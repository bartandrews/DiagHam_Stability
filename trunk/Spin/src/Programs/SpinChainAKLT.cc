#include "Hamiltonian/SpinChainAKLTHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1ChainWithSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithSzInversionSymmetries.h"
#include "HilbertSpace/Spin3_2Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Operator/SpinS2Operator.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

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


// evaluate Hilbert space dimension for a spin 1 system with a fixed total Sz
//
// totalSz = twice the total Sz value
// sz = twice the current Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension
long ComputeSpin1HilbertSpaceDimension(int totalSz, int sz, int nbrSites);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainAKLT" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 2);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 8);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "projector-normalization", "normalize the hamiltonian as a pure sum of projector");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "additional-quadratic", "coefficient in front of the additional quadratic term (0 being the pure AKLT hamiltonian)", 0.0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "export-charpolynomial", "export the hamiltonian characteristic polynomial");  
  (*ToolsGroup) += new BooleanOption  ('\n', "export-sresolvedcharpolynomial", "export the hamiltonian characteristic polynomial resolved in S quantum number");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainAKLT -h" << endl;
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
  char* BoundaryName = new char [64];
  if (Manager.GetBoolean("projector-normalization") == false)
    {
      if (Manager.GetBoolean("use-periodic") == false)
	sprintf (BoundaryName, "open");
      else
	sprintf (BoundaryName, "closed");
    }
  else
    {
      if (Manager.GetBoolean("use-periodic") == false)
	sprintf (BoundaryName, "open_projnorm");
      else
	sprintf (BoundaryName, "closed_projnorm");
    }
  if ((SpinValue & 1) == 0)
    {
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_%saklt_quadratic_%.6f_n_%d", (SpinValue / 2), BoundaryName, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_%saklt_n_%d", (SpinValue / 2), BoundaryName, NbrSpins);
	}
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz SzSym InvSym ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz SzSym ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz InvSym ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	}
    }
  else
    {
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_2_%saklt_quadratic_%.6f_n_%d", SpinValue, BoundaryName, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_2_%saklt_n_%d", SpinValue, BoundaryName, NbrSpins);
	}
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz SzSym InvSym ", BoundaryName, SpinValue, NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz SzSym ", BoundaryName, SpinValue, NbrSpins);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz InvSym ", BoundaryName, SpinValue, NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz ", BoundaryName, SpinValue, NbrSpins);
	    }
	}
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if (Manager.GetInteger("initial-sz") != 0)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;

  if ((InitalSzValue == 0) && (Manager.GetBoolean("disable-szsymmetry") == false))
    {
      for (int SzSymmetrySector = -1; SzSymmetrySector <= 1; SzSymmetrySector += 2)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
		{
		  AbstractSpinChain* Chain = 0;
		  switch (SpinValue)
		    {
		    case 2 :
		      Chain = new Spin1ChainWithSzInversionSymmetries (NbrSpins, InversionSymmetrySector, SzSymmetrySector, InitalSzValue, 1000000);
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
		      cout << "Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << ", inversion sector=" << InversionSymmetrySector << endl; 
		      SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
		      if (Manager.GetBoolean("projector-normalization"))
			{
			  if (Manager.GetBoolean("use-periodic"))
			    {
			      Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins)) / 3.0);
			    }
			  else
			    {
			      Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins - 2)) / 3.0);
			    }
			}
		      char* TmpSzString = new char[64];
		      sprintf (TmpSzString, "%d %d %d", InitalSzValue, SzSymmetrySector, InversionSymmetrySector);
		      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		      sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d_szsym_%d", OutputFileName, InitalSzValue, InversionSymmetrySector, SzSymmetrySector);
		      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					       FirstRun, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      FirstRun = false;
		      delete[] TmpSzString;
		    }
		  delete Chain;
		}
	    }
	  else
	    {
	      AbstractSpinChain* Chain = 0;
	      switch (SpinValue)
		{
		case 2 :
		  Chain = new Spin1ChainWithSzSymmetry (NbrSpins, SzSymmetrySector, InitalSzValue, 1000000);
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
		  cout << "Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << endl; 
		  SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
		  if (Manager.GetBoolean("projector-normalization"))
		    {
		      if (Manager.GetBoolean("use-periodic"))
			{
			  Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins)) / 3.0);
			}
		      else
			{
			  Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins - 2)) / 3.0);
			}
		    }
		  char* TmpSzString = new char[64];
		  sprintf (TmpSzString, "%d %d", InitalSzValue, SzSymmetrySector);
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_sz_%d_szsym_%d", OutputFileName, InitalSzValue, SzSymmetrySector);
		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete[] TmpSzString;
		}
	      delete Chain;
	    }
	}
      InitalSzValue += 2;
    }
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	{
	  for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
	    {
	      AbstractSpinChain* Chain = 0;
	      switch (SpinValue)
		{
		case 2 :
		  Chain = new Spin1ChainWithInversionSymmetry (NbrSpins, InversionSymmetrySector, InitalSzValue, 1000000);
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
		  cout << "Sz = " << InitalSzValue << ", inversion sector=" << InversionSymmetrySector << endl; 
		  SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
		  if (Manager.GetBoolean("projector-normalization"))
		    {
		      if (Manager.GetBoolean("use-periodic"))
			{
			  Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins)) / 3.0);
			}
		      else
			{
			  Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins - 2)) / 3.0);
			}
		    }
		  char* TmpSzString = new char[64];
		  if (Manager.GetBoolean("disable-szsymmetry") == false)
		    {
		      sprintf (TmpSzString, "%d 0 %d", InitalSzValue, InversionSymmetrySector);
		    }
		  else
		    {
		      sprintf (TmpSzString, "%d %d", InitalSzValue, InversionSymmetrySector);
		    }
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d", OutputFileName, InitalSzValue, InversionSymmetrySector);
		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete[] TmpSzString;
		}
	      delete Chain;
	    }
	}
      else
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
	    case 3 :
	      Chain = new Spin3_2Chain (NbrSpins, InitalSzValue, 1000000);
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
	      cout << "Sz = " << InitalSzValue << endl; 
	      SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
	      if (Manager.GetBoolean("projector-normalization"))
		{
		  if (Manager.GetBoolean("use-periodic"))
		    {
		      Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins)) / 3.0);
		    }
		  else
		    {
		      Hamiltonian.ShiftHamiltonian(((double) (2 * NbrSpins - 2)) / 3.0);
		    }
		}
	      char* TmpSzString = new char[64];
	      if (Manager.GetBoolean("disable-szsymmetry") == false)
		{
		  sprintf (TmpSzString, "%d 0", InitalSzValue);
		}
	      else
		{
		  sprintf (TmpSzString, "%d", InitalSzValue);
		}
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);


	      if ((Manager.GetBoolean("export-charpolynomial")) || (Manager.GetBoolean("export-sresolvedcharpolynomial")))
		{
		  if (Manager.GetBoolean("export-sresolvedcharpolynomial"))
		    {
#ifdef __GMP__
		      SpinS2Operator S2Operator(Chain, NbrSpins);
		      LongIntegerMatrix TmpMatrixS2(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		      LongIntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		      LongIntegerMatrix TmpMatrix2(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		      LongIntegerMatrix TmpProjector(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		      int MinSValue = 0;
		      int MaxSValue = 2 * NbrSpins;	  
		      long* SU2Degeneracy = new long [((MaxSValue - MinSValue) / 2) + 1];
		      for (int k = MinSValue; k <= MaxSValue; k += 2)
			{
			  SU2Degeneracy[(k - MinSValue) >> 1] = ComputeSpin1HilbertSpaceDimension(k, 0, NbrSpins);
			}
		      for (int k = MinSValue; k < MaxSValue; k += 2)
			{
			  SU2Degeneracy[(k - MinSValue) >> 1] -= SU2Degeneracy[1 + ((k - MinSValue) >> 1)];
			}
		      for (int s = MinSValue; s <= MaxSValue; s += 2)
			{
			  cout << "computing 2S=" << s << " sector" << endl;
			  cout << "  building projected hamiltonian" << endl;
			  TmpProjector.SetToIdentity();
			  long TmpNormalisation = 1l;
			  for (int j = MinSValue; j <= MaxSValue; j += 2)
			    {
			      if (j != s)
				{
				  long TmpShift;
				  S2Operator.GetOperator(TmpMatrixS2);
				  TmpShift = ((long) (-j * (j + 2))) >> 2;
				  for (int k = 0; k < Chain->GetHilbertSpaceDimension(); ++k)
				    {
				      TmpMatrixS2.AddToMatrixElement(k, k, TmpShift);
				    }
				  TmpProjector.Multiply(TmpMatrixS2);
				  TmpNormalisation *= ((long) ((s * (s + 2)) - (j * (j + 2)))) >> 2;
				}
			    }
			  
			  Hamiltonian.GetHamiltonian(TmpMatrix, 3.0);
			  TmpMatrix2 = TmpMatrix * TmpProjector;
			  TmpMatrix = TmpProjector * TmpMatrix2;
			  TmpMatrix /= (TmpNormalisation);
			  //			  cout << TmpMatrix << endl;
			  cout << "  simplifying projected hamiltonian" << endl;
			  cout << "    nbr null columns=" << TmpMatrix.NbrNullColumns() << endl;
			  int* TmpProportionalColums = new int [TmpMatrix.GetNbrColumn()];
			  for (int i = 0; i < TmpMatrix.GetNbrColumn(); ++i)
			    {
			      TmpProportionalColums[i] = -1;
			    }
			  int TmpNbrProportionalColums = 0;
			  for (int i = 0; i < (TmpMatrix.GetNbrColumn() - 1); ++i)
			    {
			      if (TmpProportionalColums[i] == -1)
				{
				  for (int j = i + 1; j < TmpMatrix.GetNbrColumn(); ++j)
				    {
				      if (TmpProportionalColums[j] == -1)
					{
					  if (TmpMatrix[i].IsProportional(TmpMatrix[j]) == true)
					    {
					      TmpProportionalColums[j] = i;
					      ++TmpNbrProportionalColums;
					    }
					}
				    }
				}
			    }
			  cout << "    nbr proportional columns=" << TmpNbrProportionalColums << endl;
			  cout << "  computing characteristic polynomial" << endl;
			  mpz_t TmpNormalisation2;
			  mpz_init_set_si(TmpNormalisation2, TmpNormalisation);
			  mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();
			  for (int i = 1; i <= SU2Degeneracy[(s - MinSValue) >> 1]; ++i)
			    {		  
			      for (int j = 0; j < i; ++j)
				{
				  mpz_divexact(CharacteristicPolynomial[Chain->GetHilbertSpaceDimension() - i], CharacteristicPolynomial[Chain->GetHilbertSpaceDimension() - i], TmpNormalisation2);
				}
			    }
			  cout << "  testing polynomial (checking if compatible with 2S=" << s << " with dim=" << SU2Degeneracy[(s - MinSValue) >> 1] << ")" << endl;
			  int PolynomialTestingFlag = 0;
			  for (int i = 0; i < (Chain->GetHilbertSpaceDimension() - SU2Degeneracy[(s - MinSValue) >> 1]); ++i)
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
			  
			  char* PolynomialOutputFileName = new char[strlen(OutputFileName) + 256];
			  sprintf (PolynomialOutputFileName, "%s_sz_%d_s_%d.charpol", OutputFileName, InitalSzValue, s);
			  ofstream OutputFile;
			  OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
			  OutputFile << CharacteristicPolynomial[Chain->GetHilbertSpaceDimension() - SU2Degeneracy[(s - MinSValue) >> 1]];
			  for (int i = Chain->GetHilbertSpaceDimension() - SU2Degeneracy[(s - MinSValue) >> 1] + 1; i <= Chain->GetHilbertSpaceDimension(); ++i)
			    {
			      OutputFile << "," << CharacteristicPolynomial[i];
			    }
			  OutputFile << endl;
			}
#else
		      cout << "GMP library is required" << endl;
#endif
		    }
		  else
		    {
#ifdef __GMP__
		      LongIntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		      Hamiltonian.GetHamiltonian(TmpMatrix, 3.0);
		      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomialAssumingSymmetric();
#else
		      IntegerMatrix TmpMatrix(Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
		      Hamiltonian.GetHamiltonian(TmpMatrix);
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
		}
	      
	      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
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


// evaluate Hilbert space dimension for a spin 1 system with a fixed total Sz
//
// totalSz = twice the total Sz value
// sz = twice the current Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long ComputeSpin1HilbertSpaceDimension(int totalSz, int sz, int nbrSites)
{
  if (nbrSites == 0)
    {
      if (sz == totalSz)
	{
	  return 1l;	  
	}
      else
	{
	  return 0l;	  
	}
    }
  if (((sz + 2 * nbrSites) < totalSz) || ((sz - 2 * nbrSites) > totalSz))
    {
      return 0l;
    }
  long TmpDimension = ComputeSpin1HilbertSpaceDimension(totalSz, sz - 2, nbrSites - 1);
  TmpDimension += ComputeSpin1HilbertSpaceDimension(totalSz, sz, nbrSites - 1);
  TmpDimension += ComputeSpin1HilbertSpaceDimension(totalSz, sz + 2, nbrSites - 1);
  return TmpDimension;
}
