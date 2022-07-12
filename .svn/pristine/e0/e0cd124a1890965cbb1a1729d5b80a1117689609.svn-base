#include "Hamiltonian/PairHoppingHamiltonian.h"

#include "HilbertSpace/PairHoppingP1AsSpin1Chain.h"
#include "HilbertSpace/PairHoppingP2AsSpin1Chain.h"

#include "HilbertSpace/PairHoppingP1AsSpin1ChainLong.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainLong.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Matrix/RealMatrix.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"

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
void PairHoppingComputeCharacteristicPolynomial(PairHoppingHamiltonian* hamiltonian, AbstractSpinChain* chain, char* outputFileName, AbstractArchitecture* architecture);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("PairHoppingModel" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('p', "p-value", "value that defines the filling factor p/(2p+1)", 1);
  (*SystemGroup) += new  SingleIntegerOption ('n', "nbr-spin", "number of spin 1 (should be a multiple of p)", 4);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "export-charpolynomial", "export the hamiltonian characteristic polynomial");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PairHoppingModel -h" << endl;
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
   
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  sprintf (OutputFileName, "spin_1_%s_pairhopping_p_%d_n_%d", BoundaryName, PValue, NbrSpins);
  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
    {
      sprintf (CommentLine, " %s pair hopping model with p=%d in spin 1 language with %d sites \n# InvSym ", BoundaryName, PValue, NbrSpins);
    }
  else
    {
      sprintf (CommentLine, " %s pair hopping model with p=%d in spin 1 language with %d sites \n# ", BoundaryName, PValue, NbrSpins);
    }

  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  bool FirstRun = true;

  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
    {
      // for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
      // 	{
      // 	  AbstractSpinChain* Chain = 0;
      // 	  Chain = new Spin1ChainWithSzInversionSymmetries (NbrSpins, InversionSymmetrySector, SzSymmetrySector, InitalSzValue, 1000000);
	  
      // 	  if (Chain->GetHilbertSpaceDimension() > 0)
      // 	    {
      // 	      cout << "inversion sector=" << InversionSymmetrySector << endl; 
      // 	      SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
      // 	      char* TmpSzString = new char[64];
      // 	      sprintf (TmpSzString, "%d %d %d", InitalSzValue, SzSymmetrySector, InversionSymmetrySector);
      // 	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      // 	      sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d_szsym_%d", OutputFileName, InitalSzValue, InversionSymmetrySector, SzSymmetrySector);
      // 	      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
      // 				       FirstRun, TmpEigenstateString);
      // 	      MainTaskOperation TaskOperation (&Task);
      // 	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      // 	      FirstRun = false;
      // 	      delete[] TmpSzString;
      // 	    }
      // 	  delete Chain;
      // 	}
    }
  else
    {
      AbstractSpinChain* Chain = 0;
      if (PValue == 1)
	{
	  if (NbrSpins <= 32)
	    {
	      Chain = new PairHoppingP1AsSpin1Chain (NbrSpins, Manager.GetBoolean("use-periodic"), 1000000);
	    }
	  else
	    {
	      Chain = new PairHoppingP1AsSpin1ChainLong (NbrSpins, Manager.GetBoolean("use-periodic"), 1000000);
	    }
	}
      else
	{
	  if (PValue == 2)
	    {
	      if (NbrSpins <= 32)
		{
		  Chain = new PairHoppingP2AsSpin1Chain (NbrSpins, Manager.GetBoolean("use-periodic"), 1000000);
		}
	      else
		{
		  Chain = new PairHoppingP2AsSpin1ChainLong (NbrSpins, Manager.GetBoolean("use-periodic"), 1000000);
		}		
	    }	  
	}
      if (Chain->GetLargeHilbertSpaceDimension() > 0l)
       	{
	  cout << "Hilbert space dimension = " << Chain->GetLargeHilbertSpaceDimension() << endl;
	  // for (int i = 0; i < Chain->GetHilbertSpaceDimension(); ++i)
	  //   {
	  //     Chain->PrintState(cout, i) << endl;
	  //   }
	  PairHoppingHamiltonian Hamiltonian (Chain, NbrSpins, PValue, Manager.GetBoolean("use-periodic"));

	  char* TmpString = new char[16];
	  sprintf (TmpString, "");
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpEigenstateString, "%s", OutputFileName);

	  if (Manager.GetBoolean("export-charpolynomial"))
	    {
	      PairHoppingComputeCharacteristicPolynomial(&Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture());
	    }
	  
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete[] TmpString;
	}
      delete Chain;
    }
  return 0;
}

// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// discardFourFactor = discard a global four factor used to ensure integer numbers

void PairHoppingComputeCharacteristicPolynomial(PairHoppingHamiltonian* hamiltonian, AbstractSpinChain* chain, char* outputFileName, AbstractArchitecture* architecture)
{
#ifdef __GMP__
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
  TmpMatrix = LongIntegerMatrix(TmpRawMatrix);

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

