#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "HilbertSpace/Spin1ChainWithTranslations.h"

#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"

#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
 
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleIntegerOption NbrSpinOption ('\0', "nbr-spin", "number of spins", 8);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 1);
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &NbrSpinOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type SpinChainWithTranslation -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout, "SpinChainWithTranslation");
      return 0;
    }

  int NbrSpin = NbrSpinOption.GetInteger();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  double* CouplingConstants = new double[NbrSpin + 1];
  for (int i = 0; i <= NbrSpin; ++i)
    CouplingConstants[i] = 1.0;

  for (int k = 0; k < NbrSpin; ++k)
    {
      cout << "momentum = " << k << endl;
      Spin1ChainWithTranslations Space(NbrSpin, k, 0, 10000000, 10000000);
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      SpinChainHamiltonianWithTranslations Hamiltonian(&Space, NbrSpin, 1.0);
      if (Hamiltonian.GetHilbertSpaceDimension() < 200)
	{
	  HermitianMatrix HRep2 (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep2);
	  RealSymmetricMatrix HRep (HRep2.ConvertToSymmetricMatrix());
	  if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      for (int j = 0; j < (2 * NbrEigenvalue); j += 2)
		{
		  cout << TmpTriDiag.DiagonalElement(j) << " ";
		}
	      cout << endl;
	    }
	  else
	    {
	      cout << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  int MaxNbrIterLanczos = 4000;
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	      Lanczos = new ComplexBasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);	
	    }
	  else
	    {
	      Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 0;
	  Lanczos->SetHamiltonian(&Hamiltonian);
	  Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while (Lanczos->TestConvergence() == false)      
	    {
	      ++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	    }
	  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      exit(0);
	    }
	  cout << endl;
	  cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i < NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	    }
	  cout << endl;
	}
    }
  return 0;
}

