#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "Hamiltonian/SandGlassSpinChainHamiltonian.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"

using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
 
  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 2);
  SingleIntegerOption LzMaxOption ('l', "lzmax", "twice the z projection of the total momenutm", 1);
  SingleIntegerOption NbrSandGlassOption ('p', "nbr-sandglass", "number of sandglasses", 5);
  SingleDoubleOption JInternalOption ('j', "internal-j", "coupling constant between spin into a sand glass", 1.0);
  SingleDoubleOption JExternalOption ('J', "external-j", "coupling constant between spin between sand glasses", 1.0);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrSandGlassOption;
  OptionList += &LzMaxOption;
  OptionList += &JInternalOption;
  OptionList += &JExternalOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type ExplicitMatrixExample -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrSandGlass = NbrSandGlassOption.GetInteger();
  int LzMax = LzMaxOption.GetInteger();
  double JInternal = JInternalOption.GetDouble();
  double JExternal = JExternalOption.GetDouble();
  int NbrSpin = 5 * NbrSandGlass;

  if ((LzMax & 1) != (NbrSpin & 1))
    {
      if ((NbrSpin & 1) == 0)
	LzMax &= ~0x1;
      else
	LzMax |= 0x1;
    }

  Spin1_2Chain Chain(NbrSpin, LzMax, 1000000);

  cout << Chain.GetHilbertSpaceDimension() << endl;

  SandGlassSpinChainHamiltonian Hamiltonian ((AbstractSpinChain*) &Chain, NbrSandGlass, JInternal, JExternal);
  
  if (Chain.GetHilbertSpaceDimension() <= 512)
    {
      RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
      Hamiltonian.GetHamiltonian(HRep);
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	  HRep.Householder(TmpTriDiag, 1e-7);
	  TmpTriDiag.Diagonalize();
	  TmpTriDiag.SortMatrixUpOrder();
	  for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
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
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      BasicLanczosAlgorithm Lanczos(Architecture, 1);
      int MaxNbrIterLanczos = 200; 
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = 4;
      int NbrEigenvalue = 5;
      Lanczos.SetHamiltonian(&Hamiltonian);
      Lanczos.InitializeLanczosAlgorithm();
      cout << "Run Lanczos Algorithm" << endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
      Lanczos.RunLanczosAlgorithm(4);
      RealTriDiagonalSymmetricMatrix TmpMatrix;
      while ((Precision > 1e-13) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	{
	  Lanczos.RunLanczosAlgorithm(1);
	  TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
	  Lowest =TmpMatrix.DiagonalElement(NbrEigenvalue);// Lanczos.GetGroundStateEnergy();
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest;
	}
      cout << endl;
      cout << Lowest << " " << Precision << "  Nbr of iterations = " << CurrentNbrIterLanczos << endl;
      for (int i = 0; i <= NbrEigenvalue; ++i)
	{
	  cout << TmpMatrix.DiagonalElement(i) << " ";
	}
      cout << endl;
      cout << "------------------------------------------------------------------" << endl << endl;;
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
    }
  return 0;
}

