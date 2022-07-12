#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "Hamiltonian/NDMAPSpinChainHamiltonian.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
 
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption NbrSpinOption ('\0', "nbr-spin", "number of spins", 8);
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration (including resume run)", 3000);
  BooleanOption ReorthogonalizeOption ('\n', "reorthogonalize", "use reorthogonalized version of the lanczos algorithm", false);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 1);
  SingleDoubleOption AnisotropyConstantOption ('d', "anisotropy-constant", "value of the anisotropy constant", 0.0);
  SingleDoubleOption InPlaneAnisotropyConstantOption ('e', "inplane-anisotropy-constant", "value of the in-plane anisotropy constant", 0.0);
  SingleDoubleOption PerpendicularBFieldOption ('a', "perpendicular-bfield", "value of b field in the direction perpendicular to the chain", 0.0);
  SingleDoubleOption ParallelBFieldOption ('z', "parallel-bfield", "value of b field in the direction parallel to the chain", 0.0);
  SingleIntegerOption MemoryOption ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  BooleanOption FixedMomentumOption('\n', "fixed-momentum", "evaluate for one given momentum value");
  SingleIntegerOption FixedMomentumValueOption ('\n', "momentum", "momentum value if fixed-momentum is fixed to true", 0);
  BooleanOption DiskOption ('\n', "disk", "enable disk resume capabilities", false);
  BooleanOption ResumeOption ('\n', "resume", "resume from disk datas", false);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrIterationOption;
  OptionList += &NbrSpinOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &AnisotropyConstantOption;
  OptionList += &InPlaneAnisotropyConstantOption;
  OptionList += &PerpendicularBFieldOption;
  OptionList += &ParallelBFieldOption;
  OptionList += &MemoryOption;
  OptionList += &FixedMomentumOption;
  OptionList += &FixedMomentumValueOption;
  OptionList += &ReorthogonalizeOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type NDMAPSpinChain -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }
  int NbrSpin = NbrSpinOption.GetInteger();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  double AnisotropyConstant = AnisotropyConstantOption.GetDouble();
  double InPlaneAnisotropyConstant = InPlaneAnisotropyConstantOption.GetDouble();
  double PerpendicularBField = PerpendicularBFieldOption.GetDouble();
  double ParallelBField = ParallelBFieldOption.GetDouble();
  bool FixedMomentumFlag = FixedMomentumOption.GetBoolean();
  int FixedMomentum = FixedMomentumValueOption.GetInteger();
  bool ReorthogonalizeFlag = ReorthogonalizeOption.GetBoolean();
  int Memory = MemoryOption.GetInteger() << 20;
  int MinMomentum = 0;
  int MaxMomentum = NbrSpin - 1;
  if (FixedMomentumFlag == true)
    {
      MinMomentum = FixedMomentum;
      MaxMomentum = FixedMomentum;
    }  
  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  int VectorMemory = VectorMemoryOption.GetInteger();


  for (int k = MinMomentum; k <= MaxMomentum; ++k)
    {
      cout << "momentum = " << k << endl;
      Spin1ChainWithTranslations Space(NbrSpin, k, 10000000, 10000000);
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      NDMAPSpinChainHamiltonian Hamiltonian(&Space, NbrSpin, 1.0, 1.0, ParallelBField, PerpendicularBField, AnisotropyConstant, 
					    InPlaneAnisotropyConstant, Architecture, Memory);
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
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
	      for (int j = 0; j < HRep.GetNbrRow() ; j += 2)
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
	  if ((NbrEigenvalue == 1) || (ReorthogonalizeFlag == false))
	    {
	      Lanczos = new ComplexBasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);	
	    }
	  else
	    {
	      if (DiskFlag == false)
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage(Architecture, NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);
	    }
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 0;
	  Lanczos->SetHamiltonian(&Hamiltonian);
	  if ((DiskFlag == true) && (ResumeFlag == true))
	    Lanczos->ResumeLanczosAlgorithm();
	  else
	    Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  if (ResumeFlag == false)
	    {
	      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	      CurrentNbrIterLanczos = NbrEigenvalue + 3;
	      if ((DiskFlag == true) && (CurrentNbrIterLanczos >= NbrIterLanczos))
		{
		  NbrIterLanczos = CurrentNbrIterLanczos + 1;
		}
	    }
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Lanczos->TestConvergence() == false) &&  (((DiskFlag == true) && (CurrentNbrIterLanczos < NbrIterLanczos)) ||
							    ((DiskFlag == false) && (CurrentNbrIterLanczos < MaxNbrIterLanczos))))      
	    {
	      ++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      for (int i = 0; i < NbrEigenvalue; ++i)
		cout << TmpMatrix.DiagonalElement(i) << " ";
	      cout << "precision = " << Precision << " "<< endl;
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
	  if (Lanczos->TestConvergence() == true)
	    {
	      ComplexVector* Eigenstates =  (ComplexVector*) ((Lanczos)->GetEigenstates(NbrEigenvalue));
	      for (int i = 0; i < NbrEigenvalue; ++i)
		{
		  ComplexVector TmpVector (Space.GetHilbertSpaceDimension());
		  Hamiltonian.LowLevelMultiply(Eigenstates[i], TmpVector, 0, Space.GetHilbertSpaceDimension());
		  cout << (Eigenstates[i] * TmpVector).Re << " ";
		}
	    }
	  cout << endl;
	}
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
    }
  return 0;
}

