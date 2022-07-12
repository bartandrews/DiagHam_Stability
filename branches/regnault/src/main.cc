#include "Polynomial/Polynomial.h"
#include "Complex.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Matrix/ComplexTriDiagonalHermitianMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealLowerTriangularMatrix.h"
#include "Matrix/RealUpperTriangularMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Output/MathematicaOutput.h"
#include "DMRGAlgorithm/NonPeriodicDMRGAlgorithm.h"
#include "DMRGAlgorithm/NonPeriodicAsymmetricDMRGAlgorithm.h"
#include "DMRGAlgorithm/PeriodicDMRGAlgorithm.h"
#include "DMRGAlgorithm/DMRGBlock.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1AKLTChain.h"
#include "HilbertSpace/Spin1FullAKLTChain.h"
#include "HilbertSpace/ThierryChain.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Fermions.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "HilbertSpace/DMRGPartialHilbertSpace.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/TrappedBosons.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/Mn12Hamiltonian.h"
#include "Hamiltonian/V15Hamiltonian.h"
#include "Hamiltonian/AKLTHamiltonian.h"
#include "Hamiltonian/AnisotropicSpinChainHamiltonian.h"
#include "Hamiltonian/PeriodicAnisotropicSpinChainHamiltonian.h"
#include "Hamiltonian/SpinChainHamiltonian.h"
#include "Hamiltonian/SpinChainHamiltonianWithTranslations.h"
#include "Hamiltonian/DiamondSpinChainHamiltonian.h"  
#include "Hamiltonian/OpenDiamondSpinChainHamiltonian.h"
#include "Hamiltonian/TriangleSpinChainHamiltonian.h"
#include "Hamiltonian/DoubleTriangleSpinChainHamiltonian.h"
#include "Hamiltonian/XXZSpinChainHamiltonian.h"
#include "Hamiltonian/PeriodicXXZSpinChainHamiltonian.h"
#include "Hamiltonian/PeriodicSpinChainHamiltonian.h"
#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Hamiltonian/NonPeriodicDMRGHamiltonian.h"
#include "Hamiltonian/TrappedBosonHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereCoulombHamiltonian.h"
#include "TensorProduct/TensorProductRealVector.h"
#include "TensorProduct/TensorProductStructure.h"
#include "TensorProduct/CompositeTensorProductStructure.h"
#include "TensorProduct/FullTensorProductStructure.h"
#include "TensorProduct/TensorProductIndex.h"
#include "Tensor/OneSpaceTensor.h"
#include "Tensor/TwoSpaceTensor.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Interaction/BasicInteraction.h"
#include "Interaction/DiamondInteraction.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Operator/PeriodicAnisotropicMagnetizationOperator.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


void TensorTest();
void HilbertTest();
void BlockTest();
void DMRGSpaceTest();
void DMRGFullSpaceTest();
void AsymmetricDMRG();
RealMatrix RandomRealMatrix(int nbrRow, int nbrColumn, double range);
RealSymmetricMatrix RandomRealSymmetricMatrix(int nbrRow, double range);
RealTriDiagonalSymmetricMatrix RandomRealTriDiagonalSymmetricMatrix(int nbrRow, double range);
void GroundStateScaling (int nbrSpin, int maxNbrIterLanczos, double precision, int exactLimit);
void PeriodicDMRG();
void TestHamiltonian(AbstractHamiltonian* hamiltonian);

int main(int argc, char** argv)
{

  cout.precision(14);
  int dimM = 6;
  RealTriDiagonalSymmetricMatrix TmpM (RandomRealTriDiagonalSymmetricMatrix(dimM, 10));
  RealTriDiagonalSymmetricMatrix TmpM2;
  TmpM2.Copy(TmpM);
  TmpM2.Diagonalize(50);
  TmpM2.SortMatrixUpOrder();
  cout << TmpM2.DiagonalElement(0) << " " << TmpM2.DiagonalElement(1) << " " << TmpM2.DiagonalElement(2) << " " << TmpM2.DiagonalElement(3) << " " 
       << TmpM2.DiagonalElement(4) << " " << TmpM2.DiagonalElement(5) << endl << endl;    
  cout << TmpM << endl;
  RealMatrix Q (dimM, dimM, true);
  for (int i = 0; i < dimM; ++i)
    Q(i, i) = 1.0;
  int NbrShift = 2;
  RealMatrix Qp (dimM, dimM, true);
  for (int i = 0; i < dimM; ++i)
    for (int j = 0; j < dimM; ++j)
      Qp(i, j) = 1.0;
  double* Shifts = new double [NbrShift];
  for (int i = 0; i < NbrShift; ++i)
    Shifts[i] = TmpM2.DiagonalElement(dimM - 1 - i);
//  RealUpperTriangularMatrix R (TmpM.QRFactorization(Q));
//  cout << (R * Q) << endl;
  RealTriDiagonalSymmetricMatrix LQ (TmpM.PolynomialFilterWithExactShift(Q, Shifts, NbrShift));
  cout << Q << endl;
  cout << (Qp * Q) << endl;
  cout << Qp.Multiply(Q) << endl;
  cout << LQ << endl;
  RealMatrix Q2((Matrix&) Q);
  Q2.Transpose();
  RealMatrix Q3(LQ);  
  cout << (Q * Q3 * Q2) << endl;

/*  RealLowerTriangularMatrix L (TmpM.QLFactorization(Q));
  cout << TmpM << endl << endl << Q << endl << endl << L << endl << endl;
  RealMatrix Q2 = Q * L;
  cout << Q2 << endl;
  RealMatrix Q3 (L);
  cout << Q3 << endl;
  cout << Q * Q3 << endl;
  cout << Q3 * Q << endl;
  cout << "direct................." << endl;
  RealMatrix Q4 (dimM, dimM, true);
  for (int i = 0; i < dimM; ++i)
    Q4(i, i) = 1.0;
  RealTriDiagonalSymmetricMatrix LQ (TmpM.QLConjugaison(Q4, 1.0));
  RealMatrix Q5((Matrix&) Q4);
  Q5.Transpose();
  RealMatrix Q6(LQ);  
  cout << LQ << endl << Q4 << endl << Q6 << endl << (Q4 * Q6 * Q5) << endl;*/
  return 0;

//  GroundStateScaling (18, 200, 1e-13, 50);
//  return 0;
//  HilbertTest();
//  return 0;

//  HilbertTest();
//  return 0;

//  PeriodicDMRG();
//  AsymmetricDMRG();
//  return 0;
  //  BlockTest();
  //  FermionTest();
 // DMRGFullSpaceTest();
 // return 0;

// Spin 1 AKLT program
//  for (int NbrSpin = 4; NbrSpin < 16; NbrSpin += 2)
  int NbrBlocs = 5;
  //  int NbrSpin = NbrBlocs * 3 + 1;  
  int NbrSpin = 3;//NbrBlocs * 3+ 3;  
  {
    //  int NbrSpin = 10;
    //  Spin1AKLTChain Chain(NbrSpin - 1, 256);
    //  Spin1FullAKLTChain Chain(NbrSpin - 2, 0, 10000000);
    //      Spin1Chain Chain(NbrSpin, 256);
    //    Spin1Chain Chain(NbrSpin, 0, 4);    
    //    Spin1_2Chain Chain(6, 100000);
    
//    Spin1Chain Chain(NbrSpin, 0, 1000000);
    Spin1_2Chain Chain(NbrSpin, 1000000);
    cout << Chain.GetHilbertSpaceDimension() << endl;
    double x;
    //  for (int  i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
    //  cout << Chain.FindStateIndex(Chain.ChainDescription[i]) << " " << i << endl;
    //    Chain.PrintState(cout ,i) << (*(Chain.GetQuantumNumber(i))) << " " << i << " " << Chain.FindStateIndex(Chain.ChainDescription[i]) << endl;
    
  /* SubspaceSpaceConverter Converter;
     SzQuantumNumber Q (-4);
     AbstractHilbertSpace* SubChain = Chain.ExtractSubspace(Q, Converter);
     cout << SubChain->GetHilbertSpaceDimension() << endl;
     for (int  i = 0; i < SubChain->GetHilbertSpaceDimension(); i++)
     SubChain->PrintState(cout ,i) << endl;
     RealVector V1(Chain.GetHilbertSpaceDimension(), true);
     RealVector V2(SubChain->GetHilbertSpaceDimension(), true);
     for (int i = 0; i < SubChain->GetHilbertSpaceDimension(); i++)
     {
     if (i != 0)
     V2[i - 1] = 0.0;	  
     V2[i] = 1.0;
      cout << Converter.SubspaceToSpace(V2, V1) << endl;
      }*/
    double* CouplingConstants = new double [NbrSpin + 1];
    double* CouplingConstantsX = new double [NbrSpin + 1];
    for (int i = 0; i <= NbrSpin; i++)
      {
	CouplingConstants[i] = 1.0;
	CouplingConstantsX[i] = 1.0;
      }
    PeriodicSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
//    SpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  // PeriodicSpinChainHamiltonian H2 ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  //  PeriodicXXZSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstantsX, CouplingConstants);
  //  XXZSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstantsX, CouplingConstants);
  //  SpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  //  DiamondSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrBlocs, 1.0, 1.0, 1.0);
  //  OpenDiamondSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrBlocs, 1.0, 0.8, 0.3);
  //  TriangleSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrBlocs, 1.0, 0.5, 1.0);
  //  DoubleTriangleSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrBlocs, 2.0, 1.0);
  //  cout << H << endl;
    
    // Exact diagonalization
    RealSymmetricMatrix HRep (Chain.GetHilbertSpaceDimension());
    RealTriDiagonalSymmetricMatrix TmpTriDiag (Chain.GetHilbertSpaceDimension());
    H.GetHamiltonian(HRep);
    HRep.Householder(TmpTriDiag, 1e-14);
    TmpTriDiag.Diagonalize();
    TmpTriDiag.SortMatrixUpOrder();
    for (int j = 0; j < Chain.GetHilbertSpaceDimension(); j++)
      cout << TmpTriDiag.DiagonalElement(j) << " ";
    cout << endl;
    return 0;
      // Lanczos method
      //  BasicLanczosAlgorithm Lanczos;
//        cout << H << endl;
      //  ComplexBasicLanczosAlgorithm Lanczos;
//    AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
      AbstractArchitecture* Architecture = new SMPArchitecture(2);
//      ComplexBasicLanczosAlgorithmWithEigenstates Lanczos(Architecture);
    BasicLanczosAlgorithm Lanczos(Architecture, 1);
//      ComplexBasicLanczosAlgorithmWithGroundState Lanczos(Architecture);
      int MaxNbrIterLanczos = 200; 
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = 4;
      Lanczos.SetHamiltonian(&H);
      Lanczos.InitializeLanczosAlgorithm();
      cout << "Run Lanczos Algorithm" << endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
      Lanczos.RunLanczosAlgorithm(4);
      while ((Precision > 1e-13) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	{
	  Lanczos.RunLanczosAlgorithm(1);
	  Lowest = Lanczos.GetGroundStateEnergy();
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest;
	  cout << Lowest << " " << Precision << endl;
	}
      cout << endl;
      cout << "Nbr of spins = " << NbrSpin << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
      << CurrentNbrIterLanczos << endl;
      cout << "------------------------------------------------------------------" << endl << endl;;
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
  }
  return 0;
}

void GroundStateScaling (int nbrSpin, int maxNbrIterLanczos, double precision, int exactLimit)
{
  Spin1_2Chain Chain(nbrSpin, 1000000);
  double* CouplingConstants = new double [nbrSpin];
  for (int i = 0; i < nbrSpin; i++)
    CouplingConstants[i] = 1.0;
  List<AbstractQuantumNumber*> SzList (Chain.GetQuantumNumbers());
  AbstractQuantumNumber** TmpSz;
  ListIterator<AbstractQuantumNumber*>  IterSz (SzList);
      int MaxNbrIterLanczos = 200; 
//      AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
      AbstractArchitecture* Architecture = new SMPArchitecture(2);
  while ((TmpSz = IterSz()))
    {
      SubspaceSpaceConverter Converter;
      AbstractHilbertSpace* SubChain = Chain.ExtractSubspace(**TmpSz, Converter);      
      SpinChainHamiltonian H ((AbstractSpinChain*) SubChain, nbrSpin, CouplingConstants);      
      if (SubChain->GetHilbertSpaceDimension() >= exactLimit)
	{
	  BasicLanczosAlgorithm Lanczos (Architecture, 1);
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 4;
	  Lanczos.SetHamiltonian(&H);
	  Lanczos.InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  Lanczos.RunLanczosAlgorithm(4);
	  while ((Precision > precision) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	    {
	      Lanczos.RunLanczosAlgorithm(1);
	      Lowest = Lanczos.GetGroundStateEnergy();
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest;
	      cout << Lowest << " " << Precision << endl;
	    }
	  cout << endl;
	  cout << "Nbr of spins = " << nbrSpin << " " << (**TmpSz) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	}
      else
	{
	  RealSymmetricMatrix HRep (SubChain->GetHilbertSpaceDimension());
	  H.GetHamiltonian(HRep);
	  if (SubChain->GetHilbertSpaceDimension() == 1)
	    {
	      cout << "Nbr of spins = " << nbrSpin << " " << (**TmpSz) << " " << HRep(0, 0) << " " 
		   << precision << "  Nbr of iterations = 0" << endl;
	    }
	  else
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (SubChain->GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, precision);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      double Lowest = TmpTriDiag.DiagonalElement(0);
	      cout << "Nbr of spins = " << nbrSpin << " " << (**TmpSz) << " " << Lowest << " " << precision << "  Nbr of iterations = 0" 
		   << endl;
	    }
	  
	}
      delete SubChain;
    }
}

int main2 ()
{
  int NbrSpin = 12;
//  Spin1AKLTChain Chain(NbrSpin - 1, 256);
  Spin1FullAKLTChain Chain(NbrSpin - 2, 0, 10000000);
 // Spin1Chain Chain(NbrSpin, 256);
//  Spin1Chain Chain(3, 0, 256);
  cout << Chain.GetHilbertSpaceDimension() << endl;
  double x;
//  for (int  i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
//    Chain.PrintState(cout ,i) << (*(Chain.GetQuantumNumber(i))) << endl;
				    
  /* SubspaceSpaceConverter Converter;
  SzQuantumNumber Q (-4);
  AbstractHilbertSpace* SubChain = Chain.ExtractSubspace(Q, Converter);
  cout << SubChain->GetHilbertSpaceDimension() << endl;
  for (int  i = 0; i < SubChain->GetHilbertSpaceDimension(); i++)
    SubChain->PrintState(cout ,i) << endl;
  RealVector V1(Chain.GetHilbertSpaceDimension(), true);
  RealVector V2(SubChain->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < SubChain->GetHilbertSpaceDimension(); i++)
    {
      if (i != 0)
	V2[i - 1] = 0.0;	  
      V2[i] = 1.0;
      cout << Converter.SubspaceToSpace(V2, V1) << endl;
      }*/

  double* CouplingConstants = new double [NbrSpin];
  for (int i = 0; i < NbrSpin; i++)
    CouplingConstants[i] = 1.0;
  double* CouplingConstantsX = new double [NbrSpin];
  for (int i = 0; i < NbrSpin; i++)
    CouplingConstantsX[i] = 1.0;
  //CouplingConstantsX[1] = 0.0;
  //  CouplingConstants[1] = 0.0;
  // PeriodicSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  // PeriodicSpinChainHamiltonian H2 ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  //  XXZSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstantsX, CouplingConstants);

  SpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
//      AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
  AbstractArchitecture* Architecture = new SMPArchitecture(2);

  // Lanczos method
  BasicLanczosAlgorithm Lanczos(Architecture, 1);
  int MaxNbrIterLanczos = 200; 
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 4;
  Lanczos.SetHamiltonian(&H);
  Lanczos.InitializeLanczosAlgorithm();
  Lanczos.RunLanczosAlgorithm(4);
  while ((Precision > 1e-13) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
    {
      Lanczos.RunLanczosAlgorithm(1);
      Lowest = Lanczos.GetGroundStateEnergy();
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest;
    }
  cout << "Nbr of spins = " << NbrSpin << " " << Lowest << " " << Precision << "  Nbr of iterations = " << CurrentNbrIterLanczos << endl;
  return 0;

  // Householder method
  SpinChainHamiltonian H2 ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  RealSymmetricMatrix HRep (Chain.GetHilbertSpaceDimension());
  RealTriDiagonalSymmetricMatrix TmpTriDiag (Chain.GetHilbertSpaceDimension());
  H.GetHamiltonian(HRep);
  HRep.Householder(TmpTriDiag, 1e-7);
  TmpTriDiag.Diagonalize();
  TmpTriDiag.SortMatrixUpOrder();
  for (int j = 0; j < Chain.GetHilbertSpaceDimension(); j++)
    cout << TmpTriDiag.DiagonalElement(j) << " ";
  cout << endl;
  return 0;

  List<AbstractQuantumNumber*> ListQ = Chain.GetQuantumNumbers();
  ListIterator<AbstractQuantumNumber*> IterQ (ListQ);
  AbstractQuantumNumber** TmpQ;
  SubspaceSpaceConverter Converter2;

  /*  SzQuantumNumber Q (0);
  AbstractHilbertSpace* SubChain2 = Chain.ExtractSubspace(Q, Converter2);
  cout << "Hilbert space dimension = " << SubChain2->GetHilbertSpaceDimension() << endl;
  if (SubChain2->GetHilbertSpaceDimension() > 1)
    {
      for (int  i = 0; i < SubChain2->GetHilbertSpaceDimension(); i++)
	SubChain2->PrintState(cout ,i) << (*(SubChain2->GetQuantumNumber(i))) << endl;
      H.SetHilbertSpace(SubChain2);
      cout << H << endl;
      RealVector V1(SubChain2->GetHilbertSpaceDimension());
      for (int i = 0; i < SubChain2->GetHilbertSpaceDimension(); i++)
	V1[i] = (rand() - 32767) * 0.5;
      RealSymmetricMatrix HRep2 (SubChain2->GetHilbertSpaceDimension());
      RealTriDiagonalSymmetricMatrix TmpTriDiag2 (SubChain2->GetHilbertSpaceDimension());
      //H.Hamiltonian(HRep2);
      //      HRep2.Householder(TmpTriDiag2, 1e-7);
      H.Lanczos(SubChain2->GetHilbertSpaceDimension(), TmpTriDiag2, V1);
      TmpTriDiag2.Diagonalize();
      for (int j = 0; j < SubChain2->GetHilbertSpaceDimension(); j++)
	cout << TmpTriDiag2.DiagonalElement(j) << " ";
      cout << endl;
      }*/

  double* Energies = new double [Chain.GetHilbertSpaceDimension()];
  int* SubspaceIndices = new int [Chain.GetHilbertSpaceDimension()];
  int* EigenvectorIndices = new int [Chain.GetHilbertSpaceDimension()];
  List<RealMatrix> Eigenvectors;
  List<SubspaceSpaceConverter> Converters;
  int EnergyPosition = 0;
  int SubspaceIndex = 0;
  while ((TmpQ = IterQ()))
    {
      cout << (**TmpQ) << endl;
      SubspaceSpaceConverter Converter2;
      AbstractHilbertSpace* SubChain2 = Chain.ExtractSubspace(**TmpQ, Converter2);
      cout << "Hilbert space dimension = " << SubChain2->GetHilbertSpaceDimension() << endl;
      if (SubChain2->GetHilbertSpaceDimension() > 1)
	{
	  for (int  i = 0; i < SubChain2->GetHilbertSpaceDimension(); i++)
	    SubChain2->PrintState(cout ,i) << (*(SubChain2->GetQuantumNumber(i))) << endl;
	  H.SetHilbertSpace(SubChain2);
	  cout << H << endl;	  
	  RealSymmetricMatrix HRep2 (SubChain2->GetHilbertSpaceDimension());
	  RealMatrix TmpEigen (SubChain2->GetHilbertSpaceDimension(), 
			       SubChain2->GetHilbertSpaceDimension());
	  RealTriDiagonalSymmetricMatrix TmpTriDiag2 (SubChain2->GetHilbertSpaceDimension());
	  H.GetHamiltonian(HRep2);
	  HRep2.Householder(TmpTriDiag2, 1e-7, TmpEigen);
	  TmpTriDiag2.Diagonalize(TmpEigen);
	  for (int j = 0; j < SubChain2->GetHilbertSpaceDimension(); j++)
	    {
	      Energies[EnergyPosition] = TmpTriDiag2.DiagonalElement(j);
	      EigenvectorIndices[EnergyPosition] = j;
	      SubspaceIndices[EnergyPosition] = SubspaceIndex;
	      EnergyPosition++;
	    }	
	  Converters += Converter2;
	  Eigenvectors += TmpEigen;
	  /*	  RealVector TmpEigenVector(Chain.GetHilbertSpaceDimension()); 
	  RealVector TmpEigenVector2(Chain.GetHilbertSpaceDimension()); 
	  Converter2.SubspaceToSpace(TmpEigen[2], TmpEigenVector);
	  TmpEigenVector /= TmpEigenVector.Norm();
	  H2.Multiply (TmpEigenVector, TmpEigenVector2);
	  cout << (TmpEigenVector * TmpEigenVector2) << endl;
	  cout << endl*/;
	  SubspaceIndex++;
	}
      else
	{
	  RealVector* TmpEigen = new RealVector [1];
	  TmpEigen[0] = RealVector (1);
	  TmpEigen[0][0] = 1.0;
	  RealVector TmpEigenVector(Chain.GetHilbertSpaceDimension()); 
	  RealVector TmpEigenVector2(Chain.GetHilbertSpaceDimension()); 
	  Converter2.SubspaceToSpace((*TmpEigen), TmpEigenVector);
	  TmpEigenVector /= TmpEigenVector.Norm();
	  H2.Multiply (TmpEigenVector, TmpEigenVector2);
	  Energies[EnergyPosition] = (TmpEigenVector * TmpEigenVector2);
	  EigenvectorIndices[EnergyPosition] = 0;
	  SubspaceIndices[EnergyPosition] = SubspaceIndex;
	  EnergyPosition++;
	  SubspaceIndex++;
	  Converters += Converter2;
	  Eigenvectors += RealMatrix(TmpEigen, 1);
	}
      delete SubChain2;
    }
  cout << "Energies = " << endl;
  for (int i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
    cout << Energies[i] << " ";
  cout << endl;
  /*  int ReducedDim = Chain.GetHilbertSpaceDimension() - 1;
  double TmpDouble;
  int TmpInt;
  for (int i = 0; i < ReducedDim; i++)
    for (int j = 0; j < (ReducedDim - i); j++)
      if (Energies[j] > Energies[j + 1])
	{
	  TmpDouble = Energies[j];
	  Energies[j] = Energies[j + 1];
	  Energies[j + 1] = TmpDouble;
	  TmpInt = SubspaceIndices[j];
	  SubspaceIndices[j] = SubspaceIndices[j + 1];
	  SubspaceIndices[j + 1] = TmpInt;
	  TmpInt = EigenvectorIndices[j];
	  EigenvectorIndices[j] = EigenvectorIndices[j + 1];
	  EigenvectorIndices[j + 1] = TmpInt;
	}
  cout << "Sorted energies = " << endl;
  for (int i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
    {
      cout << Energies[i] << " " << SubspaceIndices[i] << " " << EigenvectorIndices[i] << endl;
      RealVector TmpEigenVector(Chain.GetHilbertSpaceDimension()); 
      RealVector TmpEigenVector2(Chain.GetHilbertSpaceDimension()); 
      Converters[SubspaceIndices[i]].SubspaceToSpace(Eigenvectors[SubspaceIndices[i]][EigenvectorIndices[i]], TmpEigenVector);
      TmpEigenVector /= TmpEigenVector.Norm();
      H2.Multiply (TmpEigenVector, TmpEigenVector2);
      cout << (TmpEigenVector * TmpEigenVector2) << endl;
      }*/
  int TruncatedSpaceDimension = 9;//2 * NbrSpin + 1;
  RealDiagonalMatrix EnergyMatrix (Energies, TruncatedSpaceDimension);
  RealVector* TruncatedSpaceBase = new RealVector [TruncatedSpaceDimension]; 
  for (int i = 0; i < TruncatedSpaceDimension; i++)
    {
      RealVector TmpEigenVector(Chain.GetHilbertSpaceDimension());
      Converters[SubspaceIndices[i]].SubspaceToSpace(Eigenvectors[SubspaceIndices[i]][EigenvectorIndices[i]], TmpEigenVector);
      TmpEigenVector /= TmpEigenVector.Norm();
      TruncatedSpaceBase[i] = TmpEigenVector;
    } 
  RealMatrix TransformationMatrix (TruncatedSpaceBase, TruncatedSpaceDimension);
  List<Matrix*> InteractionOperators (H2.RightInteractionOperators());
  List<Matrix*> ConjugatedInteractionOperators;
  ListIterator<Matrix*> IterOperators(InteractionOperators);
  Matrix** TmpMatrix;
  while ((TmpMatrix = IterOperators()))
    {
      ConjugatedInteractionOperators += (*TmpMatrix)->Conjugate(TransformationMatrix);
    }
  /*    {
      TmpMatrix = IterOperators();
      cout << (*((RealSymmetricMatrix*)*TmpMatrix)) << endl;
      ConjugatedInteractionOperators += (*TmpMatrix)->Conjugate(TransformationMatrix);
      TmpMatrix = IterOperators();
      cout << (*((RealAntisymmetricMatrix*)*TmpMatrix)) << endl;
      ConjugatedInteractionOperators += (*TmpMatrix)->Conjugate(TransformationMatrix);
      TmpMatrix = IterOperators();
      cout << (*((RealSymmetricMatrix*)*TmpMatrix)) << endl;
      ConjugatedInteractionOperators += (*TmpMatrix)->Conjugate(TransformationMatrix);
      }*/
  cout << TransformationMatrix << endl;
  if (ConjugatedInteractionOperators.GetNbrElement() > 0)
    {
      IterOperators.DefineList(ConjugatedInteractionOperators);
      TmpMatrix = IterOperators();
      TensorProductStructure Structure (2);
      Structure.SetDimension(0, (*TmpMatrix)->GetNbrRow());
      Structure.SetDimension(1, (*TmpMatrix)->GetNbrRow());
      double* TensorCouplingConstants = new double [3];
      TensorCouplingConstants[0] = 1.0;
      TensorCouplingConstants[1] = -1.0;
      TensorCouplingConstants[2] = 1.0;
      BasicInteraction Interaction(TensorCouplingConstants, 3, 0, 1, &Structure);
      TwoSpaceTensor tensor = Interaction.Interaction(ConjugatedInteractionOperators, 
						      ConjugatedInteractionOperators);
      /*      TwoSpaceTensor tensor (Structure, *TmpMatrix, *TmpMatrix, 0);
      cout << (*((RealSymmetricMatrix*)*TmpMatrix)) << endl;
//      cout << (*((RealSymmetricMatrix*) tensor.ElementaryMatrix)) << endl;
      TmpMatrix = IterOperators();
      tensor.AddTensorProductMatrices(*TmpMatrix, *TmpMatrix, -1);
      cout << (*((RealAntisymmetricMatrix*)*TmpMatrix)) << endl;
//      cout << (*((RealSymmetricMatrix*) tensor.ElementaryMatrix)) << endl;
      TmpMatrix = IterOperators();
      tensor.AddTensorProductMatrices(*TmpMatrix, *TmpMatrix);
      cout << (*((RealSymmetricMatrix*)*TmpMatrix)) << endl;
      //      cout << (*((RealSymmetricMatrix*) tensor.ElementaryMatrix)) << endl;*/
      tensor.AddToFirstSpace(EnergyMatrix);
      tensor.AddToSecondSpace(EnergyMatrix);
      RealSymmetricMatrix TotalH (*((RealSymmetricMatrix*) (tensor.ElementaryMatrix)));
      RealTriDiagonalSymmetricMatrix TmpTotalTriDiag (TotalH.GetNbrRow());
      TotalH.Householder(TmpTotalTriDiag, 1e-7);
      TmpTotalTriDiag.Diagonalize();
      TmpTotalTriDiag.SortMatrixUpOrder();
      for (int j = 0; j < TmpTotalTriDiag.GetNbrRow(); j++)
	cout << TmpTotalTriDiag.DiagonalElement(j) << " ";
      cout << endl;
    }

  List<AbstractQuantumNumber*> TmpListQ1 (ListQ);
  List<AbstractQuantumNumber*> TmpListQ2 (ListQ);
  ListIterator<AbstractQuantumNumber*> IterListQ1(TmpListQ1);
  AbstractQuantumNumber** SumQ = new AbstractQuantumNumber* [TmpListQ1.GetNbrElement() * 
							    TmpListQ2.GetNbrElement()];
  AbstractQuantumNumber** TmpQ1;
  AbstractQuantumNumber** TmpQ2;
  int Pos = 0;
  while ((TmpQ1 = IterListQ1()))
    {
      ListIterator<AbstractQuantumNumber*> IterListQ2(TmpListQ2);
      while ((TmpQ2 = IterListQ2()))
	{
	  SumQ[Pos++] = (**TmpQ1) + (**TmpQ2);
	}      
    }
  int Lim = (TmpListQ1.GetNbrElement() * TmpListQ2.GetNbrElement());
  int NbrSubspace = 0;
  int* NbrDirectSumSubspace = new int [TmpListQ1.GetNbrElement() * 
				      TmpListQ2.GetNbrElement()];
  int** DirectSumSubspaces = new int* [TmpListQ1.GetNbrElement() * 
				      TmpListQ2.GetNbrElement()];
  int NbrSumDeleted = 0;
  for (int i = 0; i < Lim; i++)
    {
      int NbrDirectSum = 0;
      DirectSumSubspaces[NbrSubspace] = new int [2 * (TmpListQ1.GetNbrElement() * 
						      TmpListQ2.GetNbrElement() - NbrSumDeleted)];
      NbrDirectSumSubspace[NbrSubspace] = 0;
      if (SumQ[i] != 0)
	{
	  for (int j = i + 1; j < Lim; j++)
	    if ((SumQ[j] != 0) && ((*(SumQ[j])) == (*(SumQ[i]))))
	      {
		delete SumQ[j];
		SumQ[j] = 0;
		DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum] = j / TmpListQ1.GetNbrElement();
		DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1] = 
		  j - TmpListQ1.GetNbrElement() * DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum];
		NbrDirectSumSubspace[NbrSubspace]++;
		NbrDirectSum++;
	      }
	  DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum] = i / TmpListQ1.GetNbrElement();
	  DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1] = 
	    i - TmpListQ1.GetNbrElement() * DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum];
	  NbrDirectSumSubspace[NbrSubspace++]++;
	  NbrDirectSum++;
	  NbrSumDeleted += NbrDirectSum;
	  delete SumQ[i];
	}
    }
  for (int i = 0; i < NbrSubspace; i++)
    {
      cout << i << " " << NbrDirectSumSubspace[i] << "-----------------------" << endl;
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  cout << DirectSumSubspaces[i][2 * j] << " " << DirectSumSubspaces[i][2 * j + 1] << endl;
	}
    }
  /*  double* CouplingConstants = new double [NbrSpin - 1];
  for (int i = 0; i < (NbrSpin - 1); i++)
    CouplingConstants[i] = 1.0;
  SpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstants);
  cout << H << endl;
  RealSymmetricMatrix Sx (Chain.GetHilbertSpaceDimension(), true);
  RealAntisymmetricMatrix Sy (Chain.GetHilbertSpaceDimension(), true);
  RealSymmetricMatrix Sz (Chain.GetHilbertSpaceDimension(), true);
  int Pos = 1;
  Chain.Sxi(Pos, Sx);
  Chain.Syi(Pos, Sy);
  Chain.Szi(Pos, Sz);
  cout << "Sx_" << Pos << " :" << endl;
  cout << Sx << endl;
  cout << "Sy_" << Pos << " :" << endl;
  cout << Sy << endl;
  cout << "Sz_" << Pos << " :" << endl;
  cout << Sz << endl;*/
/*  AKLTHamiltonian H (Chain, 2, 2.0, 3.0);
//  cout << H << endl;
  List<Matrix*> ListMatrix = H.LeftInteractionOperators();
  ListIterator<Matrix*> IterMatrix(ListMatrix);
  Matrix** TmpMat;
  int i = 0;
  while ((TmpMat = IterMatrix()))
    {
      if (i != 1)
	cout << (*((RealSymmetricMatrix*) (*TmpMat))) << endl;
      else
	cout << (*((RealAntisymmetricMatrix*) (*TmpMat))) << endl;
      i++;
    }*/

// tensor product test
/*  int dim = 18;
  int NbrSpace = 3;
  RealVector V(dim);
  for (int i = 0; i < dim; i++)
    V[i] = 1.0  + (double) i;
  TensorProductStructure Structure(NbrSpace);
  Structure.SetDimension(0, 3);
  Structure.SetDimension(1, 3);
  Structure.SetDimension(2, 2);
  TensorProductIndex Index(Structure);
  TensorProductRealVector TV (Structure, V);
  TensorProductRealVector TV2 (Structure);
  RealSymmetricMatrix M (3);
  M(0,0) = 1.0;
  M(0,1) = 1.5;
  M(0,2) = -0.7;
  M(1,1) = 1.0;
  M(1,2) = 1.5;
  M(2,2) = 1.0;
  RealSymmetricMatrix M2 (3);
  M2(0,0) = 1.0;
  M2(0,1) = -1.0;
  M2(1,1) = 2.0;
  OneSpaceTensor T(Structure, M, 1);
  OneSpaceTensor T2(Structure, M2, 2);
  TwoSpaceTensor T3(T, T2);
  cout << (*((RealSymmetricMatrix*)(T3. ElementaryMatrix))) << endl;
  TV2.Multiply(T3, TV);
  for (int i = 0; i < 3; i++)
    {
      Index[0] = i;
      for (int j = 0; j < 3; j++)
	{
	  Index[1] = j;
	  for (int k = 0; k < 2; k++)
	    {
	      Index[2] = k;
	      cout << "(" << i << "," << j << "," << k << ")   =   " << TV2[Index] << " " << TV[Index] << endl;
	    }
	}
    }
*/
// Spin 1 AKLT program
/*
  Spin1AKLTChain Chain(2, 256);
  cout << Chain.GetHilbertSpaceDimension() << endl;
  double x;
  for (int  i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
    Chain.PrintState(cout ,i) << endl;
  AKLTHamiltonian H (Chain, 2, 2.0, 3.0);
//  cout << H << endl;
  List<Matrix*> ListMatrix = H.LeftInteractionOperators();
  ListIterator<Matrix*> IterMatrix(ListMatrix);
  Matrix** TmpMat;
  int i = 0;
  while ((TmpMat = IterMatrix()))
    {
      if (i != 1)
	cout << (*((RealSymmetricMatrix*) (*TmpMat))) << endl;
      else
	cout << (*((RealAntisymmetricMatrix*) (*TmpMat))) << endl;
      i++;
    }*/
/*
  RealSymmetricMatrix HRep (H.GetHilbertSpaceDimension());
  H.Hamiltonian(HRep);
  RealMatrix Q (H.GetHilbertSpaceDimension(), H.GetHilbertSpaceDimension(), true);
  Q(1, 0) = 1.0;
  Q(0, 1) = 1.0;
  for (int i = 2; i < H.GetHilbertSpaceDimension(); i++)
    Q(i, i) = 1.0;
  RealTriDiagonalSymmetricMatrix D (H.GetHilbertSpaceDimension());
  HRep.Householder(D, 1e-7, Q);
//  cout << Q << endl;
  D.Diagonalize(Q);
  cout << D << endl;
  RealVector V1(H.GetHilbertSpaceDimension(), true);
  RealVector V2(H.GetHilbertSpaceDimension(), true);
  RealSymmetricMatrix HRep2 (H.GetHilbertSpaceDimension());
  H.Hamiltonian(HRep2);
//  RealSymmetricMatrix HRep3 = HRep2.Conjugate(Q);
  for (int j = 0; j < H.GetHilbertSpaceDimension(); j++)
    {
      RealVector V3(H.GetHilbertSpaceDimension(), true);
      if (j != 0)
	V1[j - 1] = 0.0;
      V1[j] = 1.0;
      V3.Multiply(Q, V1);
      V2.Multiply(HRep2, V3);
      for (int i = 0; i < H.GetHilbertSpaceDimension(); i++)
	{
	  if (V3[i] != 0.0)
	    cout << V2[i] / V3[i] << " " << V2[i] << " " << V3[i] << endl;
	}
      cout << endl;
    }
*/

/*
// V15 program

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  timeval StartingTime;
  timeval EndingTime;
  double Dt;
  gettimeofday (&(TotalStartingTime), 0);

  Spin1_2Chain Chain(15, 3, 256);
  cout << Chain.GetHilbertSpaceDimension() << endl;
//  for (int  i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
//    Chain.PrintState(cout ,i) << endl;

  cout << "pre-calculation : " << endl;
  gettimeofday (&(StartingTime), 0);
  V15Hamiltonian H (Chain, 290.5, 15.9, - 22.8, 23.4, 13.8);
//  V15Hamiltonian H (Chain, 1.0, 15.9, - 22.8, 23.4, 13.8);
  double* DiagH = new double [H.GetHilbertSpaceDimension()];
  double* OffDiagH = new double [(H.GetHilbertSpaceDimension() * (H.GetHilbertSpaceDimension() - 1)) / 2];
  RealSymmetricMatrix HRep (DiagH, OffDiagH, H.GetHilbertSpaceDimension());
  H.Hamiltonian(HRep);
  MathematicaOutput Out;
//  Out << HRep;
//  cout << Out << endl;
  gettimeofday (&(EndingTime), 0);
  Dt = (double) (EndingTime.tv_sec - StartingTime.tv_sec) + ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0);
  cout << "time = " << Dt << endl;

  double* Diag = new double [Chain.GetHilbertSpaceDimension()];
  double* UpperDiag = new double [Chain.GetHilbertSpaceDimension()];
  RealTriDiagonalSymmetricMatrix D (Diag, UpperDiag, Chain.GetHilbertSpaceDimension());
  RealVector V1(Chain.GetHilbertSpaceDimension());
  for (int i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
    V1[i] = (rand() - 32767) * 0.5;
//  int NbrIter = 50;
//  int NbrIter = 1;
  int NbrIter = 0;
//  int NbrLanczosStep = 10;
  int NbrLanczosStep = Chain.GetHilbertSpaceDimension();
  for (int i = 0; i < NbrIter; i++)
    {
      cout <<  endl << " begin Lanczos step..." << endl;
      gettimeofday (&(StartingTime), 0);
      H.FullReorthogonalizedLanczos (NbrLanczosStep, D, V1);//, 10);
//      H.Lanczos (NbrLanczosStep, D, V1);
      cout << " end Lanczos step : " << endl;
      gettimeofday (&(EndingTime), 0);
      Dt = (double) (EndingTime.tv_sec - StartingTime.tv_sec) + ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
      
      cout << endl << " begin QL step..." << endl;
      gettimeofday (&(StartingTime), 0);
      
      D.Diagonalize();
      cout << " end QL step" << endl;
      gettimeofday (&(EndingTime), 0);
      Dt = (double) (EndingTime.tv_sec - StartingTime.tv_sec) + ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;

      cout << endl << " eigenvalues : " << endl;
      D.SortMatrix();
      for (int i = 0; i < NbrLanczosStep; i++)
	cout << D.DiagonalElement(i) << "  ";
      cout << endl;
      NbrLanczosStep++;
    }
  cout << endl;
  NbrIter = 0;
  NbrLanczosStep = Chain.GetHilbertSpaceDimension();
  if (NbrIter == 0)
    {
      HRep.Householder (D, 1e-7);
      D.Diagonalize();
      D.SortMatrix();
      for (int i = 0; i < NbrLanczosStep; i++)
	cout << D.DiagonalElement(i) << "  ";
      cout << endl;
    }
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + ((TotalEndingTime.tv_usec - 
									TotalStartingTime.tv_usec) / 1000000.0);
  cout << endl << endl << "total time = " << Dt << endl;

*/


/*  cout << D.DiagonalElement(NbrIter - 1) << endl << endl;
  D.Diagonalize();
  for (int i = 0; i < NbrIter; i++)
    {
      cout << D.DiagonalElement(i) << "  ";
    }
  cout << endl;
  cout << H.Lanczos(100, V1) << endl << endl;
*/

//  for (int i = 0; i < 20; i++)
//    cout << H.Lanczos(NbrIter, V1) << endl << endl;

//  H.Lanczos(NbrIter, D, V1);
//  cout << D << endl;
//  cout << H.Lanczos(NbrIter, V1) << endl;
//  D.Diagonalize();
//  cout << D << endl;

/*  int ChainLength = 2;
  Spin1AKLTChain Chain(ChainLength);
  for (int  i = 0; i < Chain.GetHilbertSpaceDimension(); i++)
    Chain.PrintState(cout ,i) << endl;
*/
/*  int dim = 27;
  double* Diag = new double [dim];
  for (int i = 0; i < dim; i++)
    Diag[i] = 0.0;
  double* OffDiag = new double [dim * (dim - 1)];
  for (int i = 0; i < dim * (dim - 1); i++)
    OffDiag[i] = 0.0;
  HermitianMatrix H(Diag, OffDiag, dim);
  Chain.AddSziSzj(H, 0, 1);
  Chain.AddSmiSpj(H, 0, 1, 0.5);
  Chain.AddSmiSpj(H, 1, 0, 0.5);
  Chain.AddSziSzj(H, 1, 1);
  Chain.AddSmiSpj(H, 1, 2, 0.5);
  Chain.AddSmiSpj(H, 2, 1, 0.5);
  MathematicaOutput Out;
  cout << H << endl;
  ComplexVector* Vecs = new ComplexVector [dim];
  for (int i = 0; i < dim; i++)
    {
      double* RComp1 = new double [dim * 2];
      for (int j = 0; j < dim * 2; j+=2)
	{
	  RComp1[j] = rand();
	  RComp1[j + 1] = 0.0;
	}
//      RComp1[2* i] = 1.0;
      Vecs[i] = ComplexVector (RComp1, dim);
    }
  Out << H;
  Out <<  Vecs[0];
  cout << Out << endl;
  ComplexMatrix Q (Vecs, dim);
  double* Diag2 = new double [dim];
  double* UpDiag2 = new double [2 * (dim - 1)];
  RealTriDiagonalSymmetricMatrix D (Diag2, UpDiag2, dim);
  H.OrthoLanczos (dim, D, Q);
  cout << Q << endl;
  cout << D << endl;
//  Polynomial P = D.CharacteristicEquation();
  D.Diagonalize(Q);
  cout << D << endl;*/
//  P.LaguerreMethodSolver(0.000000001, 0.000000001, 1000);


/*  int Ne = 2;
  int dim = 2 * (Ne + 1);
  double Step1 = 1.0;//(double) (2 * (Ne + 1));
  Step1 *= Step1;
  double* Diag = new double [dim];
  for (int i = 0; i < dim; i++)
    Diag[i] = 2.0 * Step1;
  double* OffDiag = new double [dim * (dim - 1)];
  for (int i = 0; i < dim * (dim - 1); i++)
    OffDiag[i] = 0.0;
  int k = 0;
  for (int i = 0; i < (dim - 1); i++)
    {
      OffDiag[k] = -1.0 * Step1;
      k += 2 * (dim - 1 - i);
    }
  HermitianMatrix H(Diag, OffDiag, dim);
  cout << H << endl;
  
  ComplexVector* Vecs = new ComplexVector [dim];
  for (int i = 0; i < dim; i++)
    {
      double* RComp1 = new double [dim * 2];
      for (int j = 0; j < dim * 2; j++)
	RComp1[j] = 0.0;
      RComp1[2* i] = 1.0;
      Vecs[i] = ComplexVector (RComp1, dim);
    }

  ComplexMatrix Q (Vecs, dim);
  double* Diag2 = new double [dim];
  double* UpDiag2 = new double [2 * (dim - 1)];
  RealTriDiagonalSymmetricMatrix D (Diag2, UpDiag2, dim);
  H.Lanczos (dim, D, Q);
  MathematicaOutput Out;
  D.Diagonalize(Q);
  D.SortMatrix(Q);
  Q.Resize(Ne + 1, Ne);
  H.Resize(Ne + 1, Ne + 1);
  Q.NormalizeColumns();
  HermitianMatrix QHQ = H.Conjugate(Q);
//  HermitianMatrix H2 (Diag3, OffDiag3, dim);
  H.Resize(dim, dim);
  int N = 20;
  int MaxStepNbr = 2 * (N + Ne + 1);
  SingleParticle DMRG(Ne, MaxStepNbr);
  DMRG.StoreHamiltonian(QHQ, Ne + 1);
  DMRG.StoreInteraction(Q, Ne + 1);      
  double Scale = (double) MaxStepNbr;
  Scale *= Scale;
  double test;
  for (int i = 1; i <= N; i++)
    {
      cout << "iteration nbr = " << i << " ---------------------" << endl;
      H.Resize(dim, dim);
      //      cout << H << endl;
      DMRG.InfiniteDMRGNewHamiltonian(H, QHQ, Q, 2.0 * Step1, 2.0 * Step1, Complex(-1.0 * Step1, 0.0));
      double Step2 = 1.0;//(double) (2 * (Ne + 1 +i));
      Step2 *= Step2;
//      cout << H << endl;
      H *= (Step2 / Step1);
//      cout << H << endl;
//      Out << H;
//      cout << Out << endl;
      Step1 = Step2;
      //     QHQ *= Scale / Step1;
      Q.Resize(dim, dim);
      //      cout << Q << endl;
      H.Lanczos (dim, D, Q);
      //      cout << Out << endl;
      //     cout << Q << endl;
      D.Diagonalize(Q);
      //      cout << D << endl;
      D.SortMatrix(Q);
      test = sin ((0.5 * M_PI) / (2.0 * (Ne + 1 +i) + 1.0));
      test *= 4.0 *test;
      cout << Diag2[0] << " " << test << endl;
      //      cout << D << endl;
      Q.Resize(Ne + 1, Ne);
      H.Resize(Ne + 1, Ne + 1);
      Q.NormalizeColumns();
      QHQ = H.Conjugate(Q);
      DMRG.StoreHamiltonian(QHQ, Ne + i + 1);
      DMRG.StoreInteraction(Q, Ne + i + 1);      
//      cout << H << endl;
//      cout << (M_PI * M_PI) << endl;
    }
  cout << "------------- finite size --------------" << endl;
  int i = N + Ne + 2;
  for (int j = 0; j < 2; j++)
    {
      cout << endl << "------ left to right-----------" << endl;
      for (; i < (MaxStepNbr - Ne - 2); i++)
	{
	  H.Resize(dim, dim);
	  DMRG.FiniteDMRGNewHamiltonian(H, i - 1, MaxStepNbr - i - 1, 2.0 * Scale, 2.0 * Scale, Complex(-1.0 * Scale, 0.0));      
	  Q.Resize(dim, dim);
	  H.Lanczos (dim, D, Q);
	  D.Diagonalize(Q);
	  D.SortMatrix(Q);
	  cout << Diag2[0] << "   " << test << endl;
	  Q.Resize(Ne + 1, Ne);
	  H.Resize(Ne + 1, Ne + 1);
	  Q.NormalizeColumns();
	  QHQ = H.Conjugate(Q);
	  DMRG.StoreHamiltonian(QHQ, i);
	  DMRG.StoreInteraction(Q, i);            
	}
      i--;
      cout << endl << "------ right to left -----------" << endl;
      for (; i > (Ne + 1); i--)
	{
	  H.Resize(dim, dim);
	  DMRG.FiniteDMRGNewHamiltonian(H, i - 1, MaxStepNbr - i - 1, 2.0 * Scale, 2.0 * Scale, Complex(-1.0 * Scale, 0.0));      
	  Q.Resize(dim, dim);
	  H.Lanczos (dim, D, Q);
	  D.Diagonalize(Q);
	  D.SortMatrix(Q);
	  cout << Diag2[0] << "   " << test << endl;
	  Q.Resize(Ne + 1, Ne);
	  H.Resize(Ne + 1, Ne + 1);
	  Q.NormalizeColumns();
	  QHQ = H.Conjugate(Q);
	  DMRG.StoreHamiltonian(QHQ, i);
	  DMRG.StoreInteraction(Q, i);            
	}
      i++;
    }
    cout << endl << "------ end -----------" << endl;
}*/


// RealSymmetricMatrix test
/*
  
  int dim = 4;
  double* Diag = new double [dim];
  double* Diag2 = new double [dim];
  double* UpDiag = new double [dim];
  double* OffDiag = new double [(dim * (dim - 1)) / 2];  
  RealSymmetricMatrix H (Diag, OffDiag, dim);
  RealTriDiagonalSymmetricMatrix D (Diag2, UpDiag, dim);
  RealTriDiagonalSymmetricMatrix D2 (dim);
  H.SetMatrixElement (0, 0, 1.0);
  H.SetMatrixElement (1, 1, 2.0);
  H.SetMatrixElement (2, 2, 3.0);
  H.SetMatrixElement (3, 3, 4.0);
  H.SetMatrixElement (0, 1, 5.0);
  H.SetMatrixElement (2, 0, 6.0);
  H.SetMatrixElement (0, 3, 7.0);
  H.SetMatrixElement (1, 3, 9.0);
  H.SetMatrixElement (2, 1, 8.0);
  H.SetMatrixElement (2, 3, 10.0);
  RealSymmetricMatrix H2 (dim);
  H2.SetMatrixElement (0, 0, 1.0);
  H2.SetMatrixElement (1, 1, 2.0);
  H2.SetMatrixElement (2, 2, 3.0);
  H2.SetMatrixElement (3, 3, 4.0);
  H2.SetMatrixElement (0, 1, 5.0);
  H2.SetMatrixElement (2, 0, 6.0);
  H2.SetMatrixElement (0, 3, 7.0);
  H2.SetMatrixElement (1, 3, 9.0);
  H2.SetMatrixElement (2, 1, 8.0);
  H2.SetMatrixElement (2, 3, 10.0);
//  H.Resize (3, 3);
  cout << H << endl;
  RealMatrix Q (4, 4, true);
  for (int i = 0; i < dim ; i++)
    Q(i, i) = 1.0;
  H.Householder (D, 1e-7, Q);
//  H2.Householder (D2, 1e-7, Q);
  cout << H << endl;
//  RealSymmetricMatrix H3 = H2.Conjugate(Q);
  D.Diagonalize(Q);
  cout << D << endl;
  cout << D2 << endl;
  cout << Q << endl;  
  RealVector V1(4, true);
  RealVector V2(4, true);
  RealVector V3(4, true);
  V1[1] = 1.0;
//  Q.Transpose();
  V3.Multiply(Q, V1);
  V2.Multiply(H2, V3);
  cout << V3 << endl;
//  V2.Multiply(Q, V1);
//  V3 *= D2;
  cout << V2 << endl;
  for (int i = 0; i < 4; i++)
    {
      if (V3[i] != 0.0)
	cout << V2[i] / V3[i] << endl;
    }
 //  cout << V1 << endl;
  //V2.Multiply(H, V1);
  //cout << V2 << endl;
*/ 
 
// end program

}

////////////////////////////////////////////////////////////////////////////

void TensorTest()
{
  int dim = 3;
  RealSymmetricMatrix H1 (dim);
  H1.SetMatrixElement (0, 0, 1.0);
  H1.SetMatrixElement (1, 1,-2.0);
  H1.SetMatrixElement (2, 2, 3.0);
  H1.SetMatrixElement (0, 1, -5.0);
  H1.SetMatrixElement (2, 0, 1.0);
  H1.SetMatrixElement (2, 1, 4.0);
  RealSymmetricMatrix H2 (dim);
  H2.SetMatrixElement (0, 0, 2.0);
  H2.SetMatrixElement (1, 1, 1.0);
  H2.SetMatrixElement (2, 2, 3.0);
  H2.SetMatrixElement (0, 1, 2.0);
  H2.SetMatrixElement (2, 0, -4.0);
  H2.SetMatrixElement (2, 1, 3.0);
  RealSymmetricMatrix H3 (dim, true);
  RealSymmetricMatrix H4 (dim, true);
  TensorProductStructure* Structure = new TensorProductStructure(2);
  Structure->SetDimension(0, dim);
  Structure->SetDimension(1, dim);
  TwoSpaceTensor tensor (Structure, (Matrix*) &H1, (Matrix*) &H2, 0);
  cout << H1 << endl;
  cout << H2 << endl;
  tensor.SubTensorProductMatrices((Matrix*) &H3, (Matrix*) &H4);
  cout << (*((RealSymmetricMatrix*) tensor.ElementaryMatrix)) << endl;
}

////////////////////////////////////////////////////////////////////////////

void HilbertTest()
{
//  Fermions Space(3, 256);
  int NbrSite = 4;
  int M = 7;
//  Spin1_2ChainWithTranslations Space(NbrSite, 256);
//  Spin1_2Chain Space(NbrSite, 256);
  TrappedBosons Space (NbrSite, M);
  cout << Space.GetHilbertSpaceDimension() << endl;
  for (int  i = 0; i < Space.GetHilbertSpaceDimension(); i++)
    {
//      cout << i << " " << Space.ChainDescription[i] << " " << ~(Space.ChainDescription[i]) << " " << Space.FindStateIndex(Space.ChainDescription[i]) << endl;
      Space.PrintState(cout ,i) << " " << (*(Space.GetQuantumNumber(i))) << endl;
    }
  TrappedBosonHamiltonian Hamiltonian(&Space, M);
  TestHamiltonian(&Hamiltonian);
  return;

  List<AbstractQuantumNumber*> ListQ = Space.GetQuantumNumbers();
  ListIterator<AbstractQuantumNumber*> IterQ (ListQ);
  AbstractQuantumNumber** TmpQ;
  SubspaceSpaceConverter Converter2;
  while ((TmpQ = IterQ()))
    {
      SubspaceSpaceConverter Converter2;
      AbstractHilbertSpace* Subspace = Space.ExtractSubspace(**TmpQ, Converter2);
      if (Subspace->GetHilbertSpaceDimension() > 0)
	{
	  cout << "Subspace : " << **TmpQ << endl;
	  for (int  i = 0; i < Subspace->GetHilbertSpaceDimension(); i++)
	    Subspace->PrintState(cout ,i) << " " << (*(Subspace->GetQuantumNumber(i))) << endl;
	}
      delete Subspace;
    }
//  AbstractHamiltonian* Hamiltonian = new SpinChainHamiltonianWithTranslations(&Space, NbrSite, 1.0);
}

///////////////////////////////////////////////////////////////////////////

void TestHamiltonian(AbstractHamiltonian* hamiltonian)
{
  if (hamiltonian->GetHilbertSpaceDimension() == 1)
    {
      RealSymmetricMatrix HRep (hamiltonian->GetHilbertSpaceDimension());
      hamiltonian->GetHamiltonian(HRep);
      cout << "eigenvalues : " << endl;
      cout << HRep(0, 0) << endl;
      return;
    }
  if (hamiltonian->GetHilbertSpaceDimension() > 100)
    {
      // Lanczos method
      AbstractArchitecture* Architecture = new SMPArchitecture(2);
      BasicLanczosAlgorithm Lanczos(Architecture, 1);
      int MaxNbrIterLanczos = 200; 
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = 4;
      Lanczos.SetHamiltonian(hamiltonian);
      Lanczos.InitializeLanczosAlgorithm();
      Lanczos.RunLanczosAlgorithm(4);
      while ((Precision > 1e-13) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	{
	  Lanczos.RunLanczosAlgorithm(1);
	  Lowest = Lanczos.GetGroundStateEnergy();
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest;
	  cout  << Lowest << " " << Precision << endl;
	}
    }
  else
    {
      // householder method
      RealSymmetricMatrix HRep (hamiltonian->GetHilbertSpaceDimension());
      RealTriDiagonalSymmetricMatrix TmpTriDiag (hamiltonian->GetHilbertSpaceDimension());
      hamiltonian->GetHamiltonian(HRep);
      HRep.Householder(TmpTriDiag, 1e-7);
      TmpTriDiag.Diagonalize();
      TmpTriDiag.SortMatrixUpOrder();
      cout << "eigenvalues : " << endl;
      for (int j = 0; j < hamiltonian->GetHilbertSpaceDimension(); j++)
	cout << TmpTriDiag.DiagonalElement(j) << " ";
      cout << endl;
    }
}

///////////////////////////////////////////////////////////////////////////

void DMRGFullSpaceTest()
{
// prepare base hamiltonians

  cout.precision(14);
  int NbrSpinBlock = 4;
//  Spin1Chain ChainBlock(NbrSpinBlock, 256);
  Spin1AKLTChain ChainBlock(NbrSpinBlock - 1, 256);
  
  double* CouplingConstantsBlock = new double [NbrSpinBlock];
  for (int i = 0; i < NbrSpinBlock; i++)
    CouplingConstantsBlock[i] = 1.0;
  SpinChainHamiltonian HBlock ((AbstractSpinChain*) &ChainBlock, NbrSpinBlock, CouplingConstantsBlock);
  int NbrSpinInter = 1;
  Spin1Chain ChainInter(NbrSpinInter, 256);
  double* CouplingConstantsInter = new double [NbrSpinInter];
  for (int i = 0; i < NbrSpinInter; i++)
    CouplingConstantsInter[i] = 1.0;
  SpinChainHamiltonian HInter ((AbstractSpinChain*) &ChainInter, NbrSpinInter, CouplingConstantsInter);
  double* TensorCouplingConstants = new double [3];
  TensorCouplingConstants[0] = 1.0;
  TensorCouplingConstants[1] = -1.0;
  TensorCouplingConstants[2] = 1.0;
  BasicInteraction Interaction(TensorCouplingConstants, 3);
  AbstractArchitecture* Architecture = new SMPArchitecture(2);
  BasicLanczosAlgorithmWithEigenstates LanczosAlgorithm(Architecture, 200);
  NonPeriodicDMRGAlgorithm DMRGAlgorithm(&HBlock, &HInter, &Interaction, &Interaction, &LanczosAlgorithm,
					 100);
  SzQuantumNumber GlobalQuantumNumber (0);
  DMRGAlgorithm.Constraint(&GlobalQuantumNumber);

  int NbrIter = 45;
//  int NbrIter = 0;
  double PreviousGroundState = 0.0;
  for (int i = 0; i <= NbrIter; i++)
    {
      DMRGAlgorithm.RunDMRG(i);
      cout << (i + 5) << "  " << DMRGAlgorithm.GetNbrLanczosIteration() << "  " 
	   << DMRGAlgorithm.GetTruncationError() << "  " << DMRGAlgorithm.GetGroundStateEnergy() << "  " 
	   << ((DMRGAlgorithm.GetGroundStateEnergy() - PreviousGroundState) / 2.0) << endl;
      PreviousGroundState = DMRGAlgorithm.GetGroundStateEnergy();
    }
}

///////////////////////////////////////////////////////////////////////////

void AsymmetricDMRG()
{
// prepare base hamiltonians

/*  cout.precision(14);
  int NbrSpinBlock = 4;
  Spin1AKLTChain ChainBlock(NbrSpinBlock - 1, 256);
  
  double* CouplingConstantsBlock = new double [NbrSpinBlock];
  for (int i = 0; i < NbrSpinBlock; i++)
    CouplingConstantsBlock[i] = 1.0;
  SpinChainHamiltonian HBlock ((AbstractSpinChain*) &ChainBlock, NbrSpinBlock, CouplingConstantsBlock);
  int NbrSpinInter = 1;
  Spin1Chain ChainInter(NbrSpinInter, 256);
  double* CouplingConstantsInter = new double [NbrSpinInter];
  for (int i = 0; i < NbrSpinInter; i++)
    CouplingConstantsInter[i] = 1.0;
  SpinChainHamiltonian HInter ((AbstractSpinChain*) &ChainInter, NbrSpinInter, CouplingConstantsInter);
  double* TensorCouplingConstants = new double [3];
  TensorCouplingConstants[0] = 1.0;
  TensorCouplingConstants[1] = -1.0;
  TensorCouplingConstants[2] = 1.0;
  BasicInteraction Interaction(TensorCouplingConstants, 3);
  double* TensorCouplingConstantsLeftRight = new double [3];
  TensorCouplingConstantsLeftRight[0] = 1.0;
  TensorCouplingConstantsLeftRight[1] = -1.0;
  TensorCouplingConstantsLeftRight[2] = 1.0;
  BasicInteraction InteractionLeftRight(TensorCouplingConstantsLeftRight, 3);
  BasicLanczosAlgorithmWithEigenstates LanczosAlgorithm;
  NonPeriodicAsymmetricDMRGAlgorithm DMRGAlgorithm(&HBlock, &HInter, &Interaction, &Interaction, &InteractionLeftRight, &LanczosAlgorithm,
						   100);
  SzQuantumNumber GlobalQuantumNumber (0);
  DMRGAlgorithm.Constraint(&GlobalQuantumNumber);*/

//  int NbrIter = 45;
/*  int NbrIter = 0;
  double PreviousGroundState = 0.0;
  for (int i = 0; i <= NbrIter; i++)
    {
      DMRGAlgorithm.RunDMRG(i);
      cout << (i + 5) << "  " << DMRGAlgorithm.GetNbrLanczosIteration() << "  " 
	   << DMRGAlgorithm.GetTruncationError() << "  " << DMRGAlgorithm.GetGroundStateEnergy() << "  " 
	   << ((DMRGAlgorithm.GetGroundStateEnergy() - PreviousGroundState) / 2.0) << endl;
      PreviousGroundState = DMRGAlgorithm.GetGroundStateEnergy();
    }
  return;*/

  cout.precision(14);
  int NbrBlock = 2;
  double J1 = 1.0;
  double J2 = 1.0;//0.8;
  double J3 = 1.0;//0.3;
  Spin1_2Chain ChainBlock(3 * NbrBlock, 100000);
  OpenDiamondSpinChainHamiltonian HBlock ((AbstractSpinChain*) &ChainBlock, NbrBlock - 1, J1, J2, J3);
  int NbrInter = 1;
  Spin1_2Chain ChainInter(3 * NbrInter, 256);
  OpenDiamondSpinChainHamiltonian HInter ((AbstractSpinChain*) &ChainInter, NbrInter - 1, J1, J2, J3);
  DiamondInteraction Interaction(J2, J1);
  DiamondInteraction Interaction2(J2, J1);
  DiamondInteraction Interaction3(J2, J1);//0.0, 0.0);//J2, J1);
  AbstractArchitecture* Architecture = new MonoProcessorArchitecture();//new SMPArchitecture(2);
  BasicLanczosAlgorithmWithEigenstates LanczosAlgorithm(Architecture, 200);
  NonPeriodicAsymmetricDMRGAlgorithm DMRGAlgorithm(&HBlock, &HInter, &Interaction, &Interaction2, &Interaction3, &LanczosAlgorithm,
						   100);
  SzQuantumNumber GlobalQuantumNumber (0);
  DMRGAlgorithm.Constraint(&GlobalQuantumNumber);

//  int NbrIter = 45;
  int NbrIter = 0;
  double PreviousGroundState = 0.0;
  for (int i = 0; i <= NbrIter; i++)
    {
      cout << "blabla " << i << endl;
      DMRGAlgorithm.RunDMRG(i);
      cout << (i + 5) << "  " << DMRGAlgorithm.GetNbrLanczosIteration() << "  " 
	   << DMRGAlgorithm.GetTruncationError() << "  " << DMRGAlgorithm.GetGroundStateEnergy() << "  " 
	   << ((DMRGAlgorithm.GetGroundStateEnergy() / ((double) (6 * (NbrBlock + i + 1)))) - (PreviousGroundState / ((double) (6 * (NbrBlock + i))))) << " " << ((DMRGAlgorithm.GetGroundStateEnergy() / ((double) (6 * (NbrBlock + i + 1))))) << endl;
      PreviousGroundState = DMRGAlgorithm.GetGroundStateEnergy();
    }
}

void PeriodicDMRG()
{
// prepare base hamiltonians

  cout.precision(14);
  int NbrSpinBlock = 7;
  Spin1_2Chain ChainBlock(NbrSpinBlock, 100000);
//  Spin1Chain ChainBlock(NbrSpinBlock, 100000);
  double* CouplingConstantsBlock = new double [NbrSpinBlock];
  for (int i = 0; i < NbrSpinBlock; i++)
    CouplingConstantsBlock[i] = 1.0;
  SpinChainHamiltonian HBlock ((AbstractSpinChain*) &ChainBlock, NbrSpinBlock, CouplingConstantsBlock);
  int NbrSpinInter = 1;
  Spin1_2Chain ChainInter(NbrSpinInter, 256);
  double* CouplingConstantsInter = new double [NbrSpinInter];
  for (int i = 0; i < NbrSpinInter; i++)
    CouplingConstantsInter[i] = 1.0;
  SpinChainHamiltonian HInter ((AbstractSpinChain*) &ChainInter, NbrSpinInter, CouplingConstantsInter);
  double* TensorCouplingConstants = new double [3];
  TensorCouplingConstants[0] = 1.0;
  TensorCouplingConstants[1] = -1.0;
  TensorCouplingConstants[2] = 1.0;
  BasicInteraction Interaction(TensorCouplingConstants, 3);
  AbstractArchitecture* Architecture = new MonoProcessorArchitecture();//new SMPArchitecture(2);
  BasicLanczosAlgorithmWithEigenstates LanczosAlgorithm(Architecture, 200);
  PeriodicDMRGAlgorithm DMRGAlgorithm(&HBlock, &HInter, &Interaction, &Interaction, &LanczosAlgorithm, 150);//HBlock.GetHilbertSpaceDimension());
  SzQuantumNumber GlobalQuantumNumber (0);
  DMRGAlgorithm.Constraint(&GlobalQuantumNumber);

//  int NbrIter = 45;
  int NbrIter = 18;
  double PreviousGroundState = 0.0;
  for (int i = 0; i <= NbrIter; i++)
    {
      DMRGAlgorithm.RunDMRG(i);
      cout << (2 * ((i + 1 ) * NbrSpinInter + NbrSpinBlock)) << "  " << DMRGAlgorithm.GetNbrLanczosIteration() << "  " 
	   << DMRGAlgorithm.GetTruncationError() << "  " << DMRGAlgorithm.GetGroundStateEnergy() << endl;
      PreviousGroundState = DMRGAlgorithm.GetGroundStateEnergy();
    }
}



////////////////////////////////////////////////////////////////////////////

void BlockTest()
{
  for (int n = 0; n < 10; n++)
    {
  int SrcSize = (int) (((double) rand()) / ((double) RAND_MAX) * 10.0 + 10.0);
  int DestSize = (int) (((double) rand()) / ((double) RAND_MAX) * 10.0 + 10.0);
  cout << SrcSize << " " << DestSize << endl;
  int* RowDiv = new int [SrcSize / 2 + 1];
  int Pos = 0;
  int NbrBlock = 0;
  while (Pos < SrcSize)
    {
      RowDiv[NbrBlock++] = Pos;
      Pos += (int) (((double) rand()) / ((double) RAND_MAX) * (SrcSize / 2) + 1.0);
      cout << "row=" << Pos << endl;
    }
  RowDiv[NbrBlock] = SrcSize;

  int* ColumnDiv = new int [DestSize / 2 + 1];
  Pos = 0;
  int i = 0;
  while (i < NbrBlock)
    {
      ColumnDiv[i++] = Pos;
      Pos += (int) (((double) rand()) / ((double) RAND_MAX) * (DestSize - 1 - Pos - (NbrBlock - i - 1)) + 1.0);
      cout << Pos << endl;
    }
  ColumnDiv[NbrBlock] = DestSize;

  RealSymmetricMatrix M = RandomRealSymmetricMatrix(SrcSize, 5);
  BlockDiagonalMatrix U;
  RealMatrix U2 (SrcSize, DestSize, true);
  i = 0;
  while (i < NbrBlock)
    {
      RealMatrix TmpM = RandomRealMatrix (RowDiv[i + 1] - RowDiv[i], ColumnDiv[i + 1] - ColumnDiv[i], 5);
      U += &TmpM; 
      for (int j = 0; j < (RowDiv[i + 1] - RowDiv[i]); j++)
	for (int k = 0; k < (ColumnDiv[i + 1] - ColumnDiv[i]); k++)
	  U2(RowDiv[i] + j, ColumnDiv[i] + k) = TmpM(j, k);
       i++;
    }
 
  /*
  RealMatrix M1 = RandomRealMatrix (2, 2, 5);
  M1(0,0) = 1;
  M1(0,1) = 2;
  M1(1,0) = 3;
  M1(1,1) = -4;
  RealMatrix M3 = RandomRealMatrix(2, 2, 5);
  M3(0,0) = 1;
  M3(0,1) = 2;
  M3(1,0) = 3;
  M3(1,1) = -4;
  RealMatrix M2 = RandomRealMatrix(3, 4, 5);
  M2(0,0) = 5;
  M2(0,1) = 2;
  M2(1,0) = -1;
  M2(1,1) = 3;
  M2(2,0) = -2;
  M2(2,1) = 1;
  M2(0,2) = 5;
  M2(0,3) = 2;
  M2(1,2) = -1;
  M2(1,3) = 3;
  M2(2,2) = -2;
  M2(2,3) = 1;
  BlockDiagonalMatrix U;
  U += &M1;
  U += &M2;
  U += &M3;
  RealMatrix U2 (7, 8, true);
  U2(0,0) = 1;
  U2(0,1) = 2;
  U2(1,0) = 3;
  U2(1,1) = -4;
  U2(2,2) = 5;
  U2(2,3) = 2;
  U2(3,2) = -1;
  U2(3,3) = 3;
  U2(4,2) = -2;
  U2(4,3) = 1;
  U2(2,4) = 5;
  U2(2,5) = 2;
  U2(3,4) = -1;
  U2(3,5) = 3;
  U2(4,4) = -2;
  U2(4,5) = 1;
  U2(5,6) = 1;
  U2(5,7) = 2;
  U2(6,6) = 3;
  U2(6,7) = -4;
  RealSymmetricMatrix M = RandomRealSymmetricMatrix(7, 5);
  M(0, 0) = 1;
  M(1, 1) = 2;
  M(2, 2) = 3;
  M(3, 3) = 4;
  M(4, 4) = 5;
  M(5, 5) = 6;
  M(6, 6) = 7;
  M(0, 1) = -4;
  M(0, 2) = -1;
  M(0, 3) = 2;
  M(0, 4) = 1;
  M(0, 5) = -3;
  M(0, 6) = 4;
  M(1, 2) = 2;
  M(1, 3) = 3;
  M(1, 4) = 4;
  M(1, 5) = -3;
  M(1, 6) = 5;
  M(2, 3) = 1;
  M(2, 4) = -1;
  M(2, 5) = 2;
  M(2, 6) = -3;
  M(3, 4) = 3;
  M(3, 5) = 1;
  M(3, 6) = -2;
  M(4, 5) = 3;
  M(4, 6) = -2;
  M(5, 6) = -4;*/
  //  cout << M << endl;
  RealSymmetricMatrix* Mc = (RealSymmetricMatrix*) M.Conjugate(U);
  RealSymmetricMatrix* Mc2 = (RealSymmetricMatrix*) M.Conjugate(U2);
  bool Flag = false;
  for (int j = 0; j < DestSize; j++)
    for (int k = j; k < DestSize; k++)
      if ((*Mc)(j, k) != (*Mc2)(j, k))
	Flag = true;
  if (Flag == true)
    {
      cout << "error -------------------------" << endl;
      cout << *Mc << endl;
      cout << *Mc2 << endl;
    }
  delete Mc;
  delete Mc2;
}
}

RealMatrix RandomRealMatrix(int nbrRow, int nbrColumn, double range)
{
  RealMatrix M(nbrRow, nbrColumn);
  for (int i = 0; i < nbrRow; i++)
    for (int j = 0; j < nbrColumn; j++)
      {
	M(i, j) = (range * (double) (rand() - RAND_MAX / 2)) / RAND_MAX;
      }
  return M;
}

RealSymmetricMatrix RandomRealSymmetricMatrix(int nbrRow, double range)
{
  RealSymmetricMatrix M(nbrRow);
  for (int i = 0; i < nbrRow; i++)
    for (int j = i; j < nbrRow; j++)
      {
	M(i, j) = (range * (double) (rand() - RAND_MAX / 2)) / RAND_MAX;
      }
  return M;
}

RealTriDiagonalSymmetricMatrix RandomRealTriDiagonalSymmetricMatrix(int nbrRow, double range)
{
  RealTriDiagonalSymmetricMatrix M(nbrRow);
  for (int i = 0; i < (nbrRow - 1); i++)
    {
      M.SetMatrixElement(i, i, (range * (double) (rand() - RAND_MAX / 2)) / RAND_MAX);
      M.SetMatrixElement(i, i + 1, (range * (double) (rand() - RAND_MAX / 2)) / RAND_MAX);
    }
  M.SetMatrixElement(nbrRow - 1, nbrRow - 1, (range * (double) (rand() - RAND_MAX / 2)) / RAND_MAX);
  return M;
}
