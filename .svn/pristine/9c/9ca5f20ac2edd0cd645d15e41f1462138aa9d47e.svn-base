#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "Hamiltonian/MultipleSpinChainHamiltonian.h"
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


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
 
  int NbrSpinPerChain = 2;
  if (argc >= 2)
    NbrSpinPerChain = atoi (argv[1]);
  int NbrChain = 3;
  if (argc >= 3)
    NbrChain = atoi (argv[2]);
  int NbrSpin = NbrSpinPerChain * NbrChain;

  double J1 = 1.0;
  if (argc >= 4)
    J1 = atof (argv[3]);
  double J2 = 1.0;
  if (argc >= 5)
    J2 = atof (argv[4]);
  double J3 = 1.0;
  if (argc >= 6)
    J3 = atof (argv[5]);
  double J4 = 1.0;
  if (argc >= 7)
    J4 = atof (argv[6]);
  

  Spin1_2Chain Chain(NbrSpin, 0, 1000000);

  cout << NbrSpinPerChain << " " << NbrChain << " " << J1 << " " << J2 << " " << J3 << " " << J4 << endl;
  cout << Chain.GetHilbertSpaceDimension() << endl;

  MultipleSpinChainHamiltonian Hamiltonian ((AbstractSpinChain*) &Chain, NbrSpinPerChain, NbrChain, J1, J2, J3, J4);
  
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
      AbstractArchitecture* Architecture = new SMPArchitecture(2);
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

