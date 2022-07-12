#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/FullFermionOnTorus.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "Hamiltonian/ParticleOnTorusCoulombDMRGReadyHamiltonian.h"

#include "Interaction/InternalInteraction/ParticleOnTorusCoulombInternalInteraction.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  int NbrFermions = 4;
  if (argc >= 2)
    NbrFermions = atoi (argv[1]);
  int InvNu = 3;
  int MaxMomentum = InvNu * NbrFermions;
  int L = 0;
  double GroundStateEnergy = 0.0;
/*  if (argc >= 3)
    LzTotal = atoi (argv[3]);*/
  char* OutputNameLz = "fermions_torus.dat";
  char* OutputNameL = "fermions_torus.dat";
  double XRatio = NbrFermions / 4.0;
  if (argc >= 3)
    {
      MaxMomentum = atoi (argv[2]);
    }
  int Max = MaxMomentum;
  int MomentumConstraint = 0;
  if (argc >= 4)
    MomentumConstraint = atoi (argv[3]);
//  if (argc >= 4)
//    Max = atoi (argv[3]);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
//  int Max = ((MaxMomentum - NbrFermions) * NbrFermions);
  int TotalSize = 0;
//  XRatio = 0.7;
//  XRatio = 1.0 / XRatio;
//  NbrFermions = 3;
//  XRatio = 1.0;
//  MaxMomentum = 4;
//  for (; MaxMomentum <= Max; ++MaxMomentum)
  cout << "----------------------------------------------------------------" << endl;
  cout << " Ratio = " << XRatio << endl;
  //      cout << " LzTotal = " << L << endl;
  NbrFermions = 6;
  MaxMomentum = 6;
  FullFermionOnTorus TotalSpace (NbrFermions, MaxMomentum);//, Momentum);
  List<AbstractQuantumNumber*> QuantumNumbers (TotalSpace.GetQuantumNumbers());
  for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
    {
      cout << i << " = ";
      TotalSpace.PrintState(cout, i) << endl;
    }
  ListIterator<AbstractQuantumNumber*> QuantumNumberIter (QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber;

  ParticleOnTorusCoulombDMRGReadyHamiltonian Hamiltonian(&TotalSpace, NbrFermions, MaxMomentum, XRatio);
  if (Hamiltonian.GetHilbertSpaceDimension() < 300)
    {
      	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension(), true);
//	  Hamiltonian.GetHamiltonian(HRep);
//	  cout << HRep << endl;
	  ParticleOnTorusCoulombInternalInteraction InternalInteraction(true, MaxMomentum, XRatio);
	  List<Matrix*> ListOperators (Hamiltonian.LeftInteractionOperators());
	  InternalInteraction.AddInteraction (HRep, ListOperators);
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	  HRep.Householder(TmpTriDiag, 1e-7);
	  TmpTriDiag.Diagonalize();
	  TmpTriDiag.SortMatrixUpOrder();
	  cout << "eigenvalues : " << endl;
	  for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
	    cout << TmpTriDiag.DiagonalElement(j) << " ";
	  cout << endl << endl;
//	  cout << HRep << endl;
      	  RealSymmetricMatrix HRep2 (Hamiltonian.GetHilbertSpaceDimension(), true);
	  Hamiltonian.GetHamiltonian(HRep2);
	  RealTriDiagonalSymmetricMatrix TmpTriDiag2 (Hamiltonian.GetHilbertSpaceDimension());
	  HRep2.Householder(TmpTriDiag2, 1e-7);
	  TmpTriDiag2.Diagonalize();
	  TmpTriDiag2.SortMatrixUpOrder();
	  cout << "eigenvalues : " << endl;
	  for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
	    cout << TmpTriDiag2.DiagonalElement(j) << " ";
//	  cout << endl << HRep2 << endl;
	}
/*  while ((TmpQuantumNumber = QuantumNumberIter()))
    {
      SubspaceSpaceConverter Converter;
      FermionOnTorus* Space = (FermionOnTorus*) TotalSpace.ExtractSubspace (**TmpQuantumNumber, Converter);
      if (Space == 0)
	{
	  cout << "no subspace corresponding to " << **TmpQuantumNumber << endl;
	}
      else
	{
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      cout << i << " = ";
	      Space->PrintState(cout, i) << endl;
	    }
	  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;
	  cout << **TmpQuantumNumber << endl;
	}
    }*/
  return 0;
/*  {
      ParticleOnTorusCoulombDMRGReadyHamiltonian Hamiltonian(Space, NbrFermions, MaxMomentum, XRatio);
       if (Hamiltonian.GetHilbertSpaceDimension() < 300)
	{
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
//	  cout << HRep << endl;
	  ParticleOnTorusCoulombInternalInteraction InternalInteraction(true, XRatio);
	  List<Matrix*> ListOperators (Hamiltonian.LeftInteractionOperators());
	  InternalInteraction.AddInteraction (HRep, ListOperators);
	  if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      //	  cout << "eigenvalues : " << endl;
	      for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
		{
		  cout << TmpTriDiag.DiagonalElement(j) << " ";
//		  File << (L / 2) << " " << TmpTriDiag.DiagonalElement(j) << endl;
// (TmpTriDiag.DiagonalElement(j) - GroundStateEnergy) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      cout << HRep(0, 0) << endl;;
//	      File << (L / 2) << " " << HRep(0, 0) << endl;// - GroundStateEnergy) / (4 * M_PI)) << endl;
	    }
	}
      else
	{
//	  AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
	  AbstractArchitecture* Architecture = new SMPArchitecture(2);
	  FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture);
//	  BasicLanczosAlgorithm Lanczos(Architecture);
	  int MaxNbrIterLanczos = 1000; 
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 4;
	  Lanczos.SetHamiltonian(&Hamiltonian);
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
//	      cout << Lowest << " " << Precision << endl;
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << Lowest << " " << Precision << "  Nbr of iterations = " << CurrentNbrIterLanczos << endl;
	  cout << "------------------------------------------------------------------" << endl << endl;;
//	  File << (L / 2) << " " << Lowest << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
//      XRatio += 0.01;
      cout << "----------------------------------------------------------------" << endl;
//      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrFermions) << endl;
      File << (((double) NbrFermions) / ((double) MaxMomentum)) << " " << (GroundStateEnergy / (double) NbrFermions) << endl;
    }
  File.close();
*/

  return 0;
}
