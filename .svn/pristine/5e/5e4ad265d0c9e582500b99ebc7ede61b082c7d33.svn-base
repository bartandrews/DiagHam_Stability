#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"

#include "Operator/ParticlePolarizationOperator.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  int NbrFermions = 6;
  if (argc >= 2)
    NbrFermions = atoi (argv[1]);
  int InvNu = 3;
  int MaxMomentum = (InvNu * NbrFermions) / 2;
  int L = 0;
//  double GroundStateEnergy = 0.0;
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
  if (argc >= 4)
    Max = atoi (argv[3]);
  ofstream File;
//  File.open(OutputNameLz, ios::binary | ios::out);
//  int Max = ((MaxMomentum - NbrFermions) * NbrFermions);
  int TotalSize = 0;
//  XRatio = 0.7;
//  XRatio = 1.0 / XRatio;
//  NbrFermions = 3;
//  XRatio = 1.0;
//  MaxMomentum = 3;
  double Zeeman = 0.0;
  double ZeemanMax = 1.0;
  int NbrIter = 101;
  double ZeemanInc = (ZeemanMax - Zeeman) / ((double) (NbrIter - 1));
  for (; MaxMomentum <= Max; ++MaxMomentum)
    {     
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
//      cout << " LzTotal = " << L << endl;
      FermionOnTorusWithSpin Space (NbrFermions, MaxMomentum);
      List<AbstractQuantumNumber*> QuantumNumber (Space.GetQuantumNumbers ());
      ListIterator<AbstractQuantumNumber*> IterQuantumNumber (QuantumNumber);
      AbstractQuantumNumber** TmpQuantumNumber;
      FermionOnTorusWithSpin* TmpSpace = 0;
      double* GroundStateEnergy = new double [NbrIter];
      int* GroundStateSz = new int [NbrIter];
      for (int i = 0; i < NbrIter; ++i)
	GroundStateEnergy[i] = 0.0;
      while ((TmpQuantumNumber = IterQuantumNumber()))
	{
	  (*TmpQuantumNumber)->PrintQuantumNumber(cout) << endl;
      //	  cout << (*(SzQuantumNumber*) *TmpQuantumNumber) << endl;
	  int CurrentSz = ((SzQuantumNumber*) ((*((VectorQuantumNumber*) *TmpQuantumNumber))[0]))->GetSz();
	  SubspaceSpaceConverter Dummy;
	  if (TmpSpace != 0)
	    delete TmpSpace;
	  TmpSpace = (FermionOnTorusWithSpin*) Space.ExtractSubspace(**TmpQuantumNumber, Dummy);//new FermionOnTorusWithSpin ();
/*	  for (int i = 0; i < TmpSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      cout << i << " = ";
	      TmpSpace->PrintState(cout, i) << endl;
	    }*/
	  ParticleOnTorusCoulombWithSpinHamiltonian Hamiltonian(TmpSpace, NbrFermions, MaxMomentum, XRatio, Zeeman);
	  double CurrentZeeman = Zeeman;
	  for (int ZeemanLoop = 0; ZeemanLoop < NbrIter; ++ZeemanLoop)
	    {
	      Hamiltonian.SetMagneticG(CurrentZeeman);
	      if (Hamiltonian.GetHilbertSpaceDimension() < 100)
		{
		  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
		  Hamiltonian.GetHamiltonian(HRep);
		  if (Hamiltonian.GetHilbertSpaceDimension() > 1)
		    {
		      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
		      HRep.Householder(TmpTriDiag, 1e-7);
		      TmpTriDiag.Diagonalize();
		      TmpTriDiag.SortMatrixUpOrder();
		      if ((GroundStateEnergy[ZeemanLoop] == 0.0) || (GroundStateEnergy[ZeemanLoop] > TmpTriDiag.DiagonalElement(0)))
			{
			  GroundStateEnergy[ZeemanLoop] = TmpTriDiag.DiagonalElement(0);
			  GroundStateSz[ZeemanLoop] = CurrentSz;
			}
		      for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
			{
			  cout << TmpTriDiag.DiagonalElement(j) << " ";
			}
		      cout << endl;
		    }
		  else
		    {
		      cout << HRep(0, 0) << endl;;
		      if ((GroundStateEnergy[ZeemanLoop] == 0.0) || (GroundStateEnergy[ZeemanLoop] > HRep(0, 0)))
			{
			  GroundStateEnergy[ZeemanLoop] = HRep(0, 0);
			  GroundStateSz[ZeemanLoop] = CurrentSz;
			}
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
		  if ((GroundStateEnergy[ZeemanLoop] == 0.0) || (GroundStateEnergy[ZeemanLoop] > Lowest))
		    {
		      GroundStateEnergy[ZeemanLoop] = Lowest;
		      GroundStateSz[ZeemanLoop] = CurrentSz;
		    }
		  cout << endl;
		  cout << Lowest << " " << Precision << "  Nbr of iterations = " << CurrentNbrIterLanczos << endl;
		  cout << "------------------------------------------------------------------" << endl << endl;;
		  gettimeofday (&(TotalEndingTime), 0);
		  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
		  cout << "time = " << Dt << endl;
		}
	      CurrentZeeman += ZeemanInc;
	    }
	}
      cout << "------------------------------------------------------------------" << endl ;
      cout << "------------------------------------------------------------------" << endl << endl;
      for (int i = 0; i < NbrIter; ++i)
	cout << "Zeeman = " << (Zeeman + (((double) i) * ZeemanInc)) << "  Ground State = " << GroundStateEnergy[i] 
	     << "  with Sz value " << GroundStateSz[i] << endl;
      cout << "------------------------------------------------------------------" << endl ;
      cout << "------------------------------------------------------------------" << endl << endl;
    }
  
  return 0;




 /*     for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  Space.PrintState(cout, i) << endl;
	}
      List<AbstractQuantumNumber*> QuantumNumber (Space.GetQuantumNumbers ());
      ListIterator<AbstractQuantumNumber*> IterQuantumNumber (QuantumNumber);
      AbstractQuantumNumber** TmpQuantumNumber;
      while ((TmpQuantumNumber = IterQuantumNumber()))
	{
	  SubspaceSpaceConverter Dummy;
	  FermionOnTorusWithSpin* TmpSpace = (FermionOnTorusWithSpin*) Space.ExtractSubspace(**TmpQuantumNumber, Dummy);
	  cout << "Sz =  " << (*(SzQuantumNumber*) *TmpQuantumNumber) << endl;
	  for (int i = 0; i < TmpSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      cout << i << " = ";
	      TmpSpace->PrintState(cout, i) << endl;
	    }
	}
      double TmpCoef;
      for (int m = 0; m < Space.GetHilbertSpaceDimension(); ++m)
	for (int i = 0; i < MaxMomentum; ++i)
	  for (int j = 0; j <= i; ++j)
	    for (int k = 0; k < MaxMomentum; ++k)
	      for (int l = 0; l <= k; ++l)
		{
//		  cout << "acting on " << m  << " with " << i << " "  << j << " "   << k << endl;
		  int TmpIndex = Space.AudAddAdAd(m, i, j, k ,l, TmpCoef);
		  if (TmpIndex < Space.GetHilbertSpaceDimension())
		    cout << "acting on " << m  << " with " << i << " "  << j << " "   << k << " "    << l << " gives "  << TmpIndex << " with coef " 
			 <<  TmpCoef << endl;
		}
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;
      return 0;

     ParticleOnTorusCoulombWithSpinHamiltonian Hamiltonian(TmpSpace, NbrFermions, MaxMomentum, XRatio, 100.0);
     if (Hamiltonian.GetHilbertSpaceDimension() < 100)
	{
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
	  cout << HRep << endl;
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
	  BasicLanczosAlgorithm Lanczos(Architecture);
	  int MaxNbrIterLanczos = 200; 
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
//	  RealVector& GroundState = (RealVector&) Lanczos.GetGroundState();
//	  ParticlePolarizationOperator Operator(&Space, NbrFermions);
//	  cout << Operator.MatrixElement(GroundState, GroundState);
	}
	}
      XRatio += 0.01;
      cout << "----------------------------------------------------------------" << endl;
//      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrFermions) << endl;
//      File << (((double) NbrFermions) / ((double) MaxMomentum)) << " " << (GroundStateEnergy / (double) NbrFermions) << endl;
    }
//  File.close();
*/

  return 0;
}
