#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/ParticleOnSphereCoulombHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "MathTools/ClebschGordanCoefficients.h"

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
/*  int J1 = 4;
  int J2 = 4;
  int j;
  double coef;
  ClebschGordanCoefficients Coefficients(J1, J2);
  for (int m1 = -J1; m1 <= J1; m1 += 2)
//  int m1 = 2;
    for (int m2 = -J2; m2 <= J2; m2 += 2)
      {
	Coefficients.InitializeCoefficientIterator(m1, m2);
	while (Coefficients.Iterate(j, coef))
	  {
	    Coefficients.PrintCoefficient(cout, m1, m2, j) << endl;
	  }
      }  
  return 0;*/
  int NbrBosons = 5;
  if (argc >= 2)
    NbrBosons = atoi (argv[1]);
  int InvNu = 2;
  int LzMax = InvNu * (NbrBosons - 1);
  int L = 0;
  double GroundStateEnergy = 0.0;
/*  if (argc >= 3)
    L = atoi (argv[2]);*/
/*  if (argc >= 3)
    LzTotal = atoi (argv[3]);*/
  char* OutputNameLz = "bosons_n_9_lz.dat";
  char* OutputNameL = "bosons_l.dat";
/*  if (argc >= 4)
    OutputName = argv[4];*/
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);
  LzMax = InvNu * (NbrBosons - 1);
  int Max = (LzMax * NbrBosons);
  int TotalSize = 0;
  double** Eigenvalues = new double* [2 * Max + 1];
  int* Dimensions = new int [2 * Max + 1];
  Max = 30;
  for (int  L = 0; L <= Max; L += 2)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << " LzTotal = " << L << endl;
      BosonOnSphere Space (NbrBosons, L, LzMax);
//      FermionOnSphere Space (NbrBosons, L, LzMax);
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;
      TotalSize += Space.GetHilbertSpaceDimension();
/*     for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  Space.PrintState(cout, i) << endl;
	}*/
      ParticleOnSphereCoulombHamiltonian Hamiltonian(&Space, NbrBosons, LzMax, 0);
      if (Hamiltonian.GetHilbertSpaceDimension() < 200)
	{
//	  Dimensions[L >> 1] = Hamiltonian.GetHilbertSpaceDimension();
//	  Eigenvalues[L >> 1] = new double [Hamiltonian.GetHilbertSpaceDimension()];
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
//	  cout << HRep << endl;
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
//		  Eigenvalues[L >> 1][j] = TmpTriDiag.DiagonalElement(j);
//		  cout << TmpTriDiag.DiagonalElement(j) << " ";
		  File << (L / 2) << " " << TmpTriDiag.DiagonalElement(j) << endl;
// (TmpTriDiag.DiagonalElement(j) - GroundStateEnergy) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
//	      Eigenvalues[L >> 1][0] = HRep(0, 0);
//	      cout << HRep(0, 0) << endl;
	      //	      GroundStateEnergy = HRep(0, 0);
	      File << (L / 2) << " " << HRep(0, 0) << endl;// - GroundStateEnergy) / (4 * M_PI)) << endl;
	    }
	}
      else
	{
//	  AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
	  AbstractArchitecture* Architecture = new SMPArchitecture(2);
	  int MaxNbrIterLanczos = 1000; 
//	  BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
	  FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
	  int NbrEigenvalue = 25;
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
	  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Precision > 1e-7) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
	    {
	      Lanczos.RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue);//Lanczos.GetGroundStateEnergy();
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << (L / 2) << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
//	  File << (L / 2) << " " << Lowest << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
    }
  File.close();
/*  File.open(OutputNameL, ios::binary | ios::out);
  int SpectrumSize = 1;
  double* Spectrum = new double [TotalSize];
  bool* Degeneracy = new bool [TotalSize];
  bool Flag;
  double TmpEigenvalue;
  Spectrum[0] = Eigenvalues[Max >> 1][0];
  File << (Max >> 1) << " " << (Spectrum[0] - GroundStateEnergy) << endl;
  for (int  L = ((Max >> 1) - 1); L >= 0; --L)
    {
      for (int j = 0; (j < SpectrumSize); ++j)
	Degeneracy[j] = false;
      for (int i = 0; i < Dimensions[L]; ++i)
	{
	  Flag = false;
	  TmpEigenvalue = Eigenvalues[L][i];
	  for (int j = 0; ((j < SpectrumSize) && (Flag == false)); ++j)
	    if ((Degeneracy[j] == false) && (fabs((TmpEigenvalue - Spectrum[j]) / TmpEigenvalue) < 1e-7))
	      {
		Flag = true;
		Degeneracy[j] = true;
	      }
	  if (Flag == false)
	    {	      
	      Spectrum[SpectrumSize] = TmpEigenvalue;
	      File << L << " " << (TmpEigenvalue - GroundStateEnergy) << endl;
 	      Degeneracy[SpectrumSize] = false;
	      ++SpectrumSize;
	    }
	}
       delete[] Eigenvalues[L];
    }
  delete[] Degeneracy;
  delete[] Spectrum;
  delete[] Eigenvalues;
  delete[] Dimensions;
  File.close();*/


  return 0;
}
