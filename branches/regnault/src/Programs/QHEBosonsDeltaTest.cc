#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/ParticleOnSphereDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "BitmapPicture/AbstractBitmapPicture.h"
#include "BitmapPicture/TgaFormat.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


// return RGB value corresponding to a given scalar value
//
// x = scalar to use (0.0 = blue, 1.0 = red)
// grayScale = true if a grayscale value has to be return
// return value = RGB value
RGB RGBValue (double x, bool grayScale = false);


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
  int NbrBosons = 8;
  if (argc >= 2)
    NbrBosons = atoi (argv[1]);
  int InvNu = 2;
  int LzMax = 14;
  double GroundStateEnergy = 0.0;
  int Shift = 0;
  if (argc >= 3)
    LzMax = atoi (argv[2]);
/*  if (argc >= 3)
    LzTotal = atoi (argv[3]);*/
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "bosons_delta_n_%d_2s_%d_lz.dat", NbrBosons, LzMax);
  char* OutputNameL = "bosons_l.dat";
/*  if (argc >= 4)
    OutputName = argv[4];*/
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);
//  LzMax = InvNu * (NbrBosons - 1) + Shift;
  int Max = (LzMax * NbrBosons);
  int TotalSize = 0;
  double** Eigenvalues = new double* [2 * Max + 1];
  int* Dimensions = new int [2 * Max + 1];
//  Max = 30;
  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (argc >= 4)
    L = atoi (argv[3]);
  Max = L;
  for (; L <= Max; L += 2)
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
      ParticleOnSphereDeltaHamiltonian Hamiltonian(&Space, NbrBosons, LzMax);
      if (Hamiltonian.GetHilbertSpaceDimension() < 200)
	{
//	  Dimensions[L >> 1] = Hamiltonian.GetHilbertSpaceDimension();
//	  Eigenvalues[L >> 1] = new double [Hamiltonian.GetHilbertSpaceDimension()];
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
/*	  AbstractBitmapPicture* TmpPic = Hamiltonian.GetHamiltonianColorPicture(1e-10);
	  char* TmpPicName = new char [256];
	  sprintf (TmpPicName, "pic_lz_%d.tga", L);
	  TmpPic->SavePicture(TmpPicName);*/
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
//	      for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension() ; j++)
		{
//		  Eigenvalues[L >> 1][j] = TmpTriDiag.DiagonalElement(j);
//		  cout << TmpTriDiag.DiagonalElement(j) << " ";
		  File << (L / 2) << " " << TmpTriDiag.DiagonalElement(0) << endl;
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
	  AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
//	  AbstractArchitecture* Architecture = new SMPArchitecture(2);
	  int MaxNbrIterLanczos = 4000; //2000
//	  BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
	  int NbrEigenvalue = 2;
	  FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
//	  ImplicitRestartedArnoldiWithReorthogonalizationAlgorithm Lanczos(Architecture, NbrEigenvalue, 20);//2 * NbrEigenvalue);
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos.SetHamiltonian(&Hamiltonian);
	  Lanczos.InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Lanczos.TestConvergence() == false) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))//(NbrEigenvalue + 12)))
	    {
	      Lanczos.RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);//Lanczos.GetGroundStateEnergy();
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	    }
	  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      File << "too much Lanczos iterations" << endl;
	      File.close();
	      exit(0);
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i < NbrEigenvalue; ++i)
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
	  return 0;
	  RealVector* Eigenstates = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
	  Complex Sum (0.0, 0.0);
	  Complex Sum2 (0.0, 0.0);
	  RealVector Value(2);
	  Complex TmpValue;
	  Complex TmpValue2;
	  int PicSize = 800;
	  double XInc = M_PI / ((double) PicSize);
	  double YInc = 2.0 * M_PI / ((double) PicSize);
	  double X = 0.0;
	  double Y = 0.0;
	  double** TmpArray = new double* [PicSize];
	  ParticleOnSphereFunctionBasis Basis(LzMax);
	  double MaxValue = 0.0;
/*	  for (int i = 0; i <= LzMax; ++i)
	    {
	      X = 0.0;
	      for (int x = 0; x < PicSize; ++x)
		{
		  Value[0] = X;
		  Y = 0.0;
		  for (int y = 0; y < PicSize; ++y)
		    {
		      Value[1] = Y;
		      Basis.GetFunctionValue(Value, TmpValue, i + 1);
		      Basis.GetFunctionValue(Value, TmpValue2, i);
		      Sum += TmpValue * Conj(TmpValue2) * sin (X);
		      Y += YInc;
		    }
		  X += XInc;
		}
	      Sum *= XInc * YInc;
	      cout << i << " " << Sum << endl;
	    }*/
	  int NbrPrecalculatedValues = 0;
	  Complex* PrecalculatedValues = new Complex [LzMax + 1];
	  Complex* PrecalculatedValues2 = new Complex [LzMax + 1];
	  Value[0] = 0.0;
	  Value[1] = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, LzMax);
	      ParticleOnSphereDensityDensityOperator Operator (&Space, i, LzMax, i, LzMax);
	      PrecalculatedValues[NbrPrecalculatedValues] = Operator.MatrixElement(*Eigenstates, *Eigenstates) * TmpValue * Conj(TmpValue);
	      ParticleOnSphereDensityOperator Operator2 (&Space, i);
	      PrecalculatedValues2[NbrPrecalculatedValues] = Operator2.MatrixElement(*Eigenstates, *Eigenstates);
	      cout << PrecalculatedValues2[NbrPrecalculatedValues] << endl;
	      ++NbrPrecalculatedValues;
	    }
	  ofstream File2;
	  char* OutputNameCorr = new char [256];
	  sprintf (OutputNameCorr, "bosons_delta_corr_n_%d_2s_%d.dat", NbrBosons, LzMax);
	  File2.open(OutputNameCorr, ios::binary | ios::out);
	  File2.precision(14);
	  for (int x = 0; x < PicSize; ++x)
	    {
	      Value[0] = X;
	      int Pos = 0;
	      Sum.Re = 0.0;
	      Sum.Im = 0.0;		  
	      Sum2.Re = 0.0;
	      Sum2.Im = 0.0;		  
	      for (int i = 0; i <= LzMax; ++i)
		{
		  Basis.GetFunctionValue(Value, TmpValue, i);
//		  Sum2 += PrecalculatedValues2[Pos] * (Conj(TmpValue) * TmpValue);
		  Sum -= PrecalculatedValues[Pos] * (Conj(TmpValue) * TmpValue);
		  ++Pos;
		}
//	      Sum /= Sum2;// * PrecalculatedValues2[LzMax];
	      File2 << (X * sqrt (0.5 * LzMax )) << " " << Norm(Sum)  * (16.0 * M_PI * M_PI) / ((double) (NbrBosons * NbrBosons))<< endl;
//	      File2 << (X * sqrt (0.5 * LzMax )) << " " << Norm(Sum2) << endl;
	      X += XInc;
	    }
	  File2.close();
/*	  for (int x = 0; x < PicSize; ++x)
	    {
	      Value[0] = X;
	      Y = 0.0;
	      cout << x << endl;
	      TmpArray[x] = new double [PicSize];
	      for (int y = 0; y < PicSize; ++y)
		{
		  Value[1] = Y;
		  int Pos = 0;
		  Sum.Re = 0.0;
		  Sum.Im = 0.0;		  
		  for (int i = 0; i <= LzMax; ++i)
		    {
		      Basis.GetFunctionValue(Value, TmpValue, i);
		      Sum -= PrecalculatedValues[Pos] * (Conj(TmpValue) * TmpValue);
		      //			      Sum += PrecalculatedValues[Pos] * (Conj(TmpValue) * TmpValue);
		      ++Pos;
		    }
		  TmpArray[x][y] = Norm(Sum);
		  if (TmpArray[x][y] > MaxValue)
		    MaxValue = TmpArray[x][y];
		  Y += YInc;
		}
	      X += XInc;
	    }
	  
	  TgaFormat TmpPic (PicSize, PicSize);
	  MaxValue = 1.0 / MaxValue;
	  for (int x = 0; x < PicSize; ++x)
	    for (int y = 0; y < PicSize; ++y)
	      {
		PicRGB TmpColor(RGBValue(TmpArray[y][x] * MaxValue, false));
		TmpPic.SetPixel(x, y, TmpColor);
	      }
	  TmpPic.SavePicture("pic.tga");*/
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
    }
  File.close();
  return 0;
}

RGB RGBValue (double x, bool grayScale)
{
  if (grayScale == false)
    {
      if (x >= 1.0)
	return RGB (1.0, 0.0, 0.0);
      if (x <= 0.0)
	return RGB (0.0, 0.0, 1.0);
      if (x <= 0.25)
	return RGB (0.0, 4.0 * x, 1.0);
      if (x <= 0.5)
	return RGB (0.0, 1.0, 1.0 - 4.0 * (x - 0.25));
      if (x <= 0.75)
	return RGB (4.0 * (x - 0.5), 1.0, 0.0);
      return RGB (1.0, 1.0 - 4.0 * (x - 0.75), 0.0);  
    }
      if (x >= 1.0)
	return RGB (1.0, 1.0, 1.0);
      if (x <= 0.0)
	return RGB (0.0, 0.0, 0.0);
      return RGB (x, x, x);  
  
}
