#include "Matrix/ComplexTriDiagonalHermitianMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "Hamiltonian/PeriodicAnisotropicSpinChainHamiltonian.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "Operator/PeriodicAnisotropicMagnetizationOperator.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;


void RotateCouplingConstants (double phi, double psi, double theta, double jx, double jy, double jz,
			      double& jxx, double& jyy, double& jzz, double& jxy, double& jxz, double& jyz);

int main(int argc, char** argv)
{
  cout.precision(14);
  int NbrSpin = 12;
  Spin1_2Chain Chain(NbrSpin, 1000000);
//    Spin1_2Chain Chain(NbrSpin, 0, 1000000);
  cout << Chain.GetHilbertSpaceDimension() << endl;
  double* CouplingConstants = new double [NbrSpin];
  double* CouplingConstantsX = new double [NbrSpin];
  for (int i = 0; i < NbrSpin; i++)
    {
      CouplingConstants[i] = 5.0;
      CouplingConstantsX[i] = 1.0;
    }
  double* CouplingConstantsJXX = new double [NbrSpin];
  double* CouplingConstantsJYY = new double [NbrSpin];
  double* CouplingConstantsJZZ = new double [NbrSpin];
  double* CouplingConstantsJXY = new double [NbrSpin];
  double* CouplingConstantsJXZ = new double [NbrSpin];
  double* CouplingConstantsJYZ = new double [NbrSpin];
  double* CouplingConstantsGX = new double [NbrSpin];
  double* CouplingConstantsGY = new double [NbrSpin];
  double* CouplingConstantsGZ = new double [NbrSpin];
  double B = 0.0;
  double BMax = 10.0;
  int NbrB = 10;
  if (argc >= 2)
    sscanf (argv[1],"%lf", &B);    
  if (argc >= 3)
    sscanf (argv[2],"%lf", &BMax);    
  if (argc >= 4)
    sscanf (argv[3],"%d", &NbrB);    
  double BInc;
  if (NbrB != 1)
    BInc = (BMax - B) / (double) (NbrB - 1);
  else
    BInc = 1.0;
  for (;B <= BMax; B += BInc)
    {
      cout << "B = " << B << endl;
      double BField = -0.6717 * B ;//(0.9274078)/(1.380662) * B; 
      double JX = 0.0;
      double JY = 0.0;
      double JZ = 200.0;
      for (int i = 0; i < NbrSpin;)
	{
/*	  CouplingConstantsJXX[i] = 0.0;
	  CouplingConstantsJYY[i] = 0.0;
	  CouplingConstantsJZZ[i] = 0.0;
	  CouplingConstantsJXY[i] = 0.0;
	  CouplingConstantsJXZ[i] = 0.0;
	  CouplingConstantsJYZ[i] = 0.0;
	  CouplingConstantsGX[i] = 0.0;
	  CouplingConstantsGY[i] = 1.0;
	  CouplingConstantsGZ[i] = 0.0;
	  ++i;*/
	  double Angle = 0.0;
	  double Psi = 0.0;
	  double Theta = 0.959931088;//0.955316618;
	  double Dummy = 0.0;
	  for (int j = 0; j < 6; j += 2)
	    {
	      RotateCouplingConstants(Angle, Psi, Theta, JX, JY, JZ, CouplingConstantsJXX[i], CouplingConstantsJYY[i], 
				      CouplingConstantsJZZ[i], CouplingConstantsJXY[i], 
				      CouplingConstantsJXZ[i], CouplingConstantsJYZ[i]);
	      //	  cout << CouplingConstantsJXX[i] << " " << CouplingConstantsJYY[i] << " " << CouplingConstantsJZZ[i] << " " 
	      //	       << CouplingConstantsJXY[i] << " " << CouplingConstantsJXZ[i] << " " << CouplingConstantsJYZ[i] << endl;
//	      RotateCouplingConstants(Angle, Psi, Theta, 2.0 * BField, 2.0 * BField, 2.0 * BField, Dummy, Dummy, 
//				      CouplingConstantsGZ[i], Dummy,
//				      CouplingConstantsGX[i], CouplingConstantsGY[i]);
	      CouplingConstantsGX[i] = 0.0;
	      CouplingConstantsGY[i] = 0.0;
	      CouplingConstantsGZ[i] = 2.0 * BField;
	      ++i;
	      
	      RotateCouplingConstants(Angle, Psi, Theta, JX, JY, JZ, CouplingConstantsJXX[i], CouplingConstantsJYY[i], 
				      CouplingConstantsJZZ[i], CouplingConstantsJXY[i], 
				      CouplingConstantsJXZ[i], CouplingConstantsJYZ[i]);
	      //	  cout << CouplingConstantsJXX[i] << " " << CouplingConstantsJYY[i] << " " << CouplingConstantsJZZ[i] << " " 
	      //	       << CouplingConstantsJXY[i] << " " << CouplingConstantsJXZ[i] << " " << CouplingConstantsJYZ[i] << endl;
	      RotateCouplingConstants(Angle, Psi, Theta, 0.0 * BField, 0.0 * BField, 9.0 * BField, Dummy, Dummy, 
				      CouplingConstantsGZ[i], Dummy,
				      CouplingConstantsGX[i], CouplingConstantsGY[i]);
	      
	      ++i;
	      Angle +=  2.0 * M_PI / 3.0;
	    }
	}
      PeriodicAnisotropicSpinChainHamiltonian H ((AbstractSpinChain*) &Chain, NbrSpin, CouplingConstantsJXX, CouplingConstantsJYY, 
						 CouplingConstantsJZZ, CouplingConstantsJXY, CouplingConstantsJXZ, CouplingConstantsJYZ, 
						 CouplingConstantsGX, CouplingConstantsGY, CouplingConstantsGZ);

      // Exact diagonalization
/*      RealSymmetricMatrix HRep (Chain.GetHilbertSpaceDimension());
      RealTriDiagonalSymmetricMatrix TmpTriDiag (Chain.GetHilbertSpaceDimension());
      H.GetHamiltonian(HRep);
      HRep.Householder(TmpTriDiag, 1e-14);
      TmpTriDiag.Diagonalize();
      TmpTriDiag.SortMatrixUpOrder();
      for (int j = 0; j < Chain.GetHilbertSpaceDimension(); j++)
	cout << TmpTriDiag.DiagonalElement(j) << " ";
      cout << endl;*/
//	  return 0;
      // Lanczos method
      //  BasicLanczosAlgorithm Lanczos;
//        cout << H << endl;
      //  ComplexBasicLanczosAlgorithm Lanczos;
//      AbstractArchitecture* Architecture = new MonoProcessorArchitecture;
      AbstractArchitecture* Architecture = new SMPArchitecture(2);
//      ComplexBasicLanczosAlgorithmWithEigenstates Lanczos(Architecture);
      ComplexBasicLanczosAlgorithmWithGroundState Lanczos(Architecture);
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
      cout << "Nbr of spins = " << NbrSpin << " " << Lowest << " " << Precision << "  Nbr of iterations = " << CurrentNbrIterLanczos << endl;
      ComplexVector& GroundState = (ComplexVector&) Lanczos.GetGroundState();
      PeriodicAnisotropicMagnetizationOperator Magnetization ((AbstractSpinChain*) &Chain, NbrSpin,  
							      CouplingConstantsGX, CouplingConstantsGY, CouplingConstantsGZ);
      cout << (Magnetization.MatrixElement(GroundState, GroundState) / (0.5 * NbrSpin * BField)) << endl;
      cout << "------------------------------------------------------------------" << endl << endl;;
      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
    }
  return 0;
}

void RotateCouplingConstants (double phi, double psi, double theta, double jx, double jy, double jz,
			      double& jxx, double& jyy, double& jzz, double& jxy, double& jxz, double& jyz)
{
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  double cospsi = cos(psi);
  double sinpsi = sin(psi);
  double costheta = cos(theta);
  double sintheta = sin(theta);
  double tmp;

  jxx = cosphi * cospsi * costheta - sinphi * sinpsi;
  jxy = jx * jxx;
  jxz = - jx * jxx;
  jxx *= jxx * jx;
  tmp = cosphi * sinpsi * costheta + sinphi * cospsi;
  jxx += jy * tmp * tmp;
  tmp = cosphi * sintheta;
  jxx += jz * tmp * tmp;
  
  jyy = sinphi * cospsi * costheta + cosphi * sinpsi;
  jxy *= jyy;
  jyz = - jx * jyy;
  jyy *= jyy * jx;
  tmp = cosphi * cospsi - sinphi * sinpsi * costheta;
  jyy += jy * tmp * tmp;
  tmp = sinphi * sintheta;
  jyy += jz * tmp * tmp;
  
  jzz = cosphi * sintheta;
  jzz *= jzz * jx;
  tmp = sinphi * sintheta;
  jzz += jy * tmp * tmp;
  jzz += jz * costheta * costheta;
  
  jxy += jy * (cosphi * sinpsi * costheta + sinphi * cospsi) * (sinphi * sinpsi * costheta - cosphi * cospsi)
    + jz * cosphi * sinphi * sintheta * sintheta;
  
  jxz *= cospsi * sintheta;
  jxz -= jy * (cosphi * sinpsi * costheta + sinphi * cospsi) * (sinpsi * sintheta)
    - jz * cosphi * sintheta * costheta;
  
  jyz *= (cospsi * sintheta);
  jyz += jy * (cosphi * cospsi - sinphi * sinpsi * costheta) * (sinpsi * sintheta)
    + jz * sinphi * sintheta * costheta;
}

