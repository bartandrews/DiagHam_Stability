#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ofstream;

bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);

#define PERIODIC_HAMILTONIAN_FACTOR 150.4

int main(int argc, char** argv)
{  
  cout.precision(14);

  OptionManager Manager ("PeriodicContinuum" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += MiscGroup;

  (*PotentialGroup) += new SingleIntegerOption ('H', "H-cell", "number of cells in the z direction", 21);
  (*PotentialGroup) += new SingleDoubleOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 5.65);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "barrier", "number of cells in the well barrier", 2);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "below", "number of cells between well barrier and wetting layer", 2);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "wetting", "number of cells in wetting layer", 1);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "well", "potential in the well", 1.079);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "dot", "potential in the dot", -0.4);

  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  (*MiscGroup) += new SingleIntegerOption ('n', "nbr-eigen", "number of eigenvalues", 2);
  (*MiscGroup) += new BooleanOption ('e', "eigenstate", "evaluate eigenstates", false);
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int H = ((SingleIntegerOption*) Manager["H-cell"])->GetInteger();
  double Lz = ((SingleDoubleOption*) Manager["cell-zsize"])->GetDouble();
  int UnderBarrier = ((SingleIntegerOption*) Manager["barrier"])->GetInteger();
  int BelowWettingLayer = ((SingleIntegerOption*) Manager["below"])->GetInteger();
  int WettingWidth = ((SingleIntegerOption*) Manager["wetting"])->GetInteger();
  double WellPotential = ((SingleDoubleOption*) Manager["well"])->GetDouble();
  double DotPotential = ((SingleDoubleOption*) Manager["dot"])->GetDouble();

  double Muz = ((SingleDoubleOption*) Manager["mu-z"])->GetDouble();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();

  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();   
  bool EigenstateFlag = ((BooleanOption*) Manager["eigenstate"])->GetBoolean();

  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  
  ofstream Command; Command.open ("Command.txt", ios::out | ios::app);
  Command << "\n============================ Begin =============================" << '\n';
  Command << "The program was launched at: " <<  asctime (timeinfo) << '\n';
  Manager.DisplayOption (Command, true);
  Command << endl;

  ofstream FullOption; FullOption.open ("FullOption.txt", ios::out | ios::app);
  FullOption << "\n============================ Begin =============================" << '\n';
  FullOption << "The program was launched at: " << asctime (timeinfo) << '\n';
  Manager.DisplayOption (FullOption, false);

  double SizeZ = Lz * H;
  double** RealOverlap; double** ImaginaryOverlap;
  EvaluateWaveFunctionOverlap(H, NbrStateZ, RealOverlap, ImaginaryOverlap);  

  int Dimension = NbrStateZ;
  HermitianMatrix HamiltonianRepresentation(Dimension);

  double TmpRe = 0.0, TmpIm = 0.0; Complex TmpComplex;
  int OriginZ = NbrStateZ - 1;
  for (int n1 = 0; n1 < Dimension; ++n1)
    {
      for (int n2 = 0; n2 <= n1; ++n2)
	{	 
	  TmpRe = 0.0; TmpIm = 0.0;
	  for (int k = 0; k < UnderBarrier; ++k)
	    {
	      TmpRe += (WellPotential * RealOverlap[-n1 + n2 + OriginZ][k]);
	      TmpIm += (WellPotential * ImaginaryOverlap[-n1 + n2 + OriginZ][k]);
	    }
	  for (int k = UnderBarrier + BelowWettingLayer; k < (UnderBarrier + BelowWettingLayer + WettingWidth); ++k)
	    {
	      TmpRe += (DotPotential * RealOverlap[-n1 + n2 + OriginZ][k]);
	      TmpIm += (DotPotential * ImaginaryOverlap[-n1 + n2 + OriginZ][k]);
	    }	  
	  TmpComplex = Complex(TmpRe, TmpIm);	  
	  HamiltonianRepresentation.SetMatrixElement(n1, n2, TmpComplex);
	}
      HamiltonianRepresentation.AddToMatrixElement(n1, n1, PERIODIC_HAMILTONIAN_FACTOR * double((n1 + LowImpulsionZ) * (n1 + LowImpulsionZ)) / (Muz * SizeZ * SizeZ));
    }

  RealMatrix TmpEigenvectors (Dimension, Dimension);
  RealSymmetricMatrix RealHamiltonianRepresentation = HamiltonianRepresentation.ConvertToSymmetricMatrix();
  RealTriDiagonalSymmetricMatrix TmpTriDiag (2 * Dimension);
  RealHamiltonianRepresentation.Householder(TmpTriDiag, MACHINE_PRECISION, TmpEigenvectors);
  TmpTriDiag.Diagonalize(TmpEigenvectors);
  TmpTriDiag.SortMatrixUpOrder(TmpEigenvectors);
  for (int j = 0; j < Dimension; j++)
    cout << TmpTriDiag.DiagonalElement(2 * j) << " ";
  cout << endl;
  if (EigenstateFlag)
    {
      RealVector* Eigenstates = new RealVector [NbrEigenvalue];
      double* Eigenvalues = new double [NbrEigenvalue];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  //Eigenstates[i] = TmpEigenvectors[i];
	  Eigenstates[i] = RealVector(2 * Dimension);
	  for (int j = 0; j < Dimension; ++j)
	    {
	      Eigenstates[i][2 * j] = TmpEigenvectors[2 * i][j];
	      Eigenstates[i][2 * j + 1] = TmpEigenvectors[2 * i][j + Dimension];
	    }
	  Eigenvalues[i] = TmpTriDiag.DiagonalElement(2 * i);
	}   
      ofstream OutputFile;
      OutputFile.precision(14);
      OutputFile.open("eigenvalues", ios::binary | ios::out);
      for (int i = 0; i < NbrEigenvalue; ++i)
	OutputFile << Eigenvalues[i] << " ";
      OutputFile << endl;
      OutputFile.close();

      char* TmpFileName = new char[256];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  sprintf  (TmpFileName, "eigenvector.%d", i);
	  Eigenstates[i].WriteAsciiVector(TmpFileName);
	}
      delete[] TmpFileName;
      // test the eigenvectors
      /*
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  ComplexVector TmpVector (Dimension);
	  ComplexVector TmpEigen (Dimension);
	  for (int j = 0; j < Dimension; ++j)
	    {
	      TmpEigen.Re(j) = TmpEigenvectors[2 * i][j];
	      TmpEigen.Im(j) = TmpEigenvectors[2 * i][j + Dimension];
	    }
	  TmpVector.Multiply(HamiltonianRepresentation, TmpEigen);	  
	  cout << (TmpVector * (ComplexVector&) TmpEigen) << " ";
	}
      */
      cout << endl;
    }

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Command << "The program finished at: " <<  asctime (timeinfo);
  Command << "============================== End =============================" << '\n';
  FullOption << "The program finished at: " << asctime (timeinfo);
  FullOption << "=============================== End ============================" << '\n'; 
  Command.close(); FullOption.close(); 
  
  return 1;
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);
  int Length = (nbrState - 1) * 2 + 1;
  realArray = new double* [Length];
  imaginaryArray = new double* [Length];  
  int Origin = nbrState - 1;
  for (int delta = 0; delta < Length; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      if (delta != Origin)
	{
	  Diff = 2.0 * M_PI * double (delta - Origin);
	  Tmp = Diff / nbrStep;	
	  Diff = 1.0 / Diff;	
	  for (int i = 0; i < nbrStep; ++i)
	    {
	      realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	      imaginaryArray[delta][i] = Diff * (cos(Tmp * (i + 1)) - cos(Tmp * i));
	    }
	}
      else
	for (int i = 0; i < nbrStep; ++i)
	  {
	    realArray[delta][i] = Tmp1;
	    imaginaryArray[delta][i] = 0.0;
	  }	
    }
  return true;
}
