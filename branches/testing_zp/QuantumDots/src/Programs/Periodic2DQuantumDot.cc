#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Potential/BinaryTwoDConstantCellPotential.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

const double KineticFactor = 150.4;

double Weight_E[] = {0.128917, 0.119204, 0.103178, 0.0802035, 0.0568751, 0.0372012, 0.0217138, 0.0122897, 0.00653256, 0.00326567, 0.00153112, 0.000692242};
double Weight_H[] = {3.70785e-08, 3.38324e-07, 2.70859e-06, 1.90038e-05, 0.000116127, 0.000612045, 0.0030628, 0.0113668, 0.0345713, 0.0835946, 0.153434, 0.210873};

Complex Integral(int m, int n, int i, int P);

int main(int argc, char** argv)
{
/*
  //DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE)
  char** name = new char*[1];
  name[0] = "Periodic50.txt";
  int state[1]; state[0] = 10;
  DOSSpectra spectra(1, name, state, 5e-4, -0.011, 0.18, 1e-4);
  spectra.WriteSpectra("Periodic50");
*/

  cout.precision(14);

  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode");
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  SingleIntegerOption MValueOption ('M', "M-cell", "number of cells in the x direction", 120);
  SingleDoubleOption LeftXValueOption('\n', "left-x", "proportion of left cells omitted in x direction", 0.0);
  SingleDoubleOption RightXValueOption('\n', "right-x", "proportion of right cells omitted in x direction", 0.0);
  SingleIntegerOption NValueOption ('N', "N-cell", "number of cells in the y direction", 120);
  SingleDoubleOption LeftYValueOption('\n', "left-y", "proportion of left cells omitted in y direction", 0.0);
  SingleDoubleOption RightYValueOption('\n', "right-y", "proportion of right cells omitted in y direction", 0.0);
  SingleDoubleOption CellXSizeOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 2.97);
  SingleDoubleOption CellYSizeOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 2.97);
  SingleDoubleOption XMassOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.166);
  SingleDoubleOption YMassOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.166);  
  SingleStringOption CoefficientFileNameOption('\0', "coefficients", "name of the file where interaction coeffcients are stored");

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &VerboseOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &MValueOption;
  OptionList += &LeftXValueOption;
  OptionList += &RightXValueOption;
  OptionList += &NValueOption;
  OptionList += &LeftYValueOption;
  OptionList += &RightYValueOption;
  OptionList += &CellXSizeOption;
  OptionList += &CellYSizeOption;
  OptionList += &XMassOption;
  OptionList += &YMassOption;
  OptionList += &CoefficientFileNameOption;

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
  bool VerboseFlag = VerboseOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrCellX = MValueOption.GetInteger();
  double LeftX = LeftXValueOption.GetDouble();
  double RightX = RightXValueOption.GetDouble();
  int NbrCellY = NValueOption.GetInteger();
  double LeftY = LeftYValueOption.GetDouble();
  double RightY = RightYValueOption.GetDouble();
  double Lx = CellXSizeOption.GetDouble();
  double Ly = CellYSizeOption.GetDouble();
  double Mux = XMassOption.GetDouble();
  double Muy = YMassOption.GetDouble();
  char* CoefficientFileName = CoefficientFileNameOption.GetString();

  int AbsoluteWaveVectorX = NbrCellX / 2;
  int NbrStateX = AbsoluteWaveVectorX * 2 + 1;
  int AbsoluteWaveVectorY = NbrCellY / 2;
  int NbrStateY = AbsoluteWaveVectorY * 2 + 1;
  int Dimension = NbrStateX * NbrStateY;
  int MinX = LeftX * NbrCellX; int MaxX = (1.0 - RightX) * NbrCellX;
  int MinY = LeftY * NbrCellY; int MaxY = (1.0 - RightY) * NbrCellY;

  BinaryTwoDConstantCellPotential potential = BinaryTwoDConstantCellPotential(NbrCellX, NbrCellY);

  // coupling table in X direction
  ComplexVector** ElementaryCouplingX = new ComplexVector* [NbrStateX];
  Complex integral = Complex(0.0, 0.0);
  for (int m1 = 0; m1 < NbrStateX; ++m1)
    {
      ElementaryCouplingX[m1] = new ComplexVector[NbrStateX];
      for (int m2 = 0; m2 < m1; ++m2)
	{
	  ElementaryCouplingX[m1][m2] = ComplexVector(NbrCellX);
	  ElementaryCouplingX[m2][m1] = ComplexVector(NbrCellX);
	  for (int i = 0; i < NbrCellX; ++i)
	    {
	      integral = Integral(m1 - AbsoluteWaveVectorX, m2 - AbsoluteWaveVectorX, i, NbrCellX);
	      ElementaryCouplingX[m1][m2].Re(i) = integral.Re;
	      ElementaryCouplingX[m1][m2].Im(i) = integral.Im; 
	      ElementaryCouplingX[m2][m1].Re(i) = integral.Re;
	      ElementaryCouplingX[m2][m1].Im(i) = -integral.Im; 
	    }
	}
      ElementaryCouplingX[m1][m1] = ComplexVector(NbrCellX);
      for (int i = 0; i < NbrCellX; ++i)
	{
	  integral = Integral(m1 - AbsoluteWaveVectorX, m1 - AbsoluteWaveVectorX, i, NbrCellX);
	  ElementaryCouplingX[m1][m1].Re(i) = integral.Re;
	  ElementaryCouplingX[m1][m1].Im(i) = integral.Im; 
	}    
    }

  // coupling table in Y direction
  ComplexVector** ElementaryCouplingY = new ComplexVector* [NbrStateY];
  for (int n1 = 0; n1 < NbrStateY; ++n1)
    {
      ElementaryCouplingY[n1] = new ComplexVector[NbrStateY];
      for (int n2 = 0; n2 < n1; ++n2)
	{
	  ElementaryCouplingY[n1][n2] = ComplexVector(NbrCellY);
	  ElementaryCouplingY[n2][n1] = ComplexVector(NbrCellY);
	  for (int j = 0; j < NbrCellY; ++j)
	    {
	      integral = Integral(n1 - AbsoluteWaveVectorY, n2 - AbsoluteWaveVectorY, j, NbrCellY);
	      ElementaryCouplingY[n1][n2].Re(j) = integral.Re; 
	      ElementaryCouplingY[n1][n2].Im(j) = integral.Im; 
	      ElementaryCouplingY[n2][n1].Re(j) = integral.Re; 
	      ElementaryCouplingY[n2][n1].Im(j) = -integral.Im; 
	    }
	}
      ElementaryCouplingY[n1][n1] = ComplexVector(NbrCellY);
      for (int j = 0; j < NbrCellY; ++j)
	{
	  integral = Integral(n1 - AbsoluteWaveVectorY, n1 - AbsoluteWaveVectorY, j, NbrCellY);
	  ElementaryCouplingY[n1][n1].Re(j) = integral.Re; 
	  ElementaryCouplingY[n1][n1].Im(j) = integral.Im; 
	}
    }     

  RealVector KineticX (NbrStateX);
  RealVector KineticY (NbrStateY);
  double TmpKineticX = KineticFactor / (Mux * Lx * Lx * NbrCellX * NbrCellX);
  double TmpKineticY = KineticFactor / (Muy * Ly * Ly * NbrCellY * NbrCellY); 
  for (int i = 0; i < NbrStateX; ++i)     
    KineticX[i] = TmpKineticX * (i - AbsoluteWaveVectorX) * (i - AbsoluteWaveVectorX);

  for (int j = 0; j < NbrStateY; ++j)
    KineticY[j] = TmpKineticY * (j - AbsoluteWaveVectorY) * (j - AbsoluteWaveVectorY);

  // intermediate table to accelerate the calculation of hamiltonian matrix
  ComplexVector** Intermediate = new ComplexVector* [NbrStateX];
  Complex tmp = Complex(0.0, 0.0);
  for (int m1 = 0; m1 < NbrStateX; ++m1)
    {
      Intermediate[m1] = new ComplexVector[NbrStateX];
      for (int m2 = 0; m2 < m1; ++m2)
	{
	  Intermediate[m1][m2] = ComplexVector(NbrCellY);	  
	  Intermediate[m2][m1] = ComplexVector(NbrCellY);
	  for (int j = MinY; j < MaxY; ++j)
	    {
	      tmp = Complex(0.0, 0.0);
	      for (int i = MinX; i < MaxX; ++i)		
		tmp += (potential.GetPotential(i, j) * ElementaryCouplingX[m1][m2][i]);	      		
	      Intermediate[m1][m2].Re(j) = tmp.Re;
	      Intermediate[m1][m2].Im(j) = tmp.Im;
	      Intermediate[m2][m1].Re(j) = tmp.Re;
	      Intermediate[m2][m1].Im(j) = -tmp.Im;	      
	    }
	}
      Intermediate[m1][m1] = ComplexVector(NbrCellY);
      for (int j = MinY; j < MaxY; ++j)
	{
	  tmp = Complex(0.0, 0.0);
	  for (int i = MinX; i < MaxX; ++i)
	    tmp += (potential.GetPotential(i, j) * ElementaryCouplingX[m1][m1][i]);
	  Intermediate[m1][m1].Re(j) = tmp.Re;
	  Intermediate[m1][m1].Im(j) = tmp.Im;
	}
    }

  double ShiftHamiltonian = KineticX[0] + KineticY[0];

  HermitianMatrix HamiltonianRepresentation (Dimension);
  // fill the hamiltonian matrix
  double Re = 0.0;
  double Im = 0.0;
  ComplexVector tmpCoupling (NbrCellY);
  ComplexVector tmpIntermediate (NbrCellY);
  for (int n1 = 0; n1 < NbrStateY; ++n1)
    {
      for (int n2 = 0; n2 < n1; ++n2)
	{
	  tmpCoupling = ComplexVector(ElementaryCouplingY[n2][n1]);
	  for (int m1 = 0; m1 < NbrStateX; ++m1)
	    {
	      for (int m2 = 0; m2 < NbrStateX; ++m2)
		{
		  Re = 0.0; Im = 0.0; 
		  tmpIntermediate = ComplexVector(Intermediate[m2][m1]);
		  for (int j = MinY; j < MaxY; ++j)
		    {
		      Re += (tmpCoupling.Re(j) * tmpIntermediate.Re(j) - tmpCoupling.Im(j) * tmpIntermediate.Im(j));
		      Im += (tmpCoupling.Im(j) * tmpIntermediate.Re(j) + tmpCoupling.Re(j) * tmpIntermediate.Im(j));
		    }
		  HamiltonianRepresentation.SetMatrixElement(m2 + n2 * NbrStateX, m1 + n1 * NbrStateX, Complex(Re, Im));
		}
	    }
	}
      tmpCoupling = ComplexVector(ElementaryCouplingY[n1][n1]);
      for (int m1 = 0; m1 < NbrStateX; ++m1)
	{
	  for (int m2 = 0; m2 < m1; ++m2)
	    {
	      Re = 0.0; Im = 0.0; 
	      tmpIntermediate = ComplexVector(Intermediate[m2][m1]);
	      for (int j = MinY; j < MaxY; ++j)
		{
		  Re += (tmpCoupling.Re(j) * tmpIntermediate.Re(j) - tmpCoupling.Im(j) * tmpIntermediate.Im(j));
		  Im += (tmpCoupling.Im(j) * tmpIntermediate.Re(j) + tmpCoupling.Re(j) * tmpIntermediate.Im(j));
		}
	      HamiltonianRepresentation.SetMatrixElement(m2 + n1 * NbrStateX, m1 + n1 * NbrStateX, Complex(Re, Im));	   
	    }
	  Re = KineticX[m1] + KineticY[n1] - ShiftHamiltonian;	  
	  tmpIntermediate = ComplexVector(Intermediate[m1][m1]);

	  for (int j = MinY; j < MaxY; ++j)
	    Re += (tmpCoupling.Re(j) * tmpIntermediate.Re(j) - tmpCoupling.Im(j) * tmpIntermediate.Im(j));
	  HamiltonianRepresentation.SetMatrixElement(m1 + n1 * NbrStateX, m1 + n1 * NbrStateX, Complex(Re, 0.0));      
	}
    }


  delete[] ElementaryCouplingX; delete[] ElementaryCouplingY; delete[] Intermediate;

  //cout << HamiltonianRepresentation << endl;

  // RealVector* Eigenstates = 0;
  double* Eigenvalues = new double [NbrEigenvalue];

  // find the eigenvalues (and eigenvectors if needed)

  // architecture type (i.e. 1 CPU or multi CPU)
  AbstractArchitecture* Architecture;
  if (SMPFlag == true)
    Architecture = new SMPArchitecture(2);
  else
    Architecture = new MonoProcessorArchitecture;

  // type of lanczos algorithm (with or without reorthogonalization)
  // BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
  FullReorthogonalizedComplexLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
  // store matrix in an hamiltonian class and create corresponding hilbert space

  UndescribedHilbertSpace HilbertSpace(Dimension);
  ExplicitHamiltonian Hamiltonian(&HilbertSpace, &HamiltonianRepresentation);

  // initialization of lanczos algorithm
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
  Lanczos.SetHamiltonian(&Hamiltonian);
  Lanczos.InitializeLanczosAlgorithm();
  cout << "Run Lanczos Algorithm" << endl;
  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
  RealTriDiagonalSymmetricMatrix TmpMatrix;
  
  // run Lancos algorithm up to desired precision on the n-th eigenvalues
  while (Lanczos.TestConvergence() == false)      
    {
      ++CurrentNbrIterLanczos;
      Lanczos.RunLanczosAlgorithm(1);
      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest;
      //cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
    }
  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
    {
      cout << "too much Lanczos iterations" << endl;
      exit(0);
    }

  // store eigenvalues
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      Eigenvalues[i] = TmpMatrix.DiagonalElement(i);
      cout << (Eigenvalues[i] + ShiftHamiltonian) << '\t';
    }
  cout << endl;


  // compute eigenstates
  // ComplexVector*  Eigenstates = (ComplexVector*) Lanczos.GetEigenstates(NbrEigenvalue);
  return 0;
}

Complex Integral(int m, int n, int i, int P)
{
  if (m == n)
    return Complex(1.0 / double(P), 0.0);
  else
    {
      double tmp1 = 1.0 / (2.0 * M_PI * (n - m));
      double tmp2 = 1.0 / (tmp1 * P);
      return Complex(tmp1 * (sin(tmp2 * (i + 1)) - sin(tmp2 * i)), tmp1 * (cos(tmp2 * i) - cos(tmp2 * (i + 1))));
    }
}
