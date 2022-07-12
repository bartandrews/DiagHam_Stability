#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/RealVector.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
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

const double KineticFactor = 37.60;//h bar carre reduit

//calculer l'integral des fonctions d'indice m et n
inline double integral(int m, int n, int i, int P);

double Weight_E[] = {0.128917, 0.119204, 0.103178, 0.0802035, 0.0568751, 0.0372012, 0.0217138, 0.0122897, 0.00653256, 0.00326567, 0.00153112, 0.000692242};
double Weight_H[] = {3.70785e-08, 3.38324e-07, 2.70859e-06, 1.90038e-05, 0.000116127, 0.000612045, 0.0030628, 0.0113668, 0.0345713, 0.0835946, 0.153434, 0.210873};

int main(int argc, char** argv)
{
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

  int NbrStateX = NbrCellX; int NbrStateY = NbrCellY;
  int MinX = LeftX * NbrCellX; int MaxX = (1.0 - RightX) * NbrCellX;
  int MinY = LeftY * NbrCellY; int MaxY = (1.0 - RightY) * NbrCellY;

  // enter here hilbert space dimension
  int Dimension = NbrStateX * NbrStateY;

  BinaryTwoDConstantCellPotential Potential = BinaryTwoDConstantCellPotential(NbrCellX, NbrCellY);
  //Potential.UniformWell(0.5, 1.8, Weight_E, true);
  //ofstream PotentialFile("Potential2D.txt");
  //Potential.PrintPotential(PotentialFile).flush();

  Potential.LoadPotential(CoefficientFileName);

  // coupling table in X direction
  RealVector** ElementaryCouplingX = new RealVector* [NbrStateX];
  for (int m1 = 0; m1 < NbrStateX; ++m1)
    {
      ElementaryCouplingX[m1] = new RealVector[NbrStateX];
      for (int m2 = 0; m2 < NbrStateX; ++m2)
	{
	  ElementaryCouplingX[m1][m2] = RealVector(NbrCellX);
	  for (int i = 0; i < NbrCellX; ++i)
	    ElementaryCouplingX[m1][m2][i] = integral(m1 + 1, m2 + 1, i + 1, NbrCellX);
	}
    }

  // coupling table in Y direction
  RealVector** ElementaryCouplingY = new RealVector* [NbrStateY];
  for (int n1 = 0; n1 < NbrStateY; ++n1)
    {
      ElementaryCouplingY[n1] = new RealVector[NbrStateY];
      for (int n2 = 0; n2 < NbrStateY; ++n2)
	{
	  ElementaryCouplingY[n1][n2] = RealVector(NbrCellY);
	  for (int j = 0; j < NbrCellY; ++j)
	    ElementaryCouplingY[n1][n2][j] = integral(n1 + 1, n2 + 1, j + 1, NbrCellY);
	}
    }

  RealVector KineticX (NbrStateX);
  RealVector KineticY (NbrStateY);
  double TmpKineticX = KineticFactor / (Mux * Lx * Lx * NbrCellX * NbrCellX);
  double TmpKineticY = KineticFactor / (Muy * Ly * Ly * NbrCellY * NbrCellY);
  for (int i = 0; i < NbrStateX; ++i)
    KineticX[i] = TmpKineticX * (i + 1) * (i + 1);

  for (int j = 0; j < NbrStateY; ++j)
    KineticY[j] = TmpKineticY * (j + 1) * (j + 1);

  RealVector** Intermediate = new RealVector* [NbrStateX];
  double tmp = 0.0;
  for (int m1 = 0; m1 < NbrStateX; ++m1)
    {
      Intermediate[m1] = new RealVector[NbrStateX];
      for (int m2 = 0; m2 < m1; ++m2)
	{
	  Intermediate[m1][m2] = RealVector(NbrCellY);
	  Intermediate[m2][m1] = RealVector(NbrCellY);
	  for (int j = MinY; j < MaxY; ++j)
	    {
	      tmp = 0.0;
	      for (int i = MinX; i < MaxX; ++i)
		tmp += (Potential.GetPotential(i, j) * ElementaryCouplingX[m1][m2][i]);
	      Intermediate[m1][m2][j] = tmp;
	      Intermediate[m2][m1][j] = tmp;
	    }
	}
      Intermediate[m1][m1] = RealVector(NbrCellY);
      for (int j = MinY; j < MaxY; ++j)
	{
	  tmp = 0.0;
	  for (int i = MinX; i < MaxX; ++i)
	    tmp += (Potential.GetPotential(i, j) * ElementaryCouplingX[m1][m1][i]);
	  Intermediate[m1][m1][j] = tmp;
	}
    }

  RealSymmetricMatrix HamiltonianRepresentation (Dimension);

  double ShiftHamiltonian = KineticX[NbrStateX - 1] + KineticY[NbrStateY - 1];

  for (int n1 = 0; n1 < NbrStateY; ++n1)
    for (int m1 = 0; m1 < NbrStateX; ++m1)
      {
	for (int n2 = 0; n2 < n1; ++n2)
	  {
	    for (int m2 = 0; m2 < NbrStateX; ++m2)
	      {
		tmp = 0.0;
		for (int j = MinY; j < MaxY; ++j)
		    tmp += (ElementaryCouplingY[n2][n1][j] * Intermediate[m2][m1][j]);
	 	HamiltonianRepresentation.SetMatrixElement(m2 + n2 * NbrStateX, m1 + n1 * NbrStateX, tmp);
	      }
	  }
	for (int m2 = 0; m2 < m1; ++m2)
	  {
	    tmp = 0.0;
	    for (int j = MinY; j < MaxY; ++j)
	      tmp += (ElementaryCouplingY[n1][n1][j] * Intermediate[m2][m1][j]);
	    HamiltonianRepresentation.SetMatrixElement(m2 + n1 * NbrStateX, m1 + n1 * NbrStateX, tmp);
	  }
	tmp = KineticX[m1] + KineticY[n1] - ShiftHamiltonian;
	for (int j = MinY; j < MaxY; ++j)
	  tmp += (ElementaryCouplingY[n1][n1][j] * Intermediate[m1][m1][j]);
	HamiltonianRepresentation.SetMatrixElement(m1 + n1 * NbrStateX, m1 + n1 * NbrStateX, tmp);
      }

  delete[] ElementaryCouplingX; delete[] ElementaryCouplingY; delete[] Intermediate;

  // RealVector* Eigenstates = 0;
  double* Eigenvalues = new double [NbrEigenvalue];

  // find the eigenvalues (and eigenvectors if needed)
  
  // architecture type (i.e. 1 CPU or multi CPU)
  AbstractArchitecture* Architecture;
  if (SMPFlag == true)
    Architecture = new SMPArchitecture(2);
  else
    Architecture = new MonoProcessorArchitecture;
  
  double Precision;
  double PreviousLowest;
  double Lowest;
  int CurrentNbrIterLanczos;

  // type of lanczos algorithm (with or without reorthogonalization)
  // BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
  FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);      
  // store matrix in an hamiltonian class and create corresponding hilbert space
  cout << "End of filling the matrix" << endl;      
  UndescribedHilbertSpace HilbertSpace(Dimension);
  ExplicitHamiltonian Hamiltonian(&HilbertSpace, &HamiltonianRepresentation);
  
  
  // initialization of lanczos algorithm
  Precision = 1.0;
  PreviousLowest = 1e50;
  Lowest = PreviousLowest;
  CurrentNbrIterLanczos = NbrEigenvalue + 3;
  Lanczos.SetHamiltonian(&Hamiltonian);
  Lanczos.InitializeLanczosAlgorithm();
  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
  RealTriDiagonalSymmetricMatrix TmpMatrix;
  
  // run Lancos algorithm up to desired precision on the n-th eigenvalues
  while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
    {
      Lanczos.RunLanczosAlgorithm(1);
      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
      if (VerboseFlag == true)
	cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << endl;
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest; 
    }      
  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
    {
      cout << "too much Lanczos iterations" << endl;
      exit(0);
    }
  
  // store eigenvalues      
  for (int i = 0; i < NbrEigenvalue; ++i)    
    {
      Eigenvalues[i] = TmpMatrix.DiagonalElement(i) + ShiftHamiltonian;
      cout << Eigenvalues[i] << '\t';
    }
    
  cout << endl;
      
  //compute eigenstates
  //RealVector*  Eigenstates = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);      
  
  return 0;
}


//calculer l'integral des fonctions d'indice m et n
inline double integral(int m, int n, int i, int P){
  if (m != n)
    return M_1_PI * (((sin(i * (m - n) * M_PI / P) - sin((i - 1) * (m - n) * M_PI / P)) / (m - n)) - ((sin(i * (m + n) * M_PI / P) - sin((i - 1) * (m + n) * M_PI / P)) / (m + n)));
  else
    return (1.0 / P) * (1 + (P / (2 * M_PI * n)) * (sin(2 * M_PI * n * (i - 1) / P) - sin(2 * M_PI * n * i / P)));
}
