#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/Confined3DOneParticle.h"

#include "Hamiltonian/QuantumDots3DHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


//Definir des constantes
int M = 14;//nombre de cellules dans l'axe x
int N = 14;//nombre de cellules dans l'axe y
int H = 10;//nombre de mailles en z
const int InNz1 = 1;//interface gauche du puits ou de la boite
const int InNz2 = 2;//interface droite du puits ou de la boite
const double a = 2.97;//parametre de maille en x
const double b = 2.97;//parametre de maille en y
const double c = 2.64;//parametre de maille en z
const double Me = 0.14;//masse effective
const double h = 37.60;//h bar carre reduit

//profil de potentiel
double*** Potentiel;

//1er indice: cellule d'InN; 2e et 3e: indices de fonction
double*** TransfertX;
double*** TransfertY;
double*** TransfertZ;

//calculer l'energie cinetique
inline double cinetique(int i, int j, int k);
//entrer les elements
inline void Entree();

// evaluate matrix element
//
double EvaluateMatrixElement(int i, int j);

// read interaction coefficients in a file
//
// fileName = string containing file name
// return value = tridimensionnal array where interaction coefficients are stored (null pointer if an error occurs)
double*** ReadInteractionCoefficients(char* fileName, int nbrCellX, int nbrCellY, int nbrCellZ);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  BooleanOption LanczosOption ('l', "lanczos", "enable lanczos diagonalization algorithm", false);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode");
  BooleanOption GeneratePotentialOption ('p', "potential", "generate potential instead of reading it from a file", false);
  SingleDoubleOption ElectricFieldOption ('\n', "electric-field", "electric field strength in the quantum dot");
  SingleDoubleOption InitialPotentialOption ('\n', "initial-potential", "potential value in the leftmost part of the quantum dot");
  SingleDoubleOption FinalPotentialOption ('\n', "final-potential", "potential value at the right hand side of the quantum dot");
  SingleDoubleOption ImpurityDensityOption ('\n', "impurity-density", "impurity density in the quantum dot");
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", false);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  SingleIntegerOption MemoryOption ('\n', "memory", "amount of memory that can be used for precaching (in Mb)", 500);
  SingleIntegerOption MValueOption ('M', "M-cell", "number of cells in the x direction", 10);
  SingleIntegerOption NValueOption ('N', "N-cell", "number of cells in the y direction", 10);
  SingleIntegerOption HValueOption ('H', "H-cell", "number of cells in the z direction", 10);
  SingleIntegerOption LeftSizeOption ('\n', "left-size", "size of the leftmost part in the z direction with constant null potential (in cell unit)", 0);
  SingleIntegerOption RightSizeOption ('\n', "right-size", "size of the rightmost part in the z direction with constant potential (in cell unit)", 0);
  SingleDoubleOption RightPotentialOption ('\n', "right-potential", "potential value in the rightmost part in the z direction", 0);
  SingleDoubleOption CellXSizeOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 2.97);
  SingleDoubleOption CellYSizeOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 2.97);
  SingleDoubleOption CellZSizeOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 2.64);  
  SingleDoubleOption XMassOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.14);
  SingleDoubleOption YMassOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.14);
  SingleDoubleOption ZMassOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.14);
  SingleStringOption CoefficientFileNameOption('\0', "coefficients", "name of the file where interaction coeffcients are stored", 
					       "/home/regnault/development/DMRG/DiagHam/potentiel_10_10_10_2");
  List<AbstractOption*> OptionList;
  OptionList += &LanczosOption;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &VerboseOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &EigenstateOption;
  OptionList += &MValueOption;
  OptionList += &NValueOption;
  OptionList += &HValueOption;
  OptionList += &CellXSizeOption;
  OptionList += &CellYSizeOption;
  OptionList += &CellZSizeOption;
  OptionList += &XMassOption;
  OptionList += &YMassOption;
  OptionList += &ZMassOption;
  OptionList += &CoefficientFileNameOption;
  OptionList += &GeneratePotentialOption;
  OptionList += &ElectricFieldOption;
  OptionList += &InitialPotentialOption;
  OptionList += &FinalPotentialOption;
  OptionList += &ImpurityDensityOption;
  OptionList += &LeftSizeOption;
  OptionList += &RightSizeOption;
  OptionList += &RightPotentialOption;
  OptionList += &MemoryOption;
  
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

  int Memory = MemoryOption.GetInteger();
  bool LanczosFlag = LanczosOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  bool VerboseFlag = VerboseOption.GetBoolean();
  bool EigenstateFlag = EigenstateOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  char* CoefficientFileName = CoefficientFileNameOption.GetString(); 
  M = MValueOption.GetInteger();
  N = NValueOption.GetInteger();
  H = HValueOption.GetInteger();
  double Lx = CellXSizeOption.GetDouble();
  double Ly = CellYSizeOption.GetDouble();
  double Lz = CellZSizeOption.GetDouble();
  double Mux = XMassOption.GetDouble();
  double Muy = YMassOption.GetDouble();
  double Muz = ZMassOption.GetDouble();
  int LeftSize = LeftSizeOption.GetInteger();
  int RightSize = RightSizeOption.GetInteger();
  double RightPotential = RightPotentialOption.GetDouble();

  // initialize matrix associated to the hamiltonian
  double*** InteractionCoefficients = ReadInteractionCoefficients(CoefficientFileName, M, N, H - LeftSize - RightSize);
  
  Confined3DOneParticle Space(M, N, H);
//  Confined3DOneParticle Space(3, 3, 3);
  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);
  QuantumDots3DHamiltonian Hamiltonian(&Space, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), ((double) LeftSize) * Lz, 
				       ((double) RightSize) * Lz, RightPotential,
				       Mux, Muy, Muz, M, N, H - LeftSize - RightSize, InteractionCoefficients, Memory << 20);
  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) + 
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "precalculation time = " << Dt << endl;

  RealVector* Eigenstates = 0;
  double* Eigenvalues = 0;
  if ((LanczosFlag == false) || (Space.GetHilbertSpaceDimension() < 300))
    {
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  // diagonalize the hamiltonian
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
/*	  for (int i = 0; i < Hamiltonian.GetHilbertSpaceDimension(); ++i)
	    cout << HRep(i, i) << " ";
	  cout << endl;*/
//	  cout  << HRep << endl;
//	  cout << "dim = " << Hamiltonian.GetHilbertSpaceDimension() << endl;
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Space.GetHilbertSpaceDimension());

	  cout << "start diagonalization..." << endl;
	  if (EigenstateFlag == false)
	    {
	      HRep.Householder(TmpTriDiag, MACHINE_PRECISION);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	    }
	  else
	    {
	      RealMatrix TmpEigenvectors (Space.GetHilbertSpaceDimension(), Space.GetHilbertSpaceDimension());
	      ((RealSymmetricMatrix*) Hamiltonian.GetHamiltonian())->Householder(TmpTriDiag, MACHINE_PRECISION, TmpEigenvectors);
	      TmpTriDiag.Diagonalize(TmpEigenvectors);
	      TmpTriDiag.SortMatrixUpOrder(TmpEigenvectors);
	      Eigenstates = new RealVector [NbrEigenvalue];
	      for (int i = 0; i < NbrEigenvalue; ++i)
		{
		  Eigenstates[i] = TmpEigenvectors[i];
		}
	    }

	  // store eigenvalues
	  int Max = Hamiltonian.GetHilbertSpaceDimension();
	  NbrEigenvalue = Max;
	  if (Max > NbrEigenvalue)
	    Max = NbrEigenvalue;
	  Eigenvalues = new double [NbrEigenvalue];
	  for (int j = 0; j < Max ; j++)
	    {
	      Eigenvalues[j]= TmpTriDiag.DiagonalElement(j);
	    }
	}
      else
	{
	  cout << (*(Hamiltonian.GetHamiltonian()))(0, 0) << endl;
	}
    }
  else
    {
      
      // architecture type (i.e. 1 CPU or multi CPU)
      AbstractArchitecture* Architecture;
      if (SMPFlag == true)
	Architecture = new SMPArchitecture(2);
      else
	Architecture = new MonoProcessorArchitecture;

      double HamiltonianShift = - (37.60 * ((1.0 / (Lx * Lx * Mux)) + (1.0 / (Ly * Ly * Muy)) + (1.0 / (Lz * Lz * Muz))));
      Hamiltonian.ShiftHamiltonian (HamiltonianShift);

      // type of lanczos algorithm (with or without reorthogonalization)
//	  BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
      gettimeofday (&(PrecalculationStartingTime), 0);
      FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	
      // initialization of lanczos algorithm
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = NbrEigenvalue + 3;
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
      Eigenvalues = new double [NbrEigenvalue];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  Eigenvalues[i] = TmpMatrix.DiagonalElement(i) - HamiltonianShift;
	}
      cout << endl;
      
      //compute eigenstates
      if (EigenstateFlag == true)
	Eigenstates = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
      gettimeofday (&(PrecalculationEndingTime), 0);
      Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) + 
	((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
      cout << "diagonalisation time = " << Dt << endl;
    }
  // insert here your code using the eigenvalues and the eigenvectors
  ofstream OutputFile;
  OutputFile.precision(14);
  OutputFile.open("eigenvalues", ios::binary | ios::out);
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      OutputFile << Eigenvalues[i] << " ";
      cout << Eigenvalues[i] << " ";
    }
  cout << endl;
  OutputFile << endl;
  OutputFile.close();
  if ((EigenstateFlag == true) && (Eigenstates != 0))
    {
      char* TmpFileName = new char[256];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  sprintf  (TmpFileName, "eigenvector.%d", i);
	  Eigenstates[i].WriteAsciiVector(TmpFileName);
	}
      delete[] TmpFileName;
    }
  return 0;

  // enter here hilbert space dimension
/*  int Dimension = M*N*H;
  Potentiel = new double** [M];
  for (int i = 0; i < M; ++i)
    {
      Potentiel[i] = new double* [N];
      for (int j = 0; j < N; ++j)
	Potentiel[i][j] = new double [H];
    }
  TransfertX = new double** [M];
  for (int i = 0; i < M; ++i)
    {
      TransfertX[i] = new double* [M];
      for (int j = 0; j < M; ++j)
	TransfertX[i][j] = new double [M];
    }
  TransfertY = new double** [N];
  for (int i = 0; i < N; ++i)
    {
      TransfertY[i] = new double* [N];
      for (int j = 0; j < N; ++j)
	TransfertY[i][j] = new double [N];
    }
  TransfertZ = new double** [H];
  for (int i = 0; i < H; ++i)
    {
      TransfertZ[i] = new double* [H];
      for (int j = 0; j < H; ++j)
	TransfertZ[i][j] = new double [H];
    }
  Entree();

  RealSymmetricMatrix HamiltonianRepresentation (Dimension);

  for (int i = 0; i < Dimension; ++i)
    for (int j = 0; j <= i; ++j)
      {
	// insert here your algorithm to fill the hamiltonian representation
	HamiltonianRepresentation(j, i) = EvaluateMatrixElement(i, j);
      }

  // store matrix in an hamiltonian class and create corresponding hilbert space
  UndescribedHilbertSpace HilbertSpace(Dimension);
  ExplicitHamiltonian Hamiltonian(&HilbertSpace, &HamiltonianRepresentation);
  RealVector* Eigenstates = 0;
  double* Eigenvalues = 0;

  // find the eigenvalues (and eigenvectors if needed)
  if (LanczosFlag == false)
    {
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  // diagonalize the hamiltonian
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Dimension);

	  if (EigenstateFlag == false)
	    {
	      ((RealSymmetricMatrix*) Hamiltonian.GetHamiltonian())->Householder(TmpTriDiag, MACHINE_PRECISION);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	    }
	  else
	    {
	      RealMatrix TmpEigenvectors (Dimension, Dimension);
	      ((RealSymmetricMatrix*) Hamiltonian.GetHamiltonian())->Householder(TmpTriDiag, MACHINE_PRECISION, TmpEigenvectors);
	      TmpTriDiag.Diagonalize(TmpEigenvectors);
	      TmpTriDiag.SortMatrixUpOrder(TmpEigenvectors);
	      Eigenstates = new RealVector [NbrEigenvalue];
	      for (int i = 0; i < NbrEigenvalue; ++i)
		{
		  Eigenstates[i] = TmpEigenvectors[i];
		}
	    }

	  // store eigenvalues
	  int Max = Hamiltonian.GetHilbertSpaceDimension();
	  if (Max > NbrEigenvalue)
	    Max = NbrEigenvalue;
	  Eigenvalues = new double [NbrEigenvalue];
	  for (int j = 0; j < Max ; j++)
	    {
	      Eigenvalues[j]= TmpTriDiag.DiagonalElement(j);
	    }
	}
      else
	{
	  cout << (*(Hamiltonian.GetHamiltonian()))(0, 0) << endl;
	}
    }
  else
    {
      
      // architecture type (i.e. 1 CPU or multi CPU)
      AbstractArchitecture* Architecture;
      if (SMPFlag == true)
	Architecture = new SMPArchitecture(2);
      else
	Architecture = new MonoProcessorArchitecture;


      // type of lanczos algorithm (with or without reorthogonalization)
//	  BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
      FullReorthogonalizedLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
	
      // initialization of lanczos algorithm
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = NbrEigenvalue + 3;
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
      Eigenvalues = new double [NbrEigenvalue];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  Eigenvalues[i] = TmpMatrix.DiagonalElement(i) + (h/Me)*(pow(1.0/a,2)+pow(1.0/b,2)+pow(1.0/c,2));
	}
      cout << endl;
      
      //compute eigenstates
      if (EigenstateFlag == true)
	Eigenstates = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
    }


  // insert here your code using the eigenvalues and the eigenvectors
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      cout << Eigenvalues[i] << " ";
    }
  cout << endl;
  if ((EigenstateFlag == true) && (Eigenstates != 0))
    {
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  for (int j = 0; j < Dimension; ++j)
	    cout << Eigenstates[i][j] << " ";
	}
      cout << endl;
    }
*/
  return 0;
}

//calculer l'energie cinetique
inline double cinetique(int i, int j, int k)  {
  return (h/Me)*((pow(i/(M*a),2)+pow(j/(N*b),2)+pow(k/(H*c),2)) - (pow(1.0/a,2)+pow(1.0/b,2)+pow(1.0/c,2)));
}


// read interaction coefficients in a file
//
// fileName = string containing file name
// return value = tridimensionnal array where interaction coefficients are stored (null pointer if an error occurs)

double*** ReadInteractionCoefficients(char* fileName, int nbrCellX, int nbrCellY, int nbrCellZ)
{
  ifstream InputFile;
  InputFile.open(fileName, ios::binary | ios::in);
  double*** TmpArray = new double** [nbrCellX];
  for (int i = 0; i < nbrCellX; ++i)
    {
      TmpArray[i] = new double* [nbrCellY];
      for (int j = 0; j < nbrCellY; ++j)
	TmpArray[i][j] = new double [nbrCellZ];
    }
  for (int k = 0; k < nbrCellZ; ++k)
    for (int j = 0; j < nbrCellY; ++j)
      for (int i = 0; i < nbrCellX; ++i)
	InputFile >> TmpArray[i][j][k];
  InputFile.close();

  return TmpArray;
}

//entrer les elements
inline void Entree(){

  //entrer les elements de matrice de transfert sur une chaine lineaire
  ifstream IntegralX("IntegralX.txt");
  ifstream IntegralY("IntegralY.txt");
  ifstream IntegralZ("IntegralZ.txt");

  for (int i = 0; i < M; i++){
    for (int j = 0; j < M; j++){
      for (int k = 0; k < M; k++){
	IntegralX >> TransfertX[i][j][k];
      }
    }
    IntegralX.ignore(4,'\n');
  }
  IntegralX.close();

  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      for (int k = 0; k < N; k++){
        IntegralY >> TransfertY[i][j][k];
      }
    }
    IntegralY.ignore(4,'\n');
  }
  IntegralY.close();

  for (int i = 0; i < H; i++){
    for (int j = 0; j < H; j++){
      for (int k = 0; k < H; k++){
	IntegralZ >> TransfertZ[i][j][k];
      }
    }
    IntegralZ.ignore(4,'\n');
  }
  IntegralZ.close();

  //matrice de potentiel
  ifstream Potential("Potentiel.txt");
  for (int k = 0; k < H; k++){
    for (int j = 0; j < N; j++){
      for (int i = 0; i < M; i++){
	Potential >> Potentiel[i][j][k];
      }
    Potential.ignore(4,'\n');
    }
  }
  Potential.close();
}

// evaluate matrix element
//

double EvaluateMatrixElement(int i, int j)
{

  double tmp = 0.0;

  //indices horizontales
  int m1 = (((i)%(M*N))%M);
  int n1 = (((i)%(M*N))/M);
  int p1 = ((i)/(M*N));

  if (i == j)  
    tmp = cinetique(m1 + 1, n1 + 1, p1 + 1);

  //indices verticales
  int m2 = (((j)%(M*N))%M);
  int n2 = (((j)%(M*N))/M);
  int p2 = ((j)/(M*N));

  //entrer l'element de matrice
  for (int k = 0; k < M; k ++)
    for (int l = 0; l < N; l++)
      for (int q = InNz1; q <= InNz2; q++)
	tmp += Potentiel[k][l][q]*TransfertX[k][m1][m2]*TransfertY[l][n1][n2]*TransfertZ[q][p1][p2];
	
  return tmp;
}
