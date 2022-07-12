#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

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
const double Me = 0.166;//masse effective
const double h = 37.60;//h bar carre reduit

//profil de potentiel
double** Potentiel;

//1er indice: cellule d'InN; 2e et 3e: indices de fonction
double*** TransfertX;
double*** TransfertY;

//calculer l'energie cinetique
inline double cinetique(int i, int j);
//entrer les elements
inline void Entree();
//calculer l'integral des fonctions d'indice m et n
inline double integral(int m, int n, int i, double P);


// evaluate matrix element
//
double EvaluateMatrixElement(int i, int j);

int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  BooleanOption LanczosOption ('l', "lanczos", "enable lanczos diagonalization algorithm", false);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode");
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates");
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  SingleIntegerOption MValueOption ('M', "M-cell", "number of cells in the x direction", 10);
  SingleIntegerOption NValueOption ('N', "N-cell", "number of cells in the y direction", 10);
  SingleStringOption CoefficientFileNameOption('\0', "coefficients", "name of the file where interaction coeffcients are stored", "Potentiel_E.txt");
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
  bool LanczosFlag = LanczosOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  bool VerboseFlag = VerboseOption.GetBoolean();
  bool EigenstateFlag = EigenstateOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  char* CoefficientFileName = CoefficientFileNameOption.GetString(); 
  M = MValueOption.GetInteger();
  N = NValueOption.GetInteger();

  // initialize matrix associated to the hamiltonian

  // enter here hilbert space dimension
  int Dimension = M*N;
  Potentiel = new double* [M];
  for (int i = 0; i < M; ++i)
    {
      Potentiel[i] = new double [N];
    }
  ifstream InputFile;
  InputFile.open(CoefficientFileName, ios::binary | ios::in);
  for (int j = 0; j < N; ++j)
    for (int i = 0; i < M; ++i)
      InputFile >> Potentiel[i][j];
  InputFile.close();

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
  Entree();

  RealSymmetricMatrix HamiltonianRepresentation (Dimension);

  double tmp = 0.0;
  double tmp2 = 0.0;
  int m1;
  int n1;
  int m2;
  int n2;
  
  for (int i = 0; i < Dimension; ++i)
    {
      n1 = i / M;
      m1 = i - (n1 * M);
      for (int j = 0; j < i; ++j)
	{
	  tmp = 0.0;
	  n2 = (j / M);
	  m2 = j - (n2 * M);
	  for (int k = 0; k < M; ++k)
	    {
	      tmp2 = TransfertX[k][m1][m2];
	      for (int l = 0; l < N; ++l)
		tmp += Potentiel[k][l] * tmp2 * TransfertY[l][n1][n2];
	    }
	  HamiltonianRepresentation(j, i) = tmp;
	}
      tmp = cinetique(m1 + 1, n1 + 1);
      for (int k = 0; k < M; ++k)
	{
	  tmp2 = TransfertX[k][m1][m1];
	  for (int l = 0; l < N; ++l)
	    tmp += Potentiel[k][l] * tmp2 * TransfertY[l][n1][n1];
	}
      HamiltonianRepresentation(i, i) = tmp;
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
	  Eigenvalues[i] = TmpMatrix.DiagonalElement(i) + (h/Me)*((1.0/ (a * a)) + (1.0/ (b * b)));
	}
      cout << endl;
      
      //compute eigenstates
      if (EigenstateFlag == true)
	Eigenstates = (RealVector*) Lanczos.GetEigenstates(NbrEigenvalue);
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
}

//calculer l'energie cinetique
inline double cinetique(int i, int j)  {
  return (h/Me)*( ((i/(M*a)) * (i/(M*a))) + ((j/(N*b)) * (j/(N*b))) - (1.0/ (a * a)) - (1.0/ (b * b)));
}

//entrer les elements
inline void Entree(){

//  double temp1 = M*M*a*a; double temp2 = N*N*b*b;
  
/*  for (int i = 0; i < M; i++)
    KX[i] = (i+1)*(i+1)/temp1;
  for (int i = 0; i < N; i++)
    KY[i] = (i+1)*(i+1)/temp2;*/  
  
  //entrer les elements de matrice de transfert sur une chaine lineaire
  for (int i = 0; i < M; i++){
    for (int j = 0; j < M; j++){
      for (int k = 0; k < M; k++){
	TransfertX[i][j][k] = integral(j+1,k+1,i+1,M);
      }
    }
  }
  
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      for (int k = 0; k < N; k++){
	TransfertY[i][j][k] = integral(j+1,k+1,i+1,N);
      }
    }
  }

}

  

//calculer l'integral des fonctions d'indice m et n
inline double integral(int m, int n, int i, double P){
  if (m != n)
    return (1/M_PI)*(((sin(i*(m - n)*M_PI/P) - sin((i - 1)*(m - n)*M_PI/P))/(m -n)) - ((sin(i*(m + n)*M_PI/P) - sin((i - 1)*(m + n)*M_PI/P))/(m + n)));
  else
    return (1/P)*(1 + (P/(2*M_PI*n))*(sin(2*M_PI*n*(i - 1)/P) - sin(2*M_PI*n*i/P)));
}



  // evaluate matrix element
//

double EvaluateMatrixElement(int i, int j)
{

  double tmp = 0.0;

  
  //indices horizontales
  int m1 = (i % M);
  int n1 = (i / M);

  if (i == j)  
    tmp = cinetique(m1 + 1, n1 + 1);

  //indices verticales
  int m2 = (j % M);
  int n2 = (j / M);

  //entrer l'element de matrice
  for (int k = 0; k < M; k ++)
    for (int l = 0; l < N; l++)
      tmp += Potentiel[k][l] * TransfertX[k][m1][m2] * TransfertY[l][n1][n2];
	
  return tmp;
}
