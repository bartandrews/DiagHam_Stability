#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/ThreeDOneParticle.h"

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

#include "Tools/Potential/HardBoxPyramidQuantumDotThreeDConstantCellPotential.h"

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

//Definir des constantes
int M = 14;//nombre de cellules dans l'axe x
int N = 14;//nombre de cellules dans l'axe y
int H = 10;//nombre de mailles en z

//profil de potentiel
double*** Potentiel;

//1er indice: cellule d'InN; 2e et 3e: indices de fonction
double*** TransfertX;
double*** TransfertY;
double*** TransfertZ;

// read interaction coefficients in a file
//
// fileName = string containing file name
// return value = tridimensionnal array where interaction coefficients are stored (null pointer if an error occurs)
double*** ReadInteractionCoefficients(char* fileName, int nbrCellX, int nbrCellY, int nbrCellZ);

int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  BooleanOption LanczosOption ('l', "lanczos", "enable lanczos diagonalization algorithm", true);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode", true);
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", true);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 6);
  SingleIntegerOption MemoryOption ('\n', "memory", "amount of memory that can be used for precaching (in Mb)", 1000);
  SingleIntegerOption MValueOption ('M', "M-cell", "number of cells in the x direction", 50);
  SingleIntegerOption NValueOption ('N', "N-cell", "number of cells in the y direction", 50);
  SingleIntegerOption HValueOption ('H', "H-cell", "number of cells in the z direction", 30);
  SingleIntegerOption LeftSizeOption ('\n', "left-size", "size of the leftmost part in the z direction with constant null potential (in cell unit)", 6);
  SingleIntegerOption RightSizeOption ('\n', "right-size", "size of the rightmost part in the z direction with constant potential (in cell unit)", 10);
  SingleDoubleOption CellXSizeOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 2.97);
  SingleDoubleOption CellYSizeOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 2.97);
  SingleDoubleOption CellZSizeOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 2.64);  
  SingleDoubleOption XMassOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.5045);
  SingleDoubleOption YMassOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.5045);
  SingleDoubleOption ZMassOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 1.1);
  SingleStringOption CoefficientFileNameOption('\n', "coefficients", "name of the file where interaction coeffcients are stored", 
					       "/home/regnault/development/DMRG/DiagHam/potentiel_10_10_10_2");
  BooleanOption CarrierTypeOption('c', "carrier", "carrier type, true for hole, false for electron", true);

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
  OptionList += &LeftSizeOption;
  OptionList += &RightSizeOption;
  OptionList += &MemoryOption;
  OptionList += &CarrierTypeOption; 

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
  bool Carrier = CarrierTypeOption.GetBoolean();

  HardBoxPyramidQuantumDotThreeDConstantCellPotential* potential = new HardBoxPyramidQuantumDotThreeDConstantCellPotential(M, N, H, LeftSize, RightSize, 1, 2, 1);  
  potential->LoadPotentialWithConstantField(CoefficientFileName);

  /*
  // **** PROBABILITIES ****
  double p1 = 0.0, p2 = 0.0;

  // **** Offset ****
  double VH = 0.9, VE = 1.8;

  // *** Electric field (absolute value) ****;
  double F_Dot = 0.0245, F_Barrier = 0.001;

  // *** Dot geometry ****
  int Rb = 20, Rt = 5, w = 4;
  
  if (Carrier)
    {
      // NITRIDE HOLE           
      Mux = 0.5045; Muy = 0.5045; Muz = 1.1;
      potential->SegregationPyramidDot(p1, p2, Rb, Rt, w, -F_Barrier, F_Dot, F_Dot, -F_Barrier, VH, Lz, true, "DotInput.txt");
      ofstream potenf("DotPotential.txt");
      potential->PrintPotentialWithField(potenf);
      potenf.close();
      ofstream diagram ("Diagram.txt");
      potential->PrintDiagram(diagram);
      diagram.close();
    }
  else
    {
      // NITRIDE ELECTRON           
      Mux = 0.166; Muy = 0.166; Muz = 0.184;
      potential->ReadDiagram("../h/Diagram.txt");
      potential->SegregationPyramidDot(p1, p2, Rb, Rt, w, F_Barrier, -F_Dot, -F_Dot, F_Barrier, VE, Lz, false, "DotInput.txt");
      ofstream potenf("DotPotential.txt");
      potential->PrintPotentialWithField(potenf);
      potenf.close();
    }
  */

  ThreeDOneParticle Space(M, N, H);
  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);

  //  QuantumDots3DHamiltonian(ThreeDOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput, int memory = -1);
  QuantumDots3DHamiltonian Hamiltonian(&Space, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, M, N, H - LeftSize - RightSize, potential, Memory << 20);

  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "precalculation time = " << Dt << endl;

  ofstream Input;
  Input.open("DotInput.txt", ios::out | ios::app);
  Input << "Lattice constants: X = " << Lx << ", Y = " << Ly << ", Z = " << Lz << '\n';
  Input << "Effective masses: Mx = " << Mux << ", My = " << Muy << ", Mz = " << Muz << '\n';
  Input << "Precalculation time: " << Dt << '\n';  
 
  RealVector* Eigenstates = 0;
  double* Eigenvalues = 0;
  if ((LanczosFlag == false) || (Space.GetHilbertSpaceDimension() < 300))
    {
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  // diagonalize the hamiltonian
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
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

      double HamiltonianShift = -Hamiltonian.MaxPartialDiagonalElement();
      Hamiltonian.ShiftHamiltonian (HamiltonianShift);
      cout << "Décalage:  " << HamiltonianShift << endl;
      // type of lanczos algorithm (with or without reorthogonalization)
      // BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
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
  Input << "Diagonalization time: " << Dt << '\n';
  Input << "To verify: M, N, H, under, above = " << M << ", " << N << ", " << H << ", " << LeftSize << ", " << RightSize << endl;
  Input.close();

  return 0;
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
