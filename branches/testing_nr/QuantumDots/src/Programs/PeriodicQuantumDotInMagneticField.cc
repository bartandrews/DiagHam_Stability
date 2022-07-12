#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/PeriodicThreeDOneParticle.h"

#include "Hamiltonian/PeriodicQuantumDots3DHamiltonian.h"
#include "Hamiltonian/PeriodicQuantumDots3DHamiltonianInMagneticField.h"

#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Potential/ThreeDConstantCellPotential.h"
#include "Tools/Potential/DotEmbeddedWellThreeDConstantCellPotential.h"

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


int main(int argc, char** argv)
{  
  cout.precision(14);

  OptionManager Manager ("PeriodicQuantumDotInMagneticField" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Architecture.AddOptionGroup(&Manager);
  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += LanczosGroup;
  Manager += MiscGroup;

  (*PotentialGroup) += new SingleIntegerOption ('M', "M-cell", "number of cells in the x direction", 161);
  (*PotentialGroup) += new SingleIntegerOption ('N', "N-cell", "number of cells in the y direction", 161);
  (*PotentialGroup) += new SingleIntegerOption ('H', "H-cell", "number of cells in the z direction", 21);
  (*PotentialGroup) += new SingleDoubleOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 5.65);
  (*PotentialGroup) += new SingleDoubleOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 5.65);
  (*PotentialGroup) += new SingleDoubleOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 5.65);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "barrier", "number of cells in the well barrier", 2);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "below", "number of cells between well barrier and wetting layer", 2);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "wetting", "number of cells in wetting layer", 1);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "base", "base radius in cell unit", 18);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "height", "height of dot in cell unit", 3);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "top", "top radius in cell unit", 13);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "well", "potential in the well", 1.079);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "dot", "potential in the dot", -0.4);
  (*PotentialGroup) += new SingleDoubleOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 30);

  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowx", "lower impulsion in x direction", -40);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowy", "lower impulsion in y direction", -40);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('k', "wave", "wave vector of Bloch function in Z direction (in 1/Angstrom unit)", 0.0);

  (*LanczosGroup) += new SingleIntegerOption ('n', "nbr-eigen", "number of eigenvalues", 200);
  (*LanczosGroup) += new BooleanOption ('e', "eigenstate", "evaluate eigenstates", false);
  (*LanczosGroup) += new SingleIntegerOption ('\n', "iter-max", "maximum number of lanczos iteration", 5000);
  (*LanczosGroup) += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 500);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 400);  

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new BooleanOption ('v', "verbose", "verbose mode", false);

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

  int M = ((SingleIntegerOption*) Manager["M-cell"])->GetInteger();
  int N = ((SingleIntegerOption*) Manager["N-cell"])->GetInteger();
  int H = ((SingleIntegerOption*) Manager["H-cell"])->GetInteger();
  double Lx = ((SingleDoubleOption*) Manager["cell-xsize"])->GetDouble();
  double Ly = ((SingleDoubleOption*) Manager["cell-ysize"])->GetDouble();
  double Lz = ((SingleDoubleOption*) Manager["cell-zsize"])->GetDouble();
  int UnderBarrier = ((SingleIntegerOption*) Manager["barrier"])->GetInteger();
  int BelowWettingLayer = ((SingleIntegerOption*) Manager["below"])->GetInteger();
  int WettingWidth = ((SingleIntegerOption*) Manager["wetting"])->GetInteger();
  int BaseRadius = ((SingleIntegerOption*) Manager["base"])->GetInteger();
  int DotHeight = ((SingleIntegerOption*) Manager["height"])->GetInteger();
  int TopRadius = ((SingleIntegerOption*) Manager["top"])->GetInteger();
  double WellPotential = ((SingleDoubleOption*) Manager["well"])->GetDouble();
  double DotPotential = ((SingleDoubleOption*) Manager["dot"])->GetDouble();
  double MagneticField = ((SingleDoubleOption*) Manager["magnetic"])->GetDouble();

  double Mux = ((SingleDoubleOption*) Manager["mu-x"])->GetDouble();
  double Muy = ((SingleDoubleOption*) Manager["mu-y"])->GetDouble();
  double Muz = ((SingleDoubleOption*) Manager["mu-z"])->GetDouble();
  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int LowImpulsionX = ((SingleIntegerOption*) Manager["lowx"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int LowImpulsionY = ((SingleIntegerOption*) Manager["lowy"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  double WaveVector = ((SingleDoubleOption*) Manager["wave"])->GetDouble();

  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();   
  bool EigenstateFlag = ((BooleanOption*) Manager["eigenstate"])->GetBoolean();
  int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["iter-max"])->GetInteger();
  bool DiskFlag = ((BooleanOption*) Manager["disk"])->GetBoolean();
  bool ResumeFlag = ((BooleanOption*) Manager["resume"])->GetBoolean();
  int NbrIterLanczos = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int VectorMemory = ((SingleIntegerOption*) Manager["nbr-vector"])->GetInteger();

  bool VerboseFlag = ((BooleanOption*) Manager["verbose"])->GetBoolean();

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
 
  // DotEmbeddedWellThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int underBarrier, int belowWettingLayer, int wettingWidth, int baseRadius, int dotHeight, int topRadius)
  DotEmbeddedWellThreeDConstantCellPotential* potential = new DotEmbeddedWellThreeDConstantCellPotential(M, N, H, UnderBarrier, BelowWettingLayer, WettingWidth, BaseRadius, DotHeight, TopRadius);
  //ConstructPotential(double wellPotential, double dotPotential)
  potential->ConstructPotential(WellPotential, DotPotential);
  potential->SavePotential("DotPotential.txt");
  
  // define Hilbert space
  PeriodicThreeDOneParticle* Space = new PeriodicThreeDOneParticle (NbrStateX, LowImpulsionX, NbrStateY, LowImpulsionY, NbrStateZ, LowImpulsionZ);

  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);
  
  //cout << "General space dimension: " << GeneralSpace.GetHilbertSpaceDimension() << endl;
  cout << "Sample size in cell unit: " << M << '\t' << N << '\t' << H << endl;
  cout << "Hilbert space dimensions: " << Space->GetNbrStateX() << '\t' << Space->GetNbrStateY() << '\t' << Space->GetNbrStateZ() << endl;
  cout << "Minimal impulsions:       " << Space->GetLowerImpulsionX() << '\t' << Space->GetLowerImpulsionY() << '\t' << Space->GetLowerImpulsionZ() << endl;

  PeriodicQuantumDots3DHamiltonianInMagneticField Hamiltonian(Space, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, 0.0, 0.0, MagneticField, potential, WaveVector);

  cout << endl;
  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "Precalculation time = " << Dt << endl;

  ComplexVector* Eigenstates = 0;
  double* Eigenvalues = 0;

  // type of lanczos algorithm (with or without reorthogonalization)
  AbstractLanczosAlgorithm* Lanczos;
  if (DiskFlag == false)
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);   
  else
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage(Architecture.GetArchitecture(), NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);

  cout << "----------------------------------------------------------------" << endl;

  double HamiltonianShift = -Hamiltonian.MaxPartialDiagonalElement();
  Hamiltonian.ShiftHamiltonian (HamiltonianShift);
  cout << "Hamiltonian shift =  " << HamiltonianShift << endl;
  gettimeofday (&(PrecalculationStartingTime), 0);



  cout << "Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl; 

  // initialization of lanczos algorithm
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 0;
  Lanczos->SetHamiltonian(&Hamiltonian);
  if ((DiskFlag == true) && (ResumeFlag == true))
    Lanczos->ResumeLanczosAlgorithm();
  else
    Lanczos->InitializeLanczosAlgorithm();
  cout << "------------------- Run Lanczos Algorithm ---------------------" << endl;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  if (ResumeFlag == false)
    {
      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
      CurrentNbrIterLanczos = NbrEigenvalue + 3;
      if ((DiskFlag == true) && (CurrentNbrIterLanczos >= NbrIterLanczos))
        {
          NbrIterLanczos = CurrentNbrIterLanczos + 1;
        }
    }
  RealTriDiagonalSymmetricMatrix TmpMatrix;
  while ((Lanczos->TestConvergence() == false) &&  (((DiskFlag == true) && (CurrentNbrIterLanczos < NbrIterLanczos)) ||
                                                    ((DiskFlag == false) && (CurrentNbrIterLanczos < MaxNbrIterLanczos))))
    {
      ++CurrentNbrIterLanczos;
      Lanczos->RunLanczosAlgorithm(1);
      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest;
      if (VerboseFlag)
	cout << CurrentNbrIterLanczos << "\t" <<  TmpMatrix.DiagonalElement(0) - HamiltonianShift << "\t\t" << Lowest - HamiltonianShift << "\t\t" << Precision << endl;
    }
  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
    {
      cout << "too much Lanczos iterations" << endl;
      exit(0);
    }
  cout << "------------------ Actual eigenvalues ------------------" << endl;
  // store eigenvalues
  Eigenvalues = new double [NbrEigenvalue];
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      Eigenvalues[i] = (TmpMatrix.DiagonalElement(i) - HamiltonianShift);
      cout << Eigenvalues[i] << '\t';
    }
  cout << endl;
  if (Lanczos->TestConvergence())
    {
      //compute eigenstates
      if (EigenstateFlag == true)
        Eigenstates = (ComplexVector*) Lanczos->GetEigenstates(NbrEigenvalue);
      
      // insert here your code using the eigenvalues and the eigenvectors
      if (EigenstateFlag == true)
        {
          ofstream OutputFile;
          OutputFile.precision(14);
          OutputFile.open("eigenvalues", ios::binary | ios::out);
          for (int i = 0; i < NbrEigenvalue; ++i)
            OutputFile << Eigenvalues[i] << " ";
          OutputFile << endl;
          OutputFile.close();
        }
      
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
            
      cout << "----------------- End of calculation ---------------------" << endl;      
      cout << "     ==========  CALCULATION IS FINALIZED  =========  " << endl;
      cout << "Sample size in cell unit: " << M << '\t' << N << '\t' << H << endl;
      cout << "Hilbert space dimensions: " << Space->GetNbrStateX() << '\t' << Space->GetNbrStateY() << '\t' << Space->GetNbrStateZ() << endl;
      cout << "Minimal impulsions:       " << Space->GetLowerImpulsionX() << '\t' << Space->GetLowerImpulsionY() << '\t' << Space->GetLowerImpulsionZ() << endl;
    }
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);  
  cout << endl << "Total time = " << Dt << endl;
  delete Lanczos;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Command << "The program finished at: " <<  asctime (timeinfo);
  Command << "============================== End =============================" << '\n';
  FullOption << "The program finished at: " << asctime (timeinfo);
  FullOption << "=============================== End ============================" << '\n'; 
  Command.close(); FullOption.close(); 

  return 0;
}
