#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/PlanarRotationSymmetryZPeriodicOneParticle.h"

#include "Hamiltonian/CylindricalQuantumDots3DHamiltonian.h"

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

#include "Tools/Potential/OneDConstantCellPotential.h"
#include "Tools/Potential/ThreeDConstantCylinderPotential.h"
#include "Tools/Potential/QuantumDotThreeDConstantCylinderPotential.h"

#include "Tools/Spectra/CylinderQuantumDotSpectra.h"

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
  OptionManager Manager ("CylinderQuantumDot" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* AdditionalGroup = new OptionGroup ("additional options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Architecture.AddOptionGroup(&Manager);
  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += LanczosGroup;
  Manager += AdditionalGroup;
  Manager += MiscGroup;

  (*PotentialGroup) += new SingleDoubleOption ('\n', "radius", "radius of the supercylinder (in Angstrom unit)", 1000);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "barrier", "number of cells in the well barrier", 10.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "below", "width of the layer below the wetting layer (in Angstrom unit)", 10.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "wetting", "width of the wetting layer (in Angstrom unit)", 5.0);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "nbr-dot", "number of uniformly high layer in the dot", 3);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "base", "base radius in Angstrom unit", 100.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "height", "height of dot in Angstrom unit", 17.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "top", "top radius in Anstrom unit", 74.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "above", "width of the layer above the dot layer (in Angstrom unit)", 70.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "dot", "potential in the dot", -0.4);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "well", "potential in the well", 1.079);

  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-r", "electron effective mass in plane (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('R', "R-states", "number of states in plane", 50);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('Z', "Z-states", "number of cells in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('k', "wave", "wave vector of Bloch function in Z direction (in 1/Angstrom unit)", 0.0);

  (*LanczosGroup) += new SingleIntegerOption ('n', "nbr-eigen", "number of eigenvalues", 10);
  (*LanczosGroup) += new BooleanOption ('e', "eigenstate", "evaluate eigenstates", false);
  (*LanczosGroup) += new SingleIntegerOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup) += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 500);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 400);  

  (*AdditionalGroup) += new SingleDoubleOption ('\n', "down-barrier", "width of barrier layer just below the WL (in Angstrom unit)", 0);
  (*AdditionalGroup) += new SingleDoubleOption ('\n', "up-barrier", "width of barrier layer just above  the dot (in Angstrom unit)", 0);
  (*AdditionalGroup) += new SingleDoubleOption ('\n', "potential", "potential in the barrier (in eV unit)", 0);

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

  double SuperCylinderRadius = ((SingleDoubleOption*) Manager["radius"])->GetDouble();
  double Barrier = ((SingleDoubleOption*) Manager["barrier"])->GetDouble();
  double Below = ((SingleDoubleOption*) Manager["below"])->GetDouble();
  double WettingWidth = ((SingleDoubleOption*) Manager["wetting"])->GetDouble();
  double BaseRadius = ((SingleDoubleOption*) Manager["base"])->GetDouble();
  double DotHeight = ((SingleDoubleOption*) Manager["height"])->GetDouble();
  int DotNbr = ((SingleIntegerOption*) Manager["nbr-dot"])->GetInteger();
  double TopRadius = ((SingleDoubleOption*) Manager["top"])->GetDouble();
  double Above = ((SingleDoubleOption*) Manager["above"])->GetDouble();
  double DotPotential = ((SingleDoubleOption*) Manager["dot"])->GetDouble();
  double WellPotential = ((SingleDoubleOption*) Manager["well"])->GetDouble();

  int NbrStateR = ((SingleIntegerOption*) Manager["R-states"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["Z-states"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  int NumberM = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  double Mur = ((SingleDoubleOption*) Manager["mu-r"])->GetDouble();
  double Muz = ((SingleDoubleOption*) Manager["mu-z"])->GetDouble();
  double WaveVector = ((SingleDoubleOption*) Manager["wave"])->GetDouble();

  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();   
  bool EigenstateFlag = ((BooleanOption*) Manager["eigenstate"])->GetBoolean();
  int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["iter-max"])->GetInteger();
  bool DiskFlag = ((BooleanOption*) Manager["disk"])->GetBoolean();
  bool ResumeFlag = ((BooleanOption*) Manager["resume"])->GetBoolean();
  int NbrIterLanczos = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int VectorMemory = ((SingleIntegerOption*) Manager["nbr-vector"])->GetInteger();

  double BelowBarrier = ((SingleDoubleOption*) Manager["down-barrier"])->GetDouble();  
  double AboveBarrier = ((SingleDoubleOption*) Manager["up-barrier"])->GetDouble();  
  double BarrierPotential = ((SingleDoubleOption*) Manager["potential"])->GetDouble();

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

  // QuantumDotThreeDConstantCylinderPotential(double belowHeight, double wettingWidth, int nbrCylinderDot, double dotHeight, double baseRadius, double topRadius, double aboveHeight);
  QuantumDotThreeDConstantCylinderPotential* potential = new QuantumDotThreeDConstantCylinderPotential(Below, WettingWidth, DotNbr, DotHeight, BaseRadius, TopRadius, Above, Barrier, SuperCylinderRadius);
  // void ConstructPotential(double dotPotential, double wellPotential);
  potential->ConstructPotential(DotPotential, WellPotential);
  potential->AddBarrierPotential (BelowBarrier, AboveBarrier, BarrierPotential);

  // define Hilbert space
  // PlanarRotationSymmetryZPeriodicOneParticle(int nbrStateR, int nbrStateZ, int lowerImpulsionZ);
  PlanarRotationSymmetryZPeriodicOneParticle* Space = new PlanarRotationSymmetryZPeriodicOneParticle(NumberM, NbrStateR, NbrStateZ, LowImpulsionZ);  

  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);
  
  //cout << "General space dimension: " << GeneralSpace.GetHilbertSpaceDimension() << endl;
  //cout << "Hilbert space Component: " << Space->GetNbrStateR() << '\t' << Space->GetNbrStateZ() << endl;
  //cout << "Minimal impulsions:       " << Space->GetLowerImpulsionZ() << endl;
  //cout << "Hilbert space dimension: " << Space->GetHilbertSpaceDimension() << endl;
  
  // CylindricalQuantumDots3DHamiltonian(PlanarRotationSymmetryZPeriodicOneParticle* space, double mur, double muz, double waveVector, ThreeDConstantCylinderPotential* PotentialInput);
  CylindricalQuantumDots3DHamiltonian Hamiltonian(Space, Mur, Muz, WaveVector, potential);

  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "Precalculation time = " << Dt << endl;

  ComplexVector* Eigenstates = 0;
  double* Eigenvalues = 0;

  cout << "----------------------------------------------------------------" << endl;

  double HamiltonianShift = -Hamiltonian.MaxPartialDiagonalElement();
  Hamiltonian.ShiftHamiltonian (HamiltonianShift);

  cout << "Hamiltonian shift =  " << HamiltonianShift << endl;
  gettimeofday (&(PrecalculationStartingTime), 0);

  // type of lanczos algorithm (with or without reorthogonalization)
  AbstractLanczosAlgorithm* Lanczos;
  if (DiskFlag == false)
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);   
  else
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage(Architecture.GetArchitecture(), NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);

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
          OutputFile.open("eigenvalues", ios::binary | ios::out | ios::app);	 
          for (int i = 0; i < NbrEigenvalue; ++i)
            OutputFile << Eigenvalues[i] << " ";          
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
