#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/PairXPairYPeriodic3DOneParticle.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Hamiltonian/PeriodicQuantumDots3DHamiltonian.h"
#include "Hamiltonian/XYReflexionSymmetricPeriodic3DHamiltonian.h"

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
#include <fstream.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ofstream;

bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);

int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("ExplicitPeriodic3DQuantumDots" , "0.01");
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

  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowx", "lower impulsion in x direction", -40);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowy", "lower impulsion in y direction", -40);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);

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

  double Mux = ((SingleDoubleOption*) Manager["mu-x"])->GetDouble();
  double Muy = ((SingleDoubleOption*) Manager["mu-y"])->GetDouble();
  double Muz = ((SingleDoubleOption*) Manager["mu-z"])->GetDouble();
  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int LowImpulsionX = ((SingleIntegerOption*) Manager["lowx"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int LowImpulsionY = ((SingleIntegerOption*) Manager["lowy"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();

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

  // ConstructPotential(double wellPotential, double dotPotential)
  //potential->ConstructPotential(WellPotential, DotPotential);
  potential->LoadPotential("DotPotential.txt");

  Periodic3DOneParticle* Space = new Periodic3DOneParticle(NbrStateX, LowImpulsionX, NbrStateY, LowImpulsionY, NbrStateZ, LowImpulsionZ);

  int NbrCellX = M, NbrCellY = N, NbrCellZ = H; 
  int NbrStateX = Space->GetNbrStateX(), NbrStateY = Space->GetNbrStateY(), NbrStateZ = Space->GetNbrStateZ();  
  int LowerImpulsionX = Space->GetLowerImpulsionX(), LowerImpulsionY = Space->GetLowerImpulsionY(), LowerImpulsionZ = Space->GetLowerImpulsionZ();

  double PERIODIC_HAMILTONIAN_FACTOR = 150.4;
 
  double** RealWaveFunctionOverlapX; double** ImaginaryWaveFunctionOverlapX;
  double** RealWaveFunctionOverlapY; double** ImaginaryWaveFunctionOverlapY; 
  double** RealWaveFunctionOverlapZ; double** ImaginaryWaveFunctionOverlapZ;
  double XSize = NbrCellX * Lx, YSize = NbrCellY * Ly, ZSize = NbrCellZ * Lz;
  int Dimension = NbrStateX * NbrStateY * NbrStateZ;

  if (!EvaluateWaveFunctionOverlap(NbrCellX, NbrStateX, RealWaveFunctionOverlapX, ImaginaryWaveFunctionOverlapX))
    cout << "Error in evaluation of function overlap in X direction. Stop!" << endl;  
  if (!EvaluateWaveFunctionOverlap(NbrCellY, NbrStateY, RealWaveFunctionOverlapY, ImaginaryWaveFunctionOverlapY))
    cout << "Error in evaluation of function overlap in Y direction. Stop!" << endl;
  if (!EvaluateWaveFunctionOverlap(NbrCellZ, NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (Mux * XSize * XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (Muy * YSize * YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (Muz * ZSize * ZSize);
  
  double* KineticElements = new double[Dimension];

  double FactorX = 0.0, FactorY = 0.0;
  int TotalIndex1 = 0;
  for (int i = 0; i < NbrStateX; ++i)
    {
      FactorX = double((i + LowerImpulsionX) * (i + LowerImpulsionX)) * InvXFactor;
      for (int j = 0; j < NbrStateY; ++j)
	{
	  FactorY = double((j + LowerImpulsionY) * (j + LowerImpulsionY)) * InvYFactor + FactorX;
	  for (int k = 0; k < NbrStateZ; ++k)
	    {	      
	      KineticElements[TotalIndex1] = FactorY + double((k + LowerImpulsionZ) * (k + LowerImpulsionZ)) * InvZFactor;	      
	      ++TotalIndex1;
	    }
	}
    }

  int LengthX = (NbrStateX - 1) * 2 + 1; int LengthY = (NbrStateY - 1) * 2 + 1; int LengthZ = (NbrStateZ - 1) * 2 + 1;

  double*** TmpReal = new double** [LengthX];
  double*** TmpImaginary = new double** [LengthX];

  double TmpRe, TmpIm;
  double TmpRe2, TmpIm2;
  double* TmpRealWaveFunctionOverlapX;
  double* TmpImaginaryWaveFunctionOverlapX;
  double* TmpRealWaveFunctionOverlapY;
  double* TmpImaginaryWaveFunctionOverlapY;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;

  for (int m = 0; m < LengthX; ++m)
    {
      TmpReal[m] = new double* [LengthY];
      TmpImaginary[m] = new double* [LengthY];
      TmpRealWaveFunctionOverlapX = RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [NbrCellZ];
	  TmpImaginary[m][n] = new double [NbrCellZ];
	  TmpRealWaveFunctionOverlapY = RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < NbrCellX; ++CellX)
		    {		      
		      TmpRe += potential->GetPotential(CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpRe2 - TmpImaginaryWaveFunctionOverlapX[CellX] * TmpIm2);
		      TmpIm += potential->GetPotential(CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpIm2 + TmpImaginaryWaveFunctionOverlapX[CellX] * TmpRe2);		      
		    }
		}
	      TmpRealPrecalculatedHamiltonian[CellZ] = TmpRe;  
	      TmpImaginaryPrecalculatedHamiltonian[CellZ] = TmpIm;  
	    }
	}
    }

  double*** RealPrecalculatedHamiltonian = new double** [LengthX];
  double*** ImaginaryPrecalculatedHamiltonian = new double** [LengthX];
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  for (int m = 0; m < LengthX; ++m)
    {
      RealPrecalculatedHamiltonian[m] = new double* [LengthY];      
      ImaginaryPrecalculatedHamiltonian[m] = new double* [LengthY]; 
      for (int n = 0; n < LengthY; ++n)
	{
	  RealPrecalculatedHamiltonian[m][n] = new double [LengthZ];      
	  ImaginaryPrecalculatedHamiltonian[m][n] = new double [LengthZ]; 
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];
	  for (int p = 0; p < LengthZ; ++p)
	    {
	      TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
	      TmpImaginaryWaveFunctionOverlapZ = ImaginaryWaveFunctionOverlapZ[p];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ] - TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);
		  TmpIm += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ] + TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
		}
	      RealPrecalculatedHamiltonian[m][n][p] = TmpRe;
	      ImaginaryPrecalculatedHamiltonian[m][n][p] = TmpIm;
	    }
	}
    }
  delete[] TmpReal; delete[] TmpImaginary;

  HermitianMatrix HamiltonianRepresentation (Dimension);
  
  int m1, m2, n1, n2, p1, p2;
  int IndexX, IndexY, IndexZ;
  int TmpIndex = 0; int** TotalIndex = new int* [NbrStateX];
  for (m1 = 0; m1 < NbrStateX; ++m1) 
    {
      TotalIndex[m1] = new int [NbrStateY];
      for (n1 = 0; n1 < NbrStateY; ++n1)	
	{
	  TotalIndex[m1][n1] = (m1 * NbrStateY + n1) * NbrStateZ;
	  for (p1 = 0; p1 < NbrStateZ; ++p1)
	    {	      
	      HamiltonianRepresentation.SetMatrixElement(TmpIndex, TmpIndex, KineticElements[TmpIndex]);
	      ++TmpIndex;
	    }
	}
    }

  int OriginX = NbrStateX - 1; int OriginY = NbrStateY - 1; int OriginZ = NbrStateZ - 1;
  int Index1, Index2;
  Complex Tmp;
  int* TmpTotalIndex1; int* TmpTotalIndex2;
  for (m1 = 0; m1 < NbrStateX; ++m1)
    {
      for (n1 = 0; n1 < NbrStateY; ++n1)
	{
	  Index1 = TotalIndex[m1][n1];
	  for (p1 = 0; p1 < NbrStateZ; ++p1)
	    {	      
	      for (m2 = 0; m2 < NbrStateX; ++m2)
		{		  
		  IndexX = -m1 + m2 + OriginX;		  
		  for (n2 = 0; n2 < NbrStateY; ++n2)
		    {
		      IndexY = -n1 + n2 + OriginY;
		      Index2 = TotalIndex[m2][n2];
		      for (p2 = 0; p2 < NbrStateZ; ++p2)
			{
			  IndexZ = -p1 + p2 + OriginZ;
			  Tmp = Complex(RealPrecalculatedHamiltonian[IndexX][IndexY][IndexZ], ImaginaryPrecalculatedHamiltonian[IndexX][IndexY][IndexZ]);
			  HamiltonianRepresentation.SetMatrixElement(Index1, Index2, Tmp);
			  ++Index2;
			}

		    }
		}
	      ++Index1;
	    }
	}
    }


  // RealVector* Eigenstates = 0;
  double* Eigenvalues = new double [NbrEigenvalue];

  // find the eigenvalues (and eigenvectors if needed)

  
  double Precision;
  double PreviousLowest;
  double Lowest;
  int CurrentNbrIterLanczos;

  // type of lanczos algorithm (with or without reorthogonalization)
  AbstractLanczosAlgorithm* Lanczos;
  if (DiskFlag == false)
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);   
  else
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage(Architecture.GetArchitecture(), NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);

     
  ExplicitHamiltonian Hamiltonian(Space, &HamiltonianRepresentation);
  
  
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
      Eigenvalues[i] = TmpMatrix.DiagonalElement(i);
      cout << Eigenvalues[i] << '\t';
    }

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Command << "The program finished at: " <<  asctime (timeinfo);
  Command << "============================== End =============================" << '\n';
  FullOption << "The program finished at: " << asctime (timeinfo);
  FullOption << "=============================== End ============================" << '\n'; 
  Command.close(); FullOption.close(); 

  return 0;
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
