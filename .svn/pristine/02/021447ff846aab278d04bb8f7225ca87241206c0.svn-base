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

#include "Tools/Potential/OneDConstantCellPotential.h"
#include "Tools/Potential/ThreeDConstantCylinderPotential.h"
#include "Tools/Potential/QuantumDotThreeDConstantCylinderPotential.h"

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

#define PERIODIC_HAMILTONIAN_FACTOR 150.4
#define BLOCH_FACTOR 7.619
#define PARAMAGNETIC_FACTOR         5.802e-5
#define DIAMAGNETIC_FACTOR          2.198e-10 // = e^2 * B^2 / 8m

// evaluate the wave function overlap
//
// potential = potential of the quantum dot
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap
bool EvaluateWaveFunctionOverlap(ThreeDConstantCylinderPotential* potential, int nbrState, double** &realArray, double** &imaginaryArray);

int main(int argc, char** argv)
{  
  cout.precision(14);

  OptionManager Manager ("PeriodicOneDStructure" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += MiscGroup;

  (*PotentialGroup) += new SingleDoubleOption ('\n', "barrier", "width of the well barrier (in Angstrom unit)", 10.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "below", "width of the layer below the wetting layer (in Angstrom unit)", 10.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "wetting", "width of the wetting layer (in Angstrom unit)", 5.0);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "nbr-dot", "number of uniformly high layer in the dot", 3);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "base", "base radius in Angstrom unit", 100.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "height", "height of dot in Angstrom unit", 17.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "top", "top radius in Anstrom unit", 74.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "above", "width of the layer above the dot layer (in Angstrom unit)", 70.0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "dot", "potential in the dot", -0.4);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "well", "potential in the well", 1.079);
  (*PotentialGroup) += new SingleDoubleOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 0);

  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-r", "electron effective mass in plane (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('d', "degree", "degree of the weight function", 1);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('s', "sigma", "value of sigma for gaussian weight (in Angstrom unit)", 40);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('l', "lambda", "value of lambda for second degree (2S states) weight (in Angstrom unit)", 40);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('k', "wave", "wave vector of Bloch function in Z direction (in 1/Angstrom unit)", 0.0);

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
  double MagneticField = ((SingleDoubleOption*) Manager["magnetic"])->GetDouble();

  double Mur = ((SingleDoubleOption*) Manager["mu-r"])->GetDouble();
  double Muz = ((SingleDoubleOption*) Manager["mu-z"])->GetDouble();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  int Degree = ((SingleIntegerOption*) Manager["degree"])->GetInteger();
  int NumberM = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  double Sigma = ((SingleDoubleOption*) Manager["sigma"])->GetDouble();
  double Lambda = ((SingleDoubleOption*) Manager["lambda"])->GetDouble();
  double WaveVector = ((SingleDoubleOption*) Manager["wave"])->GetDouble();

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

  QuantumDotThreeDConstantCylinderPotential* potential = new QuantumDotThreeDConstantCylinderPotential(Below, WettingWidth, DotNbr, DotHeight, BaseRadius, TopRadius, Above, Barrier, 1000.0);
  // void ConstructPotential(double dotPotential, double wellPotential);
  potential->ConstructPotential(DotPotential, WellPotential);
  if (Degree > 2)
    {
      cout << "This degree of the weight function is not taken into account: " << Degree << endl;
      return 1;
    }
  OneDConstantCellPotential* oneDPotential;
  if (Degree == 1)    
    oneDPotential = potential->GaussianReductionOneDimension(Sigma, NumberM);
  if (Degree == 2)
    oneDPotential = potential->SecondDegreeReductionOneDimension(Lambda, Sigma, NumberM);

  double** RealOverlap; double** ImaginaryOverlap;
  EvaluateWaveFunctionOverlap(potential, NbrStateZ, RealOverlap, ImaginaryOverlap);    
  
  double SizeZ = Barrier + Below + WettingWidth + DotHeight + Above;
  int Dimension = NbrStateZ;
  HermitianMatrix HamiltonianRepresentation(Dimension);
  int NumberZ = potential->GetNbrCylinderZ();

  double TmpRe = 0.0, TmpIm = 0.0; Complex TmpComplex;
  double ShiftSquareKz = BLOCH_FACTOR * WaveVector * WaveVector / (2.0 * Muz);
  double ShiftKz = BLOCH_FACTOR * WaveVector * 2.0 * M_PI/ (Muz * SizeZ);
  for (int n1 = 0; n1 < Dimension; ++n1)
    {
      for (int n2 = 0; n2 <= n1; ++n2)
	{	 
	  TmpRe = 0.0; TmpIm = 0.0;
	  for (int k = 0; k < NumberZ; ++k)
	    {
	      //cout << "n1: " << n1 << " ,n2:" << n2 << " , k:" << k << " , real: " << RealOverlap[-n1 + n2 + OriginZ][k] << " , imaginary: " << ImaginaryOverlap[-n1 + n2 + OriginZ][k] << endl;
	      TmpRe += (oneDPotential->GetPotential(k) * RealOverlap[n1 - n2][k]);
	      TmpIm += (oneDPotential->GetPotential(k) * ImaginaryOverlap[n1 - n2][k]);
	    }
	  TmpComplex = Complex(TmpRe, TmpIm);	  
	  HamiltonianRepresentation.SetMatrixElement(n1, n2, TmpComplex);
	}
      HamiltonianRepresentation.AddToMatrixElement(n1, n1, PERIODIC_HAMILTONIAN_FACTOR * double((n1 + LowImpulsionZ) * (n1 + LowImpulsionZ)) / (Muz * SizeZ * SizeZ) + ShiftSquareKz + ShiftKz * double(n1 + LowImpulsionZ));
    }

  RealMatrix TmpEigenvectors (Dimension, Dimension);
  RealSymmetricMatrix RealHamiltonianRepresentation = HamiltonianRepresentation.ConvertToSymmetricMatrix();
  RealTriDiagonalSymmetricMatrix TmpTriDiag (2 * Dimension);
  RealHamiltonianRepresentation.Householder(TmpTriDiag, MACHINE_PRECISION, TmpEigenvectors);
  TmpTriDiag.Diagonalize(TmpEigenvectors);
  TmpTriDiag.SortMatrixUpOrder(TmpEigenvectors);

  double* tmpE = new double [NbrEigenvalue];
  for (int i = 0; i < NbrEigenvalue; ++i)  
    {
      if (Degree == 1)
	{
	  if (NumberM == 0)
	    tmpE[i] = TmpTriDiag.DiagonalElement(2 * i) + DIAMAGNETIC_FACTOR * Sigma * Sigma * MagneticField * MagneticField / Mur + PERIODIC_HAMILTONIAN_FACTOR / (4.0 * M_PI * M_PI * Mur * Sigma * Sigma);
	  if (NumberM == 1)
	    tmpE[i] = TmpTriDiag.DiagonalElement(2 * i) + PARAMAGNETIC_FACTOR * MagneticField / Mur + 2.0 * DIAMAGNETIC_FACTOR * Sigma * Sigma * MagneticField * MagneticField / Mur + PERIODIC_HAMILTONIAN_FACTOR / (2.0 * M_PI * M_PI * Mur * Sigma * Sigma);
	  if (NumberM == -1)
	    tmpE[i] = TmpTriDiag.DiagonalElement(2 * i) - PARAMAGNETIC_FACTOR * MagneticField / Mur + 2.0 * DIAMAGNETIC_FACTOR * Sigma * Sigma * MagneticField * MagneticField / Mur + PERIODIC_HAMILTONIAN_FACTOR / (2.0 * M_PI * M_PI * Mur * Sigma * Sigma);
	}
      if (Degree == 2)
	{
	  double l = Lambda * Lambda; double s = Sigma * Sigma;
	  double a = 2 * l * s / (l + s);	 
	  double denominator = 2 * l * l - 2 * l * a + a * a;
	  if (NumberM == 0)
	    tmpE[i] = TmpTriDiag.DiagonalElement(2 * i) + PERIODIC_HAMILTONIAN_FACTOR * (2 * l * l + a * a) / (4.0 * M_PI * M_PI * Mur * denominator * l) + DIAMAGNETIC_FACTOR * l * MagneticField * MagneticField * (6 * l * l - 4 * l * a + a * a) / (Mur * denominator);
	}
    }

  for (int i = 0; i < NbrEigenvalue; ++i)  
    {
      cout << tmpE[i] << " ";    

      if ((NumberM != 0) && (NumberM != 1) && (NumberM != -1))
	cout << "This momentum is not taken into account" << endl;
      cout << endl;
    }
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
      OutputFile.open("eigenvalues", ios::binary | ios::out | ios::app);
      // OutputFile << Sigma << " ";

      for (int i = 0; i < NbrEigenvalue; ++i)	
	OutputFile << tmpE[i] << " ";	  	
      OutputFile << endl;
      OutputFile.close();

      char* TmpFileName = new char[256];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  sprintf  (TmpFileName, "eigenvector.%d", i);
	  Eigenstates[i].WriteAsciiVector(TmpFileName);
	}
      delete[] TmpFileName;
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
// potential = potential of the quantum dot
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool EvaluateWaveFunctionOverlap(ThreeDConstantCylinderPotential* potential, int nbrState, double** &realArray, double** &imaginaryArray)
{
  int nbrCylinder = potential->GetNbrCylinderZ();
  double* ZPosition = new double [nbrCylinder + 1];
  ZPosition[0] = 0.0;
  for (int k = 0; k < nbrCylinder; ++k)    
    ZPosition[k + 1] = ZPosition[k] + potential->GetHeight(k);      
  double ZSize = ZPosition[nbrCylinder];

  realArray = new double* [nbrState];
  imaginaryArray = new double* [nbrState];   

  realArray[0] = new double [nbrCylinder];
  imaginaryArray[0] = new double [nbrCylinder];
  for (int k = 0; k < nbrCylinder; ++k)
    {
      realArray[0][k] = (ZPosition[k + 1] - ZPosition[k]) / ZSize;
      imaginaryArray[0][k] = 0.0;     
    }

  double Diff = 0.0, Tmp = 0.0;
  for (int delta = 1; delta < nbrState; ++delta)
    {
      realArray[delta] = new double [nbrCylinder];
      imaginaryArray[delta] = new double [nbrCylinder];
      Diff = 2.0 * M_PI * double(delta);
      Tmp = Diff / ZSize;
      Diff = 1.0 / Diff;
      for (int k = 0; k < nbrCylinder; ++k)
	{
	  realArray[delta][k] = Diff * (sin(Tmp * ZPosition[k + 1]) - sin(Tmp * ZPosition[k]));
	  imaginaryArray[delta][k] = Diff * (cos(Tmp * ZPosition[k]) - cos(Tmp * ZPosition[k + 1]));     
	}
    }
  return true;
}
