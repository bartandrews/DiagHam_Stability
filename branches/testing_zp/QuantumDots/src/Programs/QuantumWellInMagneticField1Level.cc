#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Hamiltonian/QuantumWellHamiltonianInMagneticField1Level.h"


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
  SingleDoubleOption ZSizeOption ('\n', "zsize", "size in the z direction in Angstrom", 58.7);
  SingleDoubleOption MassOption ('\n', "mass", "electron effective mass (in bare electron mass unit)", 0.050);
  SingleDoubleOption BFieldOption ('\n', "bfield", "B field (in Tesla unit)", 30.90);
  SingleDoubleOption ZKineticEnergyOption ('\n', "z-kinetic", "kinetic energy (in meV) due to the z confinment", 0.0);
  SingleIntegerOption SubbandIndexOption ('\n', "subband-index", "subband index", 1);
  SingleIntegerOption LandauIndexOption ('\n', "landau-index", "Landau level index (lowest Landau level is 0)", 0);
  SingleDoubleOption BandOffsetOption ('\n', "band-offset", "band offset value (in meV unit)", 600);
  SingleStringOption CoefficientFileNameOption('\n', "coefficients", "name of the file where interaction coeffcients are stored", 
					       "/home/regnault/development/DMRG/DiagHam/potentiel_10_10_10_2");
  SingleStringOption SavePotentialFileNameOption('\n', "save-potential", "save potential description into a given file");
  SingleStringOption LoadPotentialFileNameOption('\n', "load-potential", "load potential description from a given file");
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
  OptionList += &ZSizeOption;
  OptionList += &MassOption;
  OptionList += &BFieldOption;
  OptionList += &ZKineticEnergyOption;
  OptionList += &LandauIndexOption;
  OptionList += &SubbandIndexOption;
  OptionList += &CoefficientFileNameOption;
  OptionList += &LeftSizeOption;
  OptionList += &RightSizeOption;
  OptionList += &MemoryOption;
  OptionList += &CarrierTypeOption; 
  OptionList += &BandOffsetOption;
  OptionList += &SavePotentialFileNameOption;
  OptionList += &LoadPotentialFileNameOption;

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
  double Lx = CellXSizeOption.GetDouble();
  double Ly = CellYSizeOption.GetDouble();
  double Lz = CellZSizeOption.GetDouble();
  double ZSize = ZSizeOption.GetDouble();
  double Mass = MassOption.GetDouble();
  double BField = BFieldOption.GetDouble();
  double ZKineticEnergy = ZKineticEnergyOption.GetDouble();
  int LandauIndex = LandauIndexOption.GetInteger();
  int SubbandIndex = SubbandIndexOption.GetInteger();
  int LeftSize = LeftSizeOption.GetInteger();
  int RightSize = RightSizeOption.GetInteger();
  bool Carrier = CarrierTypeOption.GetBoolean();
  double BandOffset = BandOffsetOption.GetDouble();
  char* SavePotentialFileName = SavePotentialFileNameOption.GetString();
  char* LoadPotentialFileName = LoadPotentialFileNameOption.GetString();

  timeval CurrentTime;
  gettimeofday (&(CurrentTime), 0);
  srand(CurrentTime.tv_usec);
  QuantumWellHamiltonianInMagneticField1Level Hamiltonian (1000.0, 1000.0, ZSize, Mass, BField,  ZKineticEnergy, LandauIndex, SubbandIndex, 5.87, BandOffset, 0.53, LoadPotentialFileName);
  if (SavePotentialFileName != 0)    
    Hamiltonian.SavePotential(SavePotentialFileName);
  cout << Hamiltonian.GetHilbertSpaceDimension() << endl;
  HermitianMatrix HamiltonianRepresentation;
  Hamiltonian.GetHamiltonian(HamiltonianRepresentation);
  char FileName[256];
  sprintf (FileName,"hamiltonian%f.mat", BField);  
  HamiltonianRepresentation.WriteMatrix(FileName);

  cout << "start diagonalization" << endl;
  RealDiagonalMatrix DiagonalizedHamiltonian (HamiltonianRepresentation.GetNbrRow());
  ComplexMatrix Eigenvectors(HamiltonianRepresentation.GetNbrRow(), HamiltonianRepresentation.GetNbrRow());
  HamiltonianRepresentation.Diagonalize(DiagonalizedHamiltonian, Eigenvectors);

  sprintf (FileName,"eigenvalues%f.raw", BField);
  char EigenvectorFileName[256];
  ofstream File0;
  File0.open(FileName, ios::out); 
  File0.precision(14); 
  for (int i = 0; i < HamiltonianRepresentation.GetNbrRow(); ++i)
    {
      File0 << DiagonalizedHamiltonian[i] << endl;
      sprintf (EigenvectorFileName,"eigenvalues%f.%d.vec", BField, i);
      if (EigenstateFlag == true)
	Eigenvectors[i].WriteVector(EigenvectorFileName);
      ComplexVector TmpVec (Eigenvectors[i].GetVectorDimension());
      TmpVec.Multiply(HamiltonianRepresentation, Eigenvectors[i]);
      cout << (Eigenvectors[i] * TmpVec) << endl;
    }
  File0.close();

  return 0;
}

