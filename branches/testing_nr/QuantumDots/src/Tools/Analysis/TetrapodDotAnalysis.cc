#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Potential/TetrapodThreeDConstantCellPotential.h"

#include "Tools/Spectra/PeriodicSpectra.h"


#include <iostream>
#include <stdlib.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{  
  cout.precision(14);
  OptionManager Manager ("TetrapodDotAnalysis" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* FileGroup = new OptionGroup ("file options");

  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += MiscGroup;

  (*PotentialGroup) += new SingleIntegerOption ('M', "M-cell", "number of cells in the x direction", 100);
  (*PotentialGroup) += new SingleIntegerOption ('N', "N-cell", "number of cells in the y direction", 100);
  (*PotentialGroup) += new SingleIntegerOption ('H', "H-cell", "number of cells in the z direction", 100);
  (*PotentialGroup) += new SingleDoubleOption ('c', "cell-size", "cell size in Angstrom", 5);
  (*PotentialGroup) += new SingleIntegerOption ('\n', "below", "height of the barrier just below the tetrapod (in cell unit)", 0);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "dot-radius", "radius of the spherical dot in Angstrom unit", 40);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "arm-length", "length of the four arms in Angstrom unit", 80);
  (*PotentialGroup) += new SingleDoubleOption ('\n', "arm-radius", "radius of the arm in Angstrom unit", 20);

  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction", 31);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowx", "lower impulsion in x direction", -15);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction", 31);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowy", "lower impulsion in y direction", -15);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 31);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -15);

  (*FileGroup) += new SingleStringOption ('f', "input", "the name of the input file", "eigenvector.0");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

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
  double CellSize = ((SingleDoubleOption*) Manager["cell-size"])->GetDouble();
  int BelowTetrapod = ((SingleIntegerOption*) Manager["below"])->GetInteger();
  double DotRadius = ((SingleDoubleOption*) Manager["dot-radius"])->GetDouble();
  double ArmLength = ((SingleDoubleOption*) Manager["arm-length"])->GetDouble();
  double ArmRadius = ((SingleDoubleOption*) Manager["arm-radius"])->GetDouble();

  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int LowImpulsionX = ((SingleIntegerOption*) Manager["lowx"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int LowImpulsionY = ((SingleIntegerOption*) Manager["lowy"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  
  char* FileName = ((SingleStringOption*) Manager["input"])->GetString();

  TetrapodThreeDConstantCellPotential* Potential = new TetrapodThreeDConstantCellPotential (M, N, H, BelowTetrapod, DotRadius, ArmLength, ArmRadius, CellSize);

  Potential->ConstructPotential(0.0);

  PeriodicThreeDOneParticle* Space = new PeriodicThreeDOneParticle(NbrStateX, LowImpulsionX, NbrStateY, LowImpulsionY, NbrStateZ, LowImpulsionZ);

  PeriodicSpectra Spectra (Space, FileName);

  double Sphere = 0.0, Arm = 0.0;
  
  Spectra.GetTetrapodProbability (Potential, Sphere, Arm);
  
  cout << "The probability in the sphere: " << Sphere << endl;
  cout << "The probability in the arms: " << Arm << endl;

  return 1;
}
