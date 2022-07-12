#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "HilbertSpace/PlanarRotationSymmetryZPeriodicOneParticle.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/CylinderQuantumDotSpectra.h"

#include "Tools/Potential/QuantumDotThreeDConstantCylinderPotential.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::endl;
using std::ofstream;
using std::ios;

int main(int argc, char** argv)
{  
  cout.precision(14);
  OptionManager Manager ("DotProbability" , "0.01");
  OptionGroup* PotentialGroup = new OptionGroup ("potential options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += PotentialGroup;
  Manager += HilbertSpaceGroup;
  Manager += LanczosGroup;
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

  (*HilbertSpaceGroup) += new SingleIntegerOption ('R', "R-states", "number of states in plane", 50);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('Z', "Z-states", "number of cells in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  (*HilbertSpaceGroup) += new SingleDoubleOption ('k', "wave", "wave vector of Bloch function in Z direction (in 1/Angstrom unit)", 0.0);

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new SingleStringOption ('\n', "input", "file input", ""); 

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

  int NbrStateR = ((SingleIntegerOption*) Manager["R-states"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["Z-states"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  int NumberM = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  double WaveVector = ((SingleDoubleOption*) Manager["wave"])->GetDouble();

  char* FileName = ((SingleStringOption*) Manager["input"])->GetString();

  QuantumDotThreeDConstantCylinderPotential* potential = new QuantumDotThreeDConstantCylinderPotential(Below, WettingWidth, DotNbr, DotHeight, BaseRadius, TopRadius, Above, Barrier, SuperCylinderRadius);

  PlanarRotationSymmetryZPeriodicOneParticle* Space = new PlanarRotationSymmetryZPeriodicOneParticle (NumberM, NbrStateR, NbrStateZ, LowImpulsionZ); 

  CylinderQuantumDotSpectra* spectra = new CylinderQuantumDotSpectra(Space, FileName, 0.0);
  cout << "The probability to find the particle in the dot is: " << spectra->GetDotProbability(potential) << endl;

  return 1;
}
