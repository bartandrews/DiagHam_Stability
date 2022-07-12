#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/DOSSpectra.h"
#include "Tools/Spectra/OverlapSpectra.h"
#include "Tools/Spectra/AverageSpectra.h"
#include "Tools/Spectra/TimeResolvedPLSpectra.h"
#include "Tools/Spectra/CylinderInMagneticFieldSpectra.h"
#include "Tools/Spectra/CylinderQuantumDotSpectra.h"

#include "HilbertSpace/PlanarRotationSymmetryZPeriodicOneParticle.h"

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
  BooleanOption HelpOption ('h', "help", "display this help");
  SingleStringOption InputFile('\n', "input", "name of the input file", 0);
  SingleIntegerOption NumberRValueOption ('R', "R-states", "number of states in plane", 100);
  SingleIntegerOption NumberZValueOption ('Z', "Z-states", "number of states in z direction", 21);
  SingleIntegerOption LowZOption ('\n', "lowz", "lower impulsion in z direction", -10);
  SingleIntegerOption NumberMValueOption ('m', "momentum", "quantum number of kinetic in z direction", 0);
  SingleDoubleOption MagneticFieldOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 30);
  SingleDoubleOption SizeZOption ('z', "size-z", "size of sample in Z direction (in Angstrom unit)", 110);
  SingleDoubleOption SizeROption ('r', "size-r", "size of sample in plane (in Angstrom unit)", 1000);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &InputFile;
  OptionList += &NumberRValueOption;
  OptionList += &NumberZValueOption;
  OptionList += &LowZOption;
  OptionList += &NumberMValueOption;
  OptionList += &MagneticFieldOption;
  OptionList += &SizeZOption;
  OptionList += &SizeROption;

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

  char* FileName = InputFile.GetString();
  int NbrStateR = NumberRValueOption.GetInteger();
  int NbrStateZ = NumberZValueOption.GetInteger();
  int LowZ = LowZOption.GetInteger();
  int NumberM = NumberMValueOption.GetInteger();
  double MagneticField = MagneticFieldOption.GetDouble();
  double SizeZ = SizeZOption.GetDouble();
  double SizeR = SizeROption.GetDouble();
  PlanarRotationSymmetryZPeriodicOneParticle* space = new PlanarRotationSymmetryZPeriodicOneParticle (0, NbrStateR, NbrStateZ, LowZ);
  CylinderInMagneticFieldSpectra* spectra = new CylinderInMagneticFieldSpectra(space, FileName, MagneticField);
   
  //cout << "# " << spectra->GetSquaredRadius () << endl << endl;
  
  double delta = SizeZ / 100.0; double p = 0.0;
  double shift = 0.0;
  for (double z = shift; z <= (SizeZ + shift); z += delta)
    {
      p = spectra->ZProbabilityDensity(z, SizeZ);
      cout << (z - shift) << '\t' << p << '\n';
    }

  return 1;
}
