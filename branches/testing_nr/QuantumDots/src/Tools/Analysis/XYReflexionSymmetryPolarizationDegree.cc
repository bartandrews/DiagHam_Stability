#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "HilbertSpace/PeriodicXYReflexionZPeriodicThreeDOneParticle.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/XYReflexionSymmetricPeriodicSpectra.h"




#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("XYReflexionSymmetryPolarizationDegree" , "0.01");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* FileGroup =  new OptionGroup ("File and energy options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += MiscGroup;

  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairX", "pair function in X direction", false);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairY", "pair function in Y direciton", false);

  (*FileGroup) += new SingleStringOption('\n', "file1", "name of the electron state file", "");
  (*FileGroup) += new SingleStringOption('\n', "file2", "name of the hole state file", "");  
  (*FileGroup) += new SingleDoubleOption('g', "gap", "energy of the constrained gap ", 0.8);  
  (*FileGroup) += new SingleDoubleOption('X', "sizeX", "size of the sample (in Angstrom unit)", 500);
  (*FileGroup) += new SingleDoubleOption('Y', "sizeY", "size of the sample (in Angstrom unit)", 500);

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

  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  bool PairX = ((BooleanOption*) Manager["pairX"])->GetBoolean();
  bool PairY = ((BooleanOption*) Manager["pairY"])->GetBoolean();

  char* FileName1 = ((SingleStringOption*) Manager["file1"])->GetString();
  char* FileName2 = ((SingleStringOption*) Manager["file2"])->GetString();
  double Gap = ((SingleDoubleOption*) Manager["gap"])->GetDouble();
  double SizeX = ((SingleDoubleOption*) Manager["sizeX"])->GetDouble();
  double SizeY = ((SingleDoubleOption*) Manager["sizeY"])->GetDouble();

  PeriodicXYReflexionZPeriodicThreeDOneParticle* Space = new PeriodicXYReflexionZPeriodicThreeDOneParticle (NbrStateX, PairX, NbrStateY, PairY, NbrStateZ, LowImpulsionZ); 

  XYReflexionSymmetricPeriodicSpectra spectra(Space, FileName1);  

  //void GetDerivedOverlap (XYReflexionSymmetricPeriodic3DOneParticle* space, char* fileName, double sizeX, double sizeY, double sizeZ, double &realOverlap, double &imaginaryOverlap, double &realOverlapX, double &imaginaryOverlapX, double &realOverlapY, double &imaginaryOverlapY);
  
  double real, imaginary, realX, imaginaryX, realY, imaginaryY;
  spectra.GetDerivedOverlap(Space, FileName2, SizeX, SizeY, 1.0, real, imaginary, realX, imaginaryX, realY, imaginaryY);
  //cout << "Overlap: " << real << " " << imaginary << endl;
  //cout << "Overlap of X derived: " << realX << " " << imaginaryX << endl;
  //cout << "Overlap of Y derived: " << realY << " " << imaginaryY << endl;
  double re1, re2, im1, im2;
  
  re1 = real - 38.2 * (realX - realY) / (Gap * Gap);
  re2 = real + 38.2 * (realX - realY) / (Gap * Gap);
  im1 = imaginary - 38.2 * (imaginaryX - imaginaryY) / (Gap * Gap);
  im2 = imaginary + 38.2 * (imaginaryX - imaginaryY) / (Gap * Gap);

  double tmp1 = re1 * re1 + im1 * im1;
  double tmp2 = re2 * re2 + im2 * im2;
  //cout << tmp1 << " " << tmp2 << endl;
  cout << "Polarization degree is: " << ((tmp1 - tmp2) / (tmp1 + tmp2)) << endl;

  return 1;
}
