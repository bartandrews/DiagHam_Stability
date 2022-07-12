#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/OverlapSpectra.h"
#include "Tools/Spectra/AverageSpectra.h"
#include "Tools/Spectra/DOSSpectra.h"
#include "Tools/Spectra/PeriodicSpectra.h"
#include "Tools/Spectra/XYReflexionSymmetricPeriodicSpectra.h"
#include "Tools/Spectra/HardBoxSpectra.h"

#include "HilbertSpace/Confined3DOneParticle.h"
#include "HilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/PairXPairYPeriodic3DOneParticle.h"

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
  OptionManager Manager ("PeriodicOscillatorForce" , "0.01");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* FileGroup =  new OptionGroup ("File options");
  OptionGroup* SampleGroup =  new OptionGroup ("Sample options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += MiscGroup;

  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowx", "lower impulsion in x direction", -40);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowy", "lower impulsion in y direction", -40);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);

  (*FileGroup) += new SingleStringOption('\n', "file1", "name of the file containing the first state", "eigenvector.0");
  (*FileGroup) += new SingleStringOption('\n', "file2", "name of the file containing the second state", "eigenvector.0");    
  (*FileGroup) += new SingleStringOption('\n', "output", "prefix of the output file", "Polarization.txt");
  (*FileGroup) += new SingleDoubleOption('e', "energy", "energy difference between two states (in eV unit)", 1.0);

  (*SampleGroup) += new SingleDoubleOption('X', "sizeX", "size of the sample (in Angstrom unit)", 500);
  (*SampleGroup) += new SingleDoubleOption('Y', "sizeY", "size of the sample (in Angstrom unit)", 500);
  (*SampleGroup) += new SingleDoubleOption('Z', "sizeZ", "size of the sample (in Angstrom unit)", 500);



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
  int LowImpulsionX = ((SingleIntegerOption*) Manager["lowx"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int LowImpulsionY = ((SingleIntegerOption*) Manager["lowy"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();

  char* FileName1 = ((SingleStringOption*) Manager["file1"])->GetString();
  char* FileName2 = ((SingleStringOption*) Manager["file2"])->GetString();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();
  double Energy = ((SingleDoubleOption*) Manager["energy"])->GetDouble();

  double SizeX = ((SingleDoubleOption*) Manager["sizeX"])->GetDouble();
  double SizeY = ((SingleDoubleOption*) Manager["sizeY"])->GetDouble();
  double SizeZ = ((SingleDoubleOption*) Manager["sizeZ"])->GetDouble();

   // define Hilbert space  
   Periodic3DOneParticle* Space = new Periodic3DOneParticle(NbrStateX, LowImpulsionX, NbrStateY, LowImpulsionY, NbrStateZ, LowImpulsionZ);

  PeriodicSpectra spectra(Space, FileName1);
  double ReX, ImX, ReY, ImY, ReZ, ImZ;
  spectra.GetImpulsion(FileName2, SizeX, SizeY, SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);

  ofstream OS;
  OS.open(OutputFile,ios::out | ios::app);
  if (!OS.is_open())
    {
      cout << "Error in open the file: " << OutputFile << "Exit now" << endl;
      exit(0);
    }
  OS << Energy << '\t' << ((ReX * ReX) + (ImX * ImX)) / Energy << '\t' << ((ReY * ReY) + (ImY * ImY)) / Energy << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / Energy << endl;
  OS.close();

  return 1;
}
