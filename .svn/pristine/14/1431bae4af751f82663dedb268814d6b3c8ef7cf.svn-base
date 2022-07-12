#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

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
  OptionManager Manager ("PeriodicWaveFunctionValue" , "0.01");
  OptionGroup* PositionGroup = new OptionGroup ("Position options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* FileGroup = new OptionGroup ("file options");

  Manager += PositionGroup;
  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += MiscGroup;

  
  (*PositionGroup) += new SingleDoubleOption ('\n', "begin-x", "beginning position in X direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleDoubleOption ('\n', "end-x", "ending position in X direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleIntegerOption ('x', "nbr-pointx", "number of points in X direction", 1);

  (*PositionGroup) += new SingleDoubleOption ('\n', "begin-y", "beginning position in Y direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleDoubleOption ('\n', "end-y", "ending position in Y direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleIntegerOption ('y', "nbr-pointy", "number of points in Y direction", 1);
  
  (*PositionGroup) += new SingleDoubleOption ('\n', "begin-z", "beginning position in Z direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleDoubleOption ('\n', "end-z", "ending position in Z direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleIntegerOption ('z', "nbr-pointz", "number of points in Z direction", 1);

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

  double BeginX = ((SingleDoubleOption*) Manager["begin-x"])->GetDouble();
  double EndX = ((SingleDoubleOption*) Manager["end-x"])->GetDouble();
  int NbrPointX = ((SingleIntegerOption*) Manager["nbr-pointx"])->GetInteger();

  double BeginY = ((SingleDoubleOption*) Manager["begin-y"])->GetDouble();
  double EndY = ((SingleDoubleOption*) Manager["end-y"])->GetDouble();
  int NbrPointY = ((SingleIntegerOption*) Manager["nbr-pointy"])->GetInteger();

  double BeginZ = ((SingleDoubleOption*) Manager["begin-z"])->GetDouble();
  double EndZ = ((SingleDoubleOption*) Manager["end-z"])->GetDouble();
  int NbrPointZ = ((SingleIntegerOption*) Manager["nbr-pointz"])->GetInteger();

  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int LowImpulsionX = ((SingleIntegerOption*) Manager["lowx"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int LowImpulsionY = ((SingleIntegerOption*) Manager["lowy"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  
  char* FileName = ((SingleStringOption*) Manager["input"])->GetString();

  PeriodicThreeDOneParticle* Space = new PeriodicThreeDOneParticle(NbrStateX, LowImpulsionX, NbrStateY, LowImpulsionY, NbrStateZ, LowImpulsionZ);

  PeriodicSpectra Spectra (Space, FileName);
  
  double deltaX = 0.0;
  if (NbrPointX != 0)
    deltaX = (EndX - BeginX) / (NbrPointX - 1);

  double deltaY = 0.0;
  if (NbrPointY != 0)
    deltaY = (EndY - BeginY) / (NbrPointY - 1);
  
  double deltaZ = 0.0;
  if (NbrPointZ != 0)
    deltaZ = (EndZ - BeginZ) / (NbrPointZ - 1);

  double PositionX = BeginX, PositionY = 0.0, PositionZ = 0.0;
  double Real = 0.0, Imaginary = 0.0; double Value = 0;
  for (int i = 0; i < NbrPointX; ++i)
    {
      PositionY = BeginY;
      for (int j = 0; j < NbrPointY; ++j)
	{
	  PositionZ = BeginZ;
	  for (int k = 0; k < NbrPointZ; ++k)
	    {
	      Spectra.WaveFunctionValue (PositionX, 1.0, PositionY, 1.0, PositionZ, 1.0, Real, Imaginary);
	      Value = Real * Real + Imaginary * Imaginary;
	      cout << PositionX << '\t' << PositionY << '\t' << PositionZ << '\t' << Real << '\t' << Imaginary << '\t' << Value << '\n';
	      PositionZ += deltaZ;
	    }
	  
	  PositionY += deltaY;
	}
      PositionX += deltaX;
    }
  
  return 1;
}
