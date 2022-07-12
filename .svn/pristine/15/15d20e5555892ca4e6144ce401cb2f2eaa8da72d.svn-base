#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereFuseStates" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-vectors", "file that describes states to fuse");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "padding", "number of empty one-body states to insert between two fused Hilbert spaces", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the fused vector that will be generated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  int Padding = Manager.GetInteger("padding");
  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(((SingleStringOption*) Manager["input-states"])->GetString()) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }
  
  int LeftNbrParticles = 0;
  int LeftLzMax = 0;
  int LeftTotalLz = 0;
  bool Statistics = false;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(0, 0),
						   LeftNbrParticles, LeftLzMax, LeftTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << InputVectors(0, 0) << endl;
      return -1;
    }
  int RightNbrParticles = 0;
  int RightLzMax = 0;
  int RightTotalLz = 0;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(InputVectors(0, 1),
						   RightNbrParticles, RightLzMax, RightTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from left state name " << InputVectors(0, 1) << endl;
      return -1;
    }

  NbrParticles = RightNbrParticles + LeftNbrParticles;
  LzMax = RightLzMax + LeftLzMax + Padding;
  TotalLz = 0;
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [256];
      sprintf (OutputFileName, "fermions_fused_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
    }

  FermionOnSphere* OutputBasis = new FermionOnSphere(NbrParticles, TotalLz, LzMax);

  RealVector OutputState(OutputBasis->GetHilbertSpaceDimension(), true);
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }	  

  return 0;
}
