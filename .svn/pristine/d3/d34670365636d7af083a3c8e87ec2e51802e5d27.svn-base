#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/FermionOnTorus.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <limits>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusChangeAngle" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
	
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that corresponds to the inpit state");
  (*SystemGroup) += new SingleDoubleOption ('\n', "input-angle", "angle of the initial torus (in radian)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusChangeAngle -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file should be provided. See man page for option syntax or type FQHETorusChangeAngle -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-state")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-state") << endl;
    }  

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0; 
  int TotalKy = 0;
  bool Statistics = true;
  ParticleOnTorus* Space = 0;
  double OutputRatio = 1.0;

  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						  NbrParticles, NbrFluxQuanta, TotalKy, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }
  if (Statistics == true)
    {
      Space = new FermionOnTorus(NbrParticles, NbrFluxQuanta, TotalKy);
    }
  else
    { 
      Space = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, TotalKy);
    }

  ComplexVector State;
  if (State.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }
  char* OutputFileName;
  if (Manager.GetString("output-file") == 0)
    {
      char* TmpPos = strstr(Manager.GetString("input-state"), "_ratio_");
      if (TmpPos == 0)
	{
	  cout << "can't deduce output file name from " << Manager.GetString("input-state") << " (should contain _anagle_)" << endl;
	}
      char* TmpPos2 = strstr(TmpPos, "_k");
      if (TmpPos2 == 0)
	{
	  cout << "can't deduce output file name from " << Manager.GetString("input-state") << " (should contain _k)" << endl;
	}      
      OutputFileName = new char[strlen(Manager.GetString("input-state")) + 256];
      char TmpChar = TmpPos[0];
      TmpPos[0] = '\0';
      sprintf (OutputFileName, "%s_ratio_%.6f%s", Manager.GetString("input-state"), OutputRatio, TmpPos2);
      TmpPos[0] = TmpChar;
    }
  else
    {
      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }

  ComplexVector OutputState = Space->ChangeTorusAngle(State, Manager.GetDouble("input-angle"), 0.0, 0.0, 0.0, 0, Space->GetHilbertSpaceDimension());
  
  if (OutputState.WriteVector (OutputFileName) == false)
    {
      cout << "can't write vector file " << OutputFileName << endl;
      return -1;      
    }
  return 0;
}

