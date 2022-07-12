#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereWithSpinApplyOneBodyTransformationOperation.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/ComplexMatrix.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);
    
  // some running options and help
  OptionManager Manager ("FQHESphereWithSU2SpinRotation" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file that contains the input state");
  (*SystemGroup) += new SingleDoubleOption  ('a', "spin-angle", "rotation angle along the Sy axis (in pi units)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinRotation -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, a state file should be provided. See man page for option syntax or type FQHESphereWithSU2SpinRotation -h" << endl;
      return -1;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  bool Statistics = true;
  int TotalLz = 0;
  int TotalSz = 0;

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrFluxQuanta,
							   TotalLz, TotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  RealVector InputState;
  if (InputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }

  ParticleOnSphereWithSpin* InputSpace = 0;
  if (Statistics == false)
    {
      InputSpace = new BosonOnSphereWithSU2Spin (NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
    }
  else
    {
      cout << "not yet implemented for fermions" << endl;
      return 0;
    }

  if (InputSpace->GetLargeHilbertSpaceDimension() != InputState.GetLargeVectorDimension())
    {
      cout << "the input vector file " << Manager.GetString("input-state") << "does not have the correct dimension (is " 
	   << InputState.GetLargeVectorDimension() << ", should be " << InputSpace->GetLargeHilbertSpaceDimension() << ")" << endl;
      return -1;      

    }
  
  ParticleOnSphereWithSpin* OutputSpace = 0;
  if (Statistics == false)
    {
      OutputSpace = new BosonOnSphereWithSU2Spin (NbrParticles, TotalLz, NbrFluxQuanta);
    }
  else
    {
      cout << "not yet implemented for fermions, use the old code FQHESphereFermionsWithSpinRotation" << endl;
      return 0;
    }

  double Theta = Manager.GetDouble("spin-angle");
  RealVector OutputState;
  OutputState = OutputSpace->SU2ToSU2AllSz(InputState, InputSpace);
  if (Theta != 0.0)
    {
      RealMatrix* RotationMatrices;
      RotationMatrices = new RealMatrix [NbrFluxQuanta + 1];
      RealMatrix RotationMatrix(2,2);
      RotationMatrix.SetMatrixElement(0, 0, cos(0.5 * Theta * M_PI));
      RotationMatrix.SetMatrixElement(0, 1, -sin(0.5 * Theta * M_PI));
      RotationMatrix.SetMatrixElement(1, 0, sin(0.5 * Theta * M_PI));
      RotationMatrix.SetMatrixElement(1, 1, cos(0.5 * Theta * M_PI));
      for(int k = 0; k <= NbrFluxQuanta; ++k)
	RotationMatrices[k] = RotationMatrix;
      RealVector OutputState2(OutputState.GetLargeVectorDimension(), true);
      FQHESphereWithSpinApplyOneBodyTransformationOperation Operation1(&OutputState, &OutputState2, RotationMatrices, OutputSpace);
      Operation1.ApplyOperation(Architecture.GetArchitecture());
      OutputState = OutputState2;
      if (Manager.GetBoolean("szsymmetrized-basis"))
	{
	  if (Statistics == false)
	    {
	      BosonOnSphereWithSU2SpinSzSymmetry* OutputSpace2 = new BosonOnSphereWithSU2SpinSzSymmetry (NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetBoolean("minus-szparity"));
	      OutputState2 = OutputSpace2->ConvertToNbodyBasis(OutputState, OutputSpace);	      
	      OutputState = OutputState2;
	      delete OutputSpace2;
	    }
	  else
	    {
	      cout << "not yet implemented for fermions, use the old code FQHESphereFermionsWithSpinRotation" << endl;
	      return 0;
	    }
	}
    }
  char* TmpInputSzString = new char [16];
  sprintf (TmpInputSzString, "sz_%d", TotalSz);
  char* TmpOutputSzString = new char [32];
  sprintf (TmpOutputSzString, "syrotation_%.6f", Theta);
  char* OutputFileName = ReplaceString(Manager.GetString("input-state"), TmpInputSzString, TmpOutputSzString);
  if (OutputFileName == 0)
    {
      OutputFileName = new char[256];
      if (Statistics == false)
	{
	  sprintf(OutputFileName, "bosons_su2_dummy_n_%d_2s_%d_syrotation_%.6f_lz_%d.0.vec", NbrParticles, NbrFluxQuanta, Theta, TotalLz);
	}
      else
	{
	  sprintf(OutputFileName, "fermions_su2_dummy_n_%d_2s_%d_syrotation_%.6f_lz_%d.0.vec", NbrParticles, NbrFluxQuanta, Theta, TotalLz);
	}
    }
  if (Manager.GetBoolean("szsymmetrized-basis"))
    {
      char* TmpSzSymmetryString = new char [32];
      if (Manager.GetBoolean("minus-szparity"))
	{
	  sprintf (TmpSzSymmetryString, "_su2_szsym_-1_");
	}
      else
	{
	  sprintf (TmpSzSymmetryString, "_su2_szsym_1_");
	}
      char* OutputFileName2 = ReplaceString(OutputFileName, "_su2_", TmpSzSymmetryString);
      if (OutputFileName2 == 0)
	{
	  if (Statistics == false)
	    {
	      sprintf(OutputFileName2, "bosons_%s_dummy_n_%d_2s_%d_syrotation_%.6f_lz_%d.0.vec", TmpSzSymmetryString, NbrParticles, NbrFluxQuanta, Theta, TotalLz);
	    }
	  else
	    {
	      sprintf(OutputFileName2, "fermions_%s_dummy_n_%d_2s_%d_syrotation_%.6f_lz_%d.0.vec", TmpSzSymmetryString, NbrParticles, NbrFluxQuanta, Theta, TotalLz);
	    }
	}
      delete[] OutputFileName;
      OutputFileName = OutputFileName2;
    }
  OutputState.WriteVector(OutputFileName);
  delete[] OutputFileName;
  return 0;
}
