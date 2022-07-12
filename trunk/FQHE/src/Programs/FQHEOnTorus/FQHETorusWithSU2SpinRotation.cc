#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSpinAllSzAndMagneticTranslations.h"

#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

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

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

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
  OptionManager Manager ("FQHETorusWithSU2SpinRotation" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusWithSU2SpinRotation -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, a state file should be provided. See man page for option syntax or type FQHETorusWithSU2SpinRotation -h" << endl;
      return -1;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  bool Statistics = true;
  int KxMomentum = 0;
  int KyMomentum = 0;
  int TotalSz = 0;
  double XRatio = Manager.GetDouble("ratio");

  if (FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrFluxQuanta,
							  KxMomentum, KyMomentum, TotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  ComplexVector InputState;
  if (InputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }

  ParticleOnTorusWithSpinAndMagneticTranslations* InputSpace = 0;
  if (Statistics == false)
    {
      InputSpace = new BosonOnTorusWithSpinAndMagneticTranslations (NbrParticles, TotalSz, NbrFluxQuanta, KxMomentum, KyMomentum);
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
  
  ParticleOnTorusWithSpinAndMagneticTranslations* OutputSpace = 0;
  if (Statistics == false)
    {
      OutputSpace = new BosonOnTorusWithSpinAllSzAndMagneticTranslations (NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
    }
  else
    {
      cout << "not yet implemented for fermions" << endl;
      return 0;
    }

  double Theta = Manager.GetDouble("spin-angle");
  ComplexVector OutputState;
  OutputState = OutputSpace->SU2ToSU2AllSz(InputState, InputSpace);
  if (Theta != 0.0)
    {
      ComplexMatrix* RotationMatrices;
      RotationMatrices = new ComplexMatrix [NbrFluxQuanta];
      ComplexMatrix RotationMatrix(2,2);
      RotationMatrix.SetMatrixElement(0, 0, cos(0.5 * Theta * M_PI));
      RotationMatrix.SetMatrixElement(0, 1, -sin(0.5 * Theta * M_PI));
      RotationMatrix.SetMatrixElement(1, 0, sin(0.5 * Theta * M_PI));
      RotationMatrix.SetMatrixElement(1, 1, cos(0.5 * Theta * M_PI));
      for(int k = 0; k < NbrFluxQuanta; ++k)
	RotationMatrices[k] = RotationMatrix;
      ComplexVector OutputState2(OutputState.GetLargeVectorDimension(), true);
//      OutputSpace->TransformOneBodyBasis(OutputState, OutputState2, RotationMatrices);  
      FQHESphereWithSpinApplyOneBodyTransformationOperation Operation1(&OutputState, &OutputState2, RotationMatrices, OutputSpace);
      Operation1.ApplyOperation(Architecture.GetArchitecture());
      OutputState = OutputState2;
    }
  char* TmpInputSzString = new char [16];
  sprintf (TmpInputSzString, "sz_%d", TotalSz);
  char* TmpOutputSzString = new char [32];
  sprintf (TmpOutputSzString, "syrotation_%.6f", Theta);
  char* OutputFileName = ReplaceString(Manager.GetString("input-state"), TmpInputSzString, TmpOutputSzString);
  OutputState.WriteVector(OutputFileName);
  delete[] OutputFileName;
  return 0;
}
