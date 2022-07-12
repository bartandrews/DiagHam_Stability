#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSpinAllSzAndMagneticTranslations.h"

#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("FQHETorusWithSU2SpinComputeSx" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file that contains the input state");
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusWithSU2SpinComputeSx -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, a state file should be provided. See man page for option syntax or type FQHETorusWithSU2SpinComputeSx -h" << endl;
      return -1;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  bool Statistics = true;
  int KxMomentum = 0;
  int KyMomentum = 0;
  double XRatio = Manager.GetDouble("ratio");

  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrFluxQuanta,
						  KxMomentum, KyMomentum, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  double** PseudoPotentials  = new double*[3];
  double** OneBodyPseudoPotentials  = new double*[3];
  int* NbrPseudoPotentials  = new int[3];

  NbrPseudoPotentials[0] = 0;
  NbrPseudoPotentials[1] = 0;
  NbrPseudoPotentials[2] = 0;
  OneBodyPseudoPotentials[0] = 0;
  OneBodyPseudoPotentials[1] = 0;
  OneBodyPseudoPotentials[2] = new double[NbrFluxQuanta];
  for (int i = 0; i < NbrFluxQuanta; ++i)
    OneBodyPseudoPotentials[2][i] = 1.0;
  
  ComplexVector InputState;
  if (InputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }

  ParticleOnTorusWithSpinAndMagneticTranslations* Space = new BosonOnTorusWithSpinAllSzAndMagneticTranslations (NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian(Space, NbrParticles, NbrFluxQuanta, KxMomentum, XRatio,
													     NbrPseudoPotentials[0], PseudoPotentials[0],
													     NbrPseudoPotentials[1], PseudoPotentials[1],
													     NbrPseudoPotentials[2], PseudoPotentials[2],
													     0.0, 0.0,
													     Architecture.GetArchitecture(), 0l, 0,
													     OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], 
													     OneBodyPseudoPotentials[2]);



  ComplexVector OutputState (InputState.GetLargeVectorDimension());
  VectorHamiltonianMultiplyOperation Operation1 (Hamiltonian, &InputState, &OutputState);
  Operation1.ApplyOperation(Architecture.GetArchitecture());

  
  char* OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "sx.dat");
  Complex Tmp = OutputState * InputState;
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);
  File << " # kx ky Re(<2S_x>) <2S_x>" << endl;
  File << KxMomentum << " " << KyMomentum << " " << Tmp.Re << " " << Tmp << endl;
  File.close();
  delete[] OutputFileName;
  return 0;
}
