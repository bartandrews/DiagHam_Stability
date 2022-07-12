#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd1DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Tools/FTITightBinding/TightBindingModelSSH.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardKitaevHeisenbergModel" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-unitcells", "total number of unit cells", 4);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "delta", "hopping anisotropy",0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "cylinder", "use open boundary conditions");
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation" , "don't use translation symmetry");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");

  (*ToolsGroup) += new SingleDoubleOption ('\n', "testhermitian-error", "error threshold when testing hermiticy (0 for machine accuracy)", 0.0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardKitaevHeisenbergModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrUnitCells = Manager.GetInteger("nbr-unitcells"); 
  bool CylinderFlag = Manager.GetBoolean("cylinder");
  bool NoTranslationFlag  = Manager.GetBoolean("no-translation");
  char* StatisticPrefix = new char [64];
  sprintf (StatisticPrefix, "fermions_ssh");
  
  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_x_%d_n_%d_ns_%d", StatisticPrefix, NbrUnitCells, NbrParticles, (2 * NbrUnitCells));
  
  char* FileParameterString = new char [256];
  if (CylinderFlag)
    sprintf (FileParameterString, "cylinder_d_%g", Manager.GetDouble("delta"));
  else
    sprintf (FileParameterString, "d_%g", Manager.GetDouble("delta"));
  
  char* CommentLine = new char [256];
  sprintf (CommentLine, "");

  
  char* EigenvalueOutputFile = new char [512];
  sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);
  
  
  bool FirstRunFlag = true;

  TightBindingModelSSH TightBindingModel (NbrUnitCells,  Manager.GetDouble("delta"), CylinderFlag, Architecture.GetArchitecture(), true);
  HermitianMatrix TightBindingModelMatrix = TightBindingModel.GetRealSpaceTightBindingHamiltonian();
  
  RealSymmetricMatrix DensityDensity (2 * NbrUnitCells, true); 
  int MaxMomentum = NbrUnitCells;
  if ((CylinderFlag)|| (NoTranslationFlag))
    MaxMomentum = 1;
  
  for(int Momentum = 0; Momentum < MaxMomentum ; Momentum++)
    {
  
      long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
      ParticleOnSphere* Space = 0;
      AbstractHamiltonian* Hamiltonian = 0;

      if ((CylinderFlag)|| (NoTranslationFlag))
	Space = new   FermionOnLatticeRealSpace (NbrParticles, 2*NbrUnitCells); 
      else
	Space = new  FermionOnLatticeRealSpaceAnd1DTranslation (NbrParticles, 2*NbrUnitCells, Momentum, MaxMomentum);
      
      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      
      if ((CylinderFlag)|| (NoTranslationFlag))
	Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, 2 * NbrUnitCells, TightBindingModelMatrix,DensityDensity, Architecture.GetArchitecture(), Memory);
      else
	Hamiltonian = new ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian (Space, NbrParticles, 2 * NbrUnitCells, Momentum, MaxMomentum, TightBindingModelMatrix,DensityDensity,Architecture.GetArchitecture(),  Memory);
      
      char* ContentPrefix = new char[256];
      if ((CylinderFlag)|| (NoTranslationFlag))
	{
	  sprintf (ContentPrefix, "");
	}
      else
	{
	  sprintf (ContentPrefix, "%d",Momentum);
	}
      
      char* EigenstateOutputFile;
      char* TmpExtention = new char [512];
      if ((CylinderFlag)|| (NoTranslationFlag))
	{
	  sprintf (TmpExtention, "");
	}
      else
	sprintf (TmpExtention, "_k_%d", Momentum);
      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
      
      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
      FirstRunFlag = false;
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      cout << "------------------------------------" << endl;
      delete Hamiltonian;
      delete Space;
      delete[] EigenstateOutputFile;
      delete[] ContentPrefix;
    }
  return 0;
}


