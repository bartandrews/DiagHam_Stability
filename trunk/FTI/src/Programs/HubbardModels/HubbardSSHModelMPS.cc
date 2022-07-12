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
  OptionManager Manager ("HubbardSSHModelMPS" , "0.01");
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
  //  (*SystemGroup) += new BooleanOption  ('\n', "cylinder", "use open boundary conditions");
  //  (*SystemGroup) += new BooleanOption  ('\n', "no-translation" , "don't use translation symmetry");

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSSHModelMPS -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrUnitCells = Manager.GetInteger("nbr-unitcells"); 

  char* StatisticPrefix = new char [64];
  sprintf (StatisticPrefix, "fermions_ssh_mps");
  
  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_x_%d_n_%d_ns_%d", StatisticPrefix, NbrUnitCells, NbrParticles, (2 * NbrUnitCells));
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "d_%.6f", 1.0);

  char* MPSOutputFile = new char [512];
  sprintf(MPSOutputFile, "%s_%s.0.vec", FilePrefix, FileParameterString);
  
  FermionOnLatticeRealSpace* Space = new FermionOnLatticeRealSpace (NbrParticles, 2*NbrUnitCells); 

  Complex TmpI (0.0, 1.0);
  ComplexMatrix* MPSMatrices = new ComplexMatrix[4];
  MPSMatrices[0] = ComplexMatrix(2, 2, true);
  MPSMatrices[0].SetMatrixElement(0, 1, 1.0);
  MPSMatrices[2] = ComplexMatrix(2, 2, true);
  MPSMatrices[2].SetMatrixElement(0, 0, 1.0);
  MPSMatrices[1] = ComplexMatrix(2, 2, true);
  MPSMatrices[1].SetMatrixElement(1, 1, -1.0);
  MPSMatrices[3] = ComplexMatrix(2, 2, true);
  MPSMatrices[3].SetMatrixElement(1, 0, -1.0);
  unsigned long* TmpOccupation = new unsigned long [2 * NbrUnitCells];
  ComplexVector TmpState (Space->GetHilbertSpaceDimension(), true);
  //  ComplexVector TmpState2;
  //  TmpState2.ReadVector("fermions_ssh_x_4_n_4_ns_8_d_1.0.vec");
  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
    {
      ComplexMatrix TmpMatrix (2, 2, true);
      TmpMatrix.SetToIdentity();
      Space->GetOccupationNumber(i, TmpOccupation);
      for (int j = 0; j < NbrUnitCells; ++j)
	{
	  TmpMatrix.Multiply(MPSMatrices[TmpOccupation[2 * j] + (2 * TmpOccupation[2 * j + 1])]);
	}
//       if (SqrNorm(TmpState2[i]) > 1.0e-10)
// 	{
	  for (int j = 0; j < NbrUnitCells; ++j)
	    {
	      cout << TmpOccupation[2 * j] << " " << TmpOccupation[2 * j + 1] << " | "; 
	    }
	  cout << " : " << TmpMatrix.Tr();
// 	  cout << " " << TmpState2[i];
// 	  if (SqrNorm(0.25 * TmpMatrix.Tr() - TmpState2[i]) > 1e-10)
// 	    {
// 	      cout << " error";
// 	    }
	  cout << endl;
	  //	}
      TmpState[i] = TmpMatrix.Tr();
    }
  TmpState /= TmpState.Norm();
  TmpState.WriteVector(MPSOutputFile);
  return 0;
}


