#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnDisk.h"
#include "Hamiltonian/ParticleOnDiskDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnDiskMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEBosonsDiskDelta" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("Tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

 (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
 (*SystemGroup) += new SingleIntegerOption  ('l', "maximum-momentum", "maximum single particle momentum to study", 10, true, 1);
 (*SystemGroup) += new SingleIntegerOption  ('\n', "minimum-momentum", "minimum single particle momentum to study", 1, true, 1);

  (*LanczosGroup) += new SingleIntegerOption ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup) += new SingleIntegerOption ('\n', "full-diag", 
					      "maximum Hilbert space dimension for which full diagonalization is applied", 
					      500, true, 100);

  (*LanczosGroup) += new SingleIntegerOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDiskDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["iter-max"])->GetInteger();
  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();
  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MMin = ((SingleIntegerOption*) Manager["minimum-momentum"])->GetInteger();
  int MMax = ((SingleIntegerOption*) Manager["maximum-momentum"])->GetInteger();
  if (MMax < MMin)
    MMax = MMin;
  int FullDiagonalizationLimit = ((SingleIntegerOption*) Manager["full-diag"])->GetInteger();
  bool FirstRun = true;

  char* OutputNameLz = new char [1024];
  sprintf (OutputNameLz, "bosons_disk_delta_n_%d_lz_%d.dat", NbrBosons, MMax);

  for (int  L = MMin; L <= MMax; ++L)
    {
      BosonOnDisk Space (NbrBosons, L);
      cout << "Nbr bosons = " << NbrBosons << "    L = " << L << "    Dimension = " << Space.GetHilbertSpaceDimension() << endl;
      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
      ParticleOnDiskDeltaHamiltonian* Hamiltonian = new ParticleOnDiskDeltaHamiltonian(&Space, L);
      double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) MMax));
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "bosons_disk_delta_n_%d_lz_%d", NbrBosons, L);
	}
      QHEOnDiskMainTask Task (&Manager, &Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
    }
  return 0;
}

