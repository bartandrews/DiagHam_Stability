#include "HilbertSpace/FullBosonOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEBosonsNBodyHardCore" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-nbody", "number of particle that can interact simultaneously through the n-body hard-core interaction", 2);
  (*SystemGroup) += new  SingleStringOption ('\n', "nbody-file", "file describing which n-body hard-core interactions have to be used");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);

  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsNBodyHardCore -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrNBody = ((SingleIntegerOption*) Manager["nbr-nbody"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();  
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  double* NBodyWeightFactors = 0;
  if (((SingleStringOption*) Manager["nbody-file"])->GetString() != 0)
    {
      ConfigurationParser NBodyDefinition;
      if (NBodyDefinition.Parse(((SingleStringOption*) Manager["nbody-file"])->GetString()) == false)
	{
	  NBodyDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((NBodyDefinition.GetAsSingleInteger("NbrNBody", NbrNBody) == false) || (NbrNBody < 2))
	{
	  cout << "NbrNBody is not defined or as a wrong value in " << ((SingleStringOption*) Manager["nbody-file"])->GetString() << endl;
	  return -1;
	}
      int TmpNbrNBody;
      if ((NBodyDefinition.GetAsDoubleArray("Weights", ' ', NBodyWeightFactors, TmpNbrNBody) == false) || ((TmpNbrNBody - NbrNBody) != 1))
	{
	  cout << "Weights is not defined or as a wrong value in " << ((SingleStringOption*) Manager["nbody-file"])->GetString() << endl;
	  return -1;
	}
    }

  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "bosons_hardcore_nbody_%d_n_%d_2s_%d_full.dat", NbrNBody, NbrBosons, LzMax);
  int Max = (LzMax * NbrBosons);

  FullBosonOnSphere Space (NbrBosons, LzMax);
  //  BosonOnSphere Space (NbrBosons, -22, LzMax);
  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
  AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
  if (NBodyWeightFactors == 0)
    {
      Hamiltonian = new ParticleOnSphereNBodyHardCoreHamiltonian(&Space, NbrBosons, LzMax, NbrNBody, 0.0, 
								 Architecture.GetArchitecture(), 
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName);
    }
  else
    {
      Hamiltonian = new ParticleOnSphereNBodyHardCoreHamiltonian(&Space, NbrBosons, LzMax, NbrNBody, NBodyWeightFactors, 0.0,
								 Architecture.GetArchitecture(), 
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName);
    }

  double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) LzMax));
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
    {
      EigenvectorName = new char [256];
      sprintf (EigenvectorName, "bosons_hardcore_nbody_%d_n_%d_2s_%d_full", NbrNBody, NbrBosons, LzMax);
    }
  QHEOnSphereMainTask Task (&Manager, &Space, Hamiltonian, 0, Shift, OutputNameLz, true, EigenvectorName);
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  delete Hamiltonian;
  if (EigenvectorName != 0)
    {
      delete[] EigenvectorName;
    }
}
