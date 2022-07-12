#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereCoulombDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEBosonsCoulomb" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "add-delta", "ratio betwen coulomd interaction and delta interaction (pure coulomb if 0)", 0.0);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);

  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsCoulomb -h" << endl;
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
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  double CoulombRatio = ((SingleDoubleOption*) Manager["add-delta"])->GetDouble();
  bool FirstRun = true;

  char* OutputNameLz = new char [256];
  if (CoulombRatio == 0.0)
    sprintf (OutputNameLz, "bosons_coulomb_n_%d_2s_%d_lz.dat", NbrBosons, LzMax);
  else
    sprintf (OutputNameLz, "bosons_coulomb_delta_%f_n_%d_2s_%d_lz.dat", CoulombRatio, NbrBosons, LzMax);
  int Max = (LzMax * NbrBosons);
  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != (InitialLz & 1))
	L += 1;
    }
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      BosonOnSphere Space (NbrBosons, L, LzMax);
      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian;
      if (CoulombRatio == 0.0)
	Hamiltonian = new ParticleOnSphereCoulombHamiltonian(&Space, NbrBosons, LzMax, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName);
      else
	Hamiltonian = new ParticleOnSphereCoulombDeltaHamiltonian(&Space, NbrBosons, LzMax, CoulombRatio, Architecture.GetArchitecture(), Memory);
      double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) LzMax));
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [128];
	  if (CoulombRatio == 0.0)
	    sprintf (EigenvectorName, "bosons_coulomb_n_%d_2s_%d_lz_%d", NbrBosons, LzMax, L);
	  else
	    sprintf (EigenvectorName, "bosons_coulomb_delta_n_%d_2s_%d_lz_%d", NbrBosons, LzMax, L);
	}
      QHEOnSphereMainTask Task (&Manager, &Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
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
