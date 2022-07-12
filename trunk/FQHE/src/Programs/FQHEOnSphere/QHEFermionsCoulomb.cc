#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "Hamiltonian/ParticleOnSphereCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereCoulombDeltaHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <stdio.h>
#include <climits>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("QHEFermionsCoulomb" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 10);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "add-delta", "add a delta interaction component", false);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "ratio", "ratio between coulomd interaction and delta interaction", 1.0);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 300, 
						true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
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
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  if (LONG_MAX>>20 < Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;  
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  bool DeltaFlag = ((BooleanOption*) Manager["add-delta"])->GetBoolean();
  double CoulombRatio = ((SingleIntegerOption*) Manager["ratio"])->GetInteger();
  bool FirstRun = true;


  char* OutputNameLz = new char [256];
  if (DeltaFlag == false)
    sprintf (OutputNameLz, "fermions_coulomb_n_%d_2s_%d_lz.dat", NbrFermions, LzMax);
  else
    sprintf (OutputNameLz, "fermions_coulomb_delta_%f_n_%d_2s_%d_lz.dat", CoulombRatio, NbrFermions, LzMax);

  int Max = ((LzMax - NbrFermions + 1) * NbrFermions);

  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != 0)
	L |= 1;
      else
	L &= ~0x1;
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
      cout << "----------------------------------------------------------------" << endl;
      cout << " LzTotal = " << L << endl;
      ParticleOnSphere* Space;
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  Space = new FermionOnSphere(NbrFermions, L, LzMax, MemorySpace);	  
	}
      else
	{
	  Space = new FermionOnSphereUnlimited(NbrFermions, L, LzMax, MemorySpace);	  
	}	
#else
      if (LzMax <= 31)
	{
	  Space = new FermionOnSphere(NbrFermions, L, LzMax, MemorySpace);	  
	}
      else
	{
	  Space = new FermionOnSphereUnlimited(NbrFermions, L, LzMax, MemorySpace);	  
	}	
#endif
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian;
      if (DeltaFlag == false)
	Hamiltonian = new ParticleOnSphereCoulombHamiltonian(Space, NbrFermions, LzMax, Architecture.GetArchitecture(), Memory, 
							     LoadPrecalculationFileName);
      else
	Hamiltonian = new ParticleOnSphereCoulombDeltaHamiltonian(Space, NbrFermions, LzMax, CoulombRatio, 
								  Architecture.GetArchitecture(), Memory);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [64];
	  if (DeltaFlag == false)
	    sprintf (EigenvectorName, "fermions_coulomb_n_%d_2s_%d_lz_%d", NbrFermions, LzMax, L);
	  else
	    sprintf (EigenvectorName, "fermions_coulomb_delta_%f_n_%d_2s_%d_lz_%d", CoulombRatio, NbrFermions, LzMax, L);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, 0.0, OutputNameLz, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete Space;
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
    }

  return 0;
}
