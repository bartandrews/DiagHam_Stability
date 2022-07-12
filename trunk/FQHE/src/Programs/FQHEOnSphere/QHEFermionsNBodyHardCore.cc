#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

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
  OptionManager Manager ("QHEFermionsNBodyHardCore" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);  
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  OptionGroup* LanczosGroup = Manager.GetOptionGroup("Lanczos options");
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-nbody", "number of particle that can interact simultaneously through the n-body hard-core interaction", 2);
  (*SystemGroup) += new SingleStringOption ('\n', "nbody-file", "file describing which n-body hard-core interactions have to be used");
  (*SystemGroup) += new BooleanOption ('\n', "add-impurities", "add two impurities (one at each pole)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "north-potential", "potential assosciated to the impurity at the north pole", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "south-potential", "potential assosciated to the impurity at the south pole", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL, only useful when adding impurities)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

//   (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
//   (*LanczosGroup) += new SingleIntegerOption  ('\n', "full-diag", 
// 						"maximum Hilbert space dimension for which full diagonalization is applied", 
// 						500, true, 100);

//   (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
//   (*LanczosGroup) += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
//   (*LanczosGroup) += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
//   (*LanczosGroup) += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
//   (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
//   (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
//   (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
//   (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
//   (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
// 					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
//   (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
//   (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
//   (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
//   (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
//   (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new  SingleDoubleOption ('\n', "energy-shift", "apply a given shift to all energies while doing Lanczos iterations (final values do not include this shift)" , 0.0);
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsNBodyHardCore -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrNBody = ((SingleIntegerOption*) Manager["nbr-nbody"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();  
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  bool FirstRun = true;
  double* NBodyWeightFactors = 0;
  double* PseudoPotentials = 0;
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
      int TmpNbrPseudoPotentials;
      if (NBodyDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
    }

  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_hardcore_nbody_%d_n_%d_2s_%d_lz.dat", NbrNBody, NbrFermions, LzMax);

  int Max = ((LzMax - NbrFermions + 1) * NbrFermions);

  int  L =InitialLz;
  if (L < -Max)
    L = -Max;
  else
    if (L > Max)
      L = Max;
  if ((abs(Max) & 1) != (abs(InitialLz) & 1))
    L += 1;
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
       ParticleOnSphere* Space;
//        if (((BooleanOption*) Manager["add-impurities"])->GetBoolean() == false)
// 	 {
#ifdef __64_BITS__
      if (LzMax <= 63)
        {
	  if ((SymmetrizedBasis == false) || (L != 0))
	    Space = new FermionOnSphere(NbrFermions, L, LzMax, MemorySpace);
	  else
	    Space = new FermionOnSphereSymmetricBasis(NbrFermions, LzMax, MemorySpace);	  
        }
      else
        {
          Space = new FermionOnSphereUnlimited(NbrFermions, L, LzMax, MemorySpace);
        }
#else
      if (LzMax <= 31)
        {
	  if ((SymmetrizedBasis == false) || (L != 0))
	    Space = new FermionOnSphere(NbrFermions, L, LzMax, MemorySpace);
	  else
	    Space = new FermionOnSphereSymmetricBasis(NbrFermions, LzMax, MemorySpace);	  
        }
      else
        {
          Space = new FermionOnSphereUnlimited(NbrFermions, L, LzMax, MemorySpace);
        }
#endif
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (((BooleanOption*) Manager["add-impurities"])->GetBoolean() == false)
	{
	  if (PseudoPotentials == 0)
	    if (NBodyWeightFactors == 0)
	      {
		Hamiltonian = new ParticleOnSphereNBodyHardCoreHamiltonian(Space, NbrFermions, LzMax, NbrNBody, 
									   ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), 
									   Architecture.GetArchitecture(), 
									   Memory, DiskCacheFlag,
									   LoadPrecalculationFileName);
	      }
	    else
	      {
		Hamiltonian = new ParticleOnSphereNBodyHardCoreHamiltonian(Space, NbrFermions, LzMax, NbrNBody, NBodyWeightFactors,
									   ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), 
									   Architecture.GetArchitecture(), 
									   Memory, DiskCacheFlag,
									   LoadPrecalculationFileName);
	      }
	  else
	    Hamiltonian = new ParticleOnSphereNBodyHardCoreHamiltonian(Space, NbrFermions, LzMax, NbrNBody, NBodyWeightFactors,
								       ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), 
								       PseudoPotentials,
								       Architecture.GetArchitecture(), 
								       Memory, DiskCacheFlag,
								       LoadPrecalculationFileName);
	}
      else
	{
	  if (NBodyWeightFactors == 0)
	    {
	      Hamiltonian = new ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(Space, NbrFermions, LzMax, NbrNBody, 
											  ((SingleDoubleOption*) Manager["north-potential"])->GetDouble(), 
											  ((SingleDoubleOption*) Manager["south-potential"])->GetDouble(),
											  ((SingleIntegerOption*) Manager["landau-level"])->GetInteger(),
											  ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), 
											  Architecture.GetArchitecture(), 
											  Memory, DiskCacheFlag,
											  LoadPrecalculationFileName);
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(Space, NbrFermions, LzMax, NbrNBody, NBodyWeightFactors,
											  ((SingleDoubleOption*) Manager["north-potential"])->GetDouble(), 
											  ((SingleDoubleOption*) Manager["south-potential"])->GetDouble(),
											  ((SingleIntegerOption*) Manager["landau-level"])->GetInteger(),
											  ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), 
											  Architecture.GetArchitecture(), 
											  Memory, DiskCacheFlag,
											  LoadPrecalculationFileName);
	    }
	}
      double Shift = - 0.5 * ((double) (NbrFermions * NbrFermions)) / (0.5 * ((double) LzMax)) + ((SingleDoubleOption*) Manager["energy-shift"])->GetDouble();
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [256];
	  sprintf (EigenvectorName, "fermions_hardcore_nbody_%d_n_%d_2s_%d_lz_%d", NbrNBody, NbrFermions, LzMax, L);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
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
