#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"

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


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsTwoBodyGenericHaldaneBasis" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "if non zero, override energy shift using the indicated value ", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

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
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
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


  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int InitialLz = 0;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();  
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
  bool FirstRun = true;
  double* PseudoPotentials = 0;
  double* OneBodyPotentials = 0;
  if (((SingleStringOption*) Manager["interaction-file"])->GetString() == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(((SingleStringOption*) Manager["interaction-file"])->GetString()) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      int TmpNbrPseudoPotentials;
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == false)
	{
	  cout << "Weights is not defined or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;
	}
      if (TmpNbrPseudoPotentials != (LzMax +1))
	{
	  cout << "Invalid number of pseudo-potentials" << endl;
	  return -1;	  
	}
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
    }

  char* OutputNameLz = new char [256 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf (OutputNameLz, "fermions_haldane_%s_n_%d_2s_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax);

  int* ReferenceState = 0;
  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
    {
      ReferenceState = new int[LzMax + 1];
      for (int i = 0; i <= LzMax; ++i)
	ReferenceState[i] = 0;
      if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "laughlin") == 0)
	for (int i = 0; i <= LzMax; i += 3)
	  ReferenceState[i] = 1;
      else
	if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "pfaffian") == 0)
	  for (int i = 0; i <= LzMax; i += 4)
	    {
	      ReferenceState[i] = 1;
	      ReferenceState[i + 1] = 1;
	    }
	else
	  if (strcasecmp(((SingleStringOption*) Manager["reference-state"])->GetString(), "readrezayi3") == 0)
	    for (int i = 0; i <= LzMax; i += 5)
	      {
		ReferenceState[i] = 1;
		ReferenceState[i + 1] = 1;
		ReferenceState[i + 2] = 1;
	      }
	  else
	    {
	      cout << "unknown reference state " << ((SingleStringOption*) Manager["reference-state"])->GetString() << endl;
	      return -1;
	    }
    }
  else
    {
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
	{
	  cout << "LzMax is not defined or as a wrong value" << endl;
	  return -1;
	}
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	{
	  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
	  return -1;     
	}
      if (MaxNbrLz != (LzMax + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return -1;     
	}
    }

  FermionOnSphereHaldaneBasis* Space = 0;
  if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
    Space = new FermionOnSphereHaldaneBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString(), MemorySpace);
  else
    Space = new FermionOnSphereHaldaneBasis(NbrParticles, InitialLz, LzMax, ReferenceState, MemorySpace);
  if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
    {
      Space->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
      return 0;
    }
  cout << "InitialLz = " << InitialLz << " " << Space->GetHilbertSpaceDimension() <<endl;
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
  if (OneBodyPotentials == 0)
    Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
							 ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
							 Architecture.GetArchitecture(), 
							 Memory, DiskCacheFlag,
							 LoadPrecalculationFileName);
  else
    Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials, OneBodyPotentials,
							 ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(),
							 Architecture.GetArchitecture(), 
							 Memory, DiskCacheFlag,
							 LoadPrecalculationFileName);
  double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
  if (((SingleDoubleOption*) Manager["energy-shift"])->GetDouble() != 0.0)
    Shift = ((SingleDoubleOption*) Manager["energy-shift"])->GetDouble();
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
    {
      EigenvectorName = new char [64];
      sprintf (EigenvectorName, "fermions_haldane_%s_n_%d_2s_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, InitialLz);
    }
  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, InitialLz, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  if (EigenvectorName != 0)
    {
      delete[] EigenvectorName;
    }
  delete Hamiltonian;
  if (FirstRun == true)
    FirstRun = false;

  return 0;
}
