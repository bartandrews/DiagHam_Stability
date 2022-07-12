#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include "Hamiltonian/ParticleOnSphereWithSpinGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "Options/Options.h"

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


  OptionManager Manager ("FQHESphereBosonsWithSpin" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(false, true, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);

  (*SystemGroup) += new BooleanOption ('\n', "l2-s2-only", "compose Hamiltonian only of L2 and S2 terms");
  
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "fast-disk", "use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm");
  (*LanczosGroup) += new  BooleanOption ('\n', "resume-fastdisk", "resume the fast-disk mode Lanczos algorithm from a stopped one (for example due to computer crash)");
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "start the Haldane algorithm from Laughlin state with exponent m)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "s2-memory", "amount of memory that can be allocated for fast multiplication of s2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "l2-memory", "amount of memory that can be allocated for fast multiplication of l2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsWithSpin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");

  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;  
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  bool onDiskCacheFlag = ((BooleanOption*) Manager["allow-disk-storage"])->GetBoolean();
  bool FirstRun = true;
  double** PseudoPotentials  = new double*[10];
  for (int i = 0; i < 3; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
	PseudoPotentials[i][j] = 0.0;
    };
  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;

  int NbrUp = (NbrBosons + SzTotal) >> 1;
  int NbrDown = (NbrBosons - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

  if (Manager.GetString("interaction-file") == 0)
    {
      if (!Manager.GetBoolean("l2-s2-only"))
	{
	  cout << "an interaction file has to be provided" << endl;
	  return -1;
	}
    }
  else
    {
      if (FQHESphereSU2GetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials,
					   OneBodyPotentialUpUp, OneBodyPotentialDownDown) == false)
	return -1;
    }

  char* InteractionName;
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';    
  if (Manager.GetBoolean("l2-s2-only"))
    {      
      InteractionName=new char[50];
      if (Manager.GetDouble("l2-factor") != 0.0)
	{
	  if (Manager.GetDouble("s2-factor") != 0.0)
	    sprintf(InteractionName,"l2_%g_s2_%g",Manager.GetDouble("l2-factor"), Manager.GetDouble("s2-factor"));
	  else
	    sprintf(InteractionName,"l2_%g",Manager.GetDouble("l2-factor"));
	}
      else
	if (Manager.GetDouble("s2-factor") != 0.0)
	  sprintf(InteractionName,"s2_%g",Manager.GetDouble("s2-factor"));
    }
  else
    {
      InteractionName=Manager.GetString("interaction-name");
      if (Manager.GetDouble("l2-factor") != 0.0)
	{
	  if (Manager.GetDouble("s2-factor") != 0.0)
	    sprintf(ExtraTerms,"_l2_%g_s2_%g",Manager.GetDouble("l2-factor"), Manager.GetDouble("s2-factor"));
	  else
	    sprintf(ExtraTerms,"_l2_%g",Manager.GetDouble("l2-factor"));
	}
      else
	if (Manager.GetDouble("s2-factor") != 0.0)
	  sprintf(ExtraTerms,"_s2_%g",Manager.GetDouble("s2-factor"));
    }
  char* OutputNameLz = new char [512 + strlen(InteractionName)];
  sprintf (OutputNameLz, "bosons_sphere_su2_%s%s_n_%d_2s_%d_sz_%d_lz.dat", InteractionName, ExtraTerms,
	   NbrBosons, LzMax, SzTotal);

  int Max = (LzMax * (NbrUp+NbrDown));
  cout << "maximum Lz value = " << Max << endl;

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
  if (NbrLz > 0)
    {
      if (L + (2 * (NbrLz - 1)) < Max)
	Max = L + (2 * (NbrLz - 1));
    }
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;
      ParticleOnSphereWithSpin* Space = 0;
      Space = (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(L);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	Space->PrintState(cout, i) << endl;

      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian;

      if (Manager.GetBoolean("l2-s2-only"))
	{
	  if (((SingleDoubleOption*) Manager["l2-factor"])->GetDouble() != 0.0)
	    {
	      Hamiltonian = new ParticleOnSphereWithSpinL2Hamiltonian(Space, NbrBosons, LzMax, L, 
								      Architecture.GetArchitecture(),
								      Manager.GetDouble("l2-factor"),
						      ((unsigned long)Manager.GetInteger("l2-memory")) << 20);
	      
	      
	      if (((SingleDoubleOption*) Manager["s2-factor"])->GetDouble() != 0.0)
		Hamiltonian->AddS2(L, SzTotal, Manager.GetDouble("s2-factor"),
				   ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
	    }
	  else
	    Hamiltonian = new ParticleOnSphereWithSpinS2Hamiltonian(Space, NbrBosons, LzMax, L, SzTotal,
								    Architecture.GetArchitecture(),
								      Manager.GetDouble("s2-factor"),
						      ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
	}
      else // full Hamiltonian
	{
	  Hamiltonian = new ParticleOnSphereWithSpinGenericHamiltonian(Space, NbrBosons, LzMax, PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown, NULL, 
								       Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);
	  
	  if (((SingleDoubleOption*) Manager["s2-factor"])->GetDouble() != 0.0)
	    Hamiltonian->AddS2(L, SzTotal, ((SingleDoubleOption*) Manager["s2-factor"])->GetDouble(), ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
	  if (((SingleDoubleOption*) Manager["l2-factor"])->GetDouble() != 0.0)
	    Hamiltonian->AddL2(L, SzTotal, ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), ((unsigned long)Manager.GetInteger("l2-memory")) << 20);
	}
      
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "bosons_sphere_su2_%s%s_n_%d_2s_%d_sz_%d_lz_%d",
		   InteractionName, ExtraTerms,
		   NbrBosons, LzMax, SzTotal, L);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete Space;      
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	  EigenvectorName = 0;
	}
      if (FirstRun == true)
	FirstRun = false;
      if (HaldaneBasisFlag) return 0; // only one subspace defined...
    }
  delete[] OutputNameLz;
  delete[] ExtraTerms;
  return 0;
}


