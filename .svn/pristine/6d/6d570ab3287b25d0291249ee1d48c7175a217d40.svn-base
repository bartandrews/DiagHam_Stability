#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSqueezedBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"

#include "Hamiltonian/ParticleOnSphereWithSpinGenericHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

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


  OptionManager Manager ("FQHESphereFermionsWithSpin" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, false, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

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
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");

  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
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

  int NbrUp = (NbrFermions + SzTotal) >> 1;
  int NbrDown = (NbrFermions - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

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
      double* TmpPseudoPotentials;
      bool Flag = false;
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials" << endl;
	      return -1;	  
	    }
	  for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	      PseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials"] != 0)
	  {
	    cout << "Pseudopotentials has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsUpUp" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[0][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["PseudopotentialsUpUp"] != 0)
	  {
	    cout << "PseudopotentialsUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsDownDown" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[1][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["PseudopotentialsDownDown"] != 0)
	  {
	    cout << "PseudopotentialsDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDown" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[2][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["PseudopotentialsUpDown"] != 0)
	  {
	    cout << "PseudopotentialsUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', OneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', OneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	      return -1;
	    }
	}
    }

  char* OutputNameLz = new char [512 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf (OutputNameLz, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
	   NbrFermions, LzMax, SzTotal);

  int Max = (((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown));
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
      if (HaldaneBasisFlag == false)
	{
	  Space = (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(L);
	}
      else // using Squeezed Haldane Basis!
	{
	  int* ReferenceState = 0;
	  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
	    {
	      ReferenceState = new int[LzMax + 1];
	      for (int i = 0; i <= LzMax; ++i)
		ReferenceState[i] = 0;
	      int tmpNbrParticles=0;
	      if (Manager.GetInteger("laughlin-exponent")>0)		
		for (int i = 0; i <= LzMax; i += Manager.GetInteger("laughlin-exponent"))
		  {
		    ReferenceState[i] = 3;
		    tmpNbrParticles+=2;
		  }
	      if (tmpNbrParticles != NbrFermions)
		{
		  cout << "Wrong shift for this laughlin seed state"<< endl;
		  exit(1);
		}
	    }
	  else
	    {
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrFermions) == false) || (NbrFermions <= 0))
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
	  // assign Hilbert space here - no further symmetries implemented, yet...
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		  Space = new FermionOnSphereWithSpinSqueezedBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  Space = new FermionOnSphereWithSpinSqueezedBasis(NbrFermions, L, LzMax, SzTotal, ReferenceState, MemorySpace);
		if (Manager.GetString("save-hilbert") != 0)
		  {
		    ((FermionOnSphereWithSpinSqueezedBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			 return 0;
		  }
	      }
	  
	}
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEHamiltonian* Hamiltonian;      
      Hamiltonian = new ParticleOnSphereWithSpinGenericHamiltonian(Space, NbrFermions, LzMax, PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown, NULL, 
								   Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);

      
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d",
		   ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
		   NbrFermions, LzMax, SzTotal, L);
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
  return 0;
}


