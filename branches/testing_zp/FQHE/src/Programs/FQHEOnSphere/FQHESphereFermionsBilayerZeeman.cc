#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"
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


  OptionManager Manager ("FQHESphereFermionsBilayer" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup ("system options");  
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
    
  Architecture.AddOptionGroup(&Manager);
  Manager += SystemGroup;
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 9);
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new SingleDoubleOption  ('t', "tunnelling", "tunnelling splitting Delta_SAS", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "zeeman", "Zeeman field", 0.0);

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

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);	
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);

  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "start the Haldane algorithm from Laughlin state with exponent m)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");

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
  
  int NbrFermions = Manager.GetInteger("nbr-particles");
  double DeltaSAS = Manager.GetDouble("tunnelling"); 
  double Zeeman = Manager.GetDouble("zeeman"); 
  int LzMax = Manager.GetInteger("lzmax");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool MinusParity=Manager.GetBoolean("minus-lzparity");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long)Manager.GetInteger("fast-search")) << 20;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  bool onDiskCacheFlag = Manager.GetBoolean("disk");
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
  double* OneBodyPotentialUpDown = 0;

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
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
	    cout << "Pseudopotentials has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
	    cout << "PseudopotentialsUpUp has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
	    cout << "PseudopotentialsDownDown has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
	    cout << "PseudopotentialsUpDown has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', OneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', OneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
    }

  if (DeltaSAS!=0.0)
    {
      OneBodyPotentialUpDown = new double[LzMax + 1];
      for (int i=0; i<=LzMax; ++i)
	OneBodyPotentialUpDown[i] = -1.0*DeltaSAS;
    }

  if (Zeeman!=0.0)
    {
      OneBodyPotentialUpUp = new double[LzMax + 1];
      OneBodyPotentialDownDown = new double[LzMax + 1];
      for (int i=0; i<=LzMax; ++i)
	   {
		  OneBodyPotentialUpUp[i] = -1.0*Zeeman;
		  OneBodyPotentialDownDown[i] = Zeeman;
	    }
    }

  char* OutputNameLz = new char [512 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf (OutputNameLz, "fermions_sphere_su2_%s_n_%d_2s_%d_t_%f_h_%f_lz.dat", Manager.GetString("interaction-name"), 
	   NbrFermions, LzMax, DeltaSAS, Zeeman);

  
  int Max = (((LzMax - (NbrFermions+1)/2 + 1) * ((NbrFermions+1)/2)) + ((LzMax - (NbrFermions)/2 + 1) * (NbrFermions/2)));
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
      cout << "lz="<<L<<endl;
      ParticleOnSphereWithSpin* Space = 0; 

	//^^^^^^^^^^^^^^^^^^^^^^^HILBERT SPACE^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef __64_BITS__
	  if (LzMax <= 28)
#else
	  if (LzMax <= 13)
#endif
	  {  
		if ((LzSymmetrizedBasis == false) || (L != 0))
	      Space = new FermionOnSphereWithSpinAllSz(NbrFermions, L, LzMax, MemorySpace);
	    else
	      {
			if (Manager.GetString("load-hilbert") != 0)
				Space = new FermionOnSphereWithSpinAllSzLzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
			else
		  		Space = new FermionOnSphereWithSpinAllSzLzSymmetry(NbrFermions, LzMax, MinusParity, MemorySpace);
			if (Manager.GetString("save-hilbert") != 0)
		    		((FermionOnSphereWithSpinAllSzLzSymmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
	       }
		}
	  else
#ifdef __128_BIT_LONGLONG__
	    if (LzMax <= 63)
#else
	      if (LzMax <= 31)
#endif
		cout<<"Not implemented yet"<<endl;
	/*
		{
		  if ((SymmetrizedBasis == false) || (totalLz != 0))
		    Space = new FermionOnSphereLong(NbrParticles, totalLz, LzMax, MemorySpace);
		  else
		    {
		      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
			Space = new FermionOnSphereSymmetricBasisLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		      else
			Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		      if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
			{
			  ((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			  return 0;
			}
		    }
		}
	     else
		Space = new FermionOnSphereUnlimited(NbrParticles, totalLz, LzMax, MemorySpace);
	}
		*/
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


	  
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();
      
      AbstractQHEHamiltonian* Hamiltonian;      
      Hamiltonian = new ParticleOnSphereWithSpinGenericHamiltonian(Space, NbrFermions, LzMax, PseudoPotentials,
								   OneBodyPotentialUpUp, OneBodyPotentialDownDown,
								   OneBodyPotentialUpDown, 
								   Architecture.GetArchitecture(), Memory, onDiskCacheFlag,
								   LoadPrecalculationFileName);

      
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "fermions_sphere_su2_%s_n_%d_2s_%d_t_%f_h_%f_lz_%d",
		   ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
		   NbrFermions, LzMax, DeltaSAS, Zeeman, L);
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
    }
  delete[] OutputNameLz;
  for (int i = 0; i < 3; ++i) delete[] PseudoPotentials[i];
  delete[] PseudoPotentials;
  return 0;
}


