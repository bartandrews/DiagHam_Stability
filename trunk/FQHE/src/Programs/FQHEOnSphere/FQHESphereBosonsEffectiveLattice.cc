// ***************************************************************************************************** //
// *			     Bosons on lattice near n_phi=1/2                                          * //
// ***************************************************************************************************** //

#include "config.h"

#include "Options/Options.h"

#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "Hamiltonian/ParticleOnSphereEffectiveLatticeHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstdio>
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


  OptionManager Manager ("FQHESphereBosonsEffectiveLattice" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup ("system options");  
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
    
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 9);
  //(*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  //(*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new SingleDoubleOption  ('a', "alpha", "deviation of flux-density from n_phi=1/2", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('t', "tunnelling", "tunnelling splitting Delta_SAS that couples to S_z term", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('b', "density-imbalance", "density imbalance term that couples to S_x", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);
  (*SystemGroup) += new BooleanOption ('\n', "project-l2", "add a projector onto the L2 groundstate");
  
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
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "l2-memory", "precalculation memory for L^2 operator",1000);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "l2-nbr-vectors", "number of states stored for L^2 projection",10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-storage", "additional number of vectors in RAM when using projected Lanczos", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-iter-max", "maximum number of iterations for internal lanczos",100);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "projector-precision", "define Lanczos precision for projection (0 if automatically defined by the program)", 1e-14);
  (*LanczosGroup) += new BooleanOption ('\n', "restart-projection", "allow lanczos projections to be restarted if full convergence not yet reached");
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-blockvectors", "use file that describe a set of initial vectors for the block Lanczos algorithm (syntax : InitialVectors=vec0.vec vec1.vec ...)", 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "fast-disk", "use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm");
  (*LanczosGroup) += new  BooleanOption ('\n', "resume-fastdisk", "resume the fast-disk mode Lanczos algorithm from a stopped one (for example due to computer crash)");

  

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);	
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", 1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  double Alpha = Manager.GetDouble("alpha");
  double DeltaSAS = Manager.GetDouble("tunnelling");
  int PairParity = Manager.GetInteger("pair-parity");
  if ((PairParity>=0)&&(DeltaSAS!=0.0))
    {
      cout << "Attention, pair-parity requires single particle tunnelling to vanish. Ignoring requested parity."<<endl;
      PairParity = -1;
    }
  double DensityImbalance = Manager.GetDouble("density-imbalance"); 
  int LzMax = Manager.GetInteger("lzmax");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long)Manager.GetInteger("fast-search")) << 20;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  bool onDiskCacheFlag = Manager.GetBoolean("disk");
  bool FirstRun = true;
  double** PseudoPotentials = 0;
  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double* OneBodyPotentialUpDown = 0;

  if (Manager.GetString("interaction-file") != 0)
    {
      PseudoPotentials  = new double*[4];
      for (int i = 0; i < 4; ++i)
	{
	  PseudoPotentials[i] = new double[LzMax + 1];
	  for (int j = 0; j <= LzMax; ++j)
	    PseudoPotentials[i][j] = 0.0;
	};
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
     if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsMixed", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsMixed" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[3][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["PseudopotentialsMixed"] != 0)
	  {
	    cout << "PseudopotentialsMixed has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
      OneBodyPotentialUpUp = new double[LzMax + 1];
      OneBodyPotentialDownDown = new double[LzMax + 1];
      for (int i=0; i<=LzMax; ++i)  
	{
	  OneBodyPotentialUpUp[i] = -1.0*DeltaSAS;
	  OneBodyPotentialDownDown[i] = DeltaSAS;
	}
    }

  if (DensityImbalance!=0.0)
    {
      OneBodyPotentialUpDown = new double[LzMax + 1];
      for (int i=0; i<=LzMax; ++i)  
	OneBodyPotentialUpDown[i] = -DensityImbalance;
    }


  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';
  int Offset=0;
  if (PairParity>=0)
    Offset+=sprintf(ExtraTerms,"_parity_%d",PairParity);
  if (Manager.GetBoolean("project-l2"))
    sprintf(ExtraTerms+Offset,"_Pl2");

  if (Manager.GetString("interaction-file")!=NULL)
    sprintf (OutputNameLz, "bosons_sphere_eff_su2%s_%s_a_%g_n_%d_2s_%d_t_%g_b_%g_lz.dat", ExtraTerms, Manager.GetString("interaction-name"), Alpha, NbrBosons, LzMax, DeltaSAS, DensityImbalance);
  else
    sprintf (OutputNameLz, "bosons_sphere_eff_su2%s_a_%g_n_%d_2s_%d_t_%g_b_%g_lz.dat", ExtraTerms, Alpha, NbrBosons, LzMax, DeltaSAS, DensityImbalance);

  
  int Max = LzMax * NbrBosons;
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
      double Shift = 0.0;
      cout << "lz="<<L<<endl;
      ParticleOnSphereWithSpin* Space = 0; 

      //^^^^^^^^^^^^^^^^^^^^^^^HILBERT SPACE^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      if (PairParity>=0)
	Space = new BosonOnSphereWithSpinAllSz(NbrBosons, L, LzMax, PairParity, MemorySpace);
      else
	Space = new BosonOnSphereWithSpinAllSz(NbrBosons, L, LzMax, MemorySpace);
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();
      
      AbstractQHEHamiltonian* Hamiltonian;
      if (PseudoPotentials!=0)
	for (int i=0; i<=LzMax; ++i)
	  cout << "PP["<<i<<"]="<<PseudoPotentials[0][i]<<" "<<PseudoPotentials[1][i]<<" "<<PseudoPotentials[2][i]<<" "<<PseudoPotentials[3][i]<<endl;
      Hamiltonian = new ParticleOnSphereEffectiveLatticeHamiltonian(Space, NbrBosons, LzMax, Alpha, PseudoPotentials,
								    OneBodyPotentialUpUp, OneBodyPotentialDownDown,
								    OneBodyPotentialUpDown, 
								    Architecture.GetArchitecture(), Memory, onDiskCacheFlag,
								    LoadPrecalculationFileName);
      
      Hamiltonian->ShiftHamiltonian(Shift);
      cout << "Shift="<<Shift<<endl;
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}

      int NbrProjectors = 0;
      AbstractHamiltonian** Projectors = NULL;
      if (Manager.GetBoolean("project-l2")) ++NbrProjectors;
      Projectors = new AbstractHamiltonian*[NbrProjectors];
      NbrProjectors = 0;
      if (Manager.GetBoolean("project-l2"))
	{
	  AbstractHamiltonian* L2Projector =
	    new ParticleOnSphereWithSpinL2Hamiltonian(Space, NbrBosons, LzMax, L, 
						      Architecture.GetArchitecture(), 1.0, ((long)Manager.GetInteger("l2-memory"))<<20);
	  L2Projector->ShiftHamiltonian(-0.25*(double)L*(L+2.0));
	  Projectors[NbrProjectors++]=L2Projector;
	}
      
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [512];
	  if (Manager.GetString("interaction-file")!=NULL)
	    sprintf (EigenvectorName, "bosons_sphere_eff_su2%s_%s_a_%g_n_%d_2s_%d_t_%g_b_%g_lz_%d", ExtraTerms, Manager.GetString("interaction-name"), Alpha, NbrBosons, LzMax, DeltaSAS, DensityImbalance, L);
	  else
	    sprintf (EigenvectorName, "bosons_sphere_eff_su2%s_a_%g_n_%d_2s_%d_t_%g_b_%g_lz_%d", ExtraTerms, Alpha, NbrBosons, LzMax, DeltaSAS, DensityImbalance, L);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax, Projectors, NbrProjectors);
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
  if (PseudoPotentials!=NULL)
    {
      for (int i = 0; i < 4; ++i) delete[] PseudoPotentials[i];
      delete[] PseudoPotentials;
    }
  if (DeltaSAS!=0.0)
    {
      delete[] OneBodyPotentialUpUp;
      delete[] OneBodyPotentialDownDown;
    }
  if (DensityImbalance!=0.0)
    delete[] OneBodyPotentialUpDown;
  if (ExtraTerms!=NULL)
    delete [] ExtraTerms;
  return 0;
}
