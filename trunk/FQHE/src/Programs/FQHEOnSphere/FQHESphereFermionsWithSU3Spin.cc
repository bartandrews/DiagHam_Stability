#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinZ3Symmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzZ3Symmetry.h"

#include "Hamiltonian/ParticleOnSphereWithSU3SpinGenericHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <cstring>
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

  OptionManager Manager ("FQHESphereFermionsWithSU3Spin" , "0.01");
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
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption  ('t', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "tzsymmetrized-basis", "use Tz <-> -Tz symmetrized version of the basis (only valid if total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "z3symmetrized-basis", "use Z3 symmetrized version of the basis (only valid if total-y=0 and total-tz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSU3Spin -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool TzSymmetrizedBasis = Manager.GetBoolean("tzsymmetrized-basis");
  bool Z3SymmetrizedBasis = Manager.GetBoolean("z3symmetrized-basis");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  bool onDiskCacheFlag = Manager.GetBoolean("allow-disk-storage");
  bool FirstRun = true;
  double** PseudoPotentials  = new double*[6];
  for (int i = 0; i < 6; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
	PseudoPotentials[i][j] = 0.0;
    };
  double* OneBodyPotential11 = 0;
  double* OneBodyPotential22 = 0;
  double* OneBodyPotential33 = 0;

  int NbrN1 = (2 * NbrFermions) + TotalY + (3 * TotalTz);
  int NbrN2 = (2 * NbrFermions) + TotalY - (3 * TotalTz);
  int NbrN3 = NbrFermions - TotalY;
  if ((NbrN1 < 0 ) || (NbrN2 < 0 ) || (NbrN3 < 0) || ((NbrN1 % 6) != 0) || ((NbrN2 % 6) != 0) || ((NbrN3 % 3) != 0))
    {
      cout << "These values of Tz and Y cannot be achieved with this particle number!" << endl;
      return -1;
    }
  NbrN1 /= 6;
  NbrN2 /= 6;
  NbrN3 /= 3;

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
	  for (int i = 0; i < 6; ++i)
	    for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	      PseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials"] != 0)
	  {
	    cout << "Pseudopotentials has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials11", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials11" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[0][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials11"] != 0)
	  {
	    cout << "Pseudopotentials11 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials22", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials22" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[3][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials22"] != 0)
	  {
	    cout << "Pseudopotentials22 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials33", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials33" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[5][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials33"] != 0)
	  {
	    cout << "Pseudopotentials33 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials12", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials12" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[1][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials12"] != 0)
	  {
	    cout << "Pseudopotentials12 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials13", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials13" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[2][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials13"] != 0)
	  {
	    cout << "Pseudopotentials13 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials23", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in Pseudopotentials23" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    PseudoPotentials[4][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials23"] != 0)
	  {
	    cout << "Pseudopotentials23 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential11", ' ', OneBodyPotential11, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotential11 has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential22", ' ', OneBodyPotential22, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotential22 has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotential33", ' ', OneBodyPotential33, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotential33 has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
    }

#ifdef __64_BITS__
  cout << "Using 64 bit processing" << endl;
#endif


  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputNameLz, "fermions_sphere_su3_%s_n_%d_2s_%d_tz_%d_y_%d_lz.dat", 
	   Manager.GetString("interaction-name"), 
	   NbrFermions, LzMax, TotalTz, TotalY);

  int Max = (((LzMax - NbrN1 + 1) * NbrN1) + (((LzMax - NbrN2 + 1) * NbrN2)) + (((LzMax - NbrN3 + 1) * NbrN3)));
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
      ParticleOnSphereWithSU3Spin* Space = 0;
      if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == false))
	{
#ifdef __64_BITS__
	  if (LzMax <= 20)
#else
	    if (LzMax <= 9)
#endif
	      {
		Space = new FermionOnSphereWithSU3Spin(NbrFermions, L, LzMax, TotalTz, TotalY, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	}
      else
	{
#ifdef __64_BITS__
	  if (LzMax > 20)
#else
	    if (LzMax > 9)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == false))
	    {
	      if (Manager.GetString("load-hilbert") == 0)
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(NbrFermions, L, LzMax, TotalY, Manager.GetBoolean("minus-tzparity"), MemorySpace);
		  if (Manager.GetString("save-hilbert") != 0)
		    {
		      ((FermionOnSphereWithSU3SpinTzSymmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		}
	      else
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
		}
	    }
	  else
	    if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == true))
	      {
		if (Manager.GetString("load-hilbert") == 0)
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(NbrFermions, L, LzMax, TotalTz, MemorySpace);
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereWithSU3SpinZ3Symmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
		else
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(Manager.GetString("load-hilbert"), MemorySpace);
		  }
	      }
	    else
	      if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == true))
		{
		  if (Manager.GetString("load-hilbert") == 0)
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(NbrFermions, L, LzMax, Manager.GetBoolean("minus-tzparity"), MemorySpace);
		      if (Manager.GetString("save-hilbert") != 0)
			{
			  ((FermionOnSphereWithSU3SpinTzZ3Symmetry*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			  return 0;
			}
		    }
		  else
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(Manager.GetString("load-hilbert"), MemorySpace);
		    }
		}
	}

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEHamiltonian* Hamiltonian;      
      Hamiltonian = new ParticleOnSphereWithSU3SpinGenericHamiltonian(Space, NbrFermions, LzMax, PseudoPotentials, 
								      OneBodyPotential11, OneBodyPotential22, OneBodyPotential33,
								      Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);

      
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [120];
	  sprintf (EigenvectorName, "fermions_sphere_su3_%s_n_%d_2s_%d_tz_%d_y_%d_lz_%d",
		   Manager.GetString("interaction-name"), 
		   NbrFermions, LzMax, TotalTz, TotalY, L);
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
  return 0;
}


