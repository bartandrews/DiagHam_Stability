#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "HilbertSpace/FermionOnDiskHaldaneBasis.h"
#include "HilbertSpace/FermionOnDiskLong.h"

#include "Hamiltonian/ParticleOnDiskGenericThreeBodyHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnDiskMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
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
  OptionManager Manager ("FQHEDiskFermionsThreeBodyGeneric" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "maximum-momentum", "maximum total angular momentum to study", 10, true, 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "minimum-momentum", "minimum total angular momentum to study", 1, true, 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "force-maxmomentum", "force the maximum single particle momentum to a particular value (negative from the number of particles and the state total angular momentum)", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskFermionsThreeBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int MMin = Manager.GetInteger("minimum-momentum");
  int MMax = Manager.GetInteger("maximum-momentum");
  if (MMin < (((NbrParticles - 1) * (NbrParticles)) / 2))
    MMin = (((NbrParticles - 1) * (NbrParticles)) / 2);
  if (MMax < MMin)
    MMax = MMin;
  int ForceMaxMomentum = Manager.GetInteger("force-maxmomentum");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  double* PseudoPotentials = 0;
  double* OneBodyPotentials = 0;
  double* ThreeBodyPotentials = 0;
  int TmpNbrThreeBodyPseudoPotentials = 0;

  int* ReferenceState = 0;
  if (HaldaneBasisFlag == true)
    {
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", ForceMaxMomentum) == false) || (ForceMaxMomentum <= 0))
	{
	  cout << "LzMax is not defined or as a wrong value" << endl;
	  return -1;
	}
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	{
	  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	  return -1;     
	}
      if (MaxNbrLz != (ForceMaxMomentum + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return -1;     
	}
      MMax = 0;
      for (int i = 1; i <= ForceMaxMomentum; ++i)
	MMax += i * ReferenceState[i];
      MMin = MMax;
    }

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
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentials", ' ', ThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == false)
	{
	  cout << "ThreebodyPseudopotentials are not defined or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	  return -1;
	}
      int TmpNbrPseudoPotentials;
      int TmpMax = 2 * ForceMaxMomentum;
      if (ForceMaxMomentum < 0)
	TmpMax = 2 * MMax;
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials < (TmpMax +1))
	    {	  
	      cout << "warning : not enougth pseudo-potentials, higher relative monentum pseudo potentials will be set to zero" << endl;
	      double* TmpPseudoPotentials = new double [TmpMax + 1];
	      int i = 0;
	      for (; i < TmpNbrPseudoPotentials; ++i)
		TmpPseudoPotentials[i] = PseudoPotentials[i];
	      for (; i <= TmpMax; ++i)
		TmpPseudoPotentials[i] = 0.0;	    
	      delete[] PseudoPotentials;
	      PseudoPotentials = TmpPseudoPotentials;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if ((TmpNbrPseudoPotentials < (ForceMaxMomentum +1)) || ((ForceMaxMomentum < 0) && ((TmpNbrPseudoPotentials < (MMax +1)))))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
    }

  char* OutputNameLz = new char [1024 + strlen(Manager.GetString("interaction-name"))];
  if (ForceMaxMomentum >= 0)
    sprintf (OutputNameLz, "fermions_disk_%s_n_%d_lzmax_%d_lz_%d.dat", Manager.GetString("interaction-name"), NbrParticles, ForceMaxMomentum, MMax);
  else
    sprintf (OutputNameLz, "fermions_disk_%s_n_%d_lz_%d.dat", Manager.GetString("interaction-name"), NbrParticles, MMax);
  for (int  L = MMin; L <= MMax; ++L)
    {
      ParticleOnSphere* Space = 0;
      int TmpMaxMomentum = (L - (((NbrParticles - 1) * (NbrParticles - 2)) / 2));
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (TmpMaxMomentum <= 62)
#else
	  if (TmpMaxMomentum <= 30)
#endif
	    Space = new FermionOnDisk(NbrParticles, L, TmpMaxMomentum, MemorySpace);
	  else
#ifdef __128_BIT_LONGLONG__
	    if (TmpMaxMomentum <= 126)
#else
	      if (TmpMaxMomentum <= 62)
#endif
		Space = new FermionOnDiskLong(NbrParticles, L, TmpMaxMomentum, MemorySpace);
	      else
		Space = new FermionOnDiskUnlimited(NbrParticles, L, TmpMaxMomentum, MemorySpace);
	}
      else
	{
#ifdef __64_BITS__
	  if (TmpMaxMomentum <= 62)
#else
	    if (TmpMaxMomentum <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  Space = new FermionOnDiskHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  Space = new FermionOnDiskHaldaneBasis(NbrParticles, L, TmpMaxMomentum, ReferenceState, MemorySpace);
		if (Manager.GetString("save-hilbert") != 0)
		  {
		    ((FermionOnDiskHaldaneBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		    return 0;
		  }
	      }
	}
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Hamiltonian = new ParticleOnDiskGenericThreeBodyHamiltonian(Space, NbrParticles, TmpMaxMomentum, ThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials - 1,
								  PseudoPotentials, OneBodyPotentials,
								  Architecture.GetArchitecture(), 
								  Memory, DiskCacheFlag,
								  LoadPrecalculationFileName);
      double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) MMax));
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [256];
	  if (ForceMaxMomentum >= 0)
	    sprintf (EigenvectorName, "fermions_disk_%s_n_%d_lzmax_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, ForceMaxMomentum, L);
	  else
	    sprintf (EigenvectorName, "fermions_disk_%s_n_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, L);
	}
      
      QHEOnDiskMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      delete Hamiltonian;
      delete Space;
      if (FirstRun == true)
	FirstRun = false;
    }

  return 0;
}
