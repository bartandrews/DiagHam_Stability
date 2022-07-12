#include "HilbertSpace/BosonOnDisk.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnDiskHaldaneBasisShort.h"
#include "Hamiltonian/ParticleOnDiskNBodyHardCoreHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnDiskMainTask.h"

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
  OptionManager Manager ("FQHEDiskBosonsNBodyHardCore" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("Tools options");
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
  (*SystemGroup) += new SingleIntegerOption  ('l', "maximum-momentum", "maximum single particle momentum to study", 10, true, 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "minimum-momentum", "minimum single particle momentum to study", 1, true, 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "force-maxmomentum", "force the maximum single particle momentum to a particular value (negative from the number of particles and the state total angular momentum)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-nbody", "number of particle that can interact simultaneously through the n-body hard-core interaction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskBosonsNBodyHardCore -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int NbrNBody = ((SingleIntegerOption*) Manager["nbr-nbody"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int MMin = ((SingleIntegerOption*) Manager["minimum-momentum"])->GetInteger();
  int MMax = ((SingleIntegerOption*) Manager["maximum-momentum"])->GetInteger();
  if (MMax < MMin)
    MMax = MMin;
  int ForceMaxMomentum = ((SingleIntegerOption*) Manager["force-maxmomentum"])->GetInteger();
  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  bool FirstRun = true;

  int* ReferenceState = 0;
  if (HaldaneBasisFlag == true)
    {
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrBosons) == false) || (NbrBosons <= 0))
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
	  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
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

  char* OutputNameLz = new char [1024];
  if (ForceMaxMomentum >= 0)
    if (HaldaneBasisFlag == true)
      sprintf (OutputNameLz, "bosons_disk_haldane_hardcore_nbody_%d_n_%d_lzmax_%d_lz_%d.dat", NbrNBody, NbrBosons, ForceMaxMomentum, MMax);
    else
      sprintf (OutputNameLz, "bosons_disk_hardcore_nbody_%d_n_%d_lzmax_%d_lz_%d.dat", NbrNBody, NbrBosons, ForceMaxMomentum, MMax);
  else
    if (HaldaneBasisFlag == true)
      sprintf (OutputNameLz, "bosons_disk_haldane_hardcore_nbody_%d_n_%d_lz_%d.dat", NbrNBody, NbrBosons, MMax);
    else
      sprintf (OutputNameLz, "bosons_disk_hardcore_nbody_%d_n_%d_lz_%d.dat", NbrNBody, NbrBosons, MMax);
  for (int  L = MMin; L <= MMax; ++L)
    {
      ParticleOnSphere* Space = 0;
      int TmpMaxMomentum = L;
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;
      if (HaldaneBasisFlag == false)
	{
#ifdef  __64_BITS__
	  if ((ForceMaxMomentum + NbrBosons - 1) < 63)
#else
	    if ((ForceMaxMomentum + NbrBosons - 1) < 31)	
#endif
	      Space = new BosonOnDiskShort (NbrBosons, L, ForceMaxMomentum);	  
	    else	  
	      Space = new BosonOnDisk (NbrBosons, L, ForceMaxMomentum);
	}
      else
	{
#ifdef  __64_BITS__
	  if ((ForceMaxMomentum + NbrBosons - 1) < 63)
#else
	    if ((ForceMaxMomentum + NbrBosons - 1) < 31)	
#endif
	      Space = new BosonOnDiskHaldaneBasisShort(NbrBosons, L, TmpMaxMomentum, ReferenceState);
	}
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Hamiltonian = new ParticleOnDiskNBodyHardCoreHamiltonian(Space, NbrBosons, TmpMaxMomentum, NbrNBody, Architecture.GetArchitecture(), 
							       Memory, LoadPrecalculationFileName);
      double Shift = - 0.5 * ((double) (NbrBosons * NbrBosons)) / (0.5 * ((double) MMax));
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [64];
	  if (HaldaneBasisFlag == true)
	    sprintf (EigenvectorName, "bosons_disk_haldane_hardcore_nbody_%d_n_%d_lzmax_%d_lz_%d", NbrNBody, NbrBosons, ForceMaxMomentum, L);
	  else
	    sprintf (EigenvectorName, "bosons_disk_hardcore_nbody_%d_n_%d_lzmax_%d_lz_%d", NbrNBody, NbrBosons, ForceMaxMomentum, L);
	}
      QHEOnDiskMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName);
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
