#include "HilbertSpace/ParticleOnSphereManager.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereDroplet.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

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
  OptionManager Manager ("FQHESphereL2Diagonalize" , "0.01");


  ParticleOnSphereManager ParticleManager(true, true, 1);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-fluxes1", "number of fluxes in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-particles1", "max number of particles in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-holes1", "max number of holes in a droplet", 0);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-fluxes2", "secondary condition for number of fluxes in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-particles2", "secondary condition for  max number of particles in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-holes2", "secondary condition for max number of holes in a droplet", 0);

  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "l2");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "if non zero, override energy shift using the indicated value ", -10.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of Hamiltonian", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereL2Diagonalize -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxes1 = Manager.GetInteger("nbr-fluxes1");
  int MaxNbrParticles1 = Manager.GetInteger("max-particles1");
  int MaxNbrHoles1 = Manager.GetInteger("max-holes1");

  int NbrFluxes2 = Manager.GetInteger("nbr-fluxes2");
  int MaxNbrParticles2 = Manager.GetInteger("max-particles2");
  int MaxNbrHoles2 = Manager.GetInteger("max-holes2");

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz  = Manager.GetInteger("total-lz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  bool Hermitian = !Manager.GetBoolean("no-hermitian");
  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
    sprintf (OutputNameLz, "fermions_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);
  else
    sprintf (OutputNameLz, "bosons_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);

  ParticleOnSphere* Space = new FermionOnSphereDroplet(NbrParticles, TotalLz, LzMax, NbrFluxes1, MaxNbrParticles1, MaxNbrHoles1, NbrFluxes2, MaxNbrParticles2, MaxNbrHoles2, Memory); 

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  AbstractQHEOnSphereHamiltonian* Hamiltonian = new ParticleOnSphereL2Hamiltonian(Space, NbrParticles, LzMax, TotalLz,
										  Architecture.GetArchitecture(), 
										  Manager.GetDouble("l2-factor"),
										  Memory, true, DiskCacheFlag,
										  LoadPrecalculationFileName, Hermitian);

  double Shift = Manager.GetDouble("energy-shift");
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (Manager.GetBoolean("eigenstate") == true)	
    {
      EigenvectorName = new char [64];
      if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
	sprintf (EigenvectorName, "fermions_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalLz);
      else
	sprintf (EigenvectorName, "bosons_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalLz);
    }
  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, TotalLz, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  if (EigenvectorName != 0)
    delete[] EigenvectorName;
  delete Hamiltonian;
  if (FirstRun == true)
    FirstRun = false;

  return 0;
}
