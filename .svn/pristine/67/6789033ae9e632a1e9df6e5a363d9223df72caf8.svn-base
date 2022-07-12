#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"

#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"

#include "Operator/ParticleOnSphereWithSpinSzParityOperator.h"

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

#include "LanczosAlgorithm/LanczosManager.h"

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
  OptionManager Manager ("FQHESphereWithSpinL2Diagonalize" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  ParticleOnSphereManager ParticleManager(true, true, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "l2");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of the L^2 operator ", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "if non zero, override energy shift using the indicated value ", -10.0);
  //(*SystemGroup) += new BooleanOption  ('\n', "l2-up", "calculate L2 for upspins, only");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
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
      cout << "see man page for option syntax or type FQHESphereWithSpinL2Diagonalize -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz  = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  //  bool EmulateL2Up = Manager.GetBoolean("l2-up");
  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  bool FixedSz = !(Manager.GetBoolean("all-sz"));
  if (FixedSz)
    {
      if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
	sprintf (OutputNameLz, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalSz);
      else
	sprintf (OutputNameLz, "bosons_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalSz);
    }
  else
    {
      if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
	sprintf (OutputNameLz, "fermions_sphere_su2_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);
      else
	sprintf (OutputNameLz, "bosons_sphere_su2_%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax);
    }
  
  ParticleOnSphereWithSpin* Space =  (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(TotalLz);
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
  if (Manager.GetDouble("l2-factor")!=0.0)
    {
      Hamiltonian = new ParticleOnSphereWithSpinL2Hamiltonian(Space, NbrParticles, LzMax, TotalLz,
							      Architecture.GetArchitecture(), 
							      Manager.GetDouble("l2-factor"),
							      Memory, DiskCacheFlag,
							      LoadPrecalculationFileName); // , EmulateL2Up);
      if (Manager.GetDouble("s2-factor") != 0.0)
	((AbstractQHEOnSphereWithSpinHamiltonian*) Hamiltonian)->AddS2(TotalLz, TotalSz, Manager.GetDouble("s2-factor"), Memory, FixedSz);
    }
  else
    {
      if (Manager.GetDouble("s2-factor")!=0.0)
	Hamiltonian = new ParticleOnSphereWithSpinS2Hamiltonian(Space, NbrParticles, LzMax, TotalLz, TotalSz,
								Architecture.GetArchitecture(), 
								Manager.GetDouble("s2-factor"),
								Memory, DiskCacheFlag,
								LoadPrecalculationFileName, FixedSz);
      else
	{
	  cout << "Either L2 or S2 term have to be non-zero!"<<endl;
	  exit(-1);
	}
    }
  double Shift = Manager.GetDouble("energy-shift");
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;
  if (Manager.GetBoolean("eigenstate") == true)	
    {
      EigenvectorName = new char [256 + strlen(Manager.GetString("interaction-name"))];
      if (FixedSz)
	{
	  if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
	    sprintf (EigenvectorName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalSz, TotalLz);
	  else
	    sprintf (EigenvectorName, "bosons_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalSz, TotalLz);
	}
      else
	{
	  if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
	    sprintf (EigenvectorName, "fermions_sphere_su2_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalLz);
	  else
	    sprintf (EigenvectorName, "bosons_sphere_su2_%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, TotalLz);
	}
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
