#include "HilbertSpace/FermionOnSphereWithSpinAndPairing.h"

#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "MainTask/QHEOnSphereMainTask.h"

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
  OptionManager Manager ("FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing" , "0.01");
  OptionGroup* SystemGroup  = new OptionGroup("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleDoubleOption ('\n', "charging-energy", "factor in front of the charging energy (i.e 1/(2C))", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "average-nbrparticles", "average number of particles", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "s2-memory", "amount of memory that can be allocated for fast multiplication of s2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "l2-memory", "amount of memory that can be allocated for fast multiplication of l2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LzMax = Manager.GetInteger("lzmax");
  int TotalSz = Manager.GetInteger("total-sz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;

  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double* OneBodyPotentialUpDown = 0;
  double* OneBodyPotentialPairing = 0;
  double** PseudoPotentials  = new double*[3];
  for (int i = 0; i < 3; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
    }

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereSU2GetPseudopotentialsWithPairing(Manager.GetString("interaction-file"), LzMax, PseudoPotentials, 
						      OneBodyPotentialUpUp, OneBodyPotentialDownDown, OneBodyPotentialUpDown, OneBodyPotentialPairing) == false)
	{
	  return -1;
	}
    }
  if (OneBodyPotentialUpDown != 0)
    {
      cout << "warning, OneBodyPotentialUpDown is not supported" << endl;
    }
  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputNameLz, "fermions_sphere_su2_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d_lz.dat", Manager.GetString("interaction-name"), 
	   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz);

  int MinNbrParticles = abs(TotalSz);
  int MaxNbrParticles = (2 * (LzMax + 1)) - abs(TotalSz);
  int MaxL = 0;

  for (int TmpNbrParticles = MinNbrParticles; TmpNbrParticles <= MaxNbrParticles; TmpNbrParticles += 2)
    {
      int NbrUp = (TmpNbrParticles + TotalSz) / 2;
      int NbrDown = (TmpNbrParticles - TotalSz) / 2;
      if ((NbrUp <= (LzMax + 1)) && (NbrUp >= 0) && (NbrDown <= (LzMax + 1)) && (NbrDown >= 0))
	{
	  int TmpMaxL = (((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown));
	  if (TmpMaxL > MaxL)
	    MaxL = TmpMaxL;
	}
    }

  int  L = 0;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(MaxL) & 1) != 0)
	L |= 1;
      else
	L &= ~0x1;
    }
  if (NbrLz > 0)
    {
      if (L + (2 * (NbrLz - 1)) < MaxL)
	MaxL = L + (2 * (NbrLz - 1));
    }

  for (; L <= MaxL; L += 2)
    {
      ParticleOnSphereWithSpin* Space = new FermionOnSphereWithSpinAndPairing(L, LzMax, TotalSz);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricGenericHamiltonianAndPairing (Space, LzMax, 
												   PseudoPotentials, OneBodyPotentialUpUp, 
												   OneBodyPotentialDownDown,
												   OneBodyPotentialPairing,
												   Manager.GetDouble("charging-energy"), 
												   Manager.GetDouble("average-nbrparticles"),
												   Architecture.GetArchitecture(), 
												   Memory, DiskCacheFlag,
												   LoadPrecalculationFileName);
      if (Manager.GetDouble("l2-factor") != 0.0)
	Hamiltonian->AddL2(L, TotalSz, Manager.GetDouble("l2-factor"), ((unsigned long)Manager.GetInteger("l2-memory")) << 20); 
      if (Manager.GetDouble("s2-factor") != 0.0)
	Hamiltonian->AddS2(L, TotalSz, Manager.GetDouble("s2-factor"), ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
      double Shift = - 10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "fermions_sphere_su2_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), 
		   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz, L);
	}
      
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      delete Hamiltonian;
      if (FirstRun == true)
	FirstRun = false;
      delete Space;
    }

  return 0;
}
