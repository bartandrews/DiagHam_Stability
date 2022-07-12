#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"

#include "Hamiltonian/ParticleOnSphereWithSpinGenericThreeBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

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
  OptionManager Manager ("FQHESphereBosonsWithSpinThreeBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  ParticleOnSphereManager ParticleManager(false, true, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  

  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "pairing", "add a pair hopping term to the Hamiltonian in all-sz mode", 0.0);
  // (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
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
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsThreeBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  // bool GroundFlag = Manager.GetBoolean("ground");
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int SzTotal = Manager.GetInteger("total-sz");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool AllSz = Manager.GetBoolean("all-sz");
  bool PairParity = Manager.GetInteger("pair-parity");
  bool FirstRun = true;

  int NbrUp = (NbrParticles + SzTotal) >> 1;
  int NbrDown = (NbrParticles - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  if ((NbrParticles&1)!=(SzTotal&1))
    {
      cout << "NbrParticles and SzTotal need to be of the same parity."<<endl;
      return -1;
    }

  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double** PseudoPotentials  = 0;

  double* ThreeBodyPotentials12 = 0;
  int NbrThreeBodyPseudoPotentials12 = 0;
  double* ThreeBodyPotentials32 = 0;
  int NbrThreeBodyPseudoPotentials32 = 0;

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

      bool Padding=false;
      {
	int TmpI;
	InteractionDefinition.GetAsSingleInteger("PseudopotentialPadding", TmpI);
	if (TmpI>0)
	  Padding=true;
      }
	
      double* TmpThreeBodyPotentials = 0;
      int TmpNbrThreeBodyPseudoPotentials = 0;
      NbrThreeBodyPseudoPotentials32 = 0;
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentials32", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  ThreeBodyPotentials32 = new double[TmpNbrThreeBodyPseudoPotentials];
	  NbrThreeBodyPseudoPotentials32 = TmpNbrThreeBodyPseudoPotentials;	      
	  for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
	    ThreeBodyPotentials32[j] = TmpThreeBodyPotentials[j];
	}
      else
	{
	  if (InteractionDefinition["ThreebodyPseudopotentials32"] != 0)
	    {
	      cout << "ThreebodyPseudopotentials32 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentials12", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  ThreeBodyPotentials12 = new double[TmpNbrThreeBodyPseudoPotentials];
	  NbrThreeBodyPseudoPotentials12 = TmpNbrThreeBodyPseudoPotentials;	      
	  for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
	    ThreeBodyPotentials12[j] = TmpThreeBodyPotentials[j];
	}
      else
	{
	  if (InteractionDefinition["ThreebodyPseudopotentials12"] != 0)
	    {
	      cout << "ThreebodyPseudopotentials12 has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      if ((NbrThreeBodyPseudoPotentials32 == 0) && (NbrThreeBodyPseudoPotentials12 == 0))
	{
	  cout << "ThreebodyPseudopotentials32 or ThreebodyPseudopotentials12 is not defined in in " << Manager.GetString("interaction-file") << endl;
	  return -1;	  
	}
      --NbrThreeBodyPseudoPotentials32;
      --NbrThreeBodyPseudoPotentials12;


      if ((InteractionDefinition["Pseudopotentials"] != 0) || (InteractionDefinition["PseudopotentialsUpUp"] != 0)
	  || (InteractionDefinition["PseudopotentialsDownDown"] != 0) || (InteractionDefinition["PseudopotentialsUpDown"] != 0) || (InteractionDefinition["PseudopotentialsPairTunneling"] != 0))
	{
	  PseudoPotentials  = new double*[4];
	  for (int i = 0; i < 4; ++i)
	    {
	      PseudoPotentials[i] = new double[LzMax + 1];
	      for (int j = 0; j <= LzMax; ++j)
		PseudoPotentials[i][j] = 0.0;
	    };
	  if (FQHESphereSU2GetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials,
					       OneBodyPotentialUpUp, OneBodyPotentialDownDown) == false)
	    return -1;
	}
    }

  char* OutputBaseName = new char [256 + strlen(Manager.GetString("interaction-name"))];
  if (AllSz)
    {
      if (PairParity)
	sprintf (OutputBaseName, "bosons_sphere_su2_%s_parity_%d_pp_%g_n_%d_2s_%d_lz", Manager.GetString("interaction-name"), PairParity, Manager.GetDouble("pairing"),NbrParticles, LzMax);
      else
	sprintf (OutputBaseName, "bosons_sphere_su2_%s_pp_%g_n_%d_2s_%d_lz", Manager.GetString("interaction-name"), Manager.GetDouble("pairing"),NbrParticles, LzMax);
    }
  else
    sprintf (OutputBaseName, "bosons_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz", Manager.GetString("interaction-name"), NbrParticles, LzMax, SzTotal);
  
  char* OutputNameLz = new char [strlen(OutputBaseName)+5];
  sprintf (OutputNameLz,"%s.dat",OutputBaseName);
  
  int Max = (LzMax * NbrUp) + (LzMax * NbrDown);

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
      ParticleOnSphereWithSpin* Space = (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(L);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (AllSz==false)
	{
	  if (PseudoPotentials == 0)
	    {
	      Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials32, NbrThreeBodyPseudoPotentials32, ThreeBodyPotentials12, NbrThreeBodyPseudoPotentials12,
										    Architecture.GetArchitecture(), 
										    Memory, DiskCacheFlag,
										    LoadPrecalculationFileName);
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials32, NbrThreeBodyPseudoPotentials32, ThreeBodyPotentials12, NbrThreeBodyPseudoPotentials12,
										    PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown,
										    Architecture.GetArchitecture(), 
										    Memory, DiskCacheFlag,
										    LoadPrecalculationFileName);
	    }
	}
      else
	{
	  if (PseudoPotentials == 0)
	    {
	      Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing(Space, NbrParticles, LzMax, ThreeBodyPotentials32, NbrThreeBodyPseudoPotentials32, ThreeBodyPotentials12, NbrThreeBodyPseudoPotentials12, Manager.GetDouble("pairing"), PseudoPotentials,
										    Architecture.GetArchitecture(), 
										    Memory, DiskCacheFlag,
										    LoadPrecalculationFileName);		
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing(Space, NbrParticles, LzMax, ThreeBodyPotentials32, NbrThreeBodyPseudoPotentials32, ThreeBodyPotentials12, NbrThreeBodyPseudoPotentials12, Manager.GetDouble("pairing"), 
										    PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown,
										    Architecture.GetArchitecture(), 
										    Memory, DiskCacheFlag,
										    LoadPrecalculationFileName);
	    }
	  
	}
      if (Manager.GetDouble("l2-factor") != 0.0)
	Hamiltonian->AddL2(L, SzTotal, Manager.GetDouble("l2-factor"), ((unsigned long)Manager.GetInteger("l2-memory")) << 20); 
      if (Manager.GetDouble("s2-factor") != 0.0)
	Hamiltonian->AddS2(L, SzTotal, Manager.GetDouble("s2-factor"), ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
      double Shift = - 10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [strlen(OutputBaseName)+10];
	  sprintf (EigenvectorName, "%s_%d", OutputBaseName, L);
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
  delete [] OutputBaseName;
  delete [] OutputNameLz;
  return 0;
}

