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

#include "Hamiltonian/ParticleOnSphereWithSpinGenericThreeBodyHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

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
  OptionManager Manager ("FQHESphereFermionsWithSpinThreeBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  ParticleOnSphereManager ParticleManager(true, false, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  

  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  // (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
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
  bool FirstRun = true;

  int NbrUp = (NbrParticles + SzTotal) >> 1;
  int NbrDown = (NbrParticles - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
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
	  PseudoPotentials = new double*[3];
	  for (int i = 0; i < 3; ++i)
	    {  
	      PseudoPotentials[i] = new double[LzMax + 1];
	      for (int j = 0; j <= LzMax; ++j)
		PseudoPotentials[i][j] = TmpPseudoPotentials[j];
	    }
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
	  if (PseudoPotentials == 0)
	    {
	      PseudoPotentials = new double*[3];
	      for (int i = 0; i < 3; ++i)
		PseudoPotentials[i] = new double[LzMax + 1];
	      for (int j = 0; j <= LzMax; ++j)
		PseudoPotentials[0][j] = TmpPseudoPotentials[j];
	    }
	  else
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
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDown" << endl;
	      return -1;	  
	    }
	  if (PseudoPotentials == 0)
	    {
	      cout << "Pseudopotentials and PseudopotentialsUpUp are not defined" << endl;
	      return -1;	  	      
	    }
	  else
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
	  if (PseudoPotentials == 0)
	    {
	      cout << "Pseudopotentials,PseudopotentialsUpUp and PseudopotentialsDownDown are not defined" << endl;
	      return -1;	  	      
	    }
	  else
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
	      cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
    }

  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  sprintf (OutputNameLz, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz.dat", Manager.GetString("interaction-name"), NbrParticles, LzMax, SzTotal);

  int Max = (((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown));

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

  cout << "WARNING : FQHESphereFermionsThreeBodyGeneric has not been tested, use at your own risk" << endl;
  for (; L <= Max; L += 2)
    {
      ParticleOnSphereWithSpin* Space = (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(L);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (PseudoPotentials == 0)
	Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials32, NbrThreeBodyPseudoPotentials32, ThreeBodyPotentials12, NbrThreeBodyPseudoPotentials12,
									      Architecture.GetArchitecture(), 
									      Memory, DiskCacheFlag,
									      LoadPrecalculationFileName);
      else
	Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials32, NbrThreeBodyPseudoPotentials32, ThreeBodyPotentials12, NbrThreeBodyPseudoPotentials12,
									      PseudoPotentials, 0, 0,
									      Architecture.GetArchitecture(), 
									      Memory, DiskCacheFlag,
									      LoadPrecalculationFileName);
      if (Manager.GetDouble("l2-factor") != 0.0)
	Hamiltonian->AddL2(L, SzTotal, Manager.GetDouble("l2-factor"), ((unsigned long)Manager.GetInteger("l2-memory")) << 20); 
      if (Manager.GetDouble("s2-factor") != 0.0)
	Hamiltonian->AddS2(L, SzTotal, Manager.GetDouble("s2-factor"), ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
      double Shift = - 10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, SzTotal, L);
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
