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
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
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
      cout << "see man page for option syntax or type FQHESphereFermionsThreeBodyGeneric -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();  
  bool DiskCacheFlag = ((BooleanOption*) Manager["disk-cache"])->GetBoolean();
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

  double** ThreeBodyPotentials = new double* [4];
  int* NbrThreeBodyPseudoPotentials = new int[4];

  if (((SingleStringOption*) Manager["interaction-file"])->GetString() == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(((SingleStringOption*) Manager["interaction-file"])->GetString()) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}

      double* TmpThreeBodyPotentials = 0;
      int TmpNbrThreeBodyPseudoPotentials = 0;
      NbrThreeBodyPseudoPotentials[0] = 0;
      NbrThreeBodyPseudoPotentials[1] = 0;
      NbrThreeBodyPseudoPotentials[2] = 0;
      NbrThreeBodyPseudoPotentials[3] = 0;
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentials", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  for (int i = 0; i < 4; ++i)
	    {
	      ThreeBodyPotentials[i] = new double[TmpNbrThreeBodyPseudoPotentials];
	      NbrThreeBodyPseudoPotentials[i] = TmpNbrThreeBodyPseudoPotentials;	      
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[i][j] = TmpThreeBodyPotentials[j];
	    }
	}
      else
	if (InteractionDefinition["ThreebodyPseudopotentials"] != 0)
	  {
	    cout << "ThreebodyPseudopotentials has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentialsUpUpUp", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  if (NbrThreeBodyPseudoPotentials[0] == 0)
	    {
	      NbrThreeBodyPseudoPotentials[0] = TmpNbrThreeBodyPseudoPotentials;
	      ThreeBodyPotentials[0] = new double[TmpNbrThreeBodyPseudoPotentials];
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[0][j] = TmpThreeBodyPotentials[j];	      
	    }
	  else
	    {
	      if (NbrThreeBodyPseudoPotentials[0] < TmpNbrThreeBodyPseudoPotentials)
		{
		  delete[] ThreeBodyPotentials[0];
		  ThreeBodyPotentials[0] = new double[TmpNbrThreeBodyPseudoPotentials];
		}
	      NbrThreeBodyPseudoPotentials[0] = TmpNbrThreeBodyPseudoPotentials;
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[0][j] = TmpThreeBodyPotentials[j];	      
	    }
	}
      else
	if (InteractionDefinition["ThreebodyPseudopotentialsUpUpUp"] != 0)
	  {
	    cout << "ThreebodyPseudopotentialsUpUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentialsDownDownDown", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  if (NbrThreeBodyPseudoPotentials[1] == 0)
	    {
	      NbrThreeBodyPseudoPotentials[1] = TmpNbrThreeBodyPseudoPotentials;
	      ThreeBodyPotentials[1] = new double[TmpNbrThreeBodyPseudoPotentials];
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[1][j] = TmpThreeBodyPotentials[j];	      
	    }
	  else
	    {
	      if (NbrThreeBodyPseudoPotentials[1] < TmpNbrThreeBodyPseudoPotentials)
		{
		  delete[] ThreeBodyPotentials[1];
		  ThreeBodyPotentials[1] = new double[TmpNbrThreeBodyPseudoPotentials];
		}
	      NbrThreeBodyPseudoPotentials[1] = TmpNbrThreeBodyPseudoPotentials;
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[1][j] = TmpThreeBodyPotentials[j];	      
	    }
	}
      else
	if (InteractionDefinition["ThreebodyPseudopotentialsDownDownDown"] != 0)
	  {
	    cout << "ThreebodyPseudopotentialsDownDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentialsUpUpDown", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  if (NbrThreeBodyPseudoPotentials[2] == 0)
	    {
	      NbrThreeBodyPseudoPotentials[2] = TmpNbrThreeBodyPseudoPotentials;
	      ThreeBodyPotentials[2] = new double[TmpNbrThreeBodyPseudoPotentials];
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[2][j] = TmpThreeBodyPotentials[j];	      
	    }
	  else
	    {
	      if (NbrThreeBodyPseudoPotentials[2] < TmpNbrThreeBodyPseudoPotentials)
		{
		  delete[] ThreeBodyPotentials[2];
		  ThreeBodyPotentials[2] = new double[TmpNbrThreeBodyPseudoPotentials];
		}
	      NbrThreeBodyPseudoPotentials[2] = TmpNbrThreeBodyPseudoPotentials;
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[2][j] = TmpThreeBodyPotentials[j];	      
	    }
	}
      else
	if (InteractionDefinition["ThreebodyPseudopotentialsUpUpDown"] != 0)
	  {
	    cout << "ThreebodyPseudopotentialsUpUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentialsDownDownUp", ' ', TmpThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == true)
	{
	  if (NbrThreeBodyPseudoPotentials[3] == 0)
	    {
	      NbrThreeBodyPseudoPotentials[3] = TmpNbrThreeBodyPseudoPotentials;
	      ThreeBodyPotentials[3] = new double[TmpNbrThreeBodyPseudoPotentials];
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[3][j] = TmpThreeBodyPotentials[j];	      
	    }
	  else
	    {
	      if (NbrThreeBodyPseudoPotentials[3] < TmpNbrThreeBodyPseudoPotentials)
		{
		  delete[] ThreeBodyPotentials[3];
		  ThreeBodyPotentials[3] = new double[TmpNbrThreeBodyPseudoPotentials];
		}
	      NbrThreeBodyPseudoPotentials[3] = TmpNbrThreeBodyPseudoPotentials;
	      for (int j = 0; j < TmpNbrThreeBodyPseudoPotentials; ++j)
		ThreeBodyPotentials[3][j] = TmpThreeBodyPotentials[j];	      
	    }
	}
      else
	if (InteractionDefinition["ThreebodyPseudopotentialsDownDownUp"] != 0)
	  {
	    cout << "ThreebodyPseudopotentialsDownDownUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (NbrThreeBodyPseudoPotentials[0] == 0)
	{
	  cout << "ThreebodyPseudopotentials or ThreebodyPseudopotentialsUpUpUp is not defined in in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;	  
	}
      if (NbrThreeBodyPseudoPotentials[1] == 0)
	{
	  cout << "ThreebodyPseudopotentials or ThreebodyPseudopotentialsDownDownDown is not defined in in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;	  
	}
      if (NbrThreeBodyPseudoPotentials[2] == 0)
	{
	  cout << "ThreebodyPseudopotentials or ThreebodyPseudopotentialsUpUpDown is not defined in in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;	  
	}
      if (NbrThreeBodyPseudoPotentials[3] == 0)
	{
	  cout << "ThreebodyPseudopotentials or ThreebodyPseudopotentialsDownDownUp is not defined in in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	  return -1;	  
	}
      --NbrThreeBodyPseudoPotentials[0];
      --NbrThreeBodyPseudoPotentials[1];
      --NbrThreeBodyPseudoPotentials[2];
      --NbrThreeBodyPseudoPotentials[3];

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
	    cout << "Pseudopotentials has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
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
	    cout << "PseudopotentialsUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
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
	    cout << "PseudopotentialsDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
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
	    cout << "PseudopotentialsUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', OneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', OneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	      return -1;
	    }
	}
    }

  char* OutputNameLz = new char [256 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
  sprintf (OutputNameLz, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, SzTotal);

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
  for (; L <= Max; L += 2)
    {
      ParticleOnSphereWithSpin* Space = (ParticleOnSphereWithSpin*) ParticleManager.GetHilbertSpace(L);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (PseudoPotentials == 0)
	Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials, NbrThreeBodyPseudoPotentials,
									      Architecture.GetArchitecture(), 
									      Memory, DiskCacheFlag,
									      LoadPrecalculationFileName);
      else
	Hamiltonian = new ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials, NbrThreeBodyPseudoPotentials,
									      PseudoPotentials, 0, 0,
									      Architecture.GetArchitecture(), 
									      Memory, DiskCacheFlag,
									      LoadPrecalculationFileName);
      if (((SingleDoubleOption*) Manager["l2-factor"])->GetDouble() != 0.0)
	Hamiltonian->AddL2(L, SzTotal, ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble());
      if (((SingleDoubleOption*) Manager["s2-factor"])->GetDouble() != 0.0)
	Hamiltonian->AddS2(L, SzTotal, ((SingleDoubleOption*) Manager["s2-factor"])->GetDouble());
      double Shift = - 10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, SzTotal, L);
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
