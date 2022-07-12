#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpin.h"

#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "MainTask/GenericRealMainTask.h"

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
  OptionManager Manager ("FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "if positive, fix the total number of particles", -1);
  (*SystemGroup) += new BooleanOption ('\n', "all-fixednbrparticles", "do the calculation for all the fixed particle number sector from 0 to --nbr-particles");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleDoubleOption ('\n', "charging-energy", "factor in front of the charging energy (i.e 1/(2C))", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "average-nbrparticles", "average number of particles", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "force-negativelz", "manually force to compute the negative lz sectors");
  (*SystemGroup) += new SingleStringOption ('\n', "directory", "use a specific directory for the input data instead of the current one");
  (*SystemGroup) += new BooleanOption  ('\n', "use-cylinder", "use the cylinder geometry intead of the sphere geometry");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-bosons", "use bosons instead of fermions");

  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*SystemGroup) += new BooleanOption  ('\n', "compute-sparsity", "compute sparsity of Hamiltonian matrix");
  
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-binhamiltonian", "export the hamiltonian as a binary file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereQuasiholesWithSpinTimeReversalSymmetryAndPairing -h" << endl;
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
  bool Statistics = !Manager.GetBoolean("use-bosons");
  
  int KValue = 1;
  int RValue = 2;

  bool UseCylinderFlag = Manager.GetBoolean("use-cylinder");
  double Ratio = Manager.GetDouble("aspect-ratio");
  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter > 0.0)
    {
      Ratio = (Perimeter * Perimeter) / (2.0 * M_PI * (LzMax + 1));
    }  

  char* FilePrefix = new char[512];

  if (UseCylinderFlag == true)
    {
      if (Perimeter > 0.0)	
	{
	  if (Statistics == true)
	    {
	      sprintf (FilePrefix, "fermions_cylinder_perimeter_%.6f", Perimeter);
	    }
	  else
	    {
	      sprintf (FilePrefix, "bosons_cylinder_perimeter_%.6f", Perimeter);
	    }
	}
      else
	{
	  if (Statistics == true)
	    {
	      sprintf (FilePrefix, "fermions_cylinder_ratio_%.6f", Ratio);
	    }
	  else
	    {
	      sprintf (FilePrefix, "bosons_cylinder_ratio_%.6f", Ratio);
	    }      
	}
    }
  else
    {
      if (Statistics == true)
	{
	  sprintf (FilePrefix, "fermions");
	}
      else
	{
	  sprintf (FilePrefix, "bosons");
	}
    }

  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double* OneBodyPotentialUpDown = 0;
  double* OneBodyPotentialPairing = 0;
  double** PseudoPotentials  = 0;
//   for (int i = 0; i < 3; ++i)
//     {
//       PseudoPotentials[i] = new double[LzMax + 1];
//     }

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
  if (PseudoPotentials != 0)
    {
      cout << "warning, there should be no two-body pseudopotentials. Two-body interactions are implemented through the restriction to the quasihole Hilbert space" << endl;
    }

  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  int TmpFixedNbrParticles = 0;
  if (Manager.GetInteger("nbr-particles") >= 0)
    TmpFixedNbrParticles = Manager.GetInteger("nbr-particles");
  sprintf (OutputNameLz, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_%d_2s_%d_sz_%d_lz.dat", FilePrefix, Manager.GetString("interaction-name"), 
	   Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), TmpFixedNbrParticles, LzMax, TotalSz);

  int MinNbrParticles = abs(TotalSz);
  int MaxNbrParticles = (2 * (LzMax + 1)) - abs(TotalSz);
  if (Manager.GetInteger("nbr-particles") >= 0)
    {
      MinNbrParticles = Manager.GetInteger("nbr-particles");
      MaxNbrParticles = MinNbrParticles;
      if (OneBodyPotentialPairing != 0)
	{
	  cout << "discarding superconducting coupling" << endl;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      OneBodyPotentialPairing[i] = 0.0;
	    }
	}
    }

  int MaxL = 0;
  for (int TmpNbrParticles = MinNbrParticles; TmpNbrParticles <= MaxNbrParticles; TmpNbrParticles += 2)
    {
      int NbrUp = (TmpNbrParticles + TotalSz) / 2;
      int NbrDown = (TmpNbrParticles - TotalSz) / 2;
      if ((NbrUp <= (LzMax + 1)) && (NbrUp >= 0) && (NbrDown <= (LzMax + 1)) && (NbrDown >= 0))
	{
	  // warning, this should not work for k > 1
	  if (KValue > 1)
	    cout << "please fix your code for k>1" << endl;
	  int MaxTotalLzUp = (LzMax * NbrUp) - ((KValue + RValue) * ((NbrUp - 1) * NbrUp) / KValue);
	  int MaxTotalLzDown = (LzMax * NbrDown) - ((KValue + RValue) * ((NbrDown - 1) * NbrDown) / KValue);
	  if ((MaxTotalLzUp + MaxTotalLzDown) > MaxL)
	    {
	      MaxL = (MaxTotalLzUp + MaxTotalLzDown);
	    }
	}
    }

  int  L = 0;
  L = InitialLz;
  if ((abs(MaxL) & 1) != 0)
    L |= 1;
  else
    L &= ~0x1;

  if (NbrLz > 0)
    {
      if (L + (2 * (NbrLz - 1)) < MaxL)
	MaxL = L + (2 * (NbrLz - 1));
    }

  if (Manager.GetBoolean("force-negativelz"))
    L = -MaxL;


  for (; L <= MaxL; L += 2)
    {
      int LocalNbrParticles = -1;
      int MaxLocalNbrParticles = -1;
      if (Manager.GetInteger("nbr-particles") >= 0)
	{
	  LocalNbrParticles = Manager.GetInteger("nbr-particles");
	  MaxLocalNbrParticles = LocalNbrParticles;
	  if (Manager.GetBoolean("all-fixednbrparticles") == true)
	    {
	      if ((TotalSz & 1) == 0)
		LocalNbrParticles = 0;
	      else
		LocalNbrParticles = 1;
	    }
	}
      for (; LocalNbrParticles <= MaxLocalNbrParticles; LocalNbrParticles += 2)
	{
	  QuasiholeOnSphereWithSpinAndPairing* Space = 0;
	  if (LocalNbrParticles >= 0)
	    {
	      Space = new QuasiholeOnSphereWithSpin (KValue, RValue, L, LzMax, LocalNbrParticles, TotalSz, Manager.GetString("directory"), FilePrefix);
	    }
	  else
	    {
	      Space = new QuasiholeOnSphereWithSpinAndPairing (KValue, RValue, L, LzMax, TotalSz, Manager.GetString("directory"), FilePrefix);
	    }
	  if (Space->GetLargeHilbertSpaceDimension() > 0l)
	    {
	      if (UseCylinderFlag == true)
		{
		  cout << "Ky=" << L;
		}
	      else
		{
		  cout << "Lz=" << L;
		}
	      if (Manager.GetInteger("nbr-particles") >= 0)
		cout << ", N=" << LocalNbrParticles;
	      cout << endl;

	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      AbstractHamiltonian* Hamiltonian = 0;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Hamiltonian = new ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing (Space, LzMax, OneBodyPotentialUpUp, 
													     OneBodyPotentialDownDown,
													     OneBodyPotentialPairing,
													     Manager.GetDouble("charging-energy"), 
													     Manager.GetDouble("average-nbrparticles"),
													     Architecture.GetArchitecture(), 
													     Memory, DiskCacheFlag,
													     LoadPrecalculationFileName);
	      
	      double Shift = 0.0;
	      Hamiltonian->ShiftHamiltonian(Shift);
	      char* EigenvectorName = 0;
	      if ((Manager.GetBoolean("eigenstate") == true) || (Manager.GetBoolean("all-eigenstates") == true))
		{
		  EigenvectorName = new char [256 + strlen(FilePrefix) + strlen(Manager.GetString("interaction-name"))];
		  if (LocalNbrParticles >= 0)
		    {
		      sprintf (EigenvectorName, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_%d_2s_%d_sz_%d_lz_%d", 
			       FilePrefix, Manager.GetString("interaction-name"), 
			       Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LocalNbrParticles, LzMax, TotalSz, L);
		    }
		  else
		    {
		      sprintf (EigenvectorName, "%s_su2_quasiholes_%s_cenergy_%.6f_n0_%.6f_pairing_n_0_2s_%d_sz_%d_lz_%d", 
			       FilePrefix, Manager.GetString("interaction-name"), 
			       Manager.GetDouble("charging-energy"), Manager.GetDouble("average-nbrparticles"), LzMax, TotalSz, L);
		    }
		}
	      
	      char* ContentPrefix = new char[256];
	      if (Manager.GetBoolean("all-fixednbrparticles") == true)
		{
		  sprintf (ContentPrefix, "%d %d", LocalNbrParticles, L);
		}
	      else
		{
		  sprintf (ContentPrefix, "%d", L);
		}
	      char* SubspaceLegend = new char[256];
	      if (Manager.GetBoolean("all-fixednbrparticles") == true)
		{
		  sprintf (SubspaceLegend, "N lz");
		}
	      else
		{
		  sprintf (SubspaceLegend, "lz");
		}

	      GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, Shift, OutputNameLz, FirstRun, EigenvectorName);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      if (EigenvectorName != 0)
		{
		  delete[] EigenvectorName;
		}
	      delete Hamiltonian;
	      if (FirstRun == true)
		FirstRun = false;
	    }
	  delete Space;
	}
    }
  return 0;
}
