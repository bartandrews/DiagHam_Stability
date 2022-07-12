#include "HilbertSpace/BosonOnS2xS2.h"
#include "HilbertSpace/BosonOnS2xS2HardcoreNoNearestNeighbors.h"

#include "Hamiltonian/ParticleOnCylinderxCylinderDeltaHamiltonian.h"
#include "Hamiltonian/ParticleOnCylinderxCylinderDeltaTruncatedHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/StringTools.h"

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
  OptionManager Manager ("FQHECylinderxCylinderBosonsDelta" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux1", "number of flux quanta for the first sphere", 0);
  (*SystemGroup) += new SingleIntegerOption  ('k', "nbr-flux2", "number of flux quanta for the second sphere", 0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio1", "aspect ratio of the first cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter1", "if non zero, fix the first cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio2", "aspect ratio of the second cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter2", "if non zero, fix the second cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "inital value for the z projection of the first sphere angular momentum (Lz)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-kz", "inital value for the projection of the first sphere angular momentum (Kz)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of Lz sectors to compute (0 if all possible sectors have to be computed)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-kz", "number of Kz sectors to compute (0 if all possible sectors have to be computed)", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-exclusion", "forbid any configuration with orbital multiple occupancy");
  
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
      cout << "see man page for option syntax or type FQHECylinderxCylinderBosonsDelta -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrBosons = Manager.GetInteger("nbr-particles");
  int NbrFluxQuanta1 = Manager.GetInteger("nbr-flux1");
  int NbrFluxQuanta2 = Manager.GetInteger("nbr-flux2");
  double Ratio1 = Manager.GetDouble("aspect-ratio1");
  double Perimeter1 = Manager.GetDouble("cylinder-perimeter1");
  if (Perimeter1 != 0.0)
    {
      Ratio1 = 2.0 * M_PI * (NbrFluxQuanta1 + 1) / (Perimeter1 * Perimeter1);
    }
  double Ratio2 = Manager.GetDouble("aspect-ratio2");
  double Perimeter2 = Manager.GetDouble("cylinder-perimeter2");
  if (Perimeter2 != 0.0)
    {
      Ratio2 = 2.0 * M_PI * (NbrFluxQuanta2 + 1) / (Perimeter2 * Perimeter2);
    }
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  
  bool FirstRun = true;

  char* GeometryName = new char[256];
  if (Perimeter1 > 0.0)	
    {
      sprintf (GeometryName, "cylinder_perimeter1_%.6f_perimeter2_%.6f", Perimeter1, Perimeter2);
    }
  else
    {
      sprintf (GeometryName, "cylinder_ratio1_%.6f_ratio2_%.6f", Ratio1, Ratio2);
    }
  char* OutputName = new char [256 + strlen(GeometryName)];
  sprintf (OutputName, "bosons_%s_delta_n_%d_2s1_%d_2s2_%d.dat", GeometryName, NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2);

  int MaxTotalLz = NbrFluxQuanta1 * NbrBosons;
  int MaxTotalKz = NbrFluxQuanta2 * NbrBosons;
  int MinTotalLz = MaxTotalLz & 1;
  int MinTotalKz = MaxTotalKz & 1;
  if (Manager.GetInteger("initial-lz") != 0)
    {
      MinTotalLz = Manager.GetInteger("initial-lz");
    }
  if (Manager.GetInteger("initial-kz") != 0)
    {
      MinTotalKz = Manager.GetInteger("initial-kz");
    }
  if (Manager.GetInteger("nbr-lz") > 0)
    {
      MaxTotalLz = MinTotalLz + 2 * (Manager.GetInteger("nbr-lz") - 1);
    }
  if (Manager.GetInteger("nbr-kz") > 0)
    {
      MaxTotalKz = MinTotalKz + 2 * (Manager.GetInteger("nbr-kz") - 1);
    }

  for (int TotalLz = MinTotalLz; TotalLz <= MaxTotalLz; TotalLz += 2)
    {
      for (int TotalKz = MinTotalKz; TotalKz <= MaxTotalKz; TotalKz += 2)
	{
	  ParticleOnSphere* Space = 0;
	  if (Manager.GetBoolean("use-exclusion") == false)
	    {
	      Space = new BosonOnS2xS2(NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, TotalLz, TotalKz);
	    }
	  else
	    {
	      Space = new BosonOnS2xS2HardcoreNoNearestNeighbors(NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, TotalLz, TotalKz);
	    }
     
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  //       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	  // 	Space->PrintState(cout, i);
	  
	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if (Manager.GetBoolean("use-exclusion") == false)
	    {
	      Hamiltonian = new ParticleOnCylinderxCylinderDeltaHamiltonian(Space, NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, Ratio1, Ratio2,
									    Architecture.GetArchitecture(), Memory);
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnCylinderxCylinderDeltaTruncatedHamiltonian(Space, NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, Ratio1, Ratio2,
										     Architecture.GetArchitecture(), Memory);
	    }

      
	  char* EigenvectorName = 0;
	  if (Manager.GetBoolean("eigenstate") == true)	
	    {
	      char* TmpVectorExtension = new char [64];
	      sprintf (TmpVectorExtension, "_lz_%d_kz_%d", TotalLz, TotalKz);
	      EigenvectorName = ReplaceString(OutputName, ".dat", TmpVectorExtension);
	    }
	
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", TotalLz, TotalKz);	  
	  char* SubspaceLegend = new char[256];
	  sprintf (SubspaceLegend, "Lz Kz");

	  GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, 0, OutputName, FirstRun, EigenvectorName);
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
    }
  return 0;
}
