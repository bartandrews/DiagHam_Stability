#include "HilbertSpace/BosonOnT2xT2.h"
//#include "HilbertSpace/BosonOnT2xT2Long.h"
//#include "HilbertSpace/BosonOnT2xT2HardcoreNoNearestNeighbors.h"

#include "Hamiltonian/ParticleOnT2xT2GenericTwoBodyHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("FQHET2xT2BosonsTwoBodyGeneric" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption ('\n', "ratio1", "ratio between the width in the x direction and the width in the y direction for the first torus", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "ratio2", "ratio between the width in the x direction and the width in the y direction for the second torus", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-ky1", "inital value for the ky momentum of the first torus", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-ky2", "inital value for the ky momentum of the second torus", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ky1", "number of Ky1 sectors to compute (0 if all possible sectors have to be computed)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ky2", "number of Ky2 sectors to compute (0 if all possible sectors have to be computed)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHET2xT2BosonsTwoBodyGeneric -h" << endl;
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
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  double Ratio1 = Manager.GetDouble("ratio1");
  double Ratio2 = Manager.GetDouble("ratio2");
  
  int NbrPseudoPotentials = 0;
  int* PseudoPotentialMomentum1;
  int* PseudoPotentialMomentum2;
  double* PseudoPotentials;
  char* InteractioName = 0;
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "no interaction file has been provided, assuming delta interaction" << endl;
      InteractioName = new char[8];
      sprintf(InteractioName, "delta");
      NbrPseudoPotentials = 1;
      PseudoPotentialMomentum1 = new int[NbrPseudoPotentials];
      PseudoPotentialMomentum2 = new int[NbrPseudoPotentials];
      PseudoPotentials = new double[NbrPseudoPotentials];
      PseudoPotentialMomentum1[0] = 0;
      PseudoPotentialMomentum2[0] = 0;
      PseudoPotentials[0] = 1.0;
    }
  else
    {
      MultiColumnASCIIFile InteractionFile;
      if (InteractionFile.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionFile.DumpErrors(cout);
	  return -1;
	}
      if (InteractionFile.GetNbrColumns() < 3)
	{
	  cout << "error, wrong number of columns in " << Manager.GetString("interaction-file") << endl;
	  return -1;
	}
      NbrPseudoPotentials = InteractionFile.GetNbrLines();
      if (NbrPseudoPotentials <= 0)
	{
	  cout << "error, no pseudo-potential defined in " << Manager.GetString("interaction-file") << endl;
	  return -1;
	}
      PseudoPotentialMomentum1 = InteractionFile.GetAsIntegerArray(0);
      if (PseudoPotentialMomentum1 == 0)
	{
	  InteractionFile.DumpErrors(cout);
	  return -1;	  
	}
      PseudoPotentialMomentum2 = InteractionFile.GetAsIntegerArray(1);
      if (PseudoPotentialMomentum2 == 0)
	{
	  InteractionFile.DumpErrors(cout);
	  return -1;	  
	}
      PseudoPotentials = InteractionFile.GetAsDoubleArray(2);
      if (PseudoPotentials == 0)
	{
	  InteractionFile.DumpErrors(cout);
	  return -1;	  
	}
      InteractioName = new char[strlen(Manager.GetString("interaction-file")) + 1];
      strcpy (InteractioName, Manager.GetString("interaction-file"));
   }


  bool FirstRun = true;

  char* OutputName = new char [256 + strlen(InteractioName)];
  sprintf (OutputName, "bosons_t2xt2_ratio1_%.6f_ratio2_%.6f_%s_n_%d_2s1_%d_2s2_%d.dat", Ratio1, Ratio2, InteractioName, NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2);

  int MaxTotalKy1 = NbrFluxQuanta1 - 1;
  int MaxTotalKy2 = NbrFluxQuanta2 - 1;
  int MinTotalKy1 = 0;
  int MinTotalKy2 = 0;
  if (Manager.GetInteger("initial-ky1") != 0)
    {
      MinTotalKy1 = Manager.GetInteger("initial-ky1");
    }
  if (Manager.GetInteger("initial-ky2") != 0)
    {
      MinTotalKy2 = Manager.GetInteger("initial-ky2");
    }
  if (Manager.GetInteger("nbr-ky1") > 0)
    {
      MaxTotalKy1 = MinTotalKy1 + Manager.GetInteger("nbr-ky1") - 1;
    }
  if (Manager.GetInteger("nbr-ky2") > 0)
    {
      MaxTotalKy2 = MinTotalKy2 + Manager.GetInteger("nbr-ky2") - 1;
    }

  for (int TotalKy1 = MinTotalKy1; TotalKy1 <= MaxTotalKy1; TotalKy1++)
    {
      int TmpMaxTotalKy2 = MaxTotalKy2;
      if ((NbrFluxQuanta1 == NbrFluxQuanta2) && (Manager.GetInteger("nbr-ky1") < 0) && (Manager.GetInteger("nbr-ky2") < 0)) 
	{
	  TmpMaxTotalKy2 = MinTotalKy1;
	}
      for (int TotalKy2 = MinTotalKy2; TotalKy2 <= TmpMaxTotalKy2; TotalKy2++)
	{
	  ParticleOnSphere* Space = 0;
#ifdef __128_BIT_LONGLONG__
	  if ((((NbrFluxQuanta1 + 1) * (NbrFluxQuanta2 + 1)) + NbrBosons) <= 63)
#else
	  if ((((NbrFluxQuanta1 + 1) * (NbrFluxQuanta2 + 1)) + NbrBosons) <= 31)	    
#endif
	    {
	      Space = new BosonOnT2xT2(NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, TotalKy1, TotalKy2);
	    }
	  else
	    {
	      Space = 0;
//	      Space = new BosonOnT2xT2Long(NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, TotalKy1, TotalKy2);
	    }
	  if (Space->GetHilbertSpaceDimension() > 0)
	    {
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      //       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	      // 	Space->PrintState(cout, i);
	      
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnT2xT2GenericTwoBodyHamiltonian(Space, NbrBosons, NbrFluxQuanta1, NbrFluxQuanta2, Ratio1, Ratio2,
									 NbrPseudoPotentials, PseudoPotentialMomentum1, PseudoPotentialMomentum2, 
									 PseudoPotentials, Architecture.GetArchitecture(), Memory);
	      
	      char* EigenvectorName = 0;
	      if (Manager.GetBoolean("eigenstate") == true)	
		{
		  char* TmpVectorExtension = new char [64];
		  sprintf (TmpVectorExtension, "_kz1_%d_kz2_%d", TotalKy1, TotalKy2);
		  EigenvectorName = ReplaceString(OutputName, ".dat", TmpVectorExtension);
		}
	      
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", TotalKy1, TotalKy2);	  
	      char* SubspaceLegend = new char[256];
	      sprintf (SubspaceLegend, "Ky1 Ky2");
	      
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
    }
  return 0;
}
