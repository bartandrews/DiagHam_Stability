#include "HilbertSpace/BosonOnT2xS2WithMagneticTranslationsShort.h"
#include "HilbertSpace/BosonOnT2xS2.h"
//#include "HilbertSpace/BosonOnT2xS2WithMagneticTranslations00Long.h"
//#include "HilbertSpace/BosonOnT2xS2WithMagneticTranslationsShortHardcoreNoNearestNeighbors.h"

#include "Hamiltonian/ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnT2xS2GenericTwoBodyHamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

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
  OptionManager Manager ("FQHET2xS2BosonsTwoBodyGeneric" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux1", "number of flux quanta for the torus", 0);
  (*SystemGroup) += new SingleIntegerOption  ('k', "nbr-flux2", "number of flux quanta for the sphere", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "ratio", "ratio between the width in the x direction and the width in the y direction for the torus", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total torus momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total torus momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "initial value for the Lz angular momentum on the sphere", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of Lz sectors to compute (0 if all possible sectors have to be computed)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "do not consider magnetic translation (only trivial translations along one axis)");
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the reduced Brillouin zone", false);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
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
      cout << "see man page for option syntax or type FQHET2xS2BosonsTwoBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrBosons = Manager.GetInteger("nbr-particles");
  int NbrFluxQuantumTorus = Manager.GetInteger("nbr-flux1");
  int NbrFluxQuantumSphere = Manager.GetInteger("nbr-flux2");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  double Ratio = Manager.GetDouble("ratio");
  
  int NbrPseudoPotentials = 0;
  int* PseudoPotentialMomentumTorus;
  int* PseudoPotentialMomentumSphere;
  double* PseudoPotentials;
  char* InteractioName = 0;
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "no interaction file has been provided, assuming delta interaction" << endl;
      InteractioName = new char[8];
      sprintf(InteractioName, "delta");
      NbrPseudoPotentials = 1;
      PseudoPotentialMomentumTorus = new int[NbrPseudoPotentials];
      PseudoPotentialMomentumSphere = new int[NbrPseudoPotentials];
      PseudoPotentials = new double[NbrPseudoPotentials];
      PseudoPotentialMomentumTorus[0] = 0;
      PseudoPotentialMomentumSphere[0] = 0;
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
      PseudoPotentialMomentumTorus = InteractionFile.GetAsIntegerArray(0);
      if (PseudoPotentialMomentumTorus == 0)
	{
	  InteractionFile.DumpErrors(cout);
	  return -1;	  
	}
      PseudoPotentialMomentumSphere = InteractionFile.GetAsIntegerArray(1);
      if (PseudoPotentialMomentumSphere == 0)
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

  int XMomentum = Manager.GetInteger("x-momentum");
  char* OutputName = new char [256 + strlen(InteractioName)];
  if (Manager.GetBoolean("no-translation") == false)
    {
      sprintf (OutputName, "bosons_t2xs2_ratio_%.6f_%s_n_%d_2s1_%d_2s2_%d.dat", Ratio, InteractioName, NbrBosons, NbrFluxQuantumTorus, NbrFluxQuantumSphere);
    }
  else
    {
      XMomentum = 0;
      Lanczos.SetRealAlgorithms();
      sprintf (OutputName, "bosons_kysym_t2xs2_ratio_%.6f_%s_n_%d_2s1_%d_2s2_%d.dat", Ratio, InteractioName, NbrBosons, NbrFluxQuantumTorus, NbrFluxQuantumSphere);
   }
  int YMomentum = Manager.GetInteger("y-momentum");
  int MomentumModulo = FindGCD(NbrBosons, NbrFluxQuantumTorus);
  int MaxTotalLz = NbrFluxQuantumSphere * NbrBosons;
  int MinTotalLz = MaxTotalLz & 1;
  if (Manager.GetInteger("initial-lz") != 0)
    {
      MinTotalLz = Manager.GetInteger("initial-lz");
    }
  if (Manager.GetInteger("nbr-lz") > 0)
    {
      MaxTotalLz = MinTotalLz + 2 * (Manager.GetInteger("nbr-lz") - 1);
    }
  int NbrMomenta;
  int* KxMomenta;
  int* KyMomenta;
  int* LzMomenta;
  if ((XMomentum >= 0) && (YMomentum >= 0))
    {
      NbrMomenta = ((MaxTotalLz - MinTotalLz) / 2 + 1);
      KxMomenta = new int[NbrMomenta];
      KyMomenta = new int[NbrMomenta];
      LzMomenta = new int[NbrMomenta];
      NbrMomenta = 0;
      for (int TmpLz = MinTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
	{
	  KxMomenta[NbrMomenta] = XMomentum;
	  KyMomenta[NbrMomenta] = YMomentum;
	  LzMomenta[NbrMomenta] = TmpLz;
	  NbrMomenta++;
	}
    }
  else
    {
      if (Manager.GetString("selected-points") == 0)
	{
	  if (XMomentum >= 0)
	    {
	      NbrMomenta = NbrFluxQuantumTorus * ((MaxTotalLz - MinTotalLz) / 2 + 1);
	      KxMomenta = new int[NbrMomenta];
	      KyMomenta = new int[NbrMomenta];
	      LzMomenta = new int[NbrMomenta];
	      NbrMomenta = 0;
	      for (int y = 0; y < NbrFluxQuantumTorus; ++y)
		{
		  for (int TmpLz = MinTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
		    {
		      KxMomenta[NbrMomenta] = XMomentum;
		      KyMomenta[NbrMomenta] = y;
		      LzMomenta[NbrMomenta] = TmpLz;
		      NbrMomenta++;
		    }
		}
	    }
	  else
	    {
	      if (YMomentum >= 0)
		{
		  NbrMomenta = MomentumModulo * ((MaxTotalLz - MinTotalLz) / 2 + 1);
		  KxMomenta = new int[NbrMomenta];
		  KyMomenta = new int[NbrMomenta];
		  LzMomenta = new int[NbrMomenta];
		  NbrMomenta = 0;
		  for (int x = 0; x < MomentumModulo; ++x)
		    {
		      for (int TmpLz = MinTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
			{
			  KxMomenta[NbrMomenta] = x;
			  KyMomenta[NbrMomenta] = YMomentum;
			  LzMomenta[NbrMomenta] = TmpLz;
			  NbrMomenta++;
			}
		    }
		}
	      else
		{
		  int TmpMax = MomentumModulo;
		  if ((Ratio == 1.0) && (Manager.GetBoolean("full-reducedbz") == false))
		    {
			  int TmpMax = (MomentumModulo + 1) / 2;
		    }
		  NbrMomenta = TmpMax * TmpMax * ((MaxTotalLz - MinTotalLz) / 2 + 1);
		  KxMomenta = new int[NbrMomenta];
		  KyMomenta = new int[NbrMomenta];
		  LzMomenta = new int[NbrMomenta];
		  NbrMomenta = 0;
		  for (int x = 0; x < TmpMax; ++x)
		    {
		      int y = 0;
		      if ((Ratio == 1.0) && (Manager.GetBoolean("full-reducedbz") == false))
			{
			  y = x;
			}
		      for (; y < TmpMax; ++y)
			{
			  for (int TmpLz = MinTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
			    {
			      KxMomenta[NbrMomenta] = x;
			      KyMomenta[NbrMomenta] = y;
			      LzMomenta[NbrMomenta] = TmpLz;
			      NbrMomenta++;
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  MultiColumnASCIIFile MomentumFile;
	  if (MomentumFile.Parse(Manager.GetString("selected-points")) == false)
	    {
	      MomentumFile.DumpErrors(cout);
	      return -1;
	    }
	  NbrMomenta = MomentumFile.GetNbrLines();
	  KxMomenta = MomentumFile.GetAsIntegerArray(0);
	  KyMomenta = MomentumFile.GetAsIntegerArray(1);
	}
    }


  for (int Pos = 0; Pos < NbrMomenta; ++Pos)
    {
      int TotalKx = KxMomenta[Pos];
      int TotalKy = KyMomenta[Pos];
      int TotalLz = LzMomenta[Pos];
      if (Manager.GetBoolean("no-translation") == false)
	{
	  cout << "(Kx=" << TotalKx << ", Ky=" << TotalKy << ", Lz=" << TotalLz << ")" << endl;
	  
	  ParticleOnTorusWithMagneticTranslations* Space = 0;
#ifdef __128_BIT_LONGLONG__
	  if ((((NbrFluxQuantumTorus + 1) * (NbrFluxQuantumSphere + 1)) + NbrBosons) <= 63)
#else
	    if ((((NbrFluxQuantumTorus + 1) * (NbrFluxQuantumSphere + 1)) + NbrBosons) <= 31)	    
#endif
	      {
		Space = new BosonOnT2xS2WithMagneticTranslationsShort(NbrBosons, NbrFluxQuantumTorus, TotalKx, TotalKy, NbrFluxQuantumSphere, TotalLz);
	      }
	    else
	      {
		Space = 0;
		//	      Space = new BosonOnT2xS2WithMagneticTranslationsLong(NbrBosons, NbrFluxQuantumTorus, NbrFluxQuantumSphere, TotalKy, TotalLz);
	      }
	  if (Space->GetHilbertSpaceDimension() > 0)
	    {
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian(Space, NbrBosons, NbrFluxQuantumTorus, TotalKx, Ratio, NbrFluxQuantumSphere,
												 NbrPseudoPotentials, PseudoPotentialMomentumTorus, PseudoPotentialMomentumSphere, 
												 PseudoPotentials, Architecture.GetArchitecture(), Memory);
	      
	      char* EigenvectorName = 0;
	      if (Manager.GetBoolean("eigenstate") == true)	
		{
		  char* TmpVectorExtension = new char [64];
		  sprintf (TmpVectorExtension, "_kx_%d_ky_%d_lz_%d", TotalKx, TotalKy, TotalLz);
		  EigenvectorName = ReplaceString(OutputName, ".dat", TmpVectorExtension);
		}
	      
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d %d", TotalKx, TotalKy, TotalLz);	  
	      char* SubspaceLegend = new char[256];
	      sprintf (SubspaceLegend, "Kx Ky Lz");
	      
	      GenericComplexMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, 0, OutputName, FirstRun, EigenvectorName);
	      MainTaskOperation TaskOperation (&Task);	  	 
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      delete Hamiltonian;
	      delete Space;
	      if (EigenvectorName != 0)
		{
		  delete[] EigenvectorName;
		}
	    }
	}
      else
	{
	  cout << "(Ky=" << TotalKy << ", Lz=" << TotalLz << ")" << endl;
	  
	  ParticleOnSphere* Space = 0;
#ifdef __128_BIT_LONGLONG__
	  if ((((NbrFluxQuantumTorus + 1) * (NbrFluxQuantumSphere + 1)) + NbrBosons) <= 63)
#else
	    if ((((NbrFluxQuantumTorus + 1) * (NbrFluxQuantumSphere + 1)) + NbrBosons) <= 31)	    
#endif
	      {
		Space = new BosonOnT2xS2(NbrBosons, NbrFluxQuantumTorus, NbrFluxQuantumSphere, TotalKy, TotalLz);
	      }
	    else
	      {
		Space = 0;
		//	      Space = new BosonOnT2xS2Long(NbrBosons, NbrFluxQuantumTorus, NbrFluxQuantumSphere, TotalKy, TotalLz);
	      }
	  if (Space->GetHilbertSpaceDimension() > 0)
	    {
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnT2xS2GenericTwoBodyHamiltonian(Space, NbrBosons, NbrFluxQuantumTorus, Ratio, NbrFluxQuantumSphere,
									 NbrPseudoPotentials, PseudoPotentialMomentumTorus, PseudoPotentialMomentumSphere, 
									 PseudoPotentials, Architecture.GetArchitecture(), Memory);
	      
	      char* EigenvectorName = 0;
	      if (Manager.GetBoolean("eigenstate") == true)	
		{
		  char* TmpVectorExtension = new char [64];
		  sprintf (TmpVectorExtension, "_ky_%d_lz_%d", TotalKy, TotalLz);
		  EigenvectorName = ReplaceString(OutputName, ".dat", TmpVectorExtension);
		}
	      
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", TotalKy, TotalLz);	  
	      char* SubspaceLegend = new char[256];
	      sprintf (SubspaceLegend, "Ky Lz");
	      
	      GenericRealMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, ContentPrefix, SubspaceLegend, 0, OutputName, FirstRun, EigenvectorName);
	      MainTaskOperation TaskOperation (&Task);	  	 
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      delete Hamiltonian;
	      delete Space;
	      if (EigenvectorName != 0)
		{
		  delete[] EigenvectorName;
		}
	    }
	}
      if (FirstRun == true)
	FirstRun = false;
    }       
  return 0;
}
