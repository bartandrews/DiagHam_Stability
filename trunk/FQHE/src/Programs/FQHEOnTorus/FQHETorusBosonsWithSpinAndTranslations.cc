#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSpinAllSzAndMagneticTranslations.h"

#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonianWithPairing.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"


#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);
    
  // some running options and help
  OptionManager Manager ("FQHETorusBosonsWithSpinAndTranslations" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-spin", "total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);
  (*SystemGroup) += new SingleDoubleOption   ('\n', "angle", "angle between the two fundamental cycles of the torus in pi units (0 if rectangular)", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spinup-flux", "inserted flux for particles with spin up (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spindown-flux", "inserted flux for particles with spin down (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name as it should appear in output files (if the NAME option is not used in the interaction file)", "unknown");
  (*SystemGroup) += new SingleDoubleOption ('\n', "pairing", "amplitude of the pairing term (0 if none)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the reduced Brillouin zone", false);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "apply a temporary energy shift used during the diagonalization", -1.0);
  (*SystemGroup) += new SingleDoubleOption   ('\n', "layer-separation", "for the coulomb interaction, layer separation in magnetic lengths", 0.0);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
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
      cout << "see man page for option syntax or type FQHETorusBosonsWithSpinAndTranslations -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSpin = Manager.GetInteger("total-spin");
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  double XRatio = Manager.GetDouble("ratio");
  double Angle = Manager.GetDouble("angle");
  double LayerSeparation = Manager.GetDouble("layer-separation");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();

  if ((TotalSpin & 1) != (NbrBosons & 1))
    {
      TotalSpin &= ~1;
      TotalSpin |= (NbrBosons & 1);
    }

  char* InteractionName = 0;
  double** PseudoPotentials  = new double*[3];
  double** OneBodyPseudoPotentials  = new double*[3];
  int* NbrPseudoPotentials  = new int[3];
  bool HaveCoulomb = false;
  int LandauLevel = 0;

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "warning, no interaction file has been provided. Will use coulomb interaction in the LLL" << endl;
      LandauLevel = 0;
      InteractionName = new char[128];
      sprintf (InteractionName, "coulomb");
      HaveCoulomb = true;
    }
  else
    {
      if (FQHETorusSU2GetPseudopotentials(Manager.GetString("interaction-file"), MaxMomentum, NbrPseudoPotentials, PseudoPotentials, OneBodyPseudoPotentials) == false)
	{	  
	  return -1;
	}
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  exit(-1);
	}
      if (InteractionDefinition["CoulombLandauLevel"] != NULL)
	{
	  LandauLevel = atoi(InteractionDefinition["CoulombLandauLevel"]);
	  HaveCoulomb = true;
	}
      if (InteractionDefinition["Name"] == NULL)
	{
	  if (Manager.GetString("interaction-name") != 0)
	    {
	      InteractionName = new char[strlen(Manager.GetString("interaction-name")) + 1];
	      sprintf(InteractionName, "%s", Manager.GetString("interaction-name"));
	     }
	   else
	     {
	       cout << "Attention, using unnamed interaction! Please include a line 'Name = ...'" << endl;
	       InteractionName = new char[10];
	       sprintf(InteractionName, "unnamed");
	     }
	}
      else
	{
	  InteractionName = new char[strlen(InteractionDefinition["Name"]) + 1];
	  strcpy(InteractionName, InteractionDefinition["Name"]);
	}    
    }
    
  char* OutputFileName = new char [512];
  if ((OneBodyPseudoPotentials[2] == 0) && (Manager.GetDouble("pairing") == 0.0))
    {
      if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	{
	  if (Angle == 0.0)    
	    {
	      sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_sz_%d_ratio_%f.dat", InteractionName, NbrBosons, MaxMomentum, TotalSpin, XRatio);
	    }
	  else
	    {
	      sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_sz_%d_ratio_%f_angle_%.6f.dat", InteractionName, NbrBosons, MaxMomentum, TotalSpin, XRatio, Angle);
	    }
	}
      else
	{
	  if (Angle == 0.0)    
	    {
	      sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_sz_%d_ratio_%f_fluxup_%f_fluxdown_%f.dat", InteractionName, NbrBosons, 
		       MaxMomentum, TotalSpin, XRatio, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
	    }
	  else
	    {
	      sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_sz_%d_ratio_%f_angle_%.6f_fluxup_%f_fluxdown_%f.dat", InteractionName, NbrBosons, 
		       MaxMomentum, TotalSpin, XRatio, Angle, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
	    }
	}
    }
  else
    {
      if (Manager.GetDouble("pairing") == 0.0)
	{
	  if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	    {
	      if (Angle == 0.0)    
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_ratio_%f.dat", InteractionName, NbrBosons, MaxMomentum, XRatio);
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_ratio_%f_angle_%.6f.dat", InteractionName, NbrBosons, MaxMomentum, XRatio, Angle);
		}
	    }
	  else
	    {
	      if (Angle == 0.0)    
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_ratio_%f_fluxup_%f_fluxdown_%f.dat", InteractionName, NbrBosons, 
			   MaxMomentum, XRatio, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_n_%d_2s_%d_ratio_%f_angle_%.6f_fluxup_%f_fluxdown_%f.dat", InteractionName, NbrBosons, 
			   MaxMomentum, XRatio, Angle, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
		}
	    }
	}
      else
	{
	  if ((Manager.GetDouble("spinup-flux") == 0.0) && (Manager.GetDouble("spindown-flux") == 0.0))
	    {
	      if (Angle == 0.0)    
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_pairing_%.6f_n_%d_2s_%d_ratio_%f.dat", InteractionName, Manager.GetDouble("pairing"), NbrBosons, MaxMomentum, XRatio);
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_pairing_%.6f_n_%d_2s_%d_ratio_%f_angle_%.6f.dat", InteractionName, Manager.GetDouble("pairing"), NbrBosons, MaxMomentum, XRatio, Angle);
		}
	    }
	  else
	    {
	      if (Angle == 0.0)    
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_pairing_%.6f_n_%d_2s_%d_ratio_%f_fluxup_%f_fluxdown_%f.dat", InteractionName, Manager.GetDouble("pairing"),
			   NbrBosons, MaxMomentum, XRatio, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_su2_%s_pairing_%.6f_n_%d_2s_%d_ratio_%f_angle_%.6f_fluxup_%f_fluxdown_%f.dat", InteractionName, Manager.GetDouble("pairing"),
			   NbrBosons, MaxMomentum, XRatio, Angle, Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"));
		}
	    }
	}
    }
  ofstream File;
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);

  int MomentumModulo = FindGCD(NbrBosons, MaxMomentum);
  int YMaxMomentum = (MomentumModulo - 1);

  int NbrMomenta;
  int* XMomenta;
  int* YMomenta;
  bool GenerateMomenta = false;
  if ((XMomentum >= 0) && (YMomentum >= 0))
    {
      NbrMomenta = 1;
      XMomenta = new int[1];
      YMomenta = new int[1];
      XMomenta[0] = XMomentum;
      YMomenta[0] = YMomentum;
    }
  else
    {
      if (Manager.GetString("selected-points") == 0)
	{
	  if (Manager.GetBoolean("all-points"))
	    {
	      NbrMomenta = MaxMomentum * MomentumModulo;
	      XMomenta = new int[NbrMomenta];
	      YMomenta = new int[NbrMomenta];
	      NbrMomenta = 0;
	      for (int x = 0; x < MomentumModulo; ++x)
		{
		  for (int y = 0; y < MaxMomentum; ++y)
		    {
		      XMomenta[NbrMomenta] = x;
		      YMomenta[NbrMomenta] = y;
		      NbrMomenta++;
		    }
		}
	    }
	  else
	    {
	      if (XMomentum >= 0)
		{
		  NbrMomenta = MaxMomentum;
		  XMomenta = new int[NbrMomenta];
		  YMomenta = new int[NbrMomenta];
		  NbrMomenta = 0;
		  for (int y = 0; y < MaxMomentum; ++y)
		    {
		      XMomenta[NbrMomenta] = XMomentum;
		      YMomenta[NbrMomenta] = y;
		      NbrMomenta++;
		    }
		}
	      else
		{
		  if (YMomentum >= 0)
		    {
		      NbrMomenta = MomentumModulo;
		      XMomenta = new int[NbrMomenta];
		      YMomenta = new int[NbrMomenta];
		      NbrMomenta = 0;
		      for (int x = 0; x < MomentumModulo; ++x)
			{
			  XMomenta[NbrMomenta] = x;
			  YMomenta[NbrMomenta] = YMomentum;
			  NbrMomenta++;
			}
		    }
		  else
		    {
		      int TmpMax = MomentumModulo;
		      if ((XRatio == 1.0) && (Manager.GetBoolean("full-reducedbz") == false))
			{
			  int TmpMax = (MomentumModulo + 1) / 2;
			}
		      NbrMomenta = TmpMax * TmpMax;
		      XMomenta = new int[NbrMomenta];
		      YMomenta = new int[NbrMomenta];
		      NbrMomenta = 0;
		      for (int x = 0; x < TmpMax; ++x)
			{
			  for (int y = 0; y < TmpMax; ++y)
			    {
			      XMomenta[NbrMomenta] = x;
			      YMomenta[NbrMomenta] = y;
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
	  XMomenta = MomentumFile.GetAsIntegerArray(0);
	  YMomenta = MomentumFile.GetAsIntegerArray(1);
	}
    }
  bool FirstRun = true;
  for (int Pos = 0;Pos < NbrMomenta; ++Pos)
    {
      XMomentum = XMomenta[Pos];
      YMomentum = YMomenta[Pos];
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      BosonOnTorusWithSpinAndMagneticTranslations* Space = 0;
      if ((OneBodyPseudoPotentials[2] == 0) && (Manager.GetDouble("pairing") == 0.0))
	{
//  	  if (Manager.GetString("load-hilbert") != 0)
// 	    {
	      Space = new BosonOnTorusWithSpinAndMagneticTranslations (NbrBosons, TotalSpin, MaxMomentum, XMomentum, YMomentum);
// 	    }
// 	  else
// 	    {
// 	    }
//  	  if (Manager.GetString("save-hilbert") != 0)
//  	    {	     
//  	      ((BosonOnTorusWithSpinAndMagneticTranslations*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
//  	      return 0;
//  	    }
	}
      else
	{
	  Space = new BosonOnTorusWithSpinAllSzAndMagneticTranslations (NbrBosons, MaxMomentum, XMomentum, YMomentum);
	} 

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEHamiltonian* Hamiltonian = 0;
      if (HaveCoulomb == true)
	{
	  Hamiltonian = new ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian (Space, NbrBosons, 
											      MaxMomentum, XMomentum, XRatio, LayerSeparation, 
											      Architecture.GetArchitecture(), Memory);
	}
      else
	{      
	  if (Manager.GetDouble("pairing") == 0.0)
	    {
	      if (Angle == 0.0)
		{
		  Hamiltonian = new ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian(Space, NbrBosons, MaxMomentum, XMomentum, XRatio,
												     NbrPseudoPotentials[0], PseudoPotentials[0],
												     NbrPseudoPotentials[1], PseudoPotentials[1],
												     NbrPseudoPotentials[2], PseudoPotentials[2],
												     Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
												     Architecture.GetArchitecture(), Memory, 0,
												     OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], 
												     OneBodyPseudoPotentials[2]);
		}
	      else
		{
		  Hamiltonian = new ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonian(Space, NbrBosons, MaxMomentum, XMomentum, XRatio, Angle,
													    NbrPseudoPotentials[0], PseudoPotentials[0],
													    NbrPseudoPotentials[1], PseudoPotentials[1],
													    NbrPseudoPotentials[2], PseudoPotentials[2],
													    Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
													    Architecture.GetArchitecture(), Memory, 0,
													    OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], 
													    OneBodyPseudoPotentials[2]);
		}
	    }
	  else
	    {
	      if (Angle == 0.0)
		{
		  Hamiltonian = new ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonianWithPairing(Space, NbrBosons, MaxMomentum, XMomentum, XRatio,
														NbrPseudoPotentials[0], PseudoPotentials[0],
														NbrPseudoPotentials[1], PseudoPotentials[1],
														NbrPseudoPotentials[2], PseudoPotentials[2],
														Manager.GetDouble("pairing"), 
														Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
														Architecture.GetArchitecture(), Memory, 0,
														OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], 
														OneBodyPseudoPotentials[2]);
		}
	      else
		{
		  // 	      Hamiltonian = new ParticleOnTwistedTorusWithSpinAndMagneticTranslationsGenericHamiltonianWithPairing(Space, NbrBosons, MaxMomentum, XMomentum, XRatio, Angle,
		  // 														   NbrPseudoPotentials[0], PseudoPotentials[0],
		  // 														   NbrPseudoPotentials[1], PseudoPotentials[1],
		  // 														   NbrPseudoPotentials[2], PseudoPotentials[2],
		  // 														   Manager.GetDouble("pairing"), 
		  // 														   Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
		  // 														   Architecture.GetArchitecture(), Memory, 0,
		  // 														   OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], 
		  // 														   OneBodyPseudoPotentials[2]);
		  cout << "error, twisted torus and pairing is not yet supported" << endl;
		  return -1;
		}
	    }
	}
      double Shift = Manager.GetDouble("energy-shift");
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))	
	{
	  EigenvectorName = new char [512];
	  char *TmpName = RemoveExtensionFromFileName(OutputFileName, ".dat");
	  sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	  delete [] TmpName;
	}

      FQHEOnTorusMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, YMomentum, Shift, OutputFileName, FirstRun, EigenvectorName, XMomentum);
      Task.SetKxValue(XMomentum);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      
      if (FirstRun == true)
	FirstRun = false;
      
      delete Hamiltonian;
      delete Space;
    }
  File.close();
  delete[] OutputFileName;
  return 0;
}
