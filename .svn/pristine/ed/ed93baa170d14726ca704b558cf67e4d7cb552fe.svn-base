#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnTorusWithSU3SpinAndMagneticTranslations.h"

#include "Hamiltonian/ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

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
  OptionManager Manager ("FQHETorusBosonsWithSU3SpinAndTranslations" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n1", "number of type 1 particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n2", "number of type 2 particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n3", "number of type 3 particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin1-flux", "inserted flux for particles with spin 1 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin2-flux", "inserted flux for particles with spin 2 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spin3-flux", "inserted flux for particles with spin 3 (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the reduced Brillouin zone", false);
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusBosonsWithSU3SpinAndTranslations -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  double XRatio = Manager.GetDouble("ratio");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if ((Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3")) == NbrBosons)
    {
      TotalTz = (Manager.GetInteger("nbr-n1") - Manager.GetInteger("nbr-n2"));
      TotalY = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") - (2 * Manager.GetInteger("nbr-n3")));
    }
  else
    {
      int NbrN1 = (2 * NbrBosons) + TotalY + (3 * TotalTz);
      int NbrN2 = (2 * NbrBosons) + TotalY - (3 * TotalTz);
      int NbrN3 = NbrBosons - TotalY;
      if ((NbrN1 < 0 ) || (NbrN2 < 0 ) || (NbrN3 < 0) || ((NbrN1 % 6) != 0) || ((NbrN2 % 6) != 0) || ((NbrN3 % 3) != 0))
	{
	  cout << "These values of Tz and Y cannot be achieved with this particle number!" << endl;
	  return -1;
	}
      NbrN1 /= 6;
      NbrN2 /= 6;
      NbrN3 /= 3;
    }

  char* InteractionName = 0;
  double** PseudoPotentials  = new double*[6];
  double** OneBodyPseudoPotentials  = new double*[6];
  int* NbrPseudoPotentials  = new int[6];
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHETorusSU3GetPseudopotentials(Manager.GetString("interaction-file"), MaxMomentum, NbrPseudoPotentials, PseudoPotentials, OneBodyPseudoPotentials) == false)
	return -1;
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  exit(-1);
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
  if ((Manager.GetDouble("spin1-flux") == 0.0) && (Manager.GetDouble("spin2-flux") == 0.0) && (Manager.GetDouble("spin3-flux") == 0.0))
    sprintf (OutputFileName, "bosons_torus_su3_%s_n_%d_2s_%d_tz_%d_y_%d_ratio_%f.dat", InteractionName, NbrBosons, MaxMomentum, TotalTz, TotalY, XRatio);
  else
    sprintf (OutputFileName, "bosons_torus_su3_%s_n_%d_2s_%d_tz_%d_y_%d_ratio_%f_flux1_%f_flux2_%f_flux3_%f.dat", InteractionName, 
	     NbrBosons, MaxMomentum, TotalTz, TotalY, XRatio, Manager.GetDouble("spin1-flux"), Manager.GetDouble("spin2-flux"), Manager.GetDouble("spin3-flux"));

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
  bool FirstRun = true;
  for (int Pos = 0;Pos < NbrMomenta; ++Pos)
    {
      XMomentum = XMomenta[Pos];
      YMomentum = YMomenta[Pos];
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      BosonOnTorusWithSU3SpinAndMagneticTranslations Space (NbrBosons, TotalTz, TotalY, MaxMomentum, XMomentum, YMomentum);	

      Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
      AbstractQHEHamiltonian* Hamiltonian = new ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian(&Space, NbrBosons, MaxMomentum, XMomentum, XRatio,
														    NbrPseudoPotentials, PseudoPotentials, OneBodyPseudoPotentials,
														    Manager.GetDouble("spin1-flux"), Manager.GetDouble("spin2-flux"), 
														    Manager.GetDouble("spin3-flux"), 
														    Architecture.GetArchitecture(), Memory);
      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if ( Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [512];
	  char *TmpName = RemoveExtensionFromFileName(OutputFileName, ".dat");
	  sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	  delete [] TmpName;
	}
      
      FQHEOnTorusMainTask Task (&Manager, &Space, &Lanczos, Hamiltonian, YMomentum, Shift, OutputFileName, FirstRun, EigenvectorName, XMomentum);
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
    }
  File.close();
  delete[] OutputFileName;
  return 0;
}
