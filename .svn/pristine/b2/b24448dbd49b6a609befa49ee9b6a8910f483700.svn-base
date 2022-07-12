#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslationsLong.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Hamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/FQHETorusPseudopotentialTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

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
  OptionManager Manager ("FQHETorusFermionsWithSpinAndTranslations" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "landau level index", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new SingleDoubleOption   ('d', "layerSeparation", 
					      "for bilayer simulations: layer separation in magnetic lengths", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the reduced Brillouin zone", false);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new SingleDoubleOption ('\n', "spinup-flux", "inserted flux for particles with spin up (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "spindown-flux", "inserted flux for particles with spin down (in 2pi / N_phi unit)", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "time-reversal", "Use model with time reversal symmetry", false);

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-binhamiltonian", "export the hamiltonian as a binary file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsWithSpinAndTranslations -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSpin = Manager.GetInteger("total-spin");
  int NbrIterLanczos = Manager.GetInteger("nbr-iter");
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  double XRatio = NbrFermions / 4.0;
  if (Manager.GetDouble("ratio") > 0)
    {
       XRatio = Manager.GetDouble("ratio");
    }
  double LayerSeparation = Manager.GetDouble("layerSeparation");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  int L = 0;
  double GroundStateEnergy = 0.0;
  bool HaveCoulomb = false;
  int LandauLevel = 0;
  char* InteractionName = 0;
  double** PseudoPotentials  = new double*[3];
  PseudoPotentials[0] = 0;
  PseudoPotentials[1] = 0;
  PseudoPotentials[2] = 0;
  double** OneBodyPseudoPotentials  = new double*[3];
  int* NbrPseudoPotentials  = new int[3];

  if (Manager.GetString("interaction-file") != NULL)
    {
      FQHETorusSU2GetPseudopotentials(Manager.GetString("interaction-file"), MaxMomentum, NbrPseudoPotentials, PseudoPotentials, OneBodyPseudoPotentials);
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
	  if ((InteractionDefinition["CoulombLandauLevel"] != NULL) && (PseudoPotentials[0] == 0) && (PseudoPotentials[1] == 0) && (PseudoPotentials[2] == 0))
	    {
	      InteractionName = new char[18];
	      if (LandauLevel >= 0)
		sprintf(InteractionName,"coulomb_l_%d",LandauLevel);
	      else
		sprintf(InteractionName,"graphene_l_%d",-LandauLevel);
	    }
	  else
	    {
	      cout << "Attention, using unnamed interaction! Please include a line 'Name = ...'" << endl;
	      InteractionName = new char[10];
	      sprintf(InteractionName,"unnamed");
	    }
	}
      else
	{
	  InteractionName = new char[strlen(InteractionDefinition["Name"])+1];
	  strcpy(InteractionName, InteractionDefinition["Name"]);
	}
    }
  else
    {
      OneBodyPseudoPotentials[2] = 0;
      LandauLevel = Manager.GetInteger("landau-level");
      InteractionName = new char[128];
      sprintf (InteractionName, "coulomb");
      HaveCoulomb = true;
    }


  char* OutputName = new char [512];
  if (OneBodyPseudoPotentials[2] == 0)
    {
        if (LayerSeparation == 0.0)
        {
            if (Manager.GetBoolean("time-reversal") == false)
                sprintf (OutputName, "fermions_torus_su2_%s_n_%d_2s_%d_sz_%d_ratio_%f.dat", InteractionName, NbrFermions, MaxMomentum, TotalSpin, XRatio);
            else
                sprintf (OutputName, "fermions_torus_timereversal_%s_n_%d_2s_%d_sz_%d_ratio_%f.dat", InteractionName, NbrFermions, MaxMomentum, TotalSpin, XRatio);
        }
        else
         {
            if (Manager.GetBoolean("time-reversal") == false)
                sprintf (OutputName, "fermions_torus_su2_%s_n_%d_2s_%d_d_%f_sz_%d_ratio_%f.dat", InteractionName, NbrFermions, MaxMomentum, LayerSeparation, TotalSpin, XRatio);
            else
                sprintf (OutputName, "fermions_torus_timereversal_%s_n_%d_2s_%d_d_%f_sz_%d_ratio_%f.dat", InteractionName, NbrFermions, MaxMomentum, LayerSeparation, TotalSpin, XRatio);
        }   
    }
  else
    {
      sprintf (OutputName, "fermions_torus_su2_%s_n_%d_2s_%d_ratio_%f.dat", InteractionName, NbrFermions, MaxMomentum, XRatio);
    }

  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);


  int MomentumModulo;
  if (Manager.GetBoolean("time-reversal") == false)
    {
      MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
    }
  else
    {
      MomentumModulo = FindGCD(abs(TotalSpin), MaxMomentum); 
    }

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

  bool FirstRun=true;
  for (int Pos = 0;Pos < NbrMomenta; ++Pos)
    {
      XMomentum = XMomenta[Pos];
      YMomentum = YMomenta[Pos];
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      ParticleOnTorusWithSpinAndMagneticTranslations* TotalSpace = 0;
      if (Manager.GetBoolean("time-reversal") == false)
	{
#ifdef __64_BITS__
	  if (MaxMomentum < 31)
#else
	    if (MaxMomentum < 15)
#endif
	      {
		if (OneBodyPseudoPotentials[2] == 0)
		  {
		    TotalSpace = new FermionOnTorusWithSpinAndMagneticTranslations (NbrFermions, TotalSpin, MaxMomentum, XMomentum, YMomentum);	
		  }
		else
		  {
		    TotalSpace = new FermionOnTorusWithSpinAndMagneticTranslations (NbrFermions, MaxMomentum, XMomentum, YMomentum);	
		  }
	      }
	    else
	      {
		if (OneBodyPseudoPotentials[2] == 0)
		  {
		    TotalSpace = new FermionOnTorusWithSpinAndMagneticTranslationsLong (NbrFermions, TotalSpin, MaxMomentum, XMomentum, YMomentum);	
		  }
		else
		  {
		    TotalSpace = new FermionOnTorusWithSpinAndMagneticTranslationsLong (NbrFermions, MaxMomentum, XMomentum, YMomentum);	
		  }
	      }
	}
      else
	{
	  TotalSpace = new FermionOnTorusWithSpinAndTimeReversalSymmetricMagneticTranslations (NbrFermions, TotalSpin, MaxMomentum, XMomentum, YMomentum);
	}
      cout << " Total Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
      cout << "momentum = (" << XMomentum << "," << YMomentum << ")" << endl;
      Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());	
      AbstractQHEHamiltonian* Hamiltonian = 0;
      if (HaveCoulomb == true)
	{
          if (Manager.GetBoolean("time-reversal") == false)
            Hamiltonian = new ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian (TotalSpace, NbrFermions, 
												MaxMomentum, XMomentum, XRatio, LayerSeparation, 
												Architecture.GetArchitecture(), Memory);
          else
            Hamiltonian = new ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian (TotalSpace, NbrFermions, 
														     MaxMomentum, XMomentum, XRatio, LayerSeparation, 
														     Architecture.GetArchitecture(), Memory);
	}
      else
	{
	  if (Manager.GetBoolean("time-reversal") == false)
	    Hamiltonian = new ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian (TotalSpace, NbrFermions, 
												MaxMomentum, XMomentum, XRatio,
												NbrPseudoPotentials[0], PseudoPotentials[0],
												NbrPseudoPotentials[1], PseudoPotentials[1],
												NbrPseudoPotentials[2], PseudoPotentials[2],
												Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
												Architecture.GetArchitecture(), Memory, 0, OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], OneBodyPseudoPotentials[2]);
	  else
	    Hamiltonian = new ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian (TotalSpace, NbrFermions, 
														     MaxMomentum, XMomentum, XRatio,
														     NbrPseudoPotentials[0], PseudoPotentials[0],
														     NbrPseudoPotentials[1], PseudoPotentials[1],
														     NbrPseudoPotentials[2], PseudoPotentials[2],
														     Manager.GetDouble("spinup-flux"), Manager.GetDouble("spindown-flux"),
														     Architecture.GetArchitecture(), Memory, 0, OneBodyPseudoPotentials[0], OneBodyPseudoPotentials[1], OneBodyPseudoPotentials[2]);
	}
	char* EigenvectorName = 0;
	if (Manager.GetBoolean("eigenstate"))	
	  {
	    EigenvectorName = new char [512];
	    char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	    sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	    delete [] TmpName;
	  }
	double Shift = -10.0;
	Hamiltonian->ShiftHamiltonian(Shift);      
	FQHEOnTorusMainTask Task (&Manager, TotalSpace, &Lanczos, Hamiltonian, YMomentum, Shift, OutputName, FirstRun, EigenvectorName, XMomentum);
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
	delete TotalSpace;
    }
  File.close();
  delete[] OutputName;
  return 0;
}
