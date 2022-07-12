#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/BosonOnTorusState.h"

#include "Hamiltonian/ParticleOnTorusDeltaWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstdio>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEBosonsTorusWithTranslation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

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
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleDoubleOption   ('R', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", 
					      -1);
  (*SystemGroup) += new SingleDoubleOption   ('\n', "angle", "angle between the two fundamental cycles of the torus in pi units (0 if rectangular)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('L', "landau-level", "Landau-level to be simulated", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "file describing the interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the full reduced Brillouin zone", false);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new BooleanOption  ('\n', "add-wigner", "consider the energy contribution from the Wigner crystal", false);
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleStringOption ('\n', "export-hilberttransformation", "export (in a binary file), the transformation matrix from the reduced Hilbert space to the eigentate basis");
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
      cout << "see man page for option syntax or type FQHEBosonsTorusWithTranslation -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrBosons = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  char* LoadPrecalculationFile=Manager.GetString("load-precalculation");
  int LandauLevel = 0;
  int NbrPseudopotentials = 0;
  double* Pseudopotentials = 0;
  double HaveCoulomb = false;
  char* InteractionName = 0;
  if (Manager.GetString("interaction-file") != 0)
    {
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
	  if ((InteractionDefinition["CoulombLandauLevel"] != NULL) && (InteractionDefinition["Pseudopotentials"] == NULL))
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
      InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', Pseudopotentials, NbrPseudopotentials);
    }
  else
    {
      LandauLevel = Manager.GetInteger("landau-level");
      InteractionName = new char[1];
      InteractionName[0]='\0';
      HaveCoulomb=true;
    }
  double XRatio = NbrBosons / 4.0;
  if (Manager.GetDouble("ratio") > 0)
    {
      XRatio = Manager.GetDouble("ratio");
    }

  double Angle = Manager.GetDouble("angle");

  long Memory = Manager.GetInteger("memory") << 20;

  char* OutputName = new char [512];
  char* SuffixOutputName = new char [256];
  if (Angle == 0.0)    
    {
      sprintf (SuffixOutputName, "n_%d_2s_%d_ratio_%.6f.dat", NbrBosons, MaxMomentum, XRatio);
    }  
  else
    {
      sprintf (SuffixOutputName, "n_%d_2s_%d_ratio_%.6f_angle_%.6f.dat", NbrBosons, MaxMomentum, XRatio, Angle);
    }  

  if (NbrPseudopotentials > 0)
    {
      sprintf (OutputName, "bosons_torus_%s_%s", InteractionName, SuffixOutputName);
    }
  else
    {
      if (LandauLevel > 0)
	sprintf (OutputName, "bosons_torus_coulomb_l_%d_%s", LandauLevel, SuffixOutputName);
      else
	if (LandauLevel < 0)
	  sprintf (OutputName, "bosons_torus_graphene_l_%d_%s", -LandauLevel, SuffixOutputName);
	else
	  sprintf (OutputName, "bosons_torus_coulomb_%s", SuffixOutputName);
    }

  int MomentumModulo = FindGCD(NbrBosons, MaxMomentum);
  int XMaxMomentum = (MomentumModulo - 1);
  bool GenerateMomenta = false;
  if ((XMomentum < 0)||(YMomentum < 0))
    GenerateMomenta = true;
  if (XMomentum < 0)
    XMomentum = 0;
  else
    XMaxMomentum = XMomentum;
  int YMaxMomentum = (MaxMomentum - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum;

  int NbrMomenta;
  int* XMomenta;
  int* YMomenta;
  int* Multiplicities = NULL;
  int CenterX=0, CenterY=0;

  if (GenerateMomenta == false)
    {
      NbrMomenta=1;
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
	      int Pos=0;
	      NbrMomenta = (XMaxMomentum-XMomentum+1)*(YMaxMomentum-YMomentum+1);
	      XMomenta = new int[NbrMomenta];
	      YMomenta = new int[NbrMomenta];
	      for (; XMomentum <= XMaxMomentum; ++XMomentum)
		for (int YMomentum2 = YMomentum; YMomentum2<= YMaxMomentum; ++YMomentum2)
		  {
		    XMomenta[Pos]=XMomentum;
		    YMomenta[Pos]=YMomentum2;
		    ++Pos;
		    cout << "Pos="<<Pos<<endl;
		  }
	    }
	  else // determine inequivalent states in BZ
	    {
	      if (Manager.GetBoolean("full-reducedbz"))
		{
		  int Pos=0;
		  XMaxMomentum = MomentumModulo;
		  YMaxMomentum = MomentumModulo;
		  NbrMomenta = MomentumModulo * MomentumModulo;
		  XMomenta = new int[NbrMomenta];
		  YMomenta = new int[NbrMomenta];
		  for (; XMomentum < XMaxMomentum; ++XMomentum)
		    for (int YMomentum2 = YMomentum; YMomentum2 < YMaxMomentum; ++YMomentum2)
		      {
			XMomenta[Pos] = XMomentum;
			YMomenta[Pos] = YMomentum2;
			++Pos;
		      }
		}
	      else
		{
		  CenterX=0;
		  CenterY=0;
		  if (XRatio == 1.0)
		    {
		      NbrMomenta=0;
		      for (int Kx = CenterX; Kx<=CenterX+MomentumModulo/2; ++Kx)
			for (int Ky= (Kx-CenterX)+CenterY; Ky<=CenterY+MomentumModulo/2; ++Ky)
			  {
			    ++NbrMomenta;
			  }
		      int Pos=0;
		      XMomenta = new int[NbrMomenta];
		      YMomenta = new int[NbrMomenta];
		      Multiplicities = new int[NbrMomenta];
		      for (int Kx = 0; Kx<=MomentumModulo/2; ++Kx)
			for (int Ky= Kx; Ky<=MomentumModulo/2; ++Ky, ++Pos)
			  {
			    XMomenta[Pos]=CenterX+Kx;
			    YMomenta[Pos]=CenterY+Ky;
			    if (Kx==0)
			      {
				if (Ky==0)
				  Multiplicities[Pos]=1; // BZ center
				else if (Ky==MomentumModulo/2)
				  Multiplicities[Pos]=2;
				else Multiplicities[Pos]=4;
			      }
			    else if (Kx==MomentumModulo/2)
			      {
				Multiplicities[Pos]=1; // BZ corner
			      }
			    else
			      {
				if (Ky==Kx) // diagonal ?
				  {
				    Multiplicities[Pos]=4; 
				  }
				else
				  {
				    if (Ky==MomentumModulo/2)
				      Multiplicities[Pos]=4;
				    else
				      Multiplicities[Pos]=8;
				  }
			      }
			  }
		    }
		  else // rectangular torus
		    {
		      NbrMomenta=(MomentumModulo/2+1)*(MomentumModulo/2+1);
		      int Pos=0;
		      XMomenta = new int[NbrMomenta];
		      YMomenta = new int[NbrMomenta];
		      Multiplicities = new int[NbrMomenta];
		      for (int Kx = 0; Kx<=MomentumModulo/2; ++Kx)
			for (int Ky= 0; Ky<=MomentumModulo/2; ++Ky, ++Pos)
			  {
			    XMomenta[Pos]=CenterX+Kx;
			    YMomenta[Pos]=CenterY+Ky;
			    if (Kx==0)
			      {
				if (Ky==0)
				  Multiplicities[Pos]=1; // BZ center
				else // on Gamma->X]
				  Multiplicities[Pos]=2;
			      }
			    else
			      {
				if (Ky==0)
				  Multiplicities[Pos]=2;
				else
				  {
				    if (Kx==MomentumModulo/2)
				      {
					if (Ky==MomentumModulo/2) // BZ corner?
					  Multiplicities[Pos]=1;
					else
					  Multiplicities[Pos]=2;
				      }
				    else
				      {
					if (Ky==MomentumModulo/2) // on edge?
					  Multiplicities[Pos]=2;
					else
					  Multiplicities[Pos]=4;
				      }
				  }
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
      if (Angle != 0.0)
          cout << " Angle = " << Angle << " *  Pi"<< endl;
      
      BosonOnTorusWithMagneticTranslationsShort* TotalSpace = new BosonOnTorusWithMagneticTranslationsShort(NbrBosons, MaxMomentum, XMomentum, YMomentum);
      Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

      AbstractQHEHamiltonian* Hamiltonian = 0;
      if (Angle == 0.0)
	{
	  Hamiltonian = new ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian (TotalSpace, NbrBosons, MaxMomentum, XMomentum, XRatio, 
										       HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, 
										       !Manager.GetBoolean("add-wigner"),
										       Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
	}
      else
	{
          Hamiltonian = new ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian(TotalSpace, NbrBosons, MaxMomentum, XMomentum, 
											     XRatio, Angle * M_PI, HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, 
											     !Manager.GetBoolean("add-wigner"),
											     Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
	}
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))	
	{
	  EigenvectorName = new char [512];
	  char *TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	  sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	  delete [] TmpName;
	}
      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);      
      FQHEOnTorusMainTask Task (&Manager, TotalSpace, &Lanczos, Hamiltonian, YMomentum, Shift, OutputName, FirstRun, EigenvectorName, XMomentum);
      Task.SetKxValue(XMomentum);
      if (Multiplicities != 0)
	Task.SetMultiplicity(Multiplicities[Pos]);
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

  delete[] XMomenta;
  delete[] YMomenta;
  if (Multiplicities!=0)
    delete[] Multiplicities;
  return 0;
}
