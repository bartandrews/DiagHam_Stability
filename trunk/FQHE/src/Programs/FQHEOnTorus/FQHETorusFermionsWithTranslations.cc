#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsRealHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian.h"
#include "Hamiltonian/ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <cstdio>
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
  OptionManager Manager ("FQHETorusFermionsWithTranslations" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption   ('\n', "double-gate", "assume that the Coulomb interaction is screened by a double gate (0 if no gating, otherwise provide the gate distance)", 0.0); 
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "file describing the interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the full reduced Brillouin zone", false);
  (*SystemGroup) += new BooleanOption  ('\n', "add-wigner", "consider the energy contribution from the Wigner crystal", false);
  (*SystemGroup) += new BooleanOption  ('\n', "mass-anisotropy", "use a mass anisotropy for the system");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "anisotropy", "value of the anisotropy parameter alpha (i.e. q_g^2 = alpha q_x^2 + q_y^2 / alpha)", 1.0);
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new  BooleanOption ('\n', "enable-realhamiltonian", "use a real Hamiltonian at the inversion symmetric points");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption('\n', "energy-variance", "in addition to energy expectation, also evaluate energy variance sqrt[<H^2>-<H>^2]");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsWithTranslations -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  char *LoadPrecalculationFile=Manager.GetString("load-precalculation");
  int LandauLevel=0;
  int NbrPseudopotentials=0;
  double *Pseudopotentials=NULL;
  double HaveCoulomb=false;
  char *InteractionName=NULL;
  if (Manager.GetString("interaction-file")!=NULL)
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
	      InteractionName = new char[256];
	      if (LandauLevel >= 0)
		{
		  if (Manager.GetDouble("double-gate") != 0.0)
		    {
		      sprintf(InteractionName,"coulomb_l_%d_dgate_%.6f", LandauLevel, Manager.GetDouble("double-gate"));
		    }
		  else
		    {
		      sprintf(InteractionName,"coulomb_l_%d",LandauLevel);
		    }
		}
	      else
		{
		  sprintf(InteractionName,"graphene_l_%d",-LandauLevel);
		}
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
      InteractionName = new char[128];
      if (Manager.GetDouble("double-gate") != 0.0)
	{
	  if (LandauLevel >= 0)
	    {
	      sprintf(InteractionName,"coulomb_l_%d_dgate_%.6f", LandauLevel, Manager.GetDouble("double-gate"));
	    }
	  else
	    {
	      sprintf(InteractionName,"graphene_l_%d_dgate_%.6f", -LandauLevel, Manager.GetDouble("double-gate"));
	    }
	}
      else
	{
	  if (LandauLevel >= 0)
	    {
	      sprintf(InteractionName,"coulomb_l_%d", LandauLevel);
	    }
	  else
	    {
	      sprintf(InteractionName,"graphene_l_%d", -LandauLevel);
	    }
	}
      HaveCoulomb = true;
    }
  double XRatio = NbrFermions / 4.0;
  if (Manager.GetDouble("ratio") > 0)
    {
      XRatio = Manager.GetDouble("ratio");
    }

  double Angle = Manager.GetDouble("angle");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* SuffixOutputName = new char [256];
  if (Angle == 0.0)    
    {
      if (Manager.GetBoolean("mass-anisotropy"))
	{
	  sprintf (SuffixOutputName, "anisotropy_%f_n_%d_2s_%d_ratio_%.6f.dat", Manager.GetDouble("anisotropy"),
		   NbrFermions, MaxMomentum, XRatio);
	}
      else
	{
	  sprintf (SuffixOutputName, "n_%d_2s_%d_ratio_%.6f.dat", NbrFermions, MaxMomentum, XRatio);
	}
    }  
  else
    {
      sprintf (SuffixOutputName, "n_%d_2s_%d_ratio_%.6f_angle_%.6f.dat", NbrFermions, MaxMomentum, XRatio, Angle);
    }  
  char* OutputName = new char [512 + strlen(SuffixOutputName) + strlen(InteractionName)];
  sprintf (OutputName, "fermions_torus_%s_%s", InteractionName, SuffixOutputName);

  if (Manager.GetString("eigenvalue-file") != 0)
      strcpy(OutputName, Manager.GetString("eigenvalue-file"));

  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
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
  int *XMomenta;
  int *YMomenta;
  int *Multiplicities = NULL;
  int CenterX=0, CenterY=0;

  if (GenerateMomenta == false)
    {
      NbrMomenta=1;
      XMomenta = new int[1];
      YMomenta = new int[1];
      XMomenta[0]=XMomentum;
      YMomenta[0]=YMomentum;
    }
  else
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
		{
		  for (int YMomentum2 = YMomentum; YMomentum2 < YMaxMomentum; ++YMomentum2)
		    {
		      XMomenta[Pos] = XMomentum;
		      YMomenta[Pos] = YMomentum2;
		      ++Pos;
		    }
		}
	    }
	  else
	    {
	      if (NbrFermions&1)
		{
		  CenterX=0;
		  CenterY=0;
		}
	      else
		{
		  if ((NbrFermions/MomentumModulo*MaxMomentum/MomentumModulo)&1) // p*q odd?
		    {
		      CenterX=MomentumModulo/2;
		      CenterY=MomentumModulo/2;
		    }
		  else
		    {
		      CenterX=0;
		      CenterY=0;
		    }
		}
	      if ((XRatio == 1.0) && (Manager.GetBoolean("mass-anisotropy") == false))
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
  
  bool FirstRun = true;
  bool ForceRealFlag = false;
  for (int Pos = 0; Pos < NbrMomenta; ++Pos)
    {
      XMomentum=XMomenta[Pos];
      YMomentum=YMomenta[Pos];
      
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;
      if (Angle != 0.0)
          cout << " Angle = " << Angle << " *  Pi"<< endl;
      //	FermionOnTorus TotalSpace (NbrFermions, MaxMomentum, y);
      FermionOnTorusWithMagneticTranslations *TotalSpace = new FermionOnTorusWithMagneticTranslations(NbrFermions, MaxMomentum, XMomentum, YMomentum);
      //cout << " Total Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
      //cout << "momentum = (" << XMomentum << "," << YMomentum << ")" << endl;
      Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

      AbstractQHEHamiltonian* Hamiltonian;
      if (Angle == 0.0)
	{
	  if (Manager.GetBoolean("mass-anisotropy") == false)
	    {
	      if (Manager.GetDouble("double-gate") != 0.0)
		{
		  Hamiltonian = new ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
													 XRatio, Manager.GetDouble("double-gate"), HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, !Manager.GetBoolean("add-wigner"),
													 Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
		}
	      else
		{
		  if ((Manager.GetBoolean("enable-realhamiltonian") == true) &&
		      (((XMomentum % MomentumModulo) == 0) || (((MomentumModulo & 1) == 0) && ((XMomentum % MomentumModulo) == (MomentumModulo / 2)))) &&
		      (((YMomentum % MomentumModulo) == 0) || (((MomentumModulo & 1) == 0) && ((YMomentum % MomentumModulo) == (MomentumModulo / 2)))))
		    {
		      cout << "using real hamiltonian" << endl;
		      ForceRealFlag = true;
		      Lanczos.SetRealAlgorithms();
		      Hamiltonian = new ParticleOnTorusCoulombWithMagneticTranslationsRealHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
												      XRatio, HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, !Manager.GetBoolean("add-wigner"),
												      Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
		    }
		  else
		    {		      
		      cout << "using complex hamiltonian" << endl;
		      Hamiltonian = new ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
												  XRatio, HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, !Manager.GetBoolean("add-wigner"),
												  Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
		    }
		}
	    }
	  else
	    {
	      Hamiltonian = new ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
													XRatio, Manager.GetDouble("anisotropy"), HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, !Manager.GetBoolean("add-wigner"),
													Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
	    }
	}
      else
	{
          Hamiltonian = new ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian(TotalSpace, NbrFermions, MaxMomentum, XMomentum, 
											     XRatio, Angle * M_PI, HaveCoulomb, LandauLevel, NbrPseudopotentials, Pseudopotentials, !Manager.GetBoolean("add-wigner"),
											     Architecture.GetArchitecture(), Memory, LoadPrecalculationFile);
	}
      
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate"))	
	{
	  EigenvectorName = new char [512];
	  char *TmpName=RemoveExtensionFromFileName(OutputName, ".dat");
          if (Manager.GetString("eigenstate-file") == 0)
              sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
          else
              sprintf (EigenvectorName, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), XMomentum, YMomentum);
	  delete [] TmpName;
	}

    if ( (Manager.GetString("energy-expectation") != 0 ) || (Manager.GetBoolean("energy-variance") != 0 ) )
	{

	  char* StateFileName = Manager.GetString("energy-expectation");
	  if (IsFile(StateFileName) == false)
	    {
	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
	      return -1;           
	    }
	  ComplexVector InputState;
 	  if (InputState.ReadVector(StateFileName) == false)
	    {
	      cout << "error while reading " << StateFileName << endl;
	      return -1;
	    }

	  if (InputState.GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  ComplexVector TmpState(TotalSpace->GetHilbertSpaceDimension(), true);

	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &InputState, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	
	  Complex EnergyValue = InputState * TmpState;
          cout << "<Energy>= " << EnergyValue.Re << " " << EnergyValue.Im << endl;

      if (Manager.GetBoolean("energy-variance") != 0 )
       {
   	     ComplexVector TmpState2(TotalSpace->GetHilbertSpaceDimension(), true);
	     VectorHamiltonianMultiplyOperation Operation2 (Hamiltonian, &TmpState, &TmpState2);
	     Operation2.ApplyOperation(Architecture.GetArchitecture());
	     Complex varH = InputState * TmpState2 - EnergyValue * EnergyValue;
	     cout << "(varH)^2 = " << varH.Re << " " << varH.Im << endl;
       }   

	  return 0;
	}



      double Shift = -10.0;
      Hamiltonian->ShiftHamiltonian(Shift);      
      FQHEOnTorusMainTask Task (&Manager, TotalSpace, &Lanczos, Hamiltonian, YMomentum, Shift, OutputName, FirstRun, EigenvectorName, XMomentum, 0, ForceRealFlag);
      //      Task.SetKxValue(XMomentum);
      if (Multiplicities!=0)
	Task.SetMultiplicity(Multiplicities[Pos]);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
      Lanczos.SetComplexAlgorithms();
      ForceRealFlag = false;
      delete Hamiltonian;
      delete TotalSpace;
    }

  delete [] XMomenta;
  delete [] YMomenta;
  if (Multiplicities!=0)
    delete [] Multiplicities;
  return 0;
}
