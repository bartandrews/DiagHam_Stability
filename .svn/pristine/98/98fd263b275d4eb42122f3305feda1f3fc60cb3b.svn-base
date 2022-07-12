#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/TwoDimensionalKagomeLatticeHamiltonian.h"
#include "Hamiltonian/TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian.h"

#include "HilbertSpace/AbstractSpinChain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainNewAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainNewSzSymmetryAnd2DTranslation.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinKagomeLattice" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  Manager += PrecalculationGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "cylinder", "use periodic boundary in the y direction only");
  (*SystemGroup) += new BooleanOption  ('\n', "force-negativesz", "compute negative Sz sectors");
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "do not use 2d translations");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 1);
  (*SystemGroup) += new  SingleDoubleOption ('j', "j-value", "coupling constant value for nearest neighbors", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('a', "anisotropy", "anisotropy between up and down triangles", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "break-c3", "amplitude ratio of C3 breaking term (1.0 if C3 is preserved)", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "easy-plane", "easy plane anisotropy", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new  BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx1", "first coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny1", "second coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx2", "first coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny2", "second coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "real-offset", "second coordinate in real space of the second spanning vector of the real space lattice (0 if lattice is untilted)", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "inversion-symmetry", "activate the inversion symmetry");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
//   (*SystemGroup) += new SingleIntegerOption  ('\n', "inversion-parity", "select the  inversion parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Sz symmetry)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Sz symmetry)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericPeriodicSpinChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSpins = NbrSitesX * NbrSitesY * 3;
  double JValue =  Manager.GetDouble("j-value");
  double JDownValue = JValue * Manager.GetDouble("anisotropy");
  double JEasyPlane = JValue * Manager.GetDouble("easy-plane");
  double JDownEasyPlane = JEasyPlane * Manager.GetDouble("anisotropy");
  double JC3Breaking = Manager.GetDouble("break-c3");
  bool NoTranslationFlag = Manager.GetBoolean("no-translation");
  
  int nx1 = Manager.GetInteger("nx1");
  int ny1 = Manager.GetInteger("ny1");
  int nx2 = Manager.GetInteger("nx2");
  int ny2 = Manager.GetInteger("ny2");
  
  int OffsetReal = Manager.GetInteger("real-offset");
  bool TiltedFlag = true;
  if ( ((nx1 == 0) && (ny1 == 0)) || ((nx2 == 0) && (ny2 == 0)) )
    TiltedFlag = false;
  else
    {
      if ((nx1*ny2 - nx2*ny1) != NbrSitesX * NbrSitesY)
	{
	  cout << "Boundary conditions define a lattice that has a number of sites different from NbrSiteX * NbrSiteY - should have (nx1*ny2 - nx2*ny1) = NbrSiteX * NbrSiteY " << endl;
	  return 0;
	}
      
      if ((((OffsetReal*ny2 + nx2) % NbrSitesX) != 0 || ((nx1 + OffsetReal*ny1) % NbrSitesX != 0)))
      {
	  cout << "Tilted lattice not properly defined. Should have ((offset*ny2 + nx2) % NbrSiteX) = 0 and ((nx1 + offset*ny1) % NbrSiteX = 0) to verify momentum conservation" << endl;
	  return 0;
      }
	
	
      cout << "Using tilted boundary conditions" << endl;
    }
  
  if ((NoTranslationFlag == false) && (Manager.GetBoolean("cylinder")))
  {
    cout << "Warning: 2d translations cannot be used for the cylinder geometry" << endl;
    NoTranslationFlag = true;
  }
  
    
  if (Manager.GetDouble("easy-plane") != 1.0)
    cout << "Warning: easy-plane anisotropy is not tested in this code" << endl;
  
  char* ParametersName = new char[256];
  if (Manager.GetDouble("anisotropy") == 1.0)
  {
    if (JC3Breaking == 1.0)
    {
      if (Manager.GetDouble("easy-plane") == 1.0)
	sprintf(ParametersName, "heisenberg");
      else
	sprintf(ParametersName, "jx_%.6f_jy_%.6f_jz_%.6f", JEasyPlane, JEasyPlane, JValue);
    }
    else
    {
      if (Manager.GetDouble("easy-plane") == 1.0)
	sprintf(ParametersName, "heisenberg_c3breaking_%.6f", JC3Breaking);
      else
	sprintf(ParametersName, "jx_%.6f_jy_%.6f_jz_%.6f_c3breaking_%.6f", JEasyPlane, JEasyPlane, JValue, JC3Breaking);
    }
  }
  else
  {
    if (JC3Breaking == 1.0)
    {
      if (Manager.GetDouble("easy-plane") == 1.0)
	sprintf(ParametersName, "anisotropy_jup_%.6f_jdown_%.6f", JValue, JDownValue);
      else
	sprintf(ParametersName, "jxup_%.6f_jyup_%.6f_jzup_%.6f_jxdown_%.6f_jydown_%.6f_jzdown_%.6f", JEasyPlane, JEasyPlane, JValue, JDownEasyPlane, JDownEasyPlane, JDownValue);
    }
    else
    {
      if (Manager.GetDouble("easy-plane") == 1.0)
	sprintf(ParametersName, "anisotropy_jup_%.6f_jdown_%.6f_c3breaking_%.6f", JValue, JDownValue, JC3Breaking);
      else
	sprintf(ParametersName, "jxup_%.6f_jyup_%.6f_jzup_%.6f_jxdown_%.6f_jydown_%.6f_jzdown_%.6f_c3breaking_%.6f", JEasyPlane, JEasyPlane, JValue, JDownEasyPlane, JDownEasyPlane, JDownValue, JC3Breaking);
    }
  }
  
     
  char* OutputFileName = new char [512];
  if (Manager.GetBoolean("cylinder"))
    sprintf (OutputFileName, "spin_1_2_kagome_cylinder_n_%d_x_%d_y_%d_%s", NbrSpins, NbrSitesX, NbrSitesY, ParametersName);
  else
  {
    if (NoTranslationFlag == false)
      if (TiltedFlag == false)
	sprintf (OutputFileName, "spin_1_2_kagome_n_%d_x_%d_y_%d_%s", NbrSpins, NbrSitesX, NbrSitesY, ParametersName);
      else
	sprintf (OutputFileName, "spin_1_2_kagome_n_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d_off_%d_%s", NbrSpins, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2, OffsetReal, ParametersName);
    else
      if (TiltedFlag == false)
	sprintf (OutputFileName, "spin_1_2_kagome_notranslation_n_%d_x_%d_y_%d_%s", NbrSpins, NbrSitesX, NbrSitesY, ParametersName);
      else
	sprintf (OutputFileName, "spin_1_2_kagome_notranslation_n_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d_off_%d_%s", NbrSpins, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2, OffsetReal, ParametersName);
  }
  
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
  
  
  AbstractSpinChain* Space = 0;
  
  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = NbrSpins & 1;
  if (Manager.GetBoolean("force-negativesz"))
    InitalSzValue = -MaxSzValue;
  if (Manager.GetInteger("initial-sz") > 1)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;
  int MinXMomentum = 0;
  int MaxXMomentum = NbrSitesX - 1;
  
  if (Manager.GetInteger("only-kx") >= 0)
    {
      MaxXMomentum = Manager.GetInteger("only-kx");
      MinXMomentum = MaxXMomentum;
    }
  int MinYMomentum = 0;
  int MaxYMomentum = NbrSitesY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
  {
    MaxYMomentum = Manager.GetInteger("only-ky");
    MinYMomentum = MaxYMomentum;
  }
  
  if (NoTranslationFlag)
  {
    MinXMomentum = 0;
    MaxXMomentum = 0;
    MinYMomentum = 0;	
    MaxYMomentum = 0;
  }
  
  bool SzSymmetryFlag = false;
  int MaxParity = 0;
  int MinParity = 0;
  if (Manager.GetBoolean("disable-szsymmetry") == false)
  {
    if ((InitalSzValue == 0) && (MaxSzValue == 0) && (NoTranslationFlag == false))
    {
      SzSymmetryFlag = true;
      MaxParity = 1;
      if (Manager.GetInteger("sz-parity") != 0)
      {
	MinParity = (1 - Manager.GetInteger("sz-parity")) / 2;
	MaxParity = MinParity;
      }
    }
    else
    {
      cout << "Work in the Sz = 0 sector to be able to use the Sz symmetry, and activate translation symmetry" << endl;
    }
  }
    
  bool InversionSymmetryFlag = false;
  int MaxInversionParity = 0;
  int MinInversionParity = 0;
  if (Manager.GetBoolean("inversion-symmetry"))
  {
    cout << "Error: inversion symmetry is not implemented" << endl;
//     if ((MinXMomentum == 0) && (MaxXMomentum == 0) && (MinYMomentum == 0) && (MaxYMomentum == 0))
//     {
//       InversionSymmetryFlag = true;
//       MaxInversionParity = 1;
//       if (Manager.GetInteger("inversion-parity") != 0)
//       {
// 	MinInversionParity = (1 - Manager.GetInteger("inversion-parity")) / 2;
// 	MaxInversionParity = MinInversionParity;
//       }
//     }
//     else
//     {
//       cout << "Work in the Sz = 0 sector to be able to use the Sz symmetry" << endl;
//     }
  }
    
  char* CommentLine = new char [512];
  if (NoTranslationFlag)
    sprintf (CommentLine, "spin 1/2 system with boundary conditions on the kagome lattice and %d sites in the x direction, %d sites in the y direction \n# Sz", NbrSitesX, NbrSitesY);
  else
  {
    if (SzSymmetryFlag == false)
      sprintf (CommentLine, "spin 1/2 system with boundary conditions on the kagome lattice and %d sites in the x direction, %d sites in the y direction and translations \n# Sz kx ky", NbrSitesX, NbrSitesY);
    else
      sprintf (CommentLine, "spin 1/2 system with boundary conditions on the kagome lattice and %d sites in the x direction, %d sites in the y direction, translations and Sz symmetry\n# Sz kx ky szsym", NbrSitesX, NbrSitesY);
  }
  
  
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
  {
    for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
    {
      for (int YMomentum = MinYMomentum; YMomentum <= MaxYMomentum; ++YMomentum)
	{
	  for (int parity = MinParity; parity <= MaxParity; ++parity)
	    {
// 	      for (int inversion = MinInversionParity; inversion <= MaxInversionParity; ++inversion)
// 		{
	      if (NoTranslationFlag)
	      {
		Space = new Spin1_2ChainNew (NbrSpins, InitalSzValue, 1000000);
		cout << "2Sz = " << InitalSzValue << endl; 
	      }
	      else
	      {
		if (SzSymmetryFlag == false)
		{
		  if (InversionSymmetryFlag == false)
		  {
		    Space = new Spin1_2ChainNewAnd2DTranslation(NbrSpins, InitalSzValue, XMomentum, NbrSitesX, YMomentum, NbrSitesY);
		    cout << "2Sz = " << InitalSzValue << " kx = " << XMomentum << " ky = " << YMomentum << endl; 
		  }
		  else
		  {
// 		    Space = new Spin1_2ChainNewInversionAnd2DTranslation(NbrSpins, InitalSzValue, inversion,  XMomentum, NbrSitesX, YMomentum, NbrSitesY);  
// 		    cout << "2Sz = " << InitalSzValue << " kx = " << XMomentum << " ky = " << YMomentum << " inversion Parity = " << (1 - 2*inversion) << endl; 
		  }
		}
		else
		{
		  if (Manager.GetString("load-hilbert") != 0)
		  {
		    Space = new Spin1_2ChainNewSzSymmetryAnd2DTranslation(Manager.GetString("load-hilbert"));
		    cout << "2Sz = " << InitalSzValue << " kx = " << XMomentum << " ky = " << YMomentum << " Sz Parity = " << (1 - 2*parity) << endl; 
		  }
		  else
		  {
		    if (Manager.GetBoolean("inversion-symmetry"))
		    {
		      
// 			cout << "2Sz = " << InitalSzValue << " kx = " << XMomentum << " ky = " << YMomentum << " Sz Parity = " << (1 - 2*parity) << " Inversion parity = " << (1 - 2*inversion) << endl; 
// 			Space = new Spin1_2ChainNewSzSymmetryInversionAnd2DTranslation(NbrSpins, InitalSzValue, parity, inversion, XMomentum, NbrSitesX, YMomentum, NbrSitesY);  
			
		      
		    }
		    else
		    {
		      Space = new Spin1_2ChainNewSzSymmetryAnd2DTranslation(NbrSpins, InitalSzValue, parity,  XMomentum, NbrSitesX, YMomentum, NbrSitesY);  
		      cout << "2Sz = " << InitalSzValue << " kx = " << XMomentum << " ky = " << YMomentum << " Sz Parity = " << (1 - 2*parity) << endl; 
		    }
		    if (Manager.GetString("save-hilbert") != 0)
		    {
		      Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		  }
		}
	    }
	    
	    if (Space->GetHilbertSpaceDimension() > 0)
	    {
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      char* TmpSzString = new char[64];
	      AbstractHamiltonian* Hamiltonian = 0;
	      if (NoTranslationFlag)
	      {
		Hamiltonian = new TwoDimensionalKagomeLatticeHamiltonian(Space, NbrSitesX, NbrSitesY, JValue, JDownValue, JEasyPlane, JDownEasyPlane, (!Manager.GetBoolean("cylinder")), OffsetReal);
		sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);
		sprintf (TmpSzString, "%d", InitalSzValue);
		
		GenericRealMainTask Task(&Manager, Space, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
		MainTaskOperation TaskOperation (&Task);
		TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      }
	      else
	      {
		if (JC3Breaking == 1.0)
		  Hamiltonian = new TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian(Space, XMomentum, NbrSitesX, YMomentum, NbrSitesY, JValue, JDownValue, JEasyPlane, JDownEasyPlane, OffsetReal);
		else
		  Hamiltonian = new TwoDimensionalKagomeLatticeAnd2DTranslationHamiltonian(Space, XMomentum, NbrSitesX, YMomentum, NbrSitesY, JValue, JDownValue, JC3Breaking, OffsetReal);
		if (SzSymmetryFlag == false)
		{
		  sprintf (TmpEigenstateString, "%s_sz_%d_kx_%d_ky_%d", OutputFileName, InitalSzValue, XMomentum, YMomentum);
		  sprintf (TmpSzString, "%d %d %d", InitalSzValue, XMomentum, YMomentum);
		}
		else
		  {
		    if (InversionSymmetryFlag == false)
		    {
		      sprintf (TmpEigenstateString, "%s_sz_%d_szsym_%d_kx_%d_ky_%d", OutputFileName, InitalSzValue, (1 - 2*parity), XMomentum, YMomentum);
		      sprintf (TmpSzString, "%d %d %d %d", InitalSzValue, XMomentum, YMomentum, (1 - 2*parity));
		    }
		    else
		    {
// 		      sprintf (TmpEigenstateString, "%s_sz_%d_kx_%d_ky_%d_szsym_%d_invparity_%d", OutputFileName, InitalSzValue, XMomentum, YMomentum, parity, inversion);
// 		      sprintf (TmpSzString, "%d %d %d %d %d", InitalSzValue, XMomentum, YMomentum, parity, inversion);
		    }
		  }
		Lanczos.SetComplexAlgorithms();
		GenericComplexMainTask Task(&Manager, Space, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
		MainTaskOperation TaskOperation (&Task);
		TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      }
	      
	      FirstRun = false;
	      delete Hamiltonian;
	      delete[] TmpSzString;
	      delete[] TmpEigenstateString;
	    }
	  }
	}
// 	}
    }
    }

  delete[] ParametersName;
  delete[] FullOutputFileName;
  delete[] OutputFileName;
  delete[] TmpEigenstateString;
  delete[] CommentLine;
  return 0;
}