#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergAnd2DTranslationHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardKitaevHeisenbergModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "total number of sites (if negative, guess it from the geometry file)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new BooleanOption  ('\n', "stripe", "model geometry is a stripe");
  (*SystemGroup) += new BooleanOption  ('\n', "torus", "model geometry is a torus");
  (*SystemGroup) += new BooleanOption  ('\n', "open", "model geometry is an open parallelogram");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsites-x", "number of hexagons in the x direction (periodic direction for the stripe geometry) ", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsites-y", "number of hexagons in the y direction (aperiodic direction for the stripe geometry)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);

  (*SystemGroup) += new SingleStringOption  ('\n', "geometry-file", "name of an optional file that gives the position of the different bonds and their nature");
  (*SystemGroup) += new SingleStringOption  ('\n', "geometry-name", "name of the geometry used", "unknown");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site (Hubbard) potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "isotropic-t", "isotropic spin nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "anisotropic-t", "anisotropic nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "j1", "strength of the neareast neighbor Heisenberg interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "j2", "strength of the neareast neighbor anisotropic interaction", 1.0);

  (*SystemGroup) += new BooleanOption  ('\n', "xperiodic-boundary", "use periodic boundary conditions in the x direction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-momentum", "set the momentum along the x direction (negative if all momentum sectors have to be evaluated)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-xmomentum", "number of momentum values in the x direction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "2dperiodic-boundaries", "use periodic boundary conditions in the x and y directions");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "y-momentum", "set the momentum along the y direction (negative if all momentum sectors have to be evaluated)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-ymomentum", "number of momentum values in the y direction", 2);

  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by .#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption ('\n', "testhermitian-error", "error threshold when testing hermiticy (0 for machine accuracy)", 0.0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardKitaevHeisenbergModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSites = Manager.GetInteger("nbr-sites"); 
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  bool StripeFlag = Manager.GetBoolean("stripe"); 
  bool TorusFlag = Manager.GetBoolean("torus");
  bool OpenFlag = Manager.GetBoolean("open");
  int NbrSitesX = Manager.GetInteger("nbrsites-x"); 
  int NbrSitesY = Manager.GetInteger("nbrsites-y"); 
  int* SitesA = 0;
  int* SitesB = 0;
  int* BondTypes = 0;
  int NbrBonds = 0;
  
  
 
  if ((StripeFlag == false) && (TorusFlag == false) && (OpenFlag == false) && (Manager.GetString("geometry-file") == 0))
    {
      cout << "Error. A lattice geometry has to be specified" << endl; 
      return -1;
    }
    
  if (Manager.GetString("geometry-file") != 0)
    {
      MultiColumnASCIIFile GeometryFile;
      if (GeometryFile.Parse(Manager.GetString("geometry-file")) == false)
	{
	  GeometryFile.DumpErrors(cout);
	  return -1;
	}
      if (GeometryFile.GetNbrColumns() < 3)
	{
	  cout << "Error. " << GeometryFile.GetNbrColumns() << " has a wrong number of columns" << endl;
	  return -1;
	}
      
      SitesA = GeometryFile.GetAsIntegerArray(0);
      SitesB = GeometryFile.GetAsIntegerArray(1);
      BondTypes = GeometryFile.GetAsIntegerArray(2);
      
      int LargestSiteIndex = -1; 
      NbrBonds = GeometryFile.GetNbrLines();
      for (int i = 0; i < GeometryFile.GetNbrLines(); ++i)
	{
	  if (SitesA[i] > LargestSiteIndex)
	    LargestSiteIndex = SitesA[i];
	}
      
      for (int i = 0; i < GeometryFile.GetNbrLines(); ++i)
	{
	  if (SitesB[i] > LargestSiteIndex)
	    LargestSiteIndex = SitesB[i];
	}
      NbrSites = LargestSiteIndex + 1;
    }

  if ((Manager.GetBoolean("xperiodic-boundary") == true)  && 
      (((TorusFlag == false) && ((NbrSites % Manager.GetInteger("max-xmomentum")) != 0)) ||
       ((TorusFlag == true) && ((NbrSitesX % Manager.GetInteger("max-xmomentum")) != 0))))
    {
      cout << "Error. The number of sites is not compatible with the periodicity in the x direction" << endl; 
      return -1;
    }
  
  if (OpenFlag == true)
  {
   if ((NbrSitesX == 0) || (NbrSitesY == 0))
   {
     cout << "Error. The number of sites in directions x and y must be specified in the open geometry" << endl;
     return -1; 
   }
   
//    if (NbrSites != 2 * (NbrSitesY + NbrSitesX + NbrSitesX * NbrSitesY))
//    {
//      cout << "Error. The total number of sites is not compatible with NbrSitesX and NbrSitesY in the stripe geometry" << endl;
//      return -1;
//    }
   NbrSites = 2 * (NbrSitesY + NbrSitesX + NbrSitesX * NbrSitesY);
   
   int TmpIndex = 0;
   NbrBonds = 3 * NbrSitesX * NbrSitesY + 2 * (NbrSitesX + NbrSitesY)  - 1;
   SitesA = new int [NbrBonds];
   SitesB = new int [NbrBonds];
   BondTypes = new int [NbrBonds];
   
   for (int indexX = 0; indexX < NbrSitesX + 1; ++indexX)
     {
      for (int indexY = 0; indexY < NbrSitesY + 1; ++indexY)
	{
	  int indexA = 2 * ((indexX * (NbrSitesY + 1)) + indexY);
	  int indexB = 2 * ((indexX * (NbrSitesY + 1)) + indexY) + 1;
	  int indexA10 = 2 * ((((indexX + 1)) * (NbrSitesY + 1)) + indexY);
	  int shiftedIndexX = indexX - 1;
	  int flag = 0;
	  
	  
	  int indexB1m1 = 2 * ((shiftedIndexX * (NbrSitesY + 1)) + ((indexY + 1))) + 1;
	  if (indexB1m1 < 0)
	  {
	    indexB1m1 = indexB1m1 + NbrSites;
	    flag = 1;
	  }
	  
	  if ((indexA < NbrSites) && (indexB < NbrSites) && (indexX < NbrSitesX))
	  {
	    SitesA [TmpIndex] = indexA;
	    SitesB [TmpIndex] = indexB;
	    BondTypes [TmpIndex] = 0;
	    ++TmpIndex;
	  }
	  
	  if ((indexB < NbrSites) && (indexA10 < NbrSites) && (indexX < NbrSitesX) )
	  {
	    SitesA [TmpIndex] = indexA10;
	    SitesB [TmpIndex] = indexB;
	    BondTypes [TmpIndex] = 1;
	    ++TmpIndex;
	  }
	  
	  if ((indexA < NbrSites) && (indexB1m1 < NbrSites) && (indexY < NbrSitesY))
	  {
	    SitesA [TmpIndex] = indexA;
	    SitesB [TmpIndex] = indexB1m1;
	    BondTypes [TmpIndex] = 2;
	    ++TmpIndex;
	  }
	  if (flag == 1)
	  {
	   SitesA [TmpIndex] = indexB1m1 - NbrSites + 2*NbrSitesY  + 1;
	   SitesB [TmpIndex] = indexB1m1;
	   BondTypes [TmpIndex] = 1;
	   ++ TmpIndex;
	  }
	  	  
	}
    }
   
   
  }
      
  if (StripeFlag == true)
  {
   if (Manager.GetBoolean("xperiodic-boundary") == true) 
     NbrSitesX = Manager.GetInteger("max-xmomentum");
   if (NbrSitesX == 0 || NbrSitesY == 0)
   {
     cout << "Error. The number of sites in direction x and y must be specified in stripe geometry" << endl;
     return -1; 
   }
   
//    if (NbrSites != 2 * (NbrSitesY + 1) * NbrSitesX)
//    {
//      cout << "Error. The total number of sites is not compatible with NbrSitesX and NbrSitesY in the stripe geometry" << endl;
//      return -1;
//    }
   NbrSites = 2 * (NbrSitesY + 1) * NbrSitesX;
   
   int TmpIndex = 0;
   NbrBonds = NbrSitesX * (3 * NbrSitesY + 2);
   SitesA = new int [NbrBonds];
   SitesB = new int [NbrBonds];
   BondTypes = new int [NbrBonds];
   
   for (int indexX = 0; indexX < NbrSitesX; ++indexX)
     {
      for (int indexY = 0; indexY < NbrSitesY; ++indexY)
	{
	  int indexA = 2 * ((indexX * (NbrSitesY + 1)) + indexY);
	  int indexB = 2 * ((indexX * (NbrSitesY + 1)) + indexY) + 1;
	  int indexA10 = 2 * ((((indexX + 1) % NbrSitesX) * (NbrSitesY + 1)) + indexY);
	  int shiftedIndexX = indexX - 1;
	  if (shiftedIndexX < 0)
	    {
	      shiftedIndexX += NbrSitesX;
	    }
	  int indexB1m1 = 2 * ((shiftedIndexX * (NbrSitesY + 1)) + (indexY + 1)) + 1;
	  SitesA [3*TmpIndex] = indexA;
	  SitesB [3*TmpIndex] = indexB;
	  BondTypes [3*TmpIndex] = 0;
	  
	  SitesA [3*TmpIndex + 1] = indexB;
	  SitesB [3*TmpIndex + 1] = indexA10;
	  BondTypes [3*TmpIndex + 1] = 1;
	  
	  SitesA [3*TmpIndex  + 2] = indexA;
	  SitesB [3*TmpIndex + 2] = indexB1m1;
	  BondTypes [3*TmpIndex + 2] = 2;
	  	  
	  TmpIndex += 1;
	}
    }
   int indexY = NbrSitesY;
   int TmpIndex1 = 3*TmpIndex;
   for (int indexX = 0; indexX < NbrSitesX; ++indexX)
     {
	int indexA = 2 * ((indexX * (NbrSitesY + 1)) + indexY);
	int indexB = 2 * ((indexX * (NbrSitesY + 1)) + indexY) + 1;
	int indexA10 = 2 * ((((indexX + 1) % NbrSitesX) * (NbrSitesY + 1)) + indexY);
	
	SitesA [TmpIndex1] = indexA;
	SitesB [TmpIndex1] = indexB;
	BondTypes [TmpIndex1] = 0;
  
	SitesA [TmpIndex1 + 1] = indexB;
	SitesB [TmpIndex1 + 1] = indexA10;
	BondTypes [TmpIndex1 + 1] = 1;
	
	TmpIndex1 += 2;
     }
    
  }
    
  if (TorusFlag == true)
  {
    if (Manager.GetBoolean("2dperiodic-boundaries") == true)
    {
      NbrSitesX = Manager.GetInteger("max-xmomentum");
      NbrSitesY = Manager.GetInteger("max-ymomentum");
    }
    
    if ((NbrSitesX == 0) || (NbrSitesY == 0))
    {
      cout << "Error. The number of sites in directions x and y must be specified in torus geometry" << endl;
      return -1; 
    }
    
//     if (2*NbrSitesX * NbrSitesY != NbrSites)
//     {
//      cout << "Error. The number of sites is not compatible with the periodicity of the torus" << endl;
//      return -1;
//     }
    
    NbrSites = 2*NbrSitesX * NbrSitesY;
    
    int TmpIndex = 0;
    NbrBonds = 3*NbrSitesX*NbrSitesY;
    SitesA = new int [NbrBonds];
    SitesB = new int [NbrBonds];
    BondTypes = new int [NbrBonds];
    
    for (int indexX = 0; indexX < NbrSitesX; ++indexX)
     {
      for (int indexY = 0; indexY < NbrSitesY; ++indexY)
	{
	  int indexA = 2 * ((indexX * NbrSitesY) + indexY);
	  int indexB = 2 * ((indexX * NbrSitesY) + indexY) + 1;
	  int indexA10 = 2 * ((((indexX + 1) % NbrSitesX) * NbrSitesY) + indexY);
	  int shiftedIndexX = indexX - 1;
	  if (shiftedIndexX < 0)
	    {
	      shiftedIndexX += NbrSitesX;
	    }
	  int indexB1m1 = 2 * ((shiftedIndexX * NbrSitesY) + ((indexY + 1) % NbrSitesY)) + 1;
	  SitesA [3*TmpIndex] = indexA;
	  SitesB [3*TmpIndex] = indexB;
	  BondTypes [3*TmpIndex] = 0;
	  
	  SitesA [3*TmpIndex + 1] = indexB;
	  SitesB [3*TmpIndex + 1] = indexA10;
	  BondTypes [3*TmpIndex + 1] = 1;
	  
	  SitesA [3*TmpIndex  + 2] = indexA;
	  SitesB [3*TmpIndex + 2] = indexB1m1;
	  BondTypes [3*TmpIndex + 2] = 2;
	  	  
	  TmpIndex += 1;
	}
    }
  }
    
//     for (int i = 0; i < NbrBonds; ++i)
//       cout << SitesA[i] << " " << SitesB[i] << " " << BondTypes[i] << " " << i << endl;

//   if ((StripeFlag) && ((NbrSites % 4) != 2))
//   {
//    cout << "Error: number of sites should be of the form 4n + 2 for stripe geometry " << endl; 
//    return -1;
//   }
  

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      if (Manager.GetBoolean("xperiodic-boundary") == false)
	{
	  if (Manager.GetBoolean("2dperiodic-boundaries") == false)
	    {
	      if (SzSymmetryFlag == false)
		{
		  if (GutzwillerFlag == false)
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg");
		  else
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller");
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_szsym");
		  else
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_szsym");
		}
	    }
	  else
	    {
	      if (SzSymmetryFlag == false)
		{
		  if (GutzwillerFlag == false)
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_x_%d_y_%d", (int) Manager.GetInteger("max-xmomentum"), 
			     (int) Manager.GetInteger("max-ymomentum"));
		  else
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_x_%d_y_%d", (int) Manager.GetInteger("max-xmomentum"),
			     (int) Manager.GetInteger("max-ymomentum"));
		}
	      else
		{
		  if (GutzwillerFlag == false)
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_szsym_x_%d_y_%d", (int) Manager.GetInteger("max-xmomentum"), 
			     (int) Manager.GetInteger("max-ymomentum"));
		  else
		    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_szsym_x_%d_y_%d", (int) Manager.GetInteger("max-xmomentum"),
			     (int) Manager.GetInteger("max-ymomentum"));
		}
	    }
	}
      else
	{
	  if (SzSymmetryFlag == false)
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	      else
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_szsym_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	      else
		sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_szsym_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("xperiodic-boundary") == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg");
	  else
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller");
	}
      else
	{
	  if (SzSymmetryFlag == false)
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	      else
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	    }
	  else
	    {
	      if (GutzwillerFlag == false)
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_szsym_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	      else
		sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller_szsym_x_%d", (int) Manager.GetInteger("max-xmomentum"));
	    }
	}
    }
    
  

  char* FilePrefix = new char [256];
  if (StripeFlag == true)
    {
      sprintf (FilePrefix, "%s_stripe_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    }
  else
    {
      if (TorusFlag == true)
	{
	  sprintf (FilePrefix, "%s_torus_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
	}
      else
	{
	  if (OpenFlag == true)
	  {
	    sprintf (FilePrefix, "%s_open_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
	  }
	  else
	  {
	    sprintf (FilePrefix, "%s_%s_n_%d_ns_%d", StatisticPrefix, Manager.GetString("geometry-name"), NbrParticles, NbrSites);
	  }
	}
    }
  
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t_%g_tK_%g_j1_%g_j2_%g", Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("j1"), Manager.GetDouble("j2"));

  char* CommentLine = new char [256];
  if (Manager.GetBoolean("xperiodic-boundary") == false)
    {
      if (Manager.GetBoolean("2dperiodic-boundaries") == false)
	{
	  sprintf (CommentLine, "");
	}
      else
	{
	  if (SzSymmetryFlag == false)
	    {
	      sprintf (CommentLine, "kx ky");
	    }
	  else
	    {
	      sprintf (CommentLine, "kx ky szp");
	    }
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  sprintf (CommentLine, "kx");
	}
      else
	{
	  sprintf (CommentLine, "kx szp");
	}
    }
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetDouble("u-potential") == 0.0)
	sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);
      else
	sprintf(EigenvalueOutputFile, "%s_%s_u_%f.dat", FilePrefix, FileParameterString, Manager.GetDouble("u-potential"));
    }

  bool FirstRunFlag = true;
  if ((Manager.GetBoolean("xperiodic-boundary") == false) && (Manager.GetBoolean("2dperiodic-boundaries") == false))
    {
      int SzParitySector = -1;
      int MaxSzParitySector = 1;
      if (SzSymmetryFlag == false)
	{
	  SzParitySector = 1;
	}
      else
	{
	  if (Manager.GetInteger("sz-parity") != 0)
	    {
	      SzParitySector = Manager.GetInteger("sz-parity");
	      MaxSzParitySector = SzParitySector;
	    }
	}
      for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
	{
	  ParticleOnSphereWithSpin* Space = 0;
	  AbstractHamiltonian* Hamiltonian = 0;
	  if (SzSymmetryFlag == false)
	    {
	      if (GutzwillerFlag == false)
		Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
	      else
		Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
	    }
	  else
	    {
	      bool MinusParitySector = true;
	      if (SzParitySector == 1)
		MinusParitySector = false;
	      if (GutzwillerFlag == false)
		Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, NbrSites, MinusParitySector);
	      else
		Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites, MinusParitySector);
	    }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(Space, NbrParticles, NbrSites, NbrBonds, SitesA, 
										 SitesB, BondTypes, Manager.GetDouble("isotropic-t"), 
										 Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), Manager.GetDouble("j1"), 
										 Manager.GetDouble("j2"), Architecture.GetArchitecture(), Memory);
	  
// // 	  Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(Space, NbrParticles, NbrSites, GeometryFile, Manager.GetDouble("isotropic-t"), 
// // 										 Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), Manager.GetDouble("j1"), 
// // 										 Manager.GetDouble("j2"), Architecture.GetArchitecture(), Memory);
	  
	  char* ContentPrefix = new char[256];
	  if (SzSymmetryFlag == false)
	    {
	      sprintf (ContentPrefix, "");
	    }
	  else
	    {
	      sprintf (ContentPrefix, "%d", SzParitySector);
	    }

	  char* EigenstateOutputFile;
	  if (Manager.GetString("eigenstate-file") != 0)
	    {
	      EigenstateOutputFile = new char [512];
	      sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
	    }
	  else
	    {
	      char* TmpExtention = new char [512];
	      if (SzSymmetryFlag == false)
		{
		  sprintf (TmpExtention, "");
		}
	      else
		{
		  sprintf (TmpExtention, "_szp_%d", SzParitySector);
		}
	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	    }
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete Space;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}
      return 0;
    }
  if (Manager.GetBoolean("xperiodic-boundary") == true)
    {
      int XPeriodicity = NbrSites / Manager.GetInteger("max-xmomentum");
      int MinXMomentum = 0;
      int MaxXMomentum = (NbrSites / XPeriodicity) - 1;
      if (Manager.GetInteger("x-momentum") >= 0)
	{
	  MaxXMomentum = Manager.GetInteger("x-momentum");
	  MinXMomentum = MaxXMomentum;
	}
      for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
	{
	  int SzParitySector = -1;
	  int MaxSzParitySector = 1;
	  if (SzSymmetryFlag == false)
	    {
	      SzParitySector = 1;
	    }
	  else
	    {
	      if (Manager.GetInteger("sz-parity") != 0)
		{
		  SzParitySector = Manager.GetInteger("sz-parity");
		  MaxSzParitySector = SzParitySector;
		}
	    }
	  for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
	    {
	      ParticleOnSphereWithSpin* Space = 0;
	      AbstractHamiltonian* Hamiltonian = 0;
	      if (SzSymmetryFlag == false)
		{
		  if (GutzwillerFlag == false)
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity);
		  else
		    Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity);
		}
	      else
		{
		  bool MinusParitySector = true;
		  if (SzParitySector == 1)
		     MinusParitySector = false;
		  if (GutzwillerFlag == false)
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, XPeriodicity, MinusParitySector);
		  else
		    Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, XMomentum, 
														    XPeriodicity,MinusParitySector);
		}
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      
	      Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(Space, NbrParticles, NbrSites, NbrBonds, SitesA, 
												     SitesB, BondTypes, XMomentum, XPeriodicity, 
												     Manager.GetDouble("isotropic-t"), 
												     Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), 
												     Manager.GetDouble("j1"), Manager.GetDouble("j2"), 
												     Architecture.GetArchitecture(), Memory);
	      
	      char* ContentPrefix = new char[256];
	      if (SzSymmetryFlag == false)
		{
		  sprintf (ContentPrefix, "%d", XMomentum);
		}
	      else
		{
		  sprintf (ContentPrefix, "%d %d", XMomentum, SzParitySector);
		}
	      char* EigenstateOutputFile;
	      if (Manager.GetString("eigenstate-file") != 0)
		{
		  EigenstateOutputFile = new char [512];
		  sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
		}
	      else
		{
		  char* TmpExtention = new char [512];
		  if (SzSymmetryFlag == false)
		    {
		      sprintf (TmpExtention, "_kx_%d", XMomentum);
		    }
		  else
		    {
		      sprintf (TmpExtention, "_szp_%d_kx_%d", SzParitySector, XMomentum);
		    }
		  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
		}
	      
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      FirstRunFlag = false;
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      cout << "------------------------------------" << endl;
	      delete Hamiltonian;
	      delete Space;
	      delete[] EigenstateOutputFile;
	      delete[] ContentPrefix;
	    }
	}
    }
  else
    {
      int XPeriodicity = NbrSites / Manager.GetInteger("max-xmomentum");
      int MinXMomentum = 0;
      int MaxXMomentum = (NbrSites / XPeriodicity) - 1;
      if (Manager.GetInteger("x-momentum") >= 0)
	{
	  MaxXMomentum = Manager.GetInteger("x-momentum");
	  MinXMomentum = MaxXMomentum;
	}
      for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
	{
	  int MinYMomentum = 0;
	  int MaxYMomentum = Manager.GetInteger("max-ymomentum") - 1;
	  if (Manager.GetInteger("y-momentum") >= 0)
	    {
	      MaxYMomentum = Manager.GetInteger("y-momentum");
	      MinYMomentum = MaxYMomentum;
	    }
	  for (int YMomentum = MinYMomentum; YMomentum <= MaxYMomentum; ++YMomentum)
	    {
	      int SzParitySector = -1;
	      int MaxSzParitySector = 1;
	      if (SzSymmetryFlag == false)
		{
		  SzParitySector = 1;
		}
	      else
		{
		  if (Manager.GetInteger("sz-parity") != 0)
		    {
		      SzParitySector = Manager.GetInteger("sz-parity");
		      MaxSzParitySector = SzParitySector;
		    }
		}
	      for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
		{
		  ParticleOnSphereWithSpin* Space = 0;
		  AbstractHamiltonian* Hamiltonian = 0;
		  if (SzSymmetryFlag == false)
		    {
		      cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << endl;
		      if (GutzwillerFlag == false)
			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("max-xmomentum"),
										       YMomentum, Manager.GetInteger("max-ymomentum"));
		      else
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, Manager.GetInteger("max-xmomentum"),
													      YMomentum, Manager.GetInteger("max-ymomentum"));
		    }
		  else
		    {
		      bool MinusParitySector = true;
		      if (SzParitySector == 1)
			MinusParitySector = false;
		      cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  SzParity = " << SzParitySector<< endl;
		      if (GutzwillerFlag == false)
			Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, NbrSites, MinusParitySector, 
												 XMomentum, Manager.GetInteger("max-xmomentum"),
												 YMomentum, Manager.GetInteger("max-ymomentum"));
		      else
			Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, MinusParitySector, 
															XMomentum, Manager.GetInteger("max-xmomentum"),
															YMomentum, Manager.GetInteger("max-ymomentum"));
		    }
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
// 		      Space->PrintState(cout, i) << endl;
		    }
		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  
		  Hamiltonian = new ParticleOnLatticeWithSpinKitaevHeisenbergAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites, NbrBonds, SitesA, 											SitesB, BondTypes, XMomentum, Manager.GetInteger("max-xmomentum"),
													 YMomentum, Manager.GetInteger("max-ymomentum"),
													 Manager.GetDouble("isotropic-t"), 
													 Manager.GetDouble("anisotropic-t"), Manager.GetDouble("u-potential"), 
													 Manager.GetDouble("j1"), Manager.GetDouble("j2"), 
													 Architecture.GetArchitecture(), Memory);
		  
		  char* ContentPrefix = new char[256];
		  if (SzSymmetryFlag == false)
		    {
		      sprintf (ContentPrefix, "%d %d", XMomentum, YMomentum);
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d", XMomentum, YMomentum, SzParitySector);
		    }
		  char* EigenstateOutputFile;
		  if (Manager.GetString("eigenstate-file") != 0)
		    {
		      EigenstateOutputFile = new char [512];
		      sprintf (EigenstateOutputFile, "%s", Manager.GetString("eigenstate-file"));
		    }
		  else
		    {
		      char* TmpExtention = new char [512];
		      if (SzSymmetryFlag == false)
			{
			  sprintf (TmpExtention, "_kx_%d_ky_%d", XMomentum, YMomentum);
			}
		      else
			{
			  sprintf (TmpExtention, "_szp_%d_kx_%d_ky_%d", SzParitySector, XMomentum, YMomentum);
			}
		      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
		    }
		  
		  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
		  FirstRunFlag = false;
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  cout << "------------------------------------" << endl;
		  delete Hamiltonian;
		  delete Space;
		  delete[] EigenstateOutputFile;
		  delete[] ContentPrefix;
		}
	    }
	}
    }  
  return 0;
}
