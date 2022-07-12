#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceAndC4Symmetry.h"
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceNNExclusion.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian.h"

#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian.h"
 
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLatticeFullOBC.h"
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLatticeFullOBCC4Symmetry.h"
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry.h"

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


// compute the description of the density-density interaction for the unit cell at the origin
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude
// tightBindingModel = tight binding model
void FCICheckerboardLatticeModelEdgeStatesComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
								     int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
								     bool bosonFlag, double uPotential, double vPotential, Abstract2DTightBindingModel* tightBindingModel);


int main(int argc, char** argv)
{
  OptionManager Manager ("FCICheckerboardLatticeModelEdgeStates" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx1", "first coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny1", "second coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx2", "first coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny2", "second coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "offset", "second coordinate in momentum space of the second spanning vector of the reciprocal lattice (0 if lattice is untilted)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "real-offset", "second coordinate in real space of the second spanning vector of the real space lattice (0 if lattice is untilted)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the gutzwiller projected Hilbert space");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive nearest next neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "five-body", "use a five body interaction instead of a two body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "flatband-gap", "when using the flat band model with two bands, set the one-body gap between the two bands", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "real-space", "use the real space representation when considering the system with all bands");
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bands without the translations");
  (*SystemGroup) += new BooleanOption  ('\n', "no-c4symmetry", "disable the C4 symmetry (C4 symmetry is only available when using --real-space)");
  (*SystemGroup) += new BooleanOption  ('\n', "full-obc", "use full open boundary conditions");
  (*SystemGroup) += new SingleStringOption  ('\n', "confining-potential", "ascii file describing an additional confining potential when using full OBC");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "confining-scale", "multiply the confining potential by some global scale", 1.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
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
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCICheckerboardLatticeModelEdgeStates -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  int NbrSites = 2 * NbrSitesX * NbrSitesY;
  if (Manager.GetBoolean("full-obc") == true)
    {
      NbrSites = NbrSitesX * NbrSitesY;
    }
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int nx1 = Manager.GetInteger("nx1");
  int ny1 = Manager.GetInteger("ny1");
  int nx2 = Manager.GetInteger("nx2");
  int ny2 = Manager.GetInteger("ny2");
  int Offset = Manager.GetInteger("offset");
  int OffsetReal = Manager.GetInteger("real-offset");
  if (Manager.GetBoolean("real-space") == false)
    OffsetReal = 0;
  else
    Offset = 0;
  bool TiltedFlag = true;
  if ( ((nx1 == 0) && (ny1 == 0)) || ((nx2 == 0) && (ny2 == 0)) )
    {
      TiltedFlag = false;
    }
  else
    {
      if ((nx1 * ny2 - nx2 * ny1) != (NbrSitesX * NbrSitesY))
	{
	  cout << "Boundary conditions define a lattice that has a number of sites different from NbrSiteX * NbrSiteY - should have (nx1*ny2 - nx2*ny1) = NbrSiteX * NbrSiteY " << endl;
	  return 0;
	}
      
      if ((Manager.GetBoolean("real-space") == false) && (((Offset*ny2 - ny1) % NbrSitesX) != 0 || ((nx1 - Offset*nx2) % NbrSitesX != 0)))
	{
	  cout << "Tilted lattice not properly defined. Should have ((offset*ny2 - ny1) % NbrSiteX) = 0 and ((nx1 - offset*nx2) % NbrSiteX = 0) to verify momentum conservation" << endl;
	  return 0;
	}
      
      if ( (Manager.GetBoolean("real-space") == true) && (((OffsetReal*ny2 + nx2) % NbrSitesX) != 0 || ((nx1 + OffsetReal*ny1) % NbrSitesX != 0)))
      {
	  cout << "Tilted lattice not properly defined. Should have ((offset*ny2 + nx2) % NbrSiteX) = 0 and ((nx1 + offset*ny1) % NbrSiteX = 0) to verify momentum conservation" << endl;
	  return 0;
      }
	
	
      cout << "Using tilted boundary conditions" << endl;
    }
  if (TiltedFlag == true)
    {
      if (((Manager.GetBoolean("three-body") == true) || (Manager.GetBoolean("four-body") == true) || (Manager.GetBoolean("five-body"))) && (Manager.GetBoolean("real-space") == false))
	{
	  cout << "Tilted lattice is only supported for the projected code with two body interaction" << endl;
	  return 0;
	}
    }
  bool C4SymmetryFlag = false;
  if ((Manager.GetBoolean("no-c4symmetry") == false) && (NbrSitesX == NbrSitesY) && ((NbrSitesX & 1) == 0) && (Manager.GetDouble("mu-s") == 0.0))
    {
      C4SymmetryFlag = true;
    }

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* GeometryPrefix = new char [512];
  if (Manager.GetBoolean("full-obc") == true)
    {
      sprintf(GeometryPrefix, "checkerboardlattice_fullobc");
    }
  else
    {
      sprintf(GeometryPrefix, "checkerboardlattice");
    }

  char* FilePrefix = new char [512 + strlen(GeometryPrefix) + strlen(StatisticPrefix)];
  if (Manager.GetBoolean("single-band") == false)
    {
      if (Manager.GetBoolean("real-space") == false)
	{
	  if (TiltedFlag == false)
	    sprintf (FilePrefix, "%s_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	  else
	    sprintf (FilePrefix, "%s_%s_tilted_n_%d_ns_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2);
	}
      else
	{
         if (Manager.GetBoolean("no-translation") == false)
	   {
	     if (Manager.GetBoolean ("gutzwiller") == false)
	       {
		 if (TiltedFlag == false)
		   sprintf (FilePrefix, "%s_realspace_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
		 else
		   sprintf (FilePrefix, "%s_realspace_%s_tilted_n_%d_ns_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2);
	       }
	     else
	       {
		 if (TiltedFlag == false)
		   sprintf (FilePrefix, "%s_realspace_gutzwiller_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
		 else
		   sprintf (FilePrefix, "%s_realspace_gutzwiller_%s_tilted_n_%d_ns_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2);
	       }
	   }
	 else
	   {
	     if (Manager.GetBoolean ("gutzwiller") == false)
	       {
		 if (TiltedFlag == false)
		   sprintf (FilePrefix, "%s_realspace_notranslation_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
		 else
		   sprintf (FilePrefix, "%s_realspace_notranslation_%s_tilted_n_%d_ns_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2);
	       }
	     else
	       {
		 if (TiltedFlag == false)
		   sprintf (FilePrefix, "%s_realspace_gutzwiller_notranslation_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
		 else
		   sprintf (FilePrefix, "%s_realspace_gutzwiller_notranslation_%s_tilted_n_%d_ns_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2);
	       }
	   }
	}
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	{ 
	  if (TiltedFlag == false)
	    sprintf (FilePrefix, "%s_singleband_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	  else
	    sprintf (FilePrefix, "%s_singleband_%s_tilted_n_%d_ns_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2);
	}
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    {
	      sprintf (FilePrefix, "%s_singleband_threebody_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles,  NbrSites, NbrSitesX, NbrSitesY);
	    }
	  else
	    {
	      if (Manager.GetBoolean("four-body") == true)
		{
		  sprintf (FilePrefix, "%s_singleband_fourbody_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
		}
	      else
		{
		  sprintf (FilePrefix, "%s_singleband_fivebody_%s_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, GeometryPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
		}
	    }
	}
    }
  char* CommentLine = new char [256];

  if (Manager.GetBoolean("full-obc") == true)
    {
      if (C4SymmetryFlag == true)
	{
	  sprintf (CommentLine, "eigenvalues\n# c4 ");
	}
      else
	{
	  sprintf (CommentLine, "eigenvalues\n# ");
	}
    }
  else
    {
      sprintf (CommentLine, "eigenvalues\n# kx ky ");
    }

  char* FileParameterString = new char [256];
  if (Manager.GetString("confining-potential") == 0)
    {
      sprintf (FileParameterString, "t1_%f_t2_%f_mus_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("mu-s"));
    }
  else
    {
      sprintf (FileParameterString, "confining_%f_t1_%f_t2_%f_mus_%f", Manager.GetDouble("confining-scale"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("mu-s"));
    }
  char* EigenvalueOutputFile = new char [512 + strlen(FilePrefix) + strlen(FileParameterString)];
  if (Manager.GetString("eigenvalue-file") != 0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    }
  else
    {
      if (Manager.GetBoolean("full-obc") == true)
	{
	  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s.dat", FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString);
	}
      else
	{
	  if (Manager.GetBoolean("single-band") == false)
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f.dat", FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	  else
	    {
              if (Manager.GetBoolean("flat-band") == true)
                {
		  sprintf (EigenvalueOutputFile, "%s_v_%f_%s_gx_%f_gy_%f.dat", FilePrefix, Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		}
	      else
		{
		  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f.dat", FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		}
	    }
	}
    }

//  FermionOnSquareLatticeRealSpaceNNExclusion TmpHilbert (NbrParticles, NbrSitesX, NbrSitesY);

//   for (int i = 0; i < TmpHilbert.GetHilbertSpaceDimension(); ++i)
//     TmpHilbert.PrintState(cout, i) << endl;
//   cout << "Hilbert space dimension = " << TmpHilbert.GetHilbertSpaceDimension() << endl;

//   return 0;
  
  Abstract2DTightBindingModel* TightBindingModel = 0;
  AbstractTightBindingModel* TightBindingModelOBC = 0;
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;
      if (Manager.GetBoolean("full-obc") == true)
	{
	  if (Manager.GetString("confining-potential") == 0)
	    {
	      if (C4SymmetryFlag == false)
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBC (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
											  Manager.GetDouble("mu-s"), 0, ExportOneBody);
		}
	      else
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry (NbrSitesX, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
													   0, ExportOneBody);
// 		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBCC4Symmetry (NbrSitesX, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
// 												    Manager.GetDouble("mu-s"), 0, ExportOneBody);
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile ConfiningPotentialFile;
	      if (ConfiningPotentialFile.Parse(Manager.GetString("confining-potential")) == false)
		{
		  ConfiningPotentialFile.DumpErrors(cout);
		  return -1;
		}
	      int NbrNonZeroConfiningPotentials = ConfiningPotentialFile.GetNbrLines();
	      if ((NbrNonZeroConfiningPotentials == 0) || (ConfiningPotentialFile.GetNbrColumns() < 3))
		{
		}
	      int* ConfiningPotentialXCoordinates = ConfiningPotentialFile.GetAsIntegerArray(0);
	      int* ConfiningPotentialYCoordinates = ConfiningPotentialFile.GetAsIntegerArray(1);
	      double* ConfiningPotentials = ConfiningPotentialFile.GetAsDoubleArray(2);
	      if (Manager.GetDouble("confining-scale") != 1.0)
		{
		  for (int i = 0; i < NbrNonZeroConfiningPotentials; ++i)
		    {
		      ConfiningPotentials[i] *= Manager.GetDouble("confining-scale");
		    }
		}
	      if (C4SymmetryFlag == false)
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBC (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
											  Manager.GetDouble("mu-s"), ConfiningPotentialXCoordinates, ConfiningPotentialYCoordinates,
											  ConfiningPotentials, NbrNonZeroConfiningPotentials,
											  0, ExportOneBody);
		}
	      else
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBCC4Symmetry (NbrSitesX, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
												    Manager.GetDouble("mu-s"), ConfiningPotentialXCoordinates, ConfiningPotentialYCoordinates,
												    ConfiningPotentials, NbrNonZeroConfiningPotentials,
												    0, ExportOneBody);
		}

	    }
	  TightBindingModelOBC->WriteAsciiSpectrum(EigenvalueOutputFile);
	  if (ExportOneBody == true)
	    {
	      char* BandStructureOutputFile = new char [512];
	      if (Manager.GetString("export-onebodyname") != 0)
		strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	      else
		sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	      if (Manager.GetBoolean("export-onebody") == true)
		{
		  TightBindingModelOBC->WriteBandStructure(BandStructureOutputFile);
		}
	      else
		{
		  TightBindingModelOBC->WriteBandStructureASCII(BandStructureOutputFile);
		}
	      delete[] BandStructureOutputFile;
	    }	  
	}
      else
	{
	  if (TiltedFlag == false)
	    {
	      TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
									    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 
									    Architecture.GetArchitecture(), ExportOneBody, true);
	    }
	  else
	    {
	      TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2, Offset, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), 
									    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), OffsetReal, ExportOneBody, true);
	    }
	  TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
	  double BandSpread = TightBindingModel->ComputeBandSpread(0);
	  double DirectBandGap = TightBindingModel->ComputeDirectBandGap(0);
	  cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
	  if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	    {
	      cout << "Chern number = " << TightBindingModel->ComputeChernNumber(0) << endl;
	    }
	  if (ExportOneBody == true)
	    {
	      char* BandStructureOutputFile = new char [512];
	      if (Manager.GetString("export-onebodyname") != 0)
		strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	      else
		sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	      if (Manager.GetBoolean("export-onebody") == true)
		{
		  TightBindingModel->WriteBandStructure(BandStructureOutputFile);
		}
	      else
		{
		  TightBindingModel->WriteBandStructureASCII(BandStructureOutputFile);
		}
	      delete[] BandStructureOutputFile;
	    }	  
	}
      return 0;
    }


  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  if ((Manager.GetBoolean("no-translation") == true) || (Manager.GetBoolean("full-obc") == true))
    {  
      if (C4SymmetryFlag == false)
	{
	  MaxKx = 0;
	  MaxKy = 0;
	}
      else
	{
	  MinKx %= 4;
	  MaxKx = 3;
	  MaxKy = 0;
	}
    }

  
  if (Manager.GetString("import-onebody") == 0)
    {
      if (Manager.GetBoolean("full-obc") == true)
	{
	  if (Manager.GetString("confining-potential") == 0)
	    {
	      if (C4SymmetryFlag == false)
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBC (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
											  Manager.GetDouble("mu-s"), 0, true);
		}
	      else
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBCC4Symmetry (NbrSitesX, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
												    Manager.GetDouble("mu-s"), 0, true);
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile ConfiningPotentialFile;
	      if (ConfiningPotentialFile.Parse(Manager.GetString("confining-potential")) == false)
		{
		  ConfiningPotentialFile.DumpErrors(cout);
		  return -1;
		}
	      int NbrNonZeroConfiningPotentials = ConfiningPotentialFile.GetNbrLines();
	      if ((NbrNonZeroConfiningPotentials == 0) || (ConfiningPotentialFile.GetNbrColumns() < 3))
		{
		}
	      int* ConfiningPotentialXCoordinates = ConfiningPotentialFile.GetAsIntegerArray(0);
	      int* ConfiningPotentialYCoordinates = ConfiningPotentialFile.GetAsIntegerArray(1);
	      double* ConfiningPotentials = ConfiningPotentialFile.GetAsDoubleArray(2);
	      if (Manager.GetDouble("confining-scale") != 1.0)
		{
		  for (int i = 0; i < NbrNonZeroConfiningPotentials; ++i)
		    {
		      ConfiningPotentials[i] *= Manager.GetDouble("confining-scale");
		    }
		}
	      if (C4SymmetryFlag == false)
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBC (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
											  Manager.GetDouble("mu-s"), ConfiningPotentialXCoordinates, ConfiningPotentialYCoordinates,
											  ConfiningPotentials, NbrNonZeroConfiningPotentials,
											  0, true);
		}
	      else
		{
		  TightBindingModelOBC = new TightBindingModelCheckerboardLatticeFullOBCC4Symmetry (NbrSitesX, Manager.GetDouble("t1"), Manager.GetDouble("t2"),  
												    Manager.GetDouble("mu-s"), ConfiningPotentialXCoordinates, ConfiningPotentialYCoordinates,
												    ConfiningPotentials, NbrNonZeroConfiningPotentials,
												    0, true);
		}
	    }
	  char* BandStructureOutputFile = new char [1024];
	  sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
	  TightBindingModelOBC->WriteBandStructure(BandStructureOutputFile);
	}
      else
	{
	  if (TiltedFlag == false)
	    {
	      TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, !(Manager.GetBoolean("single-band")));
	    }
	  else
	    {
	      TightBindingModel = new TightBindingModelCheckerboardLattice (NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2, Offset,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), OffsetReal, true, !(Manager.GetBoolean("single-band")));
	    }
	  
	  char* BandStructureOutputFile = new char [1024];
	  sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
	  TightBindingModel->WriteBandStructure(BandStructureOutputFile);
	}
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }
//  TightBindingModel = new TightBindingModel2DAtomicLimitLattice(NbrSitesX, NbrSitesY, 2, 1, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());  

  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  if ((Manager.GetBoolean("no-translation") == true) || (Manager.GetBoolean("full-obc") == true))
	    {  
	      if (C4SymmetryFlag == true)
		{
		  cout << "c4=" << i << " : " << endl;		  
		}
	    }
	  else
	    {
	      cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
	    }
 	  if (Manager.GetBoolean("single-band") == false)
 	    {
	      ParticleOnSphereWithSpin* Space = 0;
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      int* NbrInteractingOrbitals;
	      int** InteractingOrbitalsOrbitalIndices;
	      int** InteractingOrbitalsSpatialIndices;
	      double** InteractingOrbitalsPotentials;
	      if (Manager.GetBoolean("full-obc") == false)
		{
		  FCICheckerboardLatticeModelEdgeStatesComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices, 
										  InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
										  Manager.GetBoolean("boson"), Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
										  TightBindingModel);
		}
	      if (Manager.GetBoolean("real-space") == false)
		{
		  if (Manager.GetBoolean("boson") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 31)
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			}
		      else
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			}
		    }
		  else
		    {
		      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

		  Hamiltonian = new ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
													0, 1, 
													NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
													InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
													TightBindingModel, Manager.GetBoolean("flat-band"), 
													Manager.GetDouble("flatband-gap"),
													Architecture.GetArchitecture(), 
													Memory);

// 		  Hamiltonian = new ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,TightBindingModel,
// 											    Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"),	     
// 											    Manager.GetBoolean("flat-band"), Manager.GetDouble("flatband-gap"),
// 											    Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  ParticleOnSphere* Space = 0;
		  RealSymmetricMatrix DensityDensityInteraction(NbrSites, true);
		  if (Manager.GetBoolean("full-obc") == false)
		    {
		      for (int x = 0; x < NbrSitesX; ++x)
			{
			  for (int y = 0; y < NbrSitesY; ++y)
			    {
			      for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModel->GetNbrBands(); ++OrbitalIndex)
				{
				  for (int k = 0; k < NbrInteractingOrbitals[OrbitalIndex]; ++k)
				    {
				      DensityDensityInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex), 
										   TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndices[OrbitalIndex][2 * k], 
																		  y + InteractingOrbitalsSpatialIndices[OrbitalIndex][(2 * k) + 1], 
																		  InteractingOrbitalsOrbitalIndices[OrbitalIndex][k]), 
										   InteractingOrbitalsPotentials[OrbitalIndex][k]);
				      
				    }
				}
			    }
			}
		    }
		  else
		    {
// 		      for (int x = 0; x < (NbrSitesX - 1); ++x)
// 			{
// 			  for (int y = 0; y < (NbrSitesY - 1); ++y)
// 			    {
// 			      for (int k = 0; k < NbrInteractingOrbitals[OrbitalIndex]; ++k)
// 				{
// 				  DensityDensityInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y), 
// 									       TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y), 
// 									       InteractingOrbitalsPotentials[OrbitalIndex][k]);
				      
// 				    }
// 				}
// 			    }
// 			}
		    }
		  if (Manager.GetBoolean("boson") == true)
		    {
		      if (Manager.GetBoolean ("gutzwiller") == false)
			{
			  if ((Manager.GetBoolean("no-translation") == true) || (Manager.GetBoolean("full-obc") == true))
			    {
			      Space = new BosonOnLatticeRealSpace(NbrParticles, NbrSites);
			    }
			  else
			    {
			      Space = new BosonOnLatticeRealSpaceAnd2DTranslation(NbrParticles, NbrSites, 
										  i, NbrSitesX, j, NbrSitesY);
			    }
			}
		      else
			{
			  if ((Manager.GetBoolean("no-translation") == true) || (Manager.GetBoolean("full-obc") == true))
			    {
			      Space = new BosonOnLatticeGutzwillerProjectionRealSpace(NbrParticles, NbrSites);
			    }
			  else
			    {
			      Space = new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(NbrParticles, 
												      NbrSites, 
												      i, NbrSitesX, j, NbrSitesY);
			    }
			}
		    }
		  else
		    {
		      if (Manager.GetBoolean ("gutzwiller") == false)
			{
			  if ((Manager.GetBoolean("no-translation") == true) || ((Manager.GetBoolean("full-obc") == true)))
			    {
			      if (C4SymmetryFlag == false)
				{
				  Space = new FermionOnLatticeRealSpace(NbrParticles, NbrSites);
				}
			      else
				{
				  Space = new FermionOnSquareLatticeRealSpaceAndC4Symmetry(NbrParticles, NbrSitesX, i);
				}
			    }
			  else
			    {
			      Space = new FermionOnLatticeRealSpaceAnd2DTranslation(NbrParticles, NbrSites, 
										    i, NbrSitesX, j, NbrSitesY);
			    }
			}
		      else
			{
			  if ((Manager.GetBoolean("no-translation") == true) || (Manager.GetBoolean("full-obc") == true))
			    {
			      if (Manager.GetBoolean("full-obc") == true)
				{
				  if (C4SymmetryFlag == false)
				    {
				      Space = new FermionOnSquareLatticeRealSpaceNNExclusion(NbrParticles, NbrSitesX, NbrSitesY);
				    }
				  else
				    {
				      char* ExclusionOutputFile = new char [1024];
				      sprintf (ExclusionOutputFile, "%s_%s_obc_exclusion.dat", FilePrefix, FileParameterString);
				      ofstream ExclusionFile;
				      ExclusionFile.open(ExclusionOutputFile);
				      int TmpIndex = NbrSitesX * NbrSitesY;
				      int* NbrExcludedSites = new int[TmpIndex];
				      int** ExcludedSites = new int*[TmpIndex];
				      for (int l = 0; l < TmpIndex; ++l)
					{
					  NbrExcludedSites[l] = 0;
					  ExcludedSites[l] = 0;
					}
				      for (int x = 0; x < (NbrSitesX - 1); ++x)
					{
					  for (int y = 0; y < (NbrSitesX - 1); ++y)
					    {
					      TmpIndex = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(x, y);
					      NbrExcludedSites[TmpIndex] = 2;					  
					      ExcludedSites[TmpIndex] = new int[2];					  
					      ExcludedSites[TmpIndex][0] = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y);
					      ExcludedSites[TmpIndex][1] = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1);
					      ExclusionFile << x << " " << y << " " << (x + 1) << " " << y << " " << TmpIndex << " " << ExcludedSites[TmpIndex][0] << endl;
					      ExclusionFile << x << " " << y << " " << x << " " << (y + 1) << " " << TmpIndex << " " << ExcludedSites[TmpIndex][1] << endl;
					    }
					}
				      for (int x = 0; x < (NbrSitesX - 1); ++x)
					{
					  TmpIndex = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(x, NbrSitesY - 1);
					  NbrExcludedSites[TmpIndex] = 1;					  
					  ExcludedSites[TmpIndex] = new int[1];					  
					  ExcludedSites[TmpIndex][0] = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, NbrSitesY - 1);
					  ExclusionFile << x << " " << (NbrSitesY - 1) << " " << (x + 1) << " " << (NbrSitesY - 1) << " " << TmpIndex << " " << ExcludedSites[TmpIndex][0] << endl;
					}
				      for (int y = 0; y < (NbrSitesY - 1); ++y)
					{
					  TmpIndex = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(NbrSitesX - 1, y);
					  NbrExcludedSites[TmpIndex] = 1;					  
					  ExcludedSites[TmpIndex] = new int[1];					  
					  ExcludedSites[TmpIndex][0] = ((TightBindingModelCheckerboardLatticeFullOBC*) TightBindingModelOBC)->GetRealSpaceTightBindingLinearizedIndexSafe(NbrSitesX - 1, y + 1);
					  ExclusionFile << (NbrSitesX - 1) << " " << y << " " << (NbrSitesX - 1) << " " << (y + 1) << " " << TmpIndex << " " << ExcludedSites[TmpIndex][0] << endl;
					}

				      ExclusionFile.close();
				      Space = new FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion(NbrParticles, NbrSitesX, i, ExcludedSites, NbrExcludedSites);
				    }
				}
			      else
				{
				  Space = new FermionOnSquareLatticeRealSpaceNNExclusion(NbrParticles, TightBindingModel->GetNbrBands() * NbrSitesX, NbrSitesY);
				}
			    }
			  else
			    {
			      char* ExclusionOutputFile = new char [1024];
			      sprintf (ExclusionOutputFile, "%s_%s_exclusion.dat", FilePrefix, FileParameterString);
			      ofstream ExclusionFile;
			      ExclusionFile.open(ExclusionOutputFile);
			      int TmpIndex = TightBindingModel->GetNbrBands();
			      int* NbrExcludedSites = new int[TmpIndex];
			      int** ExcludedSites = new int*[TmpIndex];
			      TmpIndex = 0;
			      NbrExcludedSites[TmpIndex] = 4;
			      ExcludedSites[TmpIndex] = new int[NbrExcludedSites[TmpIndex]];
			      ExcludedSites[TmpIndex][0] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(0, 0, 1);				      
			      ExcludedSites[TmpIndex][1] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(0, -1, 1);				      
			      ExcludedSites[TmpIndex][2] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(-1, 0, 1);				      
			      ExcludedSites[TmpIndex][3] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(-1, -1, 1);	
			      ExclusionFile << "0 0 1 " << TmpIndex << " " << ExcludedSites[TmpIndex][0] << endl;
			      ExclusionFile << "0 -1 1 " << TmpIndex << " " << ExcludedSites[TmpIndex][1] << endl;
			      ExclusionFile << "-1 0 1 " << TmpIndex << " " << ExcludedSites[TmpIndex][2] << endl;
			      ExclusionFile << "-1 -1 1 " << TmpIndex << " " << ExcludedSites[TmpIndex][3] << endl;
			      ++TmpIndex;
			      NbrExcludedSites[TmpIndex] = 4;
			      ExcludedSites[TmpIndex] = new int[NbrExcludedSites[TmpIndex]];
			      ExcludedSites[TmpIndex][0] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(0, 0, 0);				      
			      ExcludedSites[TmpIndex][1] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(1, 0, 0);				      
			      ExcludedSites[TmpIndex][2] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(0, 1, 0);				      
			      ExcludedSites[TmpIndex][3] = TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(1, 1, 0);				      
			      ExclusionFile << "0 0 0 " << TmpIndex << " " << ExcludedSites[TmpIndex][0] << endl;
			      ExclusionFile << "1 0 0 " << TmpIndex << " " << ExcludedSites[TmpIndex][1] << endl;
			      ExclusionFile << "0 1 0 " << TmpIndex << " " << ExcludedSites[TmpIndex][2] << endl;
			      ExclusionFile << "1 1 0 " << TmpIndex << " " << ExcludedSites[TmpIndex][3] << endl;
			      ExclusionFile.close();
			      ++TmpIndex;

			      Space = new FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion(NbrParticles, NbrSites, 
												 i, NbrSitesX, j, NbrSitesY, ExcludedSites, NbrExcludedSites);			      
			      for (int i = 0; i < TmpIndex; ++i)
				{				  
				  delete[] ExcludedSites[i];
				}
			      delete[] ExcludedSites;
			      delete[] NbrExcludedSites;
			    }
			}
		    }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  HermitianMatrix TightBindingMatrix;
		  if (Manager.GetBoolean("full-obc") == true)
		    {
		      TightBindingMatrix = TightBindingModelOBC->GetRealSpaceTightBindingHamiltonian();
		    }
		  else
		    {
		      TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
		    }
		  if ((Manager.GetBoolean("no-translation") == true) || (Manager.GetBoolean("full-obc") == true))
		    {
		      if (C4SymmetryFlag == false)
			{
			  Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, NbrSites, 
										   TightBindingMatrix, DensityDensityInteraction,
										   Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  Hamiltonian = new ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian (Space, NbrParticles, NbrSites, 
												   i, 4, TightBindingMatrix, DensityDensityInteraction,
												   Architecture.GetArchitecture(), Memory);
			}
		    }
		  else
		    {
		      Hamiltonian = new ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, NbrSites, 
											       i, NbrSitesX, j, NbrSitesY,
											       TightBindingMatrix, DensityDensityInteraction,
											       Architecture.GetArchitecture(), Memory);
		    }
		}
	      char* ContentPrefix = new char[256];
	      if (Manager.GetBoolean("full-obc") == true)
		{
		  if (Manager.GetBoolean("full-obc") == true)
		    {
		      if (C4SymmetryFlag == true)
			{
			  sprintf (ContentPrefix, "%d", i);
			}
		      else
			{
			  sprintf (ContentPrefix, "");
			}
		    }
		}
	      else
		{
		  sprintf (ContentPrefix, "%d %d", i, j);
		}
	      char* EigenstateOutputFileExtension = new char [512];
	      if ((Manager.GetBoolean("real-space") == false) || ((Manager.GetBoolean("no-translation") == false) && (Manager.GetBoolean("full-obc") == false)))
		{
		  sprintf (EigenstateOutputFileExtension, "_kx_%d_ky_%d", i, j);
		}
	      else
		{
		  if (C4SymmetryFlag == true)
		    {
		      sprintf (EigenstateOutputFileExtension, "_c4_%d", i);
		    }
		  else
		    {
		      sprintf (EigenstateOutputFileExtension, "");
		    }
		}
	      char* EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", EigenstateOutputFileExtension);
	      delete[] EigenstateOutputFileExtension;

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
 	  else
 	    {
	      ParticleOnSphere* Space = 0;
	      if (Manager.GetBoolean("boson") == false)
		{
		  if ((NbrSitesX * NbrSitesY) <= 63)
		    {
		      Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  else
		    {
		      Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		}
	      else
		{
		  Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
 	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
 	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
 	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
		{ 
//  test code for the generic density-density
// 		  int* NbrInteractingOrbitals;
// 		  int** InteractingOrbitalsOrbitalIndices;
// 		  int** InteractingOrbitalsSpatialIndices;
// 		  double** InteractingOrbitalsPotentials;
// 		  FCICheckerboardLatticeModelEdgeStatesComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices, 
// 									InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
// 									Manager.GetBoolean("boson"), Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"));
// 		  Hamiltonian = new ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, 0,
// 													   NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
// 													   InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
// 													   TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
											      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
											      TightBindingModel,Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{ 
		  if (Manager.GetBoolean("three-body") == true)
		    {
		      Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
													   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
													   TightBindingModel, 		     
													   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      if (Manager.GetBoolean("four-body") == true)
			{
			  Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
													      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
													      TightBindingModel, 		     
													      Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
													      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
													      TightBindingModel, 		     
													      Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}

		    }
		}
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
              if (Manager.GetString("eigenstate-file")!=0)
                  sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
              else
		{
		  char* TmpExtention = new char [512];
		  sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
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
  return 0;
}

// compute the description of the density-density interaction for the unit cell at the origin
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude
// tightBindingModel = tight binding model

void FCICheckerboardLatticeModelEdgeStatesComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
								     int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
								     bool bosonFlag, double uPotential, double vPotential, Abstract2DTightBindingModel* tightBindingModel)
{
  nbrInteractingOrbitals = new int[2];
  interactingOrbitalsOrbitalIndices = new int*[2];
  interactingOrbitalsSpatialIndices = new int*[2];
  interactingOrbitalsPotentials = new double*[2];
  int p;
  int q;
  if (bosonFlag == false)
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 3;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 2;
	  nbrInteractingOrbitals[1] += 2;
	}
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 2];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 2];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];

      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 1;
      tightBindingModel->GetRealSpaceIndex(0, 0, p, q);
      interactingOrbitalsSpatialIndices[0][2 * Index] = p;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
      interactingOrbitalsPotentials[0][Index] = uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[0][Index] = 0;
	  tightBindingModel->GetRealSpaceIndex(1, 0, p, q);
	  interactingOrbitalsSpatialIndices[0][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[0][Index] = vPotential;
	  ++Index;	  
	  interactingOrbitalsOrbitalIndices[0][Index] = 0;
	  tightBindingModel->GetRealSpaceIndex(0, 1, p, q);
	  interactingOrbitalsSpatialIndices[0][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[0][Index] = vPotential;
	  ++Index;	  
	}
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      tightBindingModel->GetRealSpaceIndex(1, 0, p, q);
      interactingOrbitalsSpatialIndices[1][2 * Index] = p;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      tightBindingModel->GetRealSpaceIndex(0, 1, p, q);
      interactingOrbitalsSpatialIndices[1][2 * Index] = p;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      tightBindingModel->GetRealSpaceIndex(1, 1, p, q);
      interactingOrbitalsSpatialIndices[1][2 * Index] = p;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[1][Index] = 1;
	  tightBindingModel->GetRealSpaceIndex(1, 0, p, q);
	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[1][Index] = vPotential;
	  ++Index;	  
	  interactingOrbitalsOrbitalIndices[1][Index] = 1;
	  tightBindingModel->GetRealSpaceIndex(0, 1, p, q);
	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[1][Index] = vPotential;
	  ++Index;	  
	}
    }
  else
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 1;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 1;
	  nbrInteractingOrbitals[1] += 3;
	}
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 2];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 2];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];
  
      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 0;
      tightBindingModel->GetRealSpaceIndex(0, 0, p, q);
      interactingOrbitalsSpatialIndices[0][2 * Index] = p;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
      interactingOrbitalsPotentials[0][Index] = 0.5 * uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[0][Index] = 1;
	  tightBindingModel->GetRealSpaceIndex(0, 0, p, q);
	  interactingOrbitalsSpatialIndices[0][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[0][Index] = vPotential;	  
	}
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 1;
      tightBindingModel->GetRealSpaceIndex(0, 0, p, q);
      interactingOrbitalsSpatialIndices[1][2 * Index] = p;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
      interactingOrbitalsPotentials[1][Index] = 0.5 * uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  tightBindingModel->GetRealSpaceIndex(1, 0, p, q);
	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  tightBindingModel->GetRealSpaceIndex(0, 1, p, q);
	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  tightBindingModel->GetRealSpaceIndex(1, 1, p, q);
	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	}
    }
}
