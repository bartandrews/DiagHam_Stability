#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU3SpinMomentumSpace.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"

#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandFourBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandFiveBodyHamiltonian.h"

#include "Hamiltonian/ParticleOnLatticeKagomeLatticeTwoBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeThreeBandHamiltonian.h"

#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"

#include "Hamiltonian/ExplicitHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/TightBindingModelAlternativeKagomeLattice.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"

#include "MathTools/IntegerAlgebraTools.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// find the tilted lattices with an aspect ratio close to unity 
//
// nbrSiteX = number of sites in the x direction (max(nbrSiteX, nbrSiteY) sets the maximal value for |nx1|, |ny1|, |nx2| and |ny2|)
// nbrSiteY = number of sites in the y direction (max(nbrSiteX, nbrSiteY) sets the maximal value for |nx1|, |ny1|, |nx2| and |ny2|)
void FCIKagomeLatticeModelFindOptimalTiltedLattices(int nbrSiteX, int nbrSiteY);

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
void FCIKagomeLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential, Abstract2DTightBindingModel* tightBindingModel);


// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double t1, double t2, double lambda1, double lambda2);


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIKagomeLatticeModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-kx", "maximum of range of x-momentum sectors (if set, --only-kx value sets minimum; negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-ky", "maximum of range of y-momentum sectors (if set, --only-ky value sets minimum; negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new  BooleanOption  ('\n', "redundant-kx", "Calculate all kx subspaces", false);
  (*SystemGroup) += new  BooleanOption  ('\n', "redundant-ky", "Calculate all ky subspaces", false);

  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");

  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength (or on-site density-density potential for bosons)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive next nearest neighbor potential strength (or nearest neighbor density-density potential for bosons)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive next next nearest neighbor potential strength (or next nearest neighbor density-density potential for bosons)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use hard core bosons (bosonic hilbert space will be restricted to single occupations) -- real space option needs to activated");
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "five-body", "use a five body interaction instead of a two body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "real part of the nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "real part of the next nearest neighbor hopping amplitude", -0.3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l1", "imaginary part of the nearest neighbor hopping amplitude", 0.28);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l2", "imaginary part of the next nearest neighbor hopping amplitude", 0.2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx1", "first coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny1", "second coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx2", "first coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny2", "second coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "offset", "second coordinate in momentum space of the second spanning vector of the reciprocal lattice (0 if lattice is untilted or if Ny = 1)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "find-optimaltilt", "find tilted lattices with an aspect ratio close to unity", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band that has to be partially filled, should be 0 (lower band), 1 or 2 (upper band)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "two-bands", "use the two lowest energy bands", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-bands", "use the full three band model", 0);
  (*SystemGroup) += new BooleanOption ('\n', "project-threebands", "project the hamiltonian from the thre band model to the single band model");
  (*SystemGroup) += new BooleanOption  ('\n', "real-space", "use the real space representation when considering the system with all bands");
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bandswithout the translations");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice chemical potential on A site", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-highsymmetryspectrum", "only compute the one body spectrum, restricting to lines connecting the high symmetry points");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-berrycurvature", "compute the Berry curvature number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
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
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIKagomeLatticeModel -h" << endl;
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
  int Nx1 = Manager.GetInteger("nx1");
  int Ny1 = Manager.GetInteger("ny1");
  int Nx2 = Manager.GetInteger("nx2");
  int Ny2 = Manager.GetInteger("ny2");
  int Offset = Manager.GetInteger("offset");
  bool TiltedFlag = true;
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if (Manager.GetBoolean("find-optimaltilt") == true)
    {
      FCIKagomeLatticeModelFindOptimalTiltedLattices(NbrSitesX, NbrSitesY);
      return 0;
    }
  if ( ((Nx1 == 0) && (Ny1 == 0)) || ((Nx2 == 0) && (Ny2 == 0)))
    {
      TiltedFlag = false;
    }
  else
    {
      if (Manager.GetBoolean("real-space") == true)
	{
	  cout << "error, tilted lattices are not implement in real space" << endl;
	  return 0;
	}
      if (((Nx1 * Ny2) - (Nx2 * Ny1)) != NbrSitesX * NbrSitesY)
	{
	  cout << "Boundary conditions define a lattice that has a number of sites different from NbrSitesX * NbrSitesY - should have (nx1*ny2 - nx2*ny1) = nbr-sitex * nbr-sitey" << endl;
	  return 0;
	}
      if ((((Offset * Ny2) - Ny1) % NbrSitesX) != 0 || (((Nx1 - (Offset * Nx2)) % NbrSitesX) != 0))
	{
	  cout << "Tilted lattice not properly defined. Should have ((offset * ny2 - ny1) % nbr-sitex) = 0 and ((nx1 - offset*nx2) % nbr-sitex = 0) to verify momentum conservation" << endl;
	  return 0;
	}
      else
	{
	  cout << "Using tilted boundary conditions" << endl;
	}
    }


  if (Manager.GetBoolean("three-bands"))
    cout << "Warning: three bands code is only tested for bosons with on-site interaction" << endl;

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* LatticeFileName = new char [256];
  if (TiltedFlag == true)
    {
      sprintf (LatticeFileName, "kagomelatticetilted_n_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d_offset_%d", NbrParticles, NbrSitesX, NbrSitesY, Nx1, Ny1, Nx2, Ny2, Offset);
    }
  else
    {
      sprintf (LatticeFileName, "kagomelattice_n_%d_x_%d_y_%d", NbrParticles, NbrSitesX, NbrSitesY);
    }
  char* FilePrefix = new char [256 + strlen(LatticeFileName)];

  if ((Manager.GetBoolean("three-bands") == false) && (Manager.GetBoolean("two-bands") == false) && (Manager.GetBoolean("real-space") == false))
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	{ 
	  sprintf (FilePrefix, "%s_singleband_%s", StatisticPrefix, LatticeFileName);
	}
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    {
	      if (Manager.GetBoolean("flat-band") == true)
		sprintf (FilePrefix, "%s_singleband_flat_threebody_%s", StatisticPrefix, LatticeFileName);
	      else
		sprintf (FilePrefix, "%s_singleband_threebody_%s", StatisticPrefix, LatticeFileName);
	    }
	  else
	    {
	      if (Manager.GetBoolean("four-body") == true)
		{
		  if (Manager.GetBoolean("flat-band") == true)
		    sprintf (FilePrefix, "%s_singleband_flat_fourbody_%s", StatisticPrefix, LatticeFileName);
		  else
		    sprintf (FilePrefix, "%s_singleband_fourbody_%s", StatisticPrefix, LatticeFileName);
		}
	      else
		{
		  if (Manager.GetBoolean("flat-band") == true)
		    sprintf (FilePrefix, "%s_singleband_flat_fivebody_%s", StatisticPrefix, LatticeFileName);
		  else
		    sprintf (FilePrefix, "%s_singleband_fivebody_%s", StatisticPrefix, LatticeFileName);
		}
	    }
	}
    }
  else
    {
      if ((Manager.GetBoolean("three-bands") == false) && (Manager.GetBoolean("real-space") == false))
	{
	}
      else
	{
	  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	    { 
	      if (Manager.GetBoolean("real-space") == false)
		{
		  if (Manager.GetBoolean("flat-band") == false)
		    sprintf (FilePrefix, "%s_threeband_%s", StatisticPrefix, LatticeFileName);
		  else
		    sprintf (FilePrefix, "%s_threeband_flatband_%s", StatisticPrefix, LatticeFileName);
		}
	      else
		{
		  if (Manager.GetBoolean("no-translation") == false)
		    {
		      if (GutzwillerFlag == false)
			sprintf (FilePrefix, "%s_realspace_%s", StatisticPrefix, LatticeFileName);
		      else
			sprintf (FilePrefix, "%s_realspace_gutzwiller_%s", StatisticPrefix, LatticeFileName);
		    }
		  else
		    {
		      if (GutzwillerFlag == false)
			sprintf (FilePrefix, "%s_realspace_notranslation_%s", StatisticPrefix, LatticeFileName);
		      else
			sprintf (FilePrefix, "%s_realspace_gutzwiller_notranslation_%s", StatisticPrefix, LatticeFileName);
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("three-body") == true)
		{
		  if (Manager.GetBoolean("flat-band") == true)
		    sprintf (FilePrefix, "%s_threeband_flat_threebody_%s", StatisticPrefix, LatticeFileName);
		  else
		    sprintf (FilePrefix, "%s_threeband_threebody_%s", StatisticPrefix, LatticeFileName);
		}
	      else
		{
		  if (Manager.GetBoolean("four-body") == true)
		    {
		      if (Manager.GetBoolean("flat-band") == true)
			sprintf (FilePrefix, "%s_threeband_flat_fourbody_%s", StatisticPrefix, LatticeFileName);
		      else
			sprintf (FilePrefix, "%s_threeband_fourbody_%s", StatisticPrefix, LatticeFileName);
		    }
		  else
		    {
		      if (Manager.GetBoolean("flat-band") == true)
			sprintf (FilePrefix, "%s_threeband_flat_fivebody_%s", StatisticPrefix, LatticeFileName);
		      else
			sprintf (FilePrefix, "%s_threeband_fivebody_%s", StatisticPrefix, LatticeFileName);
		    }
		}
	    }
	}
    }

  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t1_%g_t2_%g_l1_%g_l2_%g", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"));

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (((Manager.GetBoolean("flat-band") == true) && (Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false)) || (GutzwillerFlag == true) )
	{
	  if ((Manager.GetDouble("v-potential") == 0.0) && (Manager.GetDouble("w-potential") == 0.0))
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenvalueOutputFile, "%s_%s_gx_%g_gy_%g.dat",FilePrefix, FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      else
		sprintf (EigenvalueOutputFile, "%s_%s_gx_%g_gy_%g_mus_%g.dat",FilePrefix, FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	  else
	    {
	      if (Manager.GetDouble("w-potential") == 0.0)
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		}
	      else
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_w_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
			     FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_w_%g_%s_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
			     FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		}
	    }
	}
      else
	{
	  if ((Manager.GetDouble("v-potential") == 0.0) && (Manager.GetDouble("w-potential") == 0.0))
	    {
	      if (Manager.GetDouble("mu-s") == 0.0)
		sprintf (EigenvalueOutputFile, "%s_u_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	      else
		sprintf (EigenvalueOutputFile, "%s_u_%g_%s_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
	    }
	  else
	    {
	      if (Manager.GetDouble("w-potential") == 0.0)
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
			     FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
			     FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));		  
		}
	      else
		{
		  if (Manager.GetDouble("mu-s") == 0.0)
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_w_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
			     FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
		  else
		    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_w_%g_%s_gx_%g_gy_%g_mus_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
			     FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"));
		}
	    }
	}
    }

  if ((Manager.GetBoolean("singleparticle-spectrum") == true) || (Manager.GetBoolean("singleparticle-highsymmetryspectrum") == true))
    {
      if (Manager.GetBoolean("singleparticle-highsymmetryspectrum") == false)
	{
	  bool ExportOneBody = false;
	  bool BlochForm = false;
	  if (Manager.GetBoolean("real-space"))
	    BlochForm = true;
	  if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true) || (Manager.GetBoolean("singleparticle-berrycurvature") == true))
	    ExportOneBody = true;
	  TightBindingModelKagomeLattice TightBindingModel(NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
							   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody, BlochForm);
	  if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	    {
	      cout << "Chern number = " << TightBindingModel.ComputeChernNumber(0) << endl;
	    }
	  if (Manager.GetBoolean("singleparticle-berrycurvature") == true)
	    {
	      cout << "Chern number = " << TightBindingModel.ComputeBerryCurvature(0, ReplaceExtensionToFileName(EigenvalueOutputFile, "dat", "berrycurvature.dat")) << endl;
	    }
	  TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
	  double BandSpread = TightBindingModel.ComputeBandSpread(0);
	  double DirectBandGap = TightBindingModel.ComputeDirectBandGap(0);
	  cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
	  if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	    {
	      char* BandStructureOutputFile = new char [512];
	      if (Manager.GetString("export-onebodyname") != 0)
		strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	      else
		sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	      if (Manager.GetBoolean("export-onebody") == true)
		{
		  TightBindingModel.WriteBandStructure(BandStructureOutputFile);
		}
	      else
		{
		  TightBindingModel.WriteBandStructureASCII(BandStructureOutputFile);
		}
	      delete[] BandStructureOutputFile;
	    }	  
	  return 0;
	}
      else
	{
	  TightBindingModelKagomeLattice TightBindingModel(2, 2,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
							   Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
							   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 
							   Architecture.GetArchitecture(), false);
	  TightBindingModel.WriteAsciiSpectrumAlongHighSymmetryPoints(EigenvalueOutputFile, NbrSitesX);
	  return 0;	  
	}
    }

  int MinKx = 0;
  int MaxKx = NbrSitesX - 1;
  if ((Manager.GetBoolean("redundant-kx")==false) && (fabs(Manager.GetDouble("gamma-x"))<1e-12)) // want to reduce zone, and no offset?
    MaxKx = NbrSitesX/2;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  if (Manager.GetInteger("max-kx") >= MinKx)
    {						
      MaxKx = Manager.GetInteger("max-kx");
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY - 1;
  
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  if (Manager.GetInteger("max-ky") >= MinKy)
    {						
      MaxKy = Manager.GetInteger("max-ky");
    }
  
  Abstract2DTightBindingModel* TightBindingModel;
  if (Manager.GetString("import-onebody") == 0)
    {
      if (TiltedFlag == false)
	{
	  TightBindingModel = new TightBindingModelKagomeLattice (NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 							      Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());
	}
      else
	{
	  TightBindingModel = new TightBindingModelAlternativeKagomeLattice (NbrSitesX, NbrSitesY, Nx1, Ny1, Nx2, Ny2, Offset, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
									     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());
	}
//       TightBindingModel = new TightBindingModel2DAtomicLimitLattice(NbrSitesX, NbrSitesY, 3, 2, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());  
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      int TmpMaxKy = MaxKy;
      if ((i == 0) && (Manager.GetBoolean("redundant-ky") == false))
	TmpMaxKy = NbrSitesY / 2;
      for (int j = MinKy; j <= TmpMaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;

	  ParticleOnSphere* Space = 0;
	  AbstractHamiltonian* Hamiltonian = 0;
	  if (Manager.GetBoolean("real-space") == false)
	    {
	      if (Manager.GetBoolean("three-bands") == false)
		{
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
		  
		  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
		    { 
		      if (TiltedFlag == false)
			{
			  Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
												TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  Hamiltonian = new ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, TightBindingModel, 
													   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
													   Manager.GetDouble("w-potential"), 
													   Manager.GetInteger("band-index"), 
													   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}
		    }
		  else
		    { 
		      if (Manager.GetBoolean("three-body") == true)
			{
			  // use unit three-body interaction by default
			  Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian
			    (Space, NbrParticles, NbrSitesX, NbrSitesY, /* three-body */ 1.0, Manager.GetDouble("u-potential"),
			     Manager.GetDouble("v-potential"),  TightBindingModel, 			 Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			}
		      else
			{
			  if (Manager.GetBoolean("four-body") == true)
			    {
			      Hamiltonian = 0;
			      // 		      Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
			      // 												    Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
			      // 												    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
			      // 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			    }
			  else
			    {
			      Hamiltonian = 0;
			      // 		      Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandFiveBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
			      // 												    Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
			      // 												    Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
			      // 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
			    }
			  
			}
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("boson") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 20)
			{
			  Space = new FermionOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			}
		      else
			{
			  //		      Space = new FermionOnSquareLatticeWithSU3SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			}
		    }
		  else
		    {
		      Space = new BosonOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  //	      return 0;
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
		  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
		    { 
		      Hamiltonian = new ParticleOnLatticeKagomeLatticeThreeBandHamiltonian((ParticleOnSphereWithSU3Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY,
											   Manager.GetDouble("u-potential"), Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"),
											   Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
											   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      Hamiltonian = 0;
		    }
		  if (Manager.GetBoolean("project-threebands") == true)
		    {
		      ParticleOnSphere* TargetSpace = 0;
		      ComplexMatrix* OneBodyBasis = ComputeSingleParticleTransformationMatrices(NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"));
		      ComplexMatrix SU3U1TransformationMatrix;
		      ComplexMatrix TransformedHRep;
		      ComplexMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
		      Hamiltonian->GetHamiltonian(HRep);
		      if (Manager.GetBoolean("boson") == false)
			{
			  //		      ComplexMatrix NBodyTransformationMatrix = ((FermionOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
			  //		      TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
			  //		      TargetSpace = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			  //		      SU3U1TransformationMatrix = ((FermionOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixSU3ToU1((FermionOnSquareLatticeMomentumSpace*) TargetSpace);
			}
		      else
			{
			  ComplexMatrix NBodyTransformationMatrix = ((BosonOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
			  TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
			  TargetSpace = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
			  SU3U1TransformationMatrix = ((BosonOnSquareLatticeWithSU3SpinMomentumSpace*) Space)->TransformationMatrixSU3ToU1((BosonOnSquareLatticeMomentumSpace*) TargetSpace, 2);
			}
		      
		      ComplexMatrix TransformedHRep2 = TransformedHRep.InvConjugate(SU3U1TransformationMatrix);
		      if (Manager.GetDouble("u-potential") != 0.0)
			TransformedHRep2 /= Manager.GetDouble("u-potential");
		      
		      RealDiagonalMatrix TmpDiag;
		      HermitianMatrix HRep2(TransformedHRep2);
		      delete Hamiltonian;
		      delete Space;
		      delete[] OneBodyBasis;
		      Hamiltonian = new ExplicitHamiltonian(TargetSpace, &HRep2);
		      Space = TargetSpace;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("no-translation") == false)
		{
		  if (Manager.GetBoolean("boson") == false)
		    {
		      Space = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(),  i, NbrSitesX, j, NbrSitesY);
		    }
		  else
		    {
		      if (GutzwillerFlag == false)
			Space = new BosonOnLatticeRealSpaceAnd2DTranslation(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), i, NbrSitesX, j, NbrSitesY);
		      else
			Space = new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), i, NbrSitesX, j, NbrSitesY);
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("boson") == false)
		    {
		      Space = new FermionOnLatticeRealSpace (NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
		    }
		  else
		    {
		      if (GutzwillerFlag == false)
			Space = new BosonOnLatticeRealSpace(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
		      else
			Space = new BosonOnLatticeGutzwillerProjectionRealSpace (NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
		    }
		}
	      
	      
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
	      
	      RealSymmetricMatrix DensityDensityInteraction(TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), true);
	      int* NbrInteractingOrbitals;
	      int** InteractingOrbitalsOrbitalIndices;
	      int** InteractingOrbitalsSpatialIndices;
	      double** InteractingOrbitalsPotentials;
	      FCIKagomeLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
							      InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
							      Manager.GetBoolean("boson"), Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), TightBindingModel);
	      for (int x = 0; x < NbrSitesX; ++x)
		{
		  for (int y = 0; y < NbrSitesY; ++y)
		    {
		      for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModel->GetNbrBands(); ++OrbitalIndex)
			{
			  for (int k = 0; k < NbrInteractingOrbitals[OrbitalIndex]; ++k)
			    {
			      // 		    cout << InteractingOrbitalsSpatialIndices[OrbitalIndex][2 * k] << " " << InteractingOrbitalsSpatialIndices[OrbitalIndex][(2 * k) + 1] << " " << InteractingOrbitalsOrbitalIndices[OrbitalIndex][k] << " " << InteractingOrbitalsPotentials[OrbitalIndex][k] << endl;
			      DensityDensityInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex), TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndices[OrbitalIndex][2 * k], y + InteractingOrbitalsSpatialIndices[OrbitalIndex][(2 * k) + 1], InteractingOrbitalsOrbitalIndices[OrbitalIndex][k]), InteractingOrbitalsPotentials[OrbitalIndex][k]);
			      
			    }
			}
		    }
		}
	      if (Manager.GetBoolean("no-translation") == false)
		Hamiltonian = new ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), i, NbrSitesX, j, NbrSitesY,
											 TightBindingMatrix, DensityDensityInteraction,
											 Architecture.GetArchitecture(), Memory);
	      else
		Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(),  TightBindingMatrix, DensityDensityInteraction,
									 Architecture.GetArchitecture(), Memory);
	    }
	  
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile;
	  if (Manager.GetString("eigenstate-file") != 0)
	    {
	      EigenstateOutputFile = new char [512];
	      sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
	    }
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
	  
	  if (Manager.GetBoolean("no-translation") == true)
	    return 0;
	}
    }
  return 0;
}


// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites

ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double t1, double t2, double lambda1, double lambda2)
{
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[nbrSitesX * nbrSitesY];
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  double NNHopping = t1;
  double NNSpinOrbit = lambda1;
  double NextNNHopping = t2;
  double NextNNSpinOrbit = lambda2;
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  int Index = (kx * nbrSitesY) + ky;
	  double KX = 0.5 * KxFactor * ((double) kx);
	  double KY = 0.5 * KyFactor * ((double) ky);
	  Complex HAB (-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HAB *= cos (KX);
	  Complex HAC(-2.0 * NNHopping, 2.0 * NNSpinOrbit);
	  HAC *= cos (KY);
	  Complex HBC(-2.0 * NNHopping, -2.0 * NNSpinOrbit);
	  HBC *= cos(KX - KY);
	  
	  Complex HAB2 (-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HAB2 *= cos (KX - 2.0 * KY);
	  Complex HAC2 (-2.0 * NextNNHopping, -2.0 * NextNNSpinOrbit);
	  HAC2 *= cos (2.0 * KX - KY);
	  Complex HBC2 (-2.0 * NextNNHopping, 2.0 * NextNNSpinOrbit);
	  HBC2 *= cos (KX + KY);
	  
	  HAB += HAB2;
	  HAC += HAC2;
	  HBC += HBC2;
	  
	  HermitianMatrix TmpOneBodyHamiltonian(3, true);	  
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	  RealDiagonalMatrix TmpDiag;	  
	  OneBodyBasis[Index] = ComplexMatrix (3, 3);
	  OneBodyBasis[Index].SetToIdentity();
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis[Index]);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis[Index]);
#endif   
	}
    }
  return OneBodyBasis;
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

void FCIKagomeLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential, Abstract2DTightBindingModel* tightBindingModel)
{
  int nbrBands = tightBindingModel->GetNbrBands();
  
  nbrInteractingOrbitals = new int[nbrBands];
  interactingOrbitalsOrbitalIndices = new int*[nbrBands];
  interactingOrbitalsSpatialIndices = new int*[nbrBands];
  interactingOrbitalsPotentials = new double*[nbrBands];
  int p;
  int q;
  if (bosonFlag == false)
    {
      nbrInteractingOrbitals[0] = 2; 
      nbrInteractingOrbitals[1] = 2;      
      nbrInteractingOrbitals[2] = 2;      
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 4; 
	  nbrInteractingOrbitals[1] += 4;      
	  nbrInteractingOrbitals[2] += 4;           
	}
      for (int i = 0; i < nbrBands; ++i)
	{
	  interactingOrbitalsOrbitalIndices[i] = new int[nbrInteractingOrbitals[i]];
	  interactingOrbitalsSpatialIndices[i] = new int[2 * nbrInteractingOrbitals[i]];
	  interactingOrbitalsPotentials[i] = new double[nbrInteractingOrbitals[i]];
	}

    int TmpIndex = 0;

  // links starting from A
    interactingOrbitalsOrbitalIndices[0][TmpIndex] = 1;
    interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 0;
    interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    interactingOrbitalsPotentials[0][TmpIndex] = uPotential;
    ++TmpIndex;
    interactingOrbitalsOrbitalIndices[0][TmpIndex] = 2;
    interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 0;
    interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    interactingOrbitalsPotentials[0][TmpIndex] = uPotential;
    ++TmpIndex;
    
    TmpIndex -= 2;

    // links starting from B
    interactingOrbitalsOrbitalIndices[1][TmpIndex] = 0;
    interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 1;
    interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    interactingOrbitalsPotentials[1][TmpIndex] = uPotential;
    ++TmpIndex;
    interactingOrbitalsOrbitalIndices[1][TmpIndex] = 2;
    interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 0;
    interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    interactingOrbitalsPotentials[1][TmpIndex] = uPotential;
    ++TmpIndex;
   
    TmpIndex -= 2;

    // links starting from C
    interactingOrbitalsOrbitalIndices[2][TmpIndex] = 0;
    interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = 0;
    interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
    interactingOrbitalsPotentials[2][TmpIndex] = uPotential;
    ++TmpIndex;
    interactingOrbitalsOrbitalIndices[2][TmpIndex] = 1;
    interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = -1;
    interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
    interactingOrbitalsPotentials[2][TmpIndex] = uPotential;
    ++TmpIndex;

    if (vPotential != 0.0)
      {
	cout << "Warning: next nearest neighbor hopping not implemented in real space" << endl;
      }


    }
  else
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 1;
      nbrInteractingOrbitals[2] = 1;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 2;
	  nbrInteractingOrbitals[1] += 2;
	  nbrInteractingOrbitals[2] += 2;
	}
	
      for (int j = 0; j < nbrBands; ++j)
      {
	interactingOrbitalsOrbitalIndices[j] = new int[nbrInteractingOrbitals[j]];
	interactingOrbitalsSpatialIndices[j] = new int[nbrInteractingOrbitals[j] * 2];
	interactingOrbitalsPotentials[j] = new double[nbrInteractingOrbitals[j]];
      }
        
      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 0;
      interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[0][Index] = 0.5 * uPotential;
      ++Index;

      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 1;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[1][Index] = 0.5 * uPotential;
      ++Index;
      
      Index = 0;
      interactingOrbitalsOrbitalIndices[2][Index] = 2;
      interactingOrbitalsSpatialIndices[2][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[2][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[2][Index] = 0.5 * uPotential;
      ++Index;

      if (vPotential != 0.0)
	{
	  int TmpIndex = 1;

	  // links starting from A
	  interactingOrbitalsOrbitalIndices[0][TmpIndex] = 1;
	  interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 0;
	  interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
	  interactingOrbitalsPotentials[0][TmpIndex] = vPotential;
	  ++TmpIndex;
	  interactingOrbitalsOrbitalIndices[0][TmpIndex] = 2;
	  interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 0;
	  interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
	  interactingOrbitalsPotentials[0][TmpIndex] = vPotential;
	  ++TmpIndex;
    
	  TmpIndex -= 2;

	  // links starting from B
	  interactingOrbitalsOrbitalIndices[1][TmpIndex] = 0;
	  interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 1;
	  interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
	  interactingOrbitalsPotentials[1][TmpIndex] = vPotential;
	  ++TmpIndex;
	  interactingOrbitalsOrbitalIndices[1][TmpIndex] = 2;
	  interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 0;
	  interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
	  interactingOrbitalsPotentials[1][TmpIndex] = vPotential;
	  ++TmpIndex;
   
	  TmpIndex -= 2;

	  // links starting from C
	  interactingOrbitalsOrbitalIndices[2][TmpIndex] = 0;
	  interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = 0;
	  interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
	  interactingOrbitalsPotentials[2][TmpIndex] = vPotential;
	  ++TmpIndex;
	  interactingOrbitalsOrbitalIndices[2][TmpIndex] = 1;
	  interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = -1;
	  interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
	  interactingOrbitalsPotentials[2][TmpIndex] = vPotential;
	  ++TmpIndex;
	}
    }
}


// find the tilted lattices with an aspect ratio close to unity 
//
// nbrSiteX = number of sites in the x direction (max(nbrSiteX, nbrSiteY) sets the maximal value for |nx1|, |ny1|, |nx2| and |ny2|)
// nbrSiteY = number of sites in the y direction (max(nbrSiteX, nbrSiteY) sets the maximal value for |nx1|, |ny1|, |nx2| and |ny2|)

void FCIKagomeLatticeModelFindOptimalTiltedLattices(int nbrSiteX, int nbrSiteY)
{
  int NbrUnitCells = nbrSiteX * nbrSiteY;
  int MaxN = nbrSiteX;
  if (MaxN < nbrSiteY)
    {
      MaxN = nbrSiteY;
    }
  int MaxNbrLattices = (2 * MaxN + 1) * (MaxN + 1) * (MaxN +1);
  double* AspectRatios = new double [MaxNbrLattices];
  int* TmpIndices = new int [MaxNbrLattices];
  int* Nx1Values = new int [MaxNbrLattices];
  int* Ny1Values = new int [MaxNbrLattices];
  int* Nx2Values = new int [MaxNbrLattices];
  int* Ny2Values = new int [MaxNbrLattices];
  int Index = 0;
  double CosAngle = cos (M_PI / 3.0);
  double SinAngle = sin (M_PI / 3.0);
  for (int Nx1 = 0; Nx1 <= MaxN; ++Nx1)
    {
      for (int Ny1 = 0; Ny1 <= MaxN; ++Ny1)
	{
	  for (int Nx2 = -MaxN; Nx2 <= MaxN; ++Nx2)
	    {
	      for (int Ny2 = -MaxN; Ny2 <= MaxN; ++Ny2)
		{
		  if (((Nx1 * Ny2) - (Nx2 * Ny1)) == NbrUnitCells)
		    {
		      double TmpAspectRatio = ((double) NbrUnitCells) * SinAngle / (((double) ((Nx1 * Nx1) + (Ny1 * Ny1))) 
										    + 2.0 * ((double) (Nx1 * Ny1)) * CosAngle);
		      if (TmpAspectRatio > 1.0)
			{
			  AspectRatios[Index] = 1.0 / TmpAspectRatio;			  
			}
		      else
			{
			  AspectRatios[Index] = TmpAspectRatio;			  
			}
		      Nx1Values[Index] = Nx1;
		      Ny1Values[Index] = Ny1;
		      Nx2Values[Index] = Nx2;
		      Ny2Values[Index] = Ny2;
		      TmpIndices[Index] = Index;
		      ++Index;
		    }
		}
	    }
	}
    }
  SortArrayDownOrdering<int>(AspectRatios, TmpIndices, Index);
  for (int i = 0; i < 100; ++i)
    {
      int Nx1 = Nx1Values[TmpIndices[i]];
      int Ny1 = Ny1Values[TmpIndices[i]];
      int Nx2 = Nx2Values[TmpIndices[i]];
      int Ny2 = Ny2Values[TmpIndices[i]];
      int Ny = FindGCD(FindGCD(NbrUnitCells, abs(Nx2)), abs(Ny2));
      int Nx = NbrUnitCells / Ny;
      double TmpAspectRatio = ((double) NbrUnitCells) * SinAngle / (((double) ((Nx1 * Nx1) + (Ny1 * Ny1))) 
								    + 2.0 * ((double) (Nx1 * Ny1)) * CosAngle);
      int Offset = 0;
      bool Flag = false;
      while ((Offset <= 20) && (Flag == false))
	{
	  if (((abs((Offset * Ny2) - Ny1) % Nx) == 0) && ((abs(Nx1 - (Offset * Nx2)) % Nx) == 0))
	    {
	      Flag = true;
	    }
	  else
	    {
	      ++Offset;
	    }
	}
      if (Flag == false)
	{
	  Offset = -1;
	  while ((Offset >= -20) && (Flag == false))
	    {
	      if (((abs((Offset * Ny2) - Ny1) % Nx) == 0) && ((abs(Nx1 - (Offset * Nx2)) % Nx) == 0))
		{
		  Flag = true;
		}
	      else
		{
		  --Offset;
		}
	    }
	}
      if (Flag == true)
	{
	  cout << "Kappa=" <<TmpAspectRatio << " Nx=" << Nx << " Ny=" << Ny << " nx1=" << Nx1 << " ny1=" << Ny1 
	       << " nx2=" << Nx2 << " ny2=" << Ny2 << "offset=" << Offset << " : " << " -x " << Nx << " -y " << Ny << " --nx1 " << Nx1 << " --ny1 " << Ny1 
	       << " --nx2 " << Nx2 << " --ny2 " << Ny2 << " --offset " << Offset << endl;
	}
    }
  delete[] AspectRatios;
  delete[] TmpIndices;
  delete[] Nx1Values;
  delete[] Ny1Values;
  delete[] Nx2Values;
  delete[] Ny2Values;
}
