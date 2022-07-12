#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
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
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian.h"
 
#include "Tools/FTITightBinding/TightBindingModel3DChargePumpedChiralHingeCubicLattice.h"
#include "Tools/FTITightBinding/Generic3DTightBindingModel.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

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
void FHI3DChargePumpedChiralHingeCubicLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
									     int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
									     bool bosonFlag, double uPotential, double vPotential,
									     Abstract3DTightBindingModel* tightBindingModel);


int main(int argc, char** argv)
{
  OptionManager Manager ("FHI3DChargePumpedChiralHingeCubicLatticeModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kz", "only evalute a given z momentum sector (negative if all kz sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the gutzwiller projected Hilbert space");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive nearest next neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "neareast neighbor hopping amplitude within the unit cell", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "neareast neighbor hopping amplitude between unit cells in the XY plane", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tz", "twice the neareast neighbor and next nearest neighbor hopping amplitude along z ", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential on even vs odd sites", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-z", "boundary condition twisting angle along z (in 2 Pi unit)", 0.0);
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
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "when using real space representation, discard all translations");
  (*SystemGroup) += new BooleanOption  ('\n', "xy-translation", "when using real space representation, only considertranslations in the x and y direction");
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
      cout << "see man page for option syntax or type FHI3DChargePumpedChiralHingeCubicLatticeModel -h" << endl;
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
  int NbrSitesZ = Manager.GetInteger("nbr-sitez"); 
  int NbrSites = 4 * NbrSitesX * NbrSitesY * NbrSitesZ;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  bool TiltedFlag = false;
  int Nx1 = 1;
  int Ny1 = 0;
  int Nx2 = 0;
  int Ny2 = 0;
  
  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }


  char* FileSystemGeometry = new char [512];
  if (TiltedFlag == false)
    {
      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_z_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, NbrSitesZ);
    }
  else
    {
      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_z_%d_tilted_nx1_%d_ny1_%d_nx2_%d_ny2_%d", NbrParticles,
	       NbrSites, NbrSitesX, NbrSitesY, NbrSitesZ, Nx1, Ny1, Nx2, Ny2);
    }
  
  char* FilePrefix = new char [512 + strlen(FileSystemGeometry)];
  if (Manager.GetBoolean("single-band") == false)
    {
      if (Manager.GetBoolean("real-space") == false)
	{
	  sprintf (FilePrefix, "%s_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
	}
      else
	{
	  if (Manager.GetBoolean("no-translation") == false)
	    {
	      if (Manager.GetBoolean("xy-translation") == false)
		{
		  if (Manager.GetBoolean ("gutzwiller") == false)
		    {
		      sprintf (FilePrefix, "%s_realspace_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
		    }
		  else
		    {
		      sprintf (FilePrefix, "%s_realspace_gutzwiller_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
		    }
		}
	      else
		{
		  if (Manager.GetBoolean ("gutzwiller") == false)
		    {
		      sprintf (FilePrefix, "%s_realspace_xytranslation_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
		    }
		  else
		    {
		      sprintf (FilePrefix, "%s_realspace_xytranslation_gutzwiller_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
		    }
		}
	    }
	  else
	    {
	     if (Manager.GetBoolean ("gutzwiller") == false)
	       {
		 sprintf (FilePrefix, "%s_realspace_notranslation_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
	       }
	     else
	       {
		 sprintf (FilePrefix, "%s_realspace_gutzwiller_notranslation_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
	       }
	   }
	}
    }
  else
    {
	sprintf (FilePrefix, "%s_singleband_3dchiralhinge_%s", StatisticPrefix, FileSystemGeometry);
    }

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky kz");
  char* FileParameterString = new char [512];
  if (Manager.GetDouble("mu-s") == 0.0)
    {
      sprintf (FileParameterString, "t1_%.3f_t2_%.3f_tz_%.3f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tz"));
    }
  else
    {
      sprintf (FileParameterString, "t1_%.3f_t2_%.3f_tz_%.3f_mus_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tz"), Manager.GetDouble("mu-s"));
    }
  
  char* FileTwistedBoundaryConditions = new char [512];
  sprintf (FileTwistedBoundaryConditions, "gx_%.4f_gy_%.4f_gz_%.4f", Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
  
  char* EigenvalueOutputFile = new char [512 + strlen(FileTwistedBoundaryConditions) + strlen(FilePrefix) + strlen(FileParameterString)];
  
  if (Manager.GetString("eigenvalue-file") != 0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetBoolean("single-band") == false)
	{
	  sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_%s.dat", FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"),
		   FileParameterString, FileTwistedBoundaryConditions);
	}
      else
	{
	  if (Manager.GetBoolean("flat-band") == true)
	    {
	      sprintf (EigenvalueOutputFile, "%s_v_%f_%s_%s.dat", FilePrefix, Manager.GetDouble("v-potential"),
		       FileParameterString, FileTwistedBoundaryConditions);
	    }
	  else
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_%s.dat", FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"),
		       FileParameterString, FileTwistedBoundaryConditions);
	    }
	}
    }
  
  Abstract3DTightBindingModel* TightBindingModel;
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;
      if (TiltedFlag == false)
	{
	  TightBindingModel = new TightBindingModel3DChargePumpedChiralHingeCubicLattice (NbrSitesX, NbrSitesY, NbrSitesZ,
											  Manager.GetDouble("t1"), Manager.GetDouble("t2"),
											  Manager.GetDouble("tz"), Manager.GetDouble("tz"),
											  Manager.GetDouble("mu-s"),
											  Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"),
											  Architecture.GetArchitecture(), ExportOneBody);
	}
      else
	{
	  cout << "tilted lattices are not implemented" << endl;
	  return 0;
	}
      TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel->ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel->ComputeDirectBandGap(1);
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
  int MinKz = 0;
  int MaxKz = NbrSitesZ - 1;
  if (Manager.GetInteger("only-kz") >= 0)
    {						
      MinKz = Manager.GetInteger("only-kz");
      MaxKz = MinKz;
    }
  if(Manager.GetBoolean("no-translation") == true)
    {  
      MaxKx = 0;
      MaxKy = 0;
      MaxKz = 0;
    }
  if (Manager.GetBoolean("xy-translation") == true)
    {
      MaxKz = 0;
    }

  
  if (Manager.GetString("import-onebody") == 0)
    {
      if (TiltedFlag == false)
	{
	  TightBindingModel = new TightBindingModel3DChargePumpedChiralHingeCubicLattice (NbrSitesX, NbrSitesY, NbrSitesZ,
											  Manager.GetDouble("t1"), Manager.GetDouble("t2"),
											  Manager.GetDouble("tz"), Manager.GetDouble("tz"),
											  Manager.GetDouble("mu-s"),
											  Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"),
											  Architecture.GetArchitecture(), true);
	}
      else
	{
	  cout << "tilted lattices are not implemented" << endl;
	  return 0;
	}
      
      char* BandStructureOutputFile = new char [64 + strlen(FilePrefix) + strlen(FileParameterString) + strlen(FileTwistedBoundaryConditions)];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString, FileTwistedBoundaryConditions);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic3DTightBindingModel(Manager.GetString("import-onebody")); 
    }

  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  for (int k = MinKz; k <= MaxKz; ++k)
	    {
	      cout << "(kx=" << i << ",ky=" << j << ",kz=" << k << ") : " << endl;
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
		  FHI3DChargePumpedChiralHingeCubicLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices, 
											  InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
											  Manager.GetBoolean("boson"), Manager.GetDouble("u-potential"),
											  Manager.GetDouble("v-potential"), TightBindingModel);
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
		      RealSymmetricMatrix DensityDensityInteraction(TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), true);
		      for (int x = 0; x < NbrSitesX; ++x)
			{
			  for (int y = 0; y < NbrSitesY; ++y)
			    {
			      for (int z = 0; z < NbrSitesZ; ++z)
				{
				  for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModel->GetNbrBands(); ++OrbitalIndex)
				    {
				      for (int l = 0; l < NbrInteractingOrbitals[OrbitalIndex]; ++l)
					{
					  DensityDensityInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y, z, OrbitalIndex), 
										       TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndices[OrbitalIndex][3 * l], 
																		      y + InteractingOrbitalsSpatialIndices[OrbitalIndex][(3 * l) + 1],
																		      z + InteractingOrbitalsSpatialIndices[OrbitalIndex][(3 * l) + 2],
																		      InteractingOrbitalsOrbitalIndices[OrbitalIndex][l]), 
										       InteractingOrbitalsPotentials[OrbitalIndex][l]);
					  
					}
				    }
				}
			    }
			}
		      if (Manager.GetBoolean("boson") == true)
			{
			  if (Manager.GetBoolean ("gutzwiller") == false)
			    {
			      if (Manager.GetBoolean("no-translation") == true)
				Space = new BosonOnLatticeRealSpace(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
			      else
				Space = new BosonOnLatticeRealSpaceAnd2DTranslation(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), 
										    i, NbrSitesX, j, NbrSitesY);
			    }
			  else
			    {
			      if (Manager.GetBoolean("no-translation") == true)
				Space = new BosonOnLatticeGutzwillerProjectionRealSpace(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
			      else
				Space = new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(NbrParticles, 
													TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), 
													i, NbrSitesX, j, NbrSitesY);
			    }
			}
		      else
			{
			  if (Manager.GetBoolean("no-translation") == true)
			    Space = new FermionOnLatticeRealSpace(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());
			  else
			    Space = new FermionOnLatticeRealSpaceAnd2DTranslation(NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), 
										  i, NbrSitesX, j, NbrSitesY);
			}
		      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		      HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
		      if (Manager.GetBoolean("no-translation") == true)
			{
			  Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), 
										   TightBindingMatrix, DensityDensityInteraction,
										   Architecture.GetArchitecture(), Memory);
			}
		      else
		    {
		      Hamiltonian = new ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), 
											       i, NbrSitesX, j, NbrSitesY,
											       TightBindingMatrix, DensityDensityInteraction,
											       Architecture.GetArchitecture(), Memory);
		    }
		    }
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, k);
		  char* EigenstateOutputFile = new char [512];
		  if ((Manager.GetBoolean("real-space") == false) || (Manager.GetBoolean("no-translation") == false))
		    {
		      char* TmpExtention = new char [512];
		      if (Manager.GetBoolean("xy-translation") == true)
			{
			  sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
			}
		      else
			{
			  sprintf (TmpExtention, "_kx_%d_ky_%d_kz_%d", i, j, k);
			}
		      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
		      delete[] TmpExtention;
		    }
		  else
		    {
		      char* TmpExtention = new char [512];
		      sprintf (TmpExtention, "");
		      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
		      delete[] TmpExtention;
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
		  //  test code for the generic density-density
		  int* NbrInteractingOrbitals;
		  int** InteractingOrbitalsOrbitalIndices;
		  int** InteractingOrbitalsSpatialIndices;
		  double** InteractingOrbitalsPotentials;
		  // FHI3DChargePumpedChiralHingeCubicLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices, 
		  // 									  InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
		  // 									  Manager.GetBoolean("boson"), Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"),
		  // 									  TightBindingModel);
// 		  Hamiltonian = new ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, 0,
// 													   NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
// 													   InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
// 													   TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  Hamiltonian = new ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
											      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
											      TightBindingModel,Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, k);
		  char* EigenstateOutputFile = new char [512];
		  if (Manager.GetString("eigenstate-file") != 0)
		    {
		      sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d_kz_%d", Manager.GetString("eigenstate-file"), i, j, k);
		    }
		  else
		    {
		      char* TmpExtention = new char [512];
		      sprintf (TmpExtention, "_kx_%d_ky_%d_kz_%d", i, j, k);
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

void FHI3DChargePumpedChiralHingeCubicLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
									     int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
									     bool bosonFlag, double uPotential, double vPotential,
									     Abstract3DTightBindingModel* tightBindingModel)
{
  nbrInteractingOrbitals = new int[4];
  interactingOrbitalsOrbitalIndices = new int*[4];
  interactingOrbitalsSpatialIndices = new int*[4];
  interactingOrbitalsPotentials = new double*[4];
  int p;
  int q;
  int r;
  if (bosonFlag == false)
    {
      nbrInteractingOrbitals[0] = 2;
      nbrInteractingOrbitals[1] = 2;
      nbrInteractingOrbitals[2] = 2;
      nbrInteractingOrbitals[3] = 2;
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 3];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 3];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];
      interactingOrbitalsOrbitalIndices[2] = new int[nbrInteractingOrbitals[2]];
      interactingOrbitalsSpatialIndices[2] = new int[nbrInteractingOrbitals[2] * 3];
      interactingOrbitalsPotentials[2] = new double[nbrInteractingOrbitals[2]];
      interactingOrbitalsOrbitalIndices[3] = new int[nbrInteractingOrbitals[3]];
      interactingOrbitalsSpatialIndices[3] = new int[nbrInteractingOrbitals[3] * 3];
      interactingOrbitalsPotentials[3] = new double[nbrInteractingOrbitals[3]];

      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 1;
      tightBindingModel->GetRealSpaceIndex(0, 0, 0, &(interactingOrbitalsSpatialIndices[0][3 * Index]));
      interactingOrbitalsPotentials[0][Index] = uPotential;
      ++Index;
      interactingOrbitalsOrbitalIndices[0][Index] = 3;
      tightBindingModel->GetRealSpaceIndex(0, 0, 0, &(interactingOrbitalsSpatialIndices[0][3 * Index]));
      interactingOrbitalsPotentials[0][Index] = uPotential;
      ++Index;
      
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 2;
      tightBindingModel->GetRealSpaceIndex(0, 0, 0, &(interactingOrbitalsSpatialIndices[1][3 * Index]));
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      tightBindingModel->GetRealSpaceIndex(1, 0, 0, &(interactingOrbitalsSpatialIndices[1][3 * Index]));
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;

      Index = 0;
      interactingOrbitalsOrbitalIndices[2][Index] = 3;
      tightBindingModel->GetRealSpaceIndex(1, 0, 0, &(interactingOrbitalsSpatialIndices[2][3 * Index]));
      interactingOrbitalsPotentials[2][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[2][Index] = 1;
      tightBindingModel->GetRealSpaceIndex(0, 1, 0, &(interactingOrbitalsSpatialIndices[2][3 * Index]));
      interactingOrbitalsPotentials[2][Index] = uPotential;		  
      ++Index;

      Index = 0;
      interactingOrbitalsOrbitalIndices[3][Index] = 2;
      tightBindingModel->GetRealSpaceIndex(0, 0, 0, &(interactingOrbitalsSpatialIndices[3][3 * Index]));
      interactingOrbitalsPotentials[3][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[3][Index] = 0;
      tightBindingModel->GetRealSpaceIndex(0, 1, 0, &(interactingOrbitalsSpatialIndices[3][3 * Index]));
      interactingOrbitalsPotentials[3][Index] = uPotential;		  
      ++Index;      
    }
  else
    {
      nbrInteractingOrbitals[0] = 0;
      nbrInteractingOrbitals[1] = 0;
      nbrInteractingOrbitals[2] = 0;
      nbrInteractingOrbitals[3] = 0;
  //     if (vPotential != 0.0)
  // 	{
  // 	  nbrInteractingOrbitals[0] += 1;
  // 	  nbrInteractingOrbitals[1] += 3;
  // 	}
      // interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      // interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 3];
      // interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      // interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      // interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 3];
      // interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];
  
  //     int Index = 0;
  //     interactingOrbitalsOrbitalIndices[0][Index] = 0;
  //     tightBindingModel->GetRealSpaceIndex(0, 0, p, q, r);
  //     interactingOrbitalsSpatialIndices[0][2 * Index] = p;
  //     interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
  //     interactingOrbitalsPotentials[0][Index] = 0.5 * uPotential;
  //     ++Index;
  //     if (vPotential != 0.0)
  // 	{
  // 	  interactingOrbitalsOrbitalIndices[0][Index] = 1;
  // 	  tightBindingModel->GetRealSpaceIndex(0, 0, p, q, r);
  // 	  interactingOrbitalsSpatialIndices[0][2 * Index] = p;
  // 	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = q;
  // 	  interactingOrbitalsPotentials[0][Index] = vPotential;	  
  // 	}
  //     Index = 0;
  //     interactingOrbitalsOrbitalIndices[1][Index] = 1;
  //     tightBindingModel->GetRealSpaceIndex(0, 0, p, q, r);
  //     interactingOrbitalsSpatialIndices[1][2 * Index] = p;
  //     interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
  //     interactingOrbitalsPotentials[1][Index] = 0.5 * uPotential;
  //     ++Index;
  //     if (vPotential != 0.0)
  // 	{
  // 	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
  // 	  tightBindingModel->GetRealSpaceIndex(1, 0, p, q, r);
  // 	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
  // 	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
  // 	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
  // 	  ++Index;
  // 	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
  // 	  tightBindingModel->GetRealSpaceIndex(0, 1, p, q, r);
  // 	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
  // 	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
  // 	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
  // 	  ++Index;
  // 	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
  // 	  tightBindingModel->GetRealSpaceIndex(1, 1, p, q, r);
  // 	  interactingOrbitalsSpatialIndices[1][2 * Index] = p;
  // 	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = q;
  // 	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
  // 	  ++Index;
  // 	}
    }
}
