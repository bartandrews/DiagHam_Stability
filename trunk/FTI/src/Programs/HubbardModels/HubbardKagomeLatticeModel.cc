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

#include "Hamiltonian/ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int GetRealSpaceTightBindingLinearizedIndexSafe(int indexX, int indexY, int indexOrbital, int nbrSitesX, int nbrSitesY);

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

void KagomeLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double NNInteraction, double NNNInteraction);


// compute the coordinate of a given point in the kagome lattice from their real space coordinates (trivial for a non-tilted lattice, same as GetRealSpaceIndex of Abstract2DTightBindingModel)
//
// i = cartesian coordinate in the x direction of the Bravais lattice
// j = cartesian coordinate in the y direction of the Bravais lattice
// p = reference on the first lattice index
// q = reference on the second lattice index
// OffsetReal = offset
void GetRealSpaceIndex (int i, int j, int& p, int& q, int offsetReal);


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardKagomeLatticeModel" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles, same as number of sites if negative", -1);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 2);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx1", "first coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny1", "second coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx2", "first coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny2", "second coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "real-offset", "second coordinate in real space of the second spanning vector of the real space lattice (0 if lattice is untilted)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
//   (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site (Hubbard) potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "jx", "strength of the neareast neighbor SxSx interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "jy", "strength of the neareast neighbor SySy interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "jz", "strength of the neareast neighbor SzSz interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('a', "anisotropy", "ratio between up and down nearest neighbor interaction, when positive", -1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "NNNjx", "strength of the next neareast neighbor SxSx interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "NNNjy", "strength of the next neareast neighbor SySy interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "NNNjz", "strength of the next neareast neighbor SzSz interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "jD", "strength of the third neareast neighbor SzSz interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "DM", "strength of the neareast neighbor Dzyaloshinskii-Moriya interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "chi", "strength of the chirality term (on triangles)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "inserted flux in the y direction", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "cylinder", "use periodic boundary conditions in one direction (y) only");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "do not use the code with 2D translations when the system is a torus");
  (*SystemGroup) += new BooleanOption  ('\n', "fixed-sz", "use the conservation of Sz");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-value", "twice the value of Sz", 0); 
//   (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
//   (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
//   (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
//   (*SystemGroup) += new SingleStringOption ('\n', "import-onebody", "import information on the tight binding model from a file");
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
      cout << "see man page for option syntax or type HubbardKitaevHeisenbergHoneycombModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  int NbrSites = 3 * NbrSitesX * NbrSitesY; 
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  if (NbrParticles < 0)
    NbrParticles = NbrSites;
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  int* SitesA = 0;
  int* SitesB = 0;
  bool CylinderFlag = Manager.GetBoolean("cylinder");
  bool NoTranslationFlag = Manager.GetBoolean("no-translation");
  if (Manager.GetBoolean("cylinder"))
    NoTranslationFlag = true;
  
  double TmpJx = Manager.GetDouble("jx");
  double TmpJy = Manager.GetDouble("jy");
  double TmpJz = Manager.GetDouble("jz");
  double AnisotropyFactor = 1.0;
  if (Manager.GetDouble("anisotropy") >= 0)
    AnisotropyFactor = Manager.GetDouble("anisotropy");
  double TmpJxDown = AnisotropyFactor * TmpJx;
  double TmpJyDown = AnisotropyFactor * TmpJy;
  double TmpJzDown = AnisotropyFactor * TmpJz;
  double TmpDM = Manager.GetDouble("DM");  
  double GammaY = 2.0 * M_PI * Manager.GetDouble("gamma-y") / ((double) NbrSitesY);  
  double TmpNNNJx = Manager.GetDouble("NNNjx");
  double TmpNNNJy = Manager.GetDouble("NNNjy");
  double TmpNNNJz = Manager.GetDouble("NNNjz");
  double TmpJD = Manager.GetDouble("jD");
  int szSector =  Manager.GetInteger("sz-value");
  
  bool ConserveSz = Manager.GetBoolean("fixed-sz");
  if (ConserveSz && ((TmpJx != TmpJy) or (TmpNNNJx != TmpNNNJy)))
  {
   cout << "Error, there is no Sz symmetry with these parameter values" << endl; 
   return 0;
  }
 
 
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
 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      if (SzSymmetryFlag == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "fermions_kagome_heisenberg");
	  else
	  {
	    if (NbrParticles != NbrSites)
	      sprintf (StatisticPrefix, "fermions_kagome_heisenberg_gutzwiller");
	    else
	      sprintf (StatisticPrefix, "spin_kagome_heisenberg");
	  }
	}
      else
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "fermions_kagome_heisenberg_szsym");
	  else
	    sprintf (StatisticPrefix, "fermions_kitaev_heisenberg_gutzwiller_szsym");
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kagome_heisenberg");
	  else
	    sprintf (StatisticPrefix, "bosons_kagome_heisenberg_gutzwiller");
	}
      else
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kagome_heisenberg_szsym");
	  else
	    sprintf (StatisticPrefix, "bosons_kagome_heisenberg_gutzwiller_szsym");
	}
    }
    
  

  char* FilePrefix = new char [256];
  if (CylinderFlag == false)
  {
    if (Manager.GetBoolean("no-translation"))
      sprintf (FilePrefix, "%s_notranslation_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    else
      sprintf (FilePrefix, "%s_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
  }
  else
    sprintf (FilePrefix, "%s_cylinder_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
  
  char* FileParameterString = new char [256];
  if (TiltedFlag == false)
    sprintf (FileParameterString, "a_%.6f_jx_%.6f_jy_%.6f_jz_%.6f_NNNjx_%.6f_NNNjy_%.6f_NNNjz_%.6f_jD_%.6f_DM_%.6f_gammay_%.6f", AnisotropyFactor, Manager.GetDouble("jx"), Manager.GetDouble("jy"), Manager.GetDouble("jz"), Manager.GetDouble("NNNjx"), Manager.GetDouble("NNNjy"), Manager.GetDouble("NNNjz"), Manager.GetDouble("jD"), Manager.GetDouble("DM"), Manager.GetDouble("gamma-y"));
  else
    sprintf (FileParameterString, "nx1_%d_ny1_%d_nx2_%d_ny2_%d_off_%d_a_%.6f_jx_%.6f_jy_%.6f_jz_%.6f_NNNjx_%.6f_NNNjy_%.6f_NNNjz_%.6f_jD_%.6f_DM_%.6f_gammay_%.6f", nx1, ny1, nx2, ny2, OffsetReal, AnisotropyFactor, Manager.GetDouble("jx"), Manager.GetDouble("jy"), Manager.GetDouble("jz"), Manager.GetDouble("NNNjx"), Manager.GetDouble("NNNjy"), Manager.GetDouble("NNNjz"), Manager.GetDouble("jD"), Manager.GetDouble("DM"), Manager.GetDouble("gamma-y"));
  
  char* CommentLine = new char [256];
  if (SzSymmetryFlag == false)
    {
      sprintf (CommentLine, "kx ky");
    }
  else
    {
      sprintf (CommentLine, "kx ky szp");
    }

  char* EigenvalueOutputFile = new char [512];
  if (ConserveSz == false)
    sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);
  else
    sprintf(EigenvalueOutputFile, "%s_%s_sz_%d.dat", FilePrefix, FileParameterString, szSector);
 
//   Abstract2DTightBindingModel* TightBindingModel;
//   if (Manager.GetBoolean("singleparticle-spectrum") == true)
//     {
//       bool ExportOneBody = false;
//       if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
// 	ExportOneBody = true;
//       TightBindingModel = new TightBindingModelKitaevHeisenbergHoneycombLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), 
// 										 Manager.GetDouble("anisotropic-phase"), Architecture.GetArchitecture(), ExportOneBody);
//       TightBindingModel->TestRealSpaceTightBindingHamiltonian();
//       cout << TightBindingModel->GetRealSpaceTightBindingNonHermitianHamiltonian() << endl;
//       TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
//       if (ExportOneBody == true)
// 	{
// 	  char* BandStructureOutputFile = new char [512];
// 	  if (Manager.GetString("export-onebodyname") != 0)
// 	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
// 	  else
// 	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
// 	  if (Manager.GetBoolean("export-onebody") == true)
// 	    {
// 	      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
// 	    }
// 	  else
// 	    {
// 	      TightBindingModel->WriteBandStructureASCII(BandStructureOutputFile);
// 	    }
// 	  delete[] BandStructureOutputFile;
// 	}	  
//       return 0;
//     }

//   if (Manager.GetString("import-onebody") == 0)
//     {
//       TightBindingModel = new TightBindingModelKitaevHeisenbergHoneycombLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), 
// 										 Manager.GetDouble("anisotropic-phase"), Architecture.GetArchitecture(), true);
//       char* BandStructureOutputFile = new char [1024];
//       sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
//       TightBindingModel->WriteBandStructure(BandStructureOutputFile);
//     }
//   else
//     {
//       TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebodydown")); 
//     }

  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
//   if (Manager.GetDouble("u-potential") != 0.0)
//     {
//       double UPotential = Manager.GetDouble("u-potential");
//       for (int i = 0; i < NbrSites; ++i)
// 	{
// 	  DensityDensityInteractionupdown.SetMatrixElement(i, i, UPotential);
// 	}
//     }
  RealSymmetricMatrix SxSxInteraction(NbrSites, true);
  RealSymmetricMatrix SySyInteraction(NbrSites, true);
  RealSymmetricMatrix SzSzInteraction(NbrSites, true);
  RealAntisymmetricMatrix SxSyInteraction(NbrSites, true);
  
  
  
  
/*  
  int* NbrInteractingOrbitals;
  int** InteractingOrbitalsOrbitalIndices;
  int** InteractingOrbitalsSpatialIndices;
  double*** InteractingOrbitalsPotentialsX;
  double*** InteractingOrbitalsPotentialsY;
  double*** InteractingOrbitalsPotentialsZ;*/
//   
//   double* TmpJx = new double[2];
//   double* TmpJy = new double[2];
//   double* TmpJz = new double[2];
//   
//   TmpJx[0] = Manager.GetDouble("jx");
//   TmpJx[1] = Manager.GetDouble("NNNjx");
//   
//   TmpJy[0] = Manager.GetDouble("jy");
//   TmpJy[1] = Manager.GetDouble("NNNjy");
//   
//   TmpJz[0] = Manager.GetDouble("jz");
//   TmpJz[1] = Manager.GetDouble("NNNjz");
//   
//   KagomeLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
// 							   InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentialsX, Manager.GetDouble("jx"), Manager.GetDouble("NNNjx"));
//   KagomeLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
// 							   InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentialsY, Manager.GetDouble("jy"), Manager.GetDouble("NNNjy"));
//   KagomeLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices,
// 							   InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentialsZ, Manager.GetDouble("jz"), Manager.GetDouble("NNNjz"));
  
  
  int p;
  int q;
  int p1;
  int q1;
  
  for (int i = 0; i < NbrSitesX; ++i)
    {
      for (int j = 0; j < NbrSitesY; ++j)
	{
	  GetRealSpaceIndex(i, j, p, q, OffsetReal);
	  
	  SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),  TmpJx);
	  SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),  TmpJy);
	  SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),  TmpJz);	  
	  SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),  TmpDM);
	  
	  
	  SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpJx * cos(0.5 * GammaY));
	  SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpJy * cos(0.5 * GammaY));
	  SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpJz);
	  SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  -TmpDM + 0.5 * (TmpJx + TmpJy) * sin(0.5 * GammaY));
	  
	  SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpJx * cos(0.5 * GammaY));
	  SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpJy * cos(0.5 * GammaY));
	  SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpJz);
	  SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p, q, 2, NbrSitesX, NbrSitesY),  TmpDM + 0.5 * (TmpJx + TmpJy) * sin(0.5 * GammaY));
	  
	  GetRealSpaceIndex(i, j + 1, p1, q1, OffsetReal);
	  SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpJxDown * cos(0.5 * GammaY));
	  SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpJyDown * cos(0.5 * GammaY));
	  SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpJzDown);
	  SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpDM + 0.5 * (TmpJx + TmpJy) * sin(0.5 * GammaY));
	  
	  if ((CylinderFlag == false) || (i < NbrSitesX - 1))
	  {
	    GetRealSpaceIndex(i + 1, j, p1, q1, OffsetReal);
	    
	    SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpJxDown);
	    SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpJyDown);
	    SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpJzDown);
	    SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  -TmpDM);
	  
	    
	    GetRealSpaceIndex(i + 1, j - 1, p1, q1, OffsetReal);
	    SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpJxDown * cos(-0.5 * GammaY));
	    SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpJyDown* cos(-0.5 * GammaY));
	    SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpJzDown);
	    SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpDM + 0.5 * (TmpJx + TmpJy) * sin(-0.5 * GammaY));
	  }
	  
	  // NNN interaction terms
	  if ((TmpNNNJx != 0.0) || (TmpNNNJy != 0.0) || (TmpNNNJz != 0.0))
	  {
// 	    cout << TmpNNNJx << " " << TmpNNNJy << " " << TmpNNNJz << endl;
	    SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 1, NbrSitesX, NbrSitesY),  TmpNNNJx * cos(-GammaY));
	    SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 1, NbrSitesX, NbrSitesY),  TmpNNNJy * cos(-GammaY));
	    SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 1, NbrSitesX, NbrSitesY),  TmpNNNJz);
	  
	    SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 2, NbrSitesX, NbrSitesY),  TmpNNNJx * cos(-0.5 * GammaY));
	    SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 2, NbrSitesX, NbrSitesY),  TmpNNNJy * cos(-0.5 * GammaY));
	    SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 2, NbrSitesX, NbrSitesY),  TmpNNNJz);
	    
	    if (((TmpNNNJx + TmpNNNJy) != 0.0) && (GammaY != 0.0))
	      {
		SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 1, NbrSitesX, NbrSitesY), 0.5 * (TmpNNNJx + TmpNNNJy) * sin( -GammaY));
		SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(i, j - 1, 2, NbrSitesX, NbrSitesY),  0.5 * (TmpNNNJx + TmpNNNJy) * sin(-0.5 * GammaY));
		
	      }
	    
	    
	    if ((CylinderFlag == false) || (i < NbrSitesX - 1))
	    {
	      GetRealSpaceIndex(i + 1, j - 1, p1, q1, OffsetReal);
	      SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 0, NbrSitesX, NbrSitesY),  TmpNNNJx * cos(GammaY));
	      SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 0, NbrSitesX, NbrSitesY),  TmpNNNJy * cos(GammaY));
	      SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 0, NbrSitesX, NbrSitesY),  TmpNNNJz);
	    
	      GetRealSpaceIndex(i + 1, j, p1, q1, OffsetReal);
	      SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpNNNJx * cos(-0.5 * GammaY));
	      SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpNNNJy * cos(-0.5 * GammaY));
	      SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY),  TmpNNNJz);
	      
	      SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpNNNJx * cos(0.5 * GammaY));
	      SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpNNNJy * cos(0.5 * GammaY));
	      SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpNNNJz);
	      
	      GetRealSpaceIndex(i + 1, j - 1, p1, q1, OffsetReal);
	      SxSxInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 2, NbrSitesX, NbrSitesY),  TmpNNNJx * cos(-0.5 * GammaY));
	      SySyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 2, NbrSitesX, NbrSitesY),  TmpNNNJy * cos(-0.5 * GammaY));
	      SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 2, NbrSitesX, NbrSitesY),  TmpNNNJz);
	      
	      if (((TmpNNNJx + TmpNNNJy) != 0.0) && (GammaY != 0.0))
	      {
		GetRealSpaceIndex(i + 1, j - 1, p1, q1, OffsetReal);
		SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 0, NbrSitesX, NbrSitesY), 0.5 * (TmpNNNJx + TmpNNNJy) * sin(-GammaY));
		GetRealSpaceIndex(i + 1, j, p1, q1, OffsetReal);
		SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 0, NbrSitesX, NbrSitesY), 0.5 * (TmpNNNJx + TmpNNNJy) * sin(-0.5 * GammaY));
		SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY), 0.5 * (TmpNNNJx + TmpNNNJy) * sin(0.5 * GammaY));
		GetRealSpaceIndex(i + 1, j - 1, p1, q1, OffsetReal);
		SxSyInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 2, NbrSitesX, NbrSitesY), 0.5 * (TmpNNNJx + TmpNNNJy) * sin(-0.5 * GammaY));
	      }
	    }	    
	  }
	 if (TmpJD != 0.0)
	  {
	     GetRealSpaceIndex(i + 1, j, p1, q1, OffsetReal);
	     SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 2, NbrSitesX, NbrSitesY),
					     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 2, NbrSitesX, NbrSitesY),  TmpJD);
	     if ((CylinderFlag == false) || (i < NbrSitesX - 1))
	     {
	      GetRealSpaceIndex(i + 1, j - 1, p1, q1, OffsetReal);
	      SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 0, NbrSitesX, NbrSitesY),
				     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1 , 0, NbrSitesX, NbrSitesY),  TmpJD);
		  
	      GetRealSpaceIndex(i, j + 1, p1, q1, OffsetReal);
	      SzSzInteraction.AddToMatrixElement(GetRealSpaceTightBindingLinearizedIndexSafe(p,q, 1, NbrSitesX, NbrSitesY),
				     GetRealSpaceTightBindingLinearizedIndexSafe(p1, q1, 1, NbrSitesX, NbrSitesY),  TmpJD);
	      }
	    }	  
	}
    }
  bool FirstRunFlag = true;

  bool TmpSzSymmetryFlag = SzSymmetryFlag;
  int MinXMomentum = 0;
  int MaxXMomentum = NbrSitesX - 1;
  
  if (Manager.GetInteger("only-kx") >= 0)
    {
      MaxXMomentum = Manager.GetInteger("only-kx");
      MinXMomentum = MaxXMomentum;
    }
  for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
    {
      int MinYMomentum = 0;
      int MaxYMomentum = NbrSitesY - 1;
      if (Manager.GetInteger("only-ky") >= 0)
	{
	  MaxYMomentum = Manager.GetInteger("only-ky");
	  MinYMomentum = MaxYMomentum;
	}
      for (int YMomentum = MinYMomentum; YMomentum <= MaxYMomentum; ++YMomentum)
	{
	    int SzParitySector = -1;
	    int MaxSzParitySector = 1;
	    if (szSector != 0)
	      TmpSzSymmetryFlag = false;
	    if (TmpSzSymmetryFlag == false)
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
		cout << "Sz parity sector = " << SzParitySector << endl;
		ParticleOnSphereWithSpin* Space = 0;
		AbstractHamiltonian* Hamiltonian = 0;
		if (TmpSzSymmetryFlag == false)
		  {
		    cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << endl;
		    if (GutzwillerFlag == false)
		    {
		      if (NoTranslationFlag)
			Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, NbrSites);
		      else
			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, NbrSitesX,
										   YMomentum, NbrSitesY);
		    }
		    else
		    {
		      if (NoTranslationFlag)
			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);		      
		      else
		      {
			if (ConserveSz == false)
			  Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, NbrSitesX,
													  YMomentum, NbrSitesY);
			else
			  Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, szSector, NbrSites, XMomentum, NbrSitesX,
													  YMomentum, NbrSitesY);
		      }
		    }
		  }
		else
		  {
		    bool MinusParitySector = true;
		    if (SzParitySector == 1)
		      MinusParitySector = false;
		    cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  SzParity = " << SzParitySector<< endl;
		    if (GutzwillerFlag == false)
		      Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, NbrSites, MinusParitySector);
		    else
		      Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites, MinusParitySector);
		  }
	      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		{
		  // 		      Space->PrintState(cout, i) << endl;
		}
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      
	      HermitianMatrix TightBindingMatrix (2 * NbrSites, true);
	      if (NoTranslationFlag)
	      {
		if ((Manager.GetDouble("DM") == 0.0) && (GammaY == 0))
		  Hamiltonian = new ParticleOnLatticeWithSpinFullRealSpaceHamiltonian(Space, NbrParticles, NbrSites,
												   TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, SxSxInteraction,
												  SySyInteraction, SzSzInteraction,
												  Architecture.GetArchitecture(), Memory);
		  else
		    Hamiltonian = new ParticleOnLatticeWithSpinFullRealSpaceHamiltonian(Space, NbrParticles, NbrSites,
												   TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, SxSxInteraction,
												  SySyInteraction, SzSzInteraction, SxSyInteraction,
												  Architecture.GetArchitecture(), Memory);
	      }
	      else
	      {
		if ((Manager.GetDouble("DM") == 0.0) && (GammaY == 0))
		  Hamiltonian = new ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites, XMomentum, NbrSitesX, YMomentum, NbrSitesY,TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, SxSxInteraction,
												  SySyInteraction, SzSzInteraction,
												  Architecture.GetArchitecture(), Memory);
		else
		  Hamiltonian = new ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites, XMomentum, NbrSitesX, YMomentum, NbrSitesY,TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, SxSxInteraction,
												  SySyInteraction, SzSzInteraction, SxSyInteraction,
												  Architecture.GetArchitecture(), Memory);
	      }
		  
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
	      char* TmpExtention = new char [512];
	      if (NoTranslationFlag)
		sprintf (TmpExtention, "");
	      else
	      {
		if (SzSymmetryFlag == false)
		  {
		    sprintf (TmpExtention, "_kx_%d_ky_%d", XMomentum, YMomentum);
		  }
		else
		  {
		    sprintf (TmpExtention, "_szp_%d_kx_%d_ky_%d", SzParitySector, XMomentum, YMomentum);
		  }
	      }

	
	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	      
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      FirstRunFlag = false;
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      cout << "------------------------------------" << endl;
	      delete Hamiltonian;
	      delete Space;
	      delete[] EigenstateOutputFile;
	      delete[] ContentPrefix;
	      if (NoTranslationFlag)
		return 0;
	      }
	    }
    }
  return 0;
}


int GetRealSpaceTightBindingLinearizedIndexSafe(int indexX, int indexY, int indexOrbital, int nbrSitesX, int nbrSitesY)
{
  if (indexX >= nbrSitesX)
    indexX -= nbrSitesX;
  if (indexX < 0)
    indexX += nbrSitesX;
  if (indexY >= nbrSitesY)
    indexY -= nbrSitesY;
  if (indexY < 0)
    indexY += nbrSitesY;
  return (indexOrbital + ((indexY  + indexX * nbrSitesY) * 3)); 
}


// compute the coordinate of a given point in the kagome lattice from their real space coordinates (trivial for a non-tilted lattice, same as GetRealSpaceIndex of Abstract2DTightBindingModel)
//
// i = cartesian coordinate in the x direction of the Bravais lattice
// j = cartesian coordinate in the y direction of the Bravais lattice
// p = reference on the first lattice index
// q = reference on the second lattice index
// OffsetReal = offset
void GetRealSpaceIndex (int i, int j, int& p, int& q, int offsetReal)
{
  p = i - offsetReal * j;
  q = j;
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

void KagomeLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double NNInteraction, double NNNInteraction)
{
  int nbrBands = 3;
  
  nbrInteractingOrbitals = new int[nbrBands];
  interactingOrbitalsOrbitalIndices = new int*[nbrBands];
  interactingOrbitalsSpatialIndices = new int*[nbrBands];
  interactingOrbitalsPotentials = new double*[nbrBands];
  nbrInteractingOrbitals[0] = 2; 
  nbrInteractingOrbitals[1] = 2;      
  nbrInteractingOrbitals[2] = 2;      
  if (NNNInteraction != 0.0)
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
  interactingOrbitalsPotentials[0][TmpIndex] = NNInteraction;
  ++TmpIndex;
  interactingOrbitalsOrbitalIndices[0][TmpIndex] = 2;
  interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 0;
  interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
  interactingOrbitalsPotentials[0][TmpIndex] = NNInteraction;
  ++TmpIndex;
  
  TmpIndex -= 2;
 // links starting from B
  interactingOrbitalsOrbitalIndices[1][TmpIndex] = 0;
  interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 1;
  interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
  interactingOrbitalsPotentials[1][TmpIndex] = NNInteraction;
  ++TmpIndex;
  interactingOrbitalsOrbitalIndices[1][TmpIndex] = 2;
  interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 0;
  interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
  interactingOrbitalsPotentials[1][TmpIndex] = NNInteraction;
  ++TmpIndex;
   
  TmpIndex -= 2;

  // links starting from C
  interactingOrbitalsOrbitalIndices[2][TmpIndex] = 0;
  interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = 0;
  interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
  interactingOrbitalsPotentials[2][TmpIndex] = NNInteraction;
  ++TmpIndex;
  interactingOrbitalsOrbitalIndices[2][TmpIndex] = 1;
  interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = -1;
  interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
  interactingOrbitalsPotentials[2][TmpIndex] = NNInteraction;
  ++TmpIndex;
  
  
  if (NNNInteraction != 0.0)
  {
      TmpIndex = 2;
      interactingOrbitalsOrbitalIndices[0][TmpIndex] = 1;
      interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = -1;
      interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 1;
      interactingOrbitalsPotentials[0][TmpIndex] = NNNInteraction; //c
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[0][TmpIndex] = 1;
      interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 0;
      interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
      interactingOrbitalsPotentials[0][TmpIndex] = NNNInteraction; //c
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[0][TmpIndex] = 2;
      interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = 1;
      interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
      interactingOrbitalsPotentials[0][TmpIndex] = NNNInteraction;
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[0][TmpIndex] = 2;
      interactingOrbitalsSpatialIndices[0][TmpIndex * 2] = -1;
      interactingOrbitalsSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
      interactingOrbitalsPotentials[0][TmpIndex] = NNNInteraction;
      
      TmpIndex = 2;
      interactingOrbitalsOrbitalIndices[1][TmpIndex] = 0;
      interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 1;
      interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
      interactingOrbitalsPotentials[1][TmpIndex] = NNNInteraction;
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[1][TmpIndex] = 0;
      interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 0;
      interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 1;
      interactingOrbitalsPotentials[1][TmpIndex] = NNNInteraction;
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[1][TmpIndex] = 2;
      interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 1;
      interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
      interactingOrbitalsPotentials[1][TmpIndex] = NNNInteraction; //c
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[1][TmpIndex] = 2;
      interactingOrbitalsSpatialIndices[1][TmpIndex * 2] = 0;
      interactingOrbitalsSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
      interactingOrbitalsPotentials[1][TmpIndex] = NNNInteraction; //c
      ++TmpIndex;
      
      TmpIndex = 2;
      interactingOrbitalsOrbitalIndices[2][TmpIndex] = 0;
      interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = 1;
      interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
      interactingOrbitalsPotentials[2][TmpIndex] = NNNInteraction; //c
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[2][TmpIndex] = 0;
      interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = -1;
      interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
      interactingOrbitalsPotentials[2][TmpIndex] = NNNInteraction; //c
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[2][TmpIndex] = 1;
      interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = 0;
      interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
      interactingOrbitalsPotentials[2][TmpIndex] = NNNInteraction;
      ++TmpIndex;
      interactingOrbitalsOrbitalIndices[2][TmpIndex] = 1;
      interactingOrbitalsSpatialIndices[2][TmpIndex * 2] = -1;
      interactingOrbitalsSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
      interactingOrbitalsPotentials[2][TmpIndex] = NNNInteraction;
      ++TmpIndex;
  }
  
}