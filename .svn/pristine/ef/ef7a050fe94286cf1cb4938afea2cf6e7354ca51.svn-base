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
#include "Hamiltonian/ParticleOnLatticeWithSpinFullRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelKitaevHeisenbergHoneycombLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

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
  OptionManager Manager ("HubbardExtendedKitaevHeisenbergHoneycombModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 2);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new BooleanOption  ('\n', "spin", "use spin model (equivalent to gutzwiller projection and NbrParticles = 2 * NbrSitesX * NbrSitesY)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site (Hubbard) potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "isotropic-t", "isotropic spin nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "anisotropic-t", "anisotropic nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "anisotropic-phase", "phase (in pi units) of the anisotropic nearest neighbor hopping amplitude", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "j1", "strength of the neareast neighbor Heisenberg interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "jK", "strength of the neareast neighbor Kitaev interaction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "jG", "strength of the neareast neighbor Heisenberg interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "j3", "strength of the neareast neighbor Kitaev interaction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "3-spin", "strength of the three-spin interaction term", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hx", "amplitude of magnetic field", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hy", "amplitude of magnetic field", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hz", "amplitude of magnetic field", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "do not use translation symmetry to compute the spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption ('\n', "import-onebody", "import information on the tight binding model from a file");
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

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  int NbrSites = 2 * NbrSitesX * NbrSitesY; 
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  bool SpinFlag = Manager.GetBoolean("spin");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  bool TranslationFlag = !Manager.GetBoolean("no-translation");
  int* SitesA = 0;
  int* SitesB = 0;
  int* BondTypes = 0;
  int NbrBonds = 0;
  
  if (SpinFlag == true)
  {
    GutzwillerFlag = true;
    NbrParticles = NbrSites;  
  }
  
  if (Manager.GetDouble("j3") != 0.0)
  {
    cout << "Error. j3 terms are not implemented." << endl;
    return -1;
  }
 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      if (SpinFlag == false)
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
	  sprintf (StatisticPrefix, "spin_kitaev_heisenberg");
	else
	  sprintf (StatisticPrefix, "spin_kitaev_heisenberg_szsym");
	  
      }
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg");
	  else
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller");
	}
      else
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_szsym");
	  else
	    sprintf (StatisticPrefix, "bosons_kitaev_heisenberg_gutzwiller_szsym");
	}
    }
    
  

  char* FilePrefix = new char [256];
  char* FileParameterString = new char [256];
  
  if (SpinFlag == false)
  {
    sprintf (FilePrefix, "%s_honeycomb_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    sprintf (FileParameterString, "t_%.6f_tK_%.6f_tKphi_%.6f_j1_%.6f_j2_%.6f_jG_%.6f_j3_%.6f_hx_%.6f_hy_%.6f_hz_%.6f", Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), Manager.GetDouble("anisotropic-phase"), Manager.GetDouble("j1"), Manager.GetDouble("jK"), Manager.GetDouble("jG"), Manager.GetDouble("j3"), Manager.GetDouble("hx"), Manager.GetDouble("hy"), Manager.GetDouble("hz"));
  }
  else
  {
    sprintf (FilePrefix, "%s_honeycomb_x_%d_y_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrSites);
    sprintf (FileParameterString, "j1_%.6f_jK_%.6f_jG_%.6f_j3_%.6f_hx_%.6f_hy_%.6f_hz_%.6f", Manager.GetDouble("j1"), Manager.GetDouble("jK"), Manager.GetDouble("jG"), Manager.GetDouble("j3"), Manager.GetDouble("hx"), Manager.GetDouble("hy"), Manager.GetDouble("hz"));
  }
  
  

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
  if (Manager.GetDouble("u-potential") == 0.0)
    sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);
  else
    sprintf(EigenvalueOutputFile, "%s_%s_u_%f.dat", FilePrefix, FileParameterString, Manager.GetDouble("u-potential"));

  Abstract2DTightBindingModel* TightBindingModel;
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModel = new TightBindingModelKitaevHeisenbergHoneycombLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), 
										 Manager.GetDouble("anisotropic-phase"), Manager.GetDouble("hx"), Manager.GetDouble("hy"), Manager.GetDouble("hz"), Architecture.GetArchitecture(), ExportOneBody);
//       TightBindingModel->TestRealSpaceTightBindingHamiltonian();
//       cout << TightBindingModel->GetRealSpaceTightBindingNonHermitianHamiltonian() << endl;
      TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
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

  if (Manager.GetString("import-onebody") == 0)
    {
      TightBindingModel = new TightBindingModelKitaevHeisenbergHoneycombLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("isotropic-t"), Manager.GetDouble("anisotropic-t"), 
										 Manager.GetDouble("anisotropic-phase"), Manager.GetDouble("hx"), Manager.GetDouble("hy"), Manager.GetDouble("hz"), Architecture.GetArchitecture(), true);
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebodydown")); 
    }

  RealSymmetricMatrix DensityDensityInteractionupup((TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) / 2, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown((TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) / 2, true);
  RealSymmetricMatrix DensityDensityInteractionupdown((TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) / 2, true);
  if (Manager.GetDouble("u-potential") != 0.0)
    {
      double UPotential = Manager.GetDouble("u-potential");
      for (int i = 0; i < NbrSites; ++i)
	{
	  DensityDensityInteractionupdown.SetMatrixElement(i, i, UPotential);
	}
    }
  RealSymmetricMatrix SxSxInteraction(NbrSites, true);
  RealSymmetricMatrix SySyInteraction(NbrSites, true);
  RealSymmetricMatrix SzSzInteraction(NbrSites, true);
  
  RealSymmetricMatrix SxSyInteraction(NbrSites, true);
  RealSymmetricMatrix SySzInteraction(NbrSites, true);
  RealSymmetricMatrix SxSzInteraction(NbrSites, true);
  
  // beware : the factor 4.0 is for compatibility reasons with the other Kitaev-Heisenberg code
//   double TmpJ1 = 4.0 * Manager.GetDouble("j1");
//   double TmpJ2 = 4.0 * Manager.GetDouble("j2");
  // Compatibility with old code: TmpJ1 = (TmpJ1OldCode - TmpJ2OldCode) * 4; TmpJK = TmpJ2OldCode * 8
  
  double TmpJ1 = Manager.GetDouble("j1");
  double TmpJK = Manager.GetDouble("jK");
  double TmpJG = Manager.GetDouble("jG");
  for (int i = 0; i < NbrSitesX; ++i)
    {
      for (int j = 0; j < NbrSitesY; ++j)
	{
	  SxSxInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,  TmpJ1);
	  SySyInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,  TmpJ1 + TmpJK);
	  SzSzInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,  TmpJ1);
// 	  SxSzInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1, 					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,  TmpJG);
					     
	  SxSxInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i - 1, j + 1, 2) >> 1,  TmpJ1);
	  SySyInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i - 1, j + 1, 2) >> 1,  TmpJ1);
	  SzSzInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i - 1, j + 1, 2) >> 1,  TmpJ1 + TmpJK);
// 	  SxSyInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 0) >> 1,					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i - 1, j + 1, 2) >> 1,  TmpJG);

	  SxSxInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i + 1, j, 0) >> 1,  TmpJ1 + TmpJK);
	  SySyInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i + 1, j, 0) >> 1,  TmpJ1);
	  SzSzInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1,
					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i + 1, j, 0) >> 1,  TmpJ1);
	  SySzInteraction.AddToMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, 2) >> 1, 					     TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(i + 1, j, 0) >> 1,  TmpJG);
	}
    }
  bool FirstRunFlag = true;

  int MinXMomentum = 0;
  int MaxXMomentum = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {
      MaxXMomentum = Manager.GetInteger("only-kx");
      MinXMomentum = MaxXMomentum;
    }
  if (TranslationFlag == false)
  {
    MinXMomentum = 0;
    MaxXMomentum = 0;
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
      if (TranslationFlag == false)
      {
	MinYMomentum = 0;
	MaxYMomentum = 0;
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
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, NbrSitesX,
										   YMomentum, NbrSitesY);
		  else
		  {
		    if (TranslationFlag)
		      Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, XMomentum, NbrSitesX,
													  YMomentum, NbrSitesY);
		    else
		      Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
		  }
		}
	      else
		{
		  bool MinusParitySector = true;
		  if (SzParitySector == 1)
		    MinusParitySector = false;
		  cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  SzParity = " << SzParitySector<< endl;
		  if (GutzwillerFlag == false)
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, NbrSites, MinusParitySector, 
											     XMomentum, NbrSitesX,
											     YMomentum, NbrSitesY);
		  else
		    Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, MinusParitySector, 
														    XMomentum, NbrSitesX,
														    YMomentum, NbrSitesY);
		}
	      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		{
		  // 		      Space->PrintState(cout, i) << endl;
		}
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      
	      HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
	      if (TranslationFlag)
		Hamiltonian = new ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites,XMomentum, NbrSitesX,
												  YMomentum, NbrSitesY, TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, SxSxInteraction,
												  SySyInteraction, SzSzInteraction, SxSyInteraction,
												  SySzInteraction, SxSzInteraction,
												  Architecture.GetArchitecture(), Memory);
	      else
		Hamiltonian = new ParticleOnLatticeWithSpinFullRealSpaceHamiltonian(Space, NbrParticles, NbrSites,TightBindingMatrix,
												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
												  DensityDensityInteractionupdown, SxSxInteraction,
												  SySyInteraction, SzSzInteraction, SxSyInteraction,
												  SySzInteraction, SxSzInteraction,
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
	      char* TmpExtention = new char [512];
	      if (TranslationFlag)
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
	      else
		sprintf (TmpExtention, "");
		
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
	    }
	}
    }  
  return 0;
}
