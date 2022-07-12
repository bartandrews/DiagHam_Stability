#include "Options/Options.h"

//#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
//#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
//#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
//#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"

#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"


#include "Hamiltonian/ParticleOnLatticeHalfContinuousHofstadterModelSingleBandHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeOFLNOrbitalTriangularLatticeTwoBandHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticePyrochloreSlabLatticeSingleBandFiveBodyHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHalfContinuousHofstadterModel.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"


#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"

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
  OptionManager Manager ("FCIHalfContinuousHofstadterModel" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "use-inversion", "only compute the momentum sectors that are not related by the inversion symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");

  (*SystemGroup) += new BooleanOption  ('\n', "two-bands", "use the two lowest energy bands", false);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);

//  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive next to nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");

  (*SystemGroup) += new SingleDoubleOption  ('\n', "laser", "strength of laser", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "flux", "flux thread", 1);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "cutOFF", "number of reciprocal lattice points", 20);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleDoubleOption ('\n', "band-flattening", "flattening factor applied to each band, each band is rescale with respect to its average value", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "band-shifting", "shift each band by this amount times the band index", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-dispersion", "use a model without dispersion and a contant gap of 10 between the two lowest bands");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest energy band");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIHalfContinuousHofstadterModel -h" << endl;
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
  int Flux = Manager.GetInteger("flux");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int BandIndex = 0;
  double LaserStrength = Manager.GetDouble("laser");
  
  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
    
    char* FilePrefix = new char [256];

    sprintf (FilePrefix, "%s_singleband_halfcontinuoushofstadtermodel_nq_%ld_n_%d_x_%d_y_%d", StatisticPrefix, Manager.GetInteger("cutOFF") , NbrParticles ,NbrSitesX, NbrSitesY);
    
/*    if (Manager.GetBoolean("two-bands") == false)
      {
	if (Manager.GetBoolean("three-body") == false)
	  { 
	    if (Manager.GetBoolean("four-body") == false)
	      { 
		sprintf (FilePrefix, "%s_singleband_oflnorbitaltriangularlattice_s_%ld_c_%d_nq_%ld_n_%d_x_%d_y_%d", StatisticPrefix, Manager.GetInteger("nbr-spin"), ChernNumber,Manager.GetInteger("cutOFF") , NbrParticles, NbrSitesX, NbrSitesY);
	      }
	    else
	      {
		sprintf (FilePrefix, "%s_singleband_fourbody_oflnorbitaltriangularlattice_s_%ld_c_%d_nq_%ld_n_%d_x_%d_y_%d", StatisticPrefix,  Manager.GetInteger("nbr-spin"), ChernNumber, Manager.GetInteger("cutOFF"), NbrParticles, NbrSitesX, NbrSitesY);
	      }
	  }
	else
	  {
	    sprintf (FilePrefix, "%s_singleband_threebody_oflnorbitaltriangularlattice_s_%ld_c_%d_nq_%ld_n_%d_x_%d_y_%d", StatisticPrefix,  Manager.GetInteger("nbr-spin"), ChernNumber, Manager.GetInteger("cutOFF"), NbrParticles, NbrSitesX, NbrSitesY);
	  }
      }
    else
      {
	if ((Manager.GetDouble("band-flattening") != 1.0) || (Manager.GetDouble("band-shifting") != 0.0))
	  {
	    sprintf (FilePrefix, "%s_twoband_flattening_%.6f_shifting_%.6f_oflnorbitaltriangularlattice_s_%ld_c_%d_nq_%ld_n_%d_x_%d_y_%d", StatisticPrefix, Manager.GetDouble("band-flattening"), Manager.GetDouble("band-shifting"), Manager.GetInteger("nbr-spin"), ChernNumber,Manager.GetInteger("cutOFF") , NbrParticles, NbrSitesX, NbrSitesY);
	  }
	else
	  {
	    sprintf (FilePrefix, "%s_twoband_oflnorbitaltriangularlattice_s_%ld_c_%d_nq_%ld_n_%d_x_%d_y_%d", StatisticPrefix, Manager.GetInteger("nbr-spin"), ChernNumber,Manager.GetInteger("cutOFF") , NbrParticles, NbrSitesX, NbrSitesY);
	  }
      }
*/

    char* FileParameterString = new char [256];
    sprintf (FileParameterString, "las_%g_phi_%d_gx_%g_gy_%g", LaserStrength, Flux, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    char* CommentLine = new char [256];
    sprintf (CommentLine, "eigenvalues\n# kx ky ");
    char* EigenvalueOutputFile = new char [512];
    if (Manager.GetString("eigenvalue-file")!=0)
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    else
      {
	sprintf (EigenvalueOutputFile, "%s_%s.dat",FilePrefix, FileParameterString);
/*
	if (Manager.GetBoolean("flat-band") == true)
	  { 
	    if (Manager.GetDouble("v-potential") == 0.0)
	      sprintf (EigenvalueOutputFile, "%s_%s.dat",FilePrefix, FileParameterString);
	    else
	      {
		sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString);
	      }
	  }
	else
	  {
	    if (Manager.GetDouble("v-potential") == 0.0)
	      sprintf (EigenvalueOutputFile, "%s_u_%g_%s.dat",FilePrefix, Manager.GetDouble("u-potential"), FileParameterString);
	  else
	    sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString);
	  }
*/
      }
    
  
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;
      
      TightBindingModelHalfContinuousHofstadterModel TightBindingModel(LaserStrength, Flux, NbrSitesX, NbrSitesY, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(),
								      Manager.GetInteger("cutOFF"),ExportOneBody);
      
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)      
	{
	  cout << "Chern number = " << TightBindingModel.ComputeChernNumber(0) << endl;
	}
      
      TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel.ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel.ComputeDirectBandGap(0);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
      if (ExportOneBody == true)
	{
	  char* BandStructureOutputFile = new char [512];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
/*	  if (Manager.GetBoolean("export-onebody") == true)
	    {
	      TightBindingModel.WriteBandStructure(BandStructureOutputFile);
	    }
	  else
	    {
	      TightBindingModel.WriteBandStructureASCII(BandStructureOutputFile);
	    }*/
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
  
  TightBindingModelHalfContinuousHofstadterModel* TightBindingModel;
  if (Manager.GetString("import-onebody") == 0)
    {
      TightBindingModel = new  TightBindingModelHalfContinuousHofstadterModel (LaserStrength, Flux, NbrSitesX, NbrSitesY, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), Manager.GetInteger("cutOFF"));
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new TightBindingModelHalfContinuousHofstadterModel(Manager.GetString("import-onebody")); 
    }

  if ((Manager.GetDouble("band-flattening") != 1.0) || (Manager.GetDouble("band-shifting") != 0.0))
    TightBindingModel->FlattenBands(Manager.GetDouble("band-flattening"), Manager.GetDouble("band-shifting"));


//  TightBindingModel2DAtomicLimitLattice  TightBindingModel(NbrSitesX, NbrSitesY,  1 , 0, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());

  int NbrMomentumSectors = 0;
  if (Manager.GetBoolean("use-inversion") == false)
    {
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      ++NbrMomentumSectors;
	    }
	}
    }
  else
    {
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      if ((i <= ((NbrSitesX - i) % NbrSitesX)) && (j <= ((NbrSitesY - j) % NbrSitesY)))
		{
		  ++NbrMomentumSectors;
		}
	    }
	}
    }
  int* KxMomentumSectors = new int [NbrMomentumSectors];
  int* KyMomentumSectors = new int [NbrMomentumSectors];
  NbrMomentumSectors = 0;
  if (Manager.GetBoolean("use-inversion") == false)
    {
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      KxMomentumSectors[NbrMomentumSectors] = i;
	      KyMomentumSectors[NbrMomentumSectors] = j;
	      ++NbrMomentumSectors;
	    }
	}
    }
  else
    {
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      if ((i <= ((NbrSitesX - i) % NbrSitesX)) && (j <= ((NbrSitesY - j) % NbrSitesY)))
		{
		  KxMomentumSectors[NbrMomentumSectors] = i;
		  KyMomentumSectors[NbrMomentumSectors] = j;
		  ++NbrMomentumSectors;
		}
	    }
	}
    }

  
  bool FirstRunFlag = true;
  for (int i = 0; i < NbrMomentumSectors; ++i)
    {
      cout << "(kx=" << KxMomentumSectors[i] << ",ky=" << KyMomentumSectors[i] << ") : " << endl;
      
      ParticleOnSphere* Space = 0;
      if (Manager.GetBoolean("two-bands") == false)
	{
	  if (Manager.GetBoolean("boson") == false)
	    {
	      if ((NbrSitesX * NbrSitesY) <= 63)
		{
		  Space = 0;
		  // Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space =0;
		  //Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	    }
	  else
	    {
	      Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomentumSectors[i], KyMomentumSectors[i]);
	    }
	}
      else
	{
	  Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomentumSectors[i], KyMomentumSectors[i]);
	}
	  
      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
      
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
      AbstractQHEHamiltonian* Hamiltonian = 0;
      if (Manager.GetBoolean("three-body") == false)
	{ 
	  if (Manager.GetBoolean("four-body") == false)
	    { 
	      if (Manager.GetBoolean("two-bands") == false)
		{
		  Hamiltonian = new ParticleOnLatticeHalfContinuousHofstadterModelSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetInteger("cutOFF") , Manager.GetDouble("u-potential"), TightBindingModel, Manager.GetBoolean("flat-band") , BandIndex, Architecture.GetArchitecture(), Memory);
		}
	      else
		{
//		  Hamiltonian = new ParticleOnLatticeOFLNOrbitalTriangularLatticeTwoBandHamiltonian( (ParticleOnSphereWithSpin*) Space, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetInteger("nbr-spin"), Manager.GetInteger("cutOFF") , Manager.GetDouble("u-potential"), TightBindingModel, Manager.GetBoolean("flat-band"), Manager.GetBoolean("no-dispersion") , Architecture.GetArchitecture(), Memory);
		}
	    }
	  else
	    {
	      //Hamiltonian = new ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), 0.0,													TightBindingModel, Manager.GetInteger("nbr-layers") - 1, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	    }
	}
      else
	{ 
	  //Hamiltonian = new ParticleOnLatticeNOrbitalSquareLatticeSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), 0.0,	TightBindingModel,0, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	}
      
      char* ContentPrefix = new char[256];
      sprintf (ContentPrefix, "%d %d", KxMomentumSectors[i], KyMomentumSectors[i]);
      char* EigenstateOutputFile = new char [512];
      if (Manager.GetString("eigenstate-file")!=0)
	sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), KxMomentumSectors[i], KyMomentumSectors[i]);
      else
	{
	  char* TmpExtention = new char [512];
	  sprintf (TmpExtention, "_kx_%d_ky_%d", KxMomentumSectors[i], KyMomentumSectors[i]);
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


