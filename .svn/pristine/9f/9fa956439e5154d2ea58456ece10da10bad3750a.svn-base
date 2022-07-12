#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpaceLong.h"

#include "Hamiltonian/ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeenericIntraUnitCellInteractionSingleBandThreeBodyHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"

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


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIAtomicLimitLatticeModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitas per unit cell", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive intra-orbital potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive inter-orbital (but same unit cell) potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four body interaction instead of a two body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "five-body", "use a five body interaction instead of a two body interaction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band that has to be partially filled", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice chemical potential on A site", -1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
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
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIAtomicLimitLatticeModel -h" << endl;
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
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int NbrOrbitals = Manager.GetInteger("nbr-orbitals");
  double* ChemicalPotentials = new double [NbrOrbitals];
  for (int i = 0; i < NbrOrbitals; ++i)
    ChemicalPotentials[i] = (double) i;

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

  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
    { 
      sprintf (FilePrefix, "%s_singleband_atomiclimit_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, (int) Manager.GetInteger("band-index"), NbrParticles, NbrSitesX, NbrSitesY);
    }
  else
    {
      if (Manager.GetBoolean("three-body") == true)
	{
	  sprintf (FilePrefix, "%s_singleband_threebody_atomiclimit_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	}
      else
	{
	  if (Manager.GetBoolean("four-body") == true)
	    {
	      sprintf (FilePrefix, "%s_singleband_fourbody_atomiclimit_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	    }
	  else
	    {
	      sprintf (FilePrefix, "%s_singleband_fivebody_atomiclimit_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
	    }
	}
    }
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "norb_%d", NbrOrbitals);
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetBoolean("flat-band") == true)
	{
	  if (Manager.GetDouble("v-potential") == 0.0)
	    {
	      sprintf (EigenvalueOutputFile, "%s_%s_gx_%g_gy_%g.dat",FilePrefix, FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	  else
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	}
      else
	{
	  if (Manager.GetDouble("v-potential") == 0.0)
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	  else
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%g_v_%g_%s_gx_%g_gy_%g.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	    }
	}
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {      
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModel2DAtomicLimitLattice TightBindingModel(NbrSitesX, NbrSitesY, NbrOrbitals, ChemicalPotentials, 
							      Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
      TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel.ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel.ComputeDirectBandGap(0);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
      if (ExportOneBody == true)
	{
	  char* BandStructureOutputFile = new char [1024];
	  if (Manager.GetString("export-onebodyname") != 0)
	    strcpy(BandStructureOutputFile, Manager.GetString("export-onebodyname"));
	  else
	    sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
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
    
  Abstract2DTightBindingModel* TightBindingModel;
  if (Manager.GetString("import-onebody") == 0)
    {
      TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, NbrOrbitals, ChemicalPotentials,  
								     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());
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
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;

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
	      if (((NbrSitesX * NbrSitesY) + ((NbrParticles * (NbrParticles - 1)) >> 1)) <= 63)
		{
		  Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space = new BosonOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	    }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (Manager.GetBoolean("five-body") == false))
	    { 
	      Hamiltonian = new ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
												      Manager.GetInteger("band-index"), Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
												      TightBindingModel, Manager.GetBoolean("flat-band"), 
												      Architecture.GetArchitecture(), Memory);
	    }
	  else
	    { 
	      if (Manager.GetBoolean("three-body") == true)
		{
		  Hamiltonian = 0;
//  		  Hamiltonian = new ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
// 											       Manager.GetDouble("u-potential"), 0.0, TightBindingModel, 		     
// 											       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  if (Manager.GetBoolean("four-body") == true)
		    {
		      Hamiltonian = 0;
// 		      Hamiltonian = new ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
// 												    Manager.GetDouble("u-potential"), &TightBindingModel, 		     
// 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      Hamiltonian = 0;
// 		      Hamiltonian = new ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandFiveBodyHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
// 												    Manager.GetDouble("u-potential"),  &TightBindingModel, 		     
// 												    Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		  
		}
	    }

	  double EnergyShift = -10.0;
	  Hamiltonian->ShiftHamiltonian(EnergyShift);
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
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, EnergyShift,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
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
  delete TightBindingModel;
  return 0;
}


