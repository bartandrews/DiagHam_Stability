#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
#include "Tools/FTITightBinding/TightBindingModelAlternativeKagomeLattice.h"
#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

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
  OptionManager Manager ("FCIAlternativeKagomeLatticeModel" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "offset", "second coordinate in momentum space of the second spanning vector of the reciprocal lattice (0 if lattice is untilted or if Ny = 1)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength (or on-site density-density potential for bosons)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive next nearest neighbor potential strength (or nearest neighbor density-density potential for bosons)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive next next nearest neighbor potential strength (or next nearest neighbor density-density potential for bosons)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "3bw-potential", "repulsive three-body nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "3bs-potential", "repulsive three-body next-to-nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", -0.3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l1", "Rashba coupling between nearest neighbor sites", 0.28);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l2", "Rashba coupling between next nearest neighbor sites", 0.2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice chemical potential on A1 site", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band that has to be partially filled, should be 0 (lower band), 1 or 2 (upper band)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest energy band");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model. The n-body interaction strength with largest n is set to unity");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption('\n', "shift", "shift energy by +1.0 to help convergence");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEAlternativeKagomeLatticeModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbr-sitex"); 
  int NbrSiteY = Manager.GetInteger("nbr-sitey"); 
  int BandIndex = Manager.GetInteger("band-index");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int nx1 = Manager.GetInteger("nx1");
  int ny1 = Manager.GetInteger("ny1");
  int nx2 = Manager.GetInteger("nx2");
  int ny2 = Manager.GetInteger("ny2");
  int offset = Manager.GetInteger("offset");
  bool TiltedFlag = true;
  if ( ((nx1 == 0) && (ny1 == 0)) || ((nx2 == 0) && (ny2 == 0)) )
    TiltedFlag = false;
  else
    {
      if ((nx1*ny2 - nx2*ny1) != NbrSiteX * NbrSiteY)
	{
	  cout << "Boundary conditions define a lattice that has a number of sites different from NbrSiteX * NbrSiteY - should have (nx1*ny2 - nx2*ny1) = NbrSiteX * NbrSiteY " << endl;
	  return 0;
	}
      if ((((offset*ny2 - ny1) % NbrSiteX) != 0) || (((nx1 - offset*nx2) % NbrSiteX) != 0))
	{
	  cout << "Tilted lattice not properly defined. Should have ((offset*ny2 - ny1) % NbrSiteX) = 0 and ((nx1 - offset*nx2) % NbrSiteX = 0) to verify momentum conservation" << endl;
	  return 0;
	}
      else
	cout << "Using tilted boundary conditions" << endl;
    }
  
  
  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    sprintf (StatisticPrefix, "fermions");
  else
    sprintf (StatisticPrefix, "bosons");
  
  char* FilePrefix = new char[512];
  int lenFilePrefix=0;
  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false) && (TiltedFlag == false))
    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_kagome_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, BandIndex, NbrParticles, NbrSiteX, NbrSiteY);
  else
    {
      if (TiltedFlag == true)
	lenFilePrefix += sprintf (FilePrefix, "%s_singleband_kagomelatticetilted_n_%d_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d", StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2);
      else
	{
        if (Manager.GetBoolean("three-body") == true)
	  lenFilePrefix += sprintf (FilePrefix, "%s_singleband_threebody_kagome_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, BandIndex, NbrParticles, NbrSiteX, NbrSiteY);
        else
	  lenFilePrefix += sprintf (FilePrefix, "%s_singleband_fourbody_kagome_band_%d_n_%d_x_%d_y_%d", StatisticPrefix, BandIndex, NbrParticles, NbrSiteX, NbrSiteY);
	}
    }
 
  if (Manager.GetDouble("w-potential") == 0.0)
    {
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_u_%f_v_%f", Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"));
    }
  else
    {
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_u_%f_v_%f_w_%f", Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"));
    }
  if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true))
    lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_3bw_%f", Manager.GetDouble("3bw-potential"));
  if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true))
    lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_3bs_%f", Manager.GetDouble("3bs-potential"));
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t1_%g_t2_%g_l1_%g_l2_%g", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"));
  lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_%s_gx_%f_gy_%f", FileParameterString,
			   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file") != 0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    sprintf(EigenvalueOutputFile, "%s.dat", FilePrefix);
  
  Abstract2DTightBindingModel *TightBindingModel;
  
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;
      if (TiltedFlag == false)
	{
	  TightBindingModel = new TightBindingModelAlternativeKagomeLattice (NbrSiteX, NbrSiteY,
									     Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
									     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
	}
      else
	{
	  TightBindingModel = new TightBindingModelKagomeLattice (NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2, offset, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
								  Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
	}
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	cout << "Chern number = " << TightBindingModel->ComputeChernNumber(BandIndex) << endl;
      TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel->ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel->ComputeDirectBandGap(0);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;
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
  int MaxKx = NbrSiteX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSiteY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
  
  if (Manager.GetString("import-onebody") == 0)
    {
      if (TiltedFlag == false)
	{
	  TightBindingModel = new TightBindingModelAlternativeKagomeLattice (NbrSiteX, NbrSiteY,
									     Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
									     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());
	}
      else
	{
	  TightBindingModel = new TightBindingModelAlternativeKagomeLattice (NbrSiteX, NbrSiteY, nx1, ny1, nx2, ny2, offset, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
									     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());
	}
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
	      if ((NbrSiteX * NbrSiteY) <= 63)
		Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
	      else
		Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, i, j);
            }
	  else
            {
	      Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
            }
	  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  AbstractQHEHamiltonian* Hamiltonian = 0;
	  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
            {
	      Hamiltonian = new ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, TightBindingModel,
											       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
											       Manager.GetDouble("w-potential"), 
											       BandIndex, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
            }
	  else
            {
	      if (Manager.GetBoolean("three-body") == true)
                {
		  Hamiltonian = new ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, TightBindingModel,
													    Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("3bw-potential"), Manager.GetDouble("3bs-potential"),
													    BandIndex, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
                }
	      else
                {
		  cout << "NotImplemented: four-body interaction is not supported at the moment." <<endl;
		  return 1;
                }
            }
	  
	  if (Manager.GetBoolean("shift"))
	    Hamiltonian->ShiftHamiltonian(1.0);
	  
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  if (Manager.GetString("eigenstate-file")!=0)
	    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
	  else
	    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", FilePrefix, i, j);
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
  return 0;
}
