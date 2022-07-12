#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian.h"
#include "Tools/FTITightBinding/TightBindingModelTwoOrbitalSquareLattice.h"
#include "Tools/FTITightBinding/TightBindingModelCylinderTwoOrbitalSquareLattice.h"
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
  OptionManager Manager ("FCISquareLatticeTwoOrbitalModel" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site Hubbard potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "uab-potential", "repulsive on-site Hubbard potential strength - ab for bosons", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive two-body nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive three-body nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "s-potential", "repulsive three-body next-to-nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "imag part of inter-orbital hopping between nearest neighbors along the x direction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "inter-orbital hopping between nearest neighbors along the y direction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t3", "intra-orbital hopping between nearest neighbors", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "folding", "folding factor for the momenta along sigma_x and sigma_y", 1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "cylinder-geometry", "only apply periodiuc boundary conditions along the x-direction");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-entspectrum", "compute the real space entanglement spectrum for a filled band");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "entspec-minenergy", "minimum energy that defines the filled band", -7.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "entspec-maxenergy", "maximum energy that defines the filled band", -0.2);
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model. The n-body interaction strength with largest n is set to unity");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
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
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCISquareLatticeTwoOrbitalModel -h" << endl;
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
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* SystemName = new char [512];
  if (Manager.GetBoolean("cylinder-geometry") == false)
    {
      sprintf (SystemName, "twoorbitals");
    }
  else
    {
      sprintf (SystemName, "cylinder_twoorbitals");
    }

  char* FilePrefix = new char [512 +strlen(SystemName)];
  int lenFilePrefix=0;
  if (Manager.GetBoolean("single-band") == false)
    {
      cout << "NotImplemented: only the single-band case is supported at the moment." <<endl;
      return 1;
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
	lenFilePrefix += sprintf (FilePrefix, "%s_singleband_%s_n_%d_x_%d_y_%d", StatisticPrefix, SystemName, NbrParticles, NbrSiteX, NbrSiteY);
      else
      {
	  if (Manager.GetBoolean("three-body") == true)
              lenFilePrefix += sprintf (FilePrefix, "%s_singleband_threebody_%s_n_%d_x_%d_y_%d", StatisticPrefix, SystemName, NbrParticles, NbrSiteX, NbrSiteY);
          else
              lenFilePrefix += sprintf (FilePrefix, "%s_singleband_fourbody_%s_n_%d_x_%d_y_%d", StatisticPrefix, SystemName, NbrParticles, NbrSiteX, NbrSiteY);
      }
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetBoolean("flat-band") == false)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_w_%f", Manager.GetDouble("w-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetDouble("s-potential") != 0.0)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_s_%f", Manager.GetDouble("s-potential"));
      //     if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) || Manager.GetBoolean("flat-band") == false)
      lenFilePrefix += sprintf (FilePrefix + lenFilePrefix, "_u_%f_uab_%f", Manager.GetDouble("u-potential"), Manager.GetDouble("uab-potential") );
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_v_%f", Manager.GetDouble("v-potential"));
      if (Manager.GetInteger("folding") == 1)
	{
	  lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_t1_%f_t2_%f_t3_%f_mus_%f_gx_%f_gy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
      else
	{
	  lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_t1_%f_t2_%f_t3_%f_mus_%f_f_%ld_gx_%f_gy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), Manager.GetInteger("folding"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
    }
  if(Manager.GetBoolean("flat-band")==true)
    {
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_flatband");
    }
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      if (Manager.GetBoolean("cylinder-geometry") == false)
	{
	  bool ExportOneBody = false;
	  if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	    ExportOneBody = true;
	  TightBindingModelTwoOrbitalSquareLattice TightBindingModel(NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), 
								     Manager.GetInteger("folding"), Manager.GetDouble("mu-s"), 
								     Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
	  if (Manager.GetBoolean("singleparticle-chernnumber") == true)
	    cout << "Chern number = " << TightBindingModel.ComputeChernNumber(0) << endl;
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
	}
      else
	{
	  bool ExportOneBody = false;
	  if (Manager.GetBoolean("singleparticle-entspectrum") == true)
	    {
	      ExportOneBody = true;
	    }
	  TightBindingModelCylinderTwoOrbitalSquareLattice TightBindingModel(NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), 
									     Manager.GetInteger("folding"), Manager.GetDouble("mu-s"), 
									     Manager.GetDouble("gamma-x"), false, Architecture.GetArchitecture(), ExportOneBody);
	  TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
	  if (Manager.GetBoolean("singleparticle-entspectrum") == true)
	    {
	      double** DensityEigenvalues = 0;
	      int* NbrDensityEigenvalues = 0;
	      double MinEnergy = Manager.GetDouble("entspec-minenergy");
	      double MaxEnergy = Manager.GetDouble("entspec-maxenergy");
	      char* EntanglementSpectrumOutputFile;
	      if (Manager.GetString("eigenvalue-file") != 0)
		{
		  EntanglementSpectrumOutputFile = new char [strlen(Manager.GetString("eigenvalue-file")) + 16];
		  sprintf(EntanglementSpectrumOutputFile, "%s.full.ent", Manager.GetString("eigenvalue-file"));
		}
	      else
		{
		  EntanglementSpectrumOutputFile = new char [strlen(FilePrefix) + 16];
		  sprintf (EntanglementSpectrumOutputFile, "%s.full.ent", FilePrefix);
		}
	      TightBindingModel.ComputeOneBodyRealSpaceEntanglementSpectrum(EntanglementSpectrumOutputFile, MinEnergy, MaxEnergy, NbrSiteY / 2);
	    }
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

  TightBindingModelTwoOrbitalSquareLattice TightBindingModel(NbrSiteX, NbrSiteY,
          Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetInteger("folding"), Manager.GetDouble("mu-s"), 
          Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());

  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
 	  if (Manager.GetBoolean("single-band") == false)
 	    {
 	    }
 	  else
 	    {
	      ParticleOnSphere* Space = 0;
	      if (Manager.GetBoolean("boson") == false)
		{
		  if ((NbrSiteX * NbrSiteY) <= 63)
		    {
		      Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		    }
		  else
		    {
		      Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		    }
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
                  Hamiltonian = new ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, &TightBindingModel,
												  Manager.GetDouble("u-potential"), Manager.GetDouble("uab-potential"), Manager.GetDouble("v-potential"), 
												  Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
              else
		{
		  if (Manager.GetBoolean("three-body") == true)
		    {
                      Hamiltonian = new ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, &TightBindingModel,
													       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
													       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
                  else
		    {
                      Hamiltonian = new ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, &TightBindingModel,
													      Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
													      Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		}
	      
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
    }
  return 0;
}

