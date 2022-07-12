#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"

#include "Hamiltonian/ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian.h"
//#include "Hamiltonian/ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian.h"
#include "Tools/FTITightBinding/TightBindingModelChern2TwoOrbitalSquareLattice.h"
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
  OptionManager Manager ("FCIChern2SquareLatticeTwoOrbitalModel" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site Hubbard potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive two-body nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive three-body nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "s-potential", "repulsive three-body next-to-nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "imag part of inter-orbital hopping between nearest neighbors along the x direction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "inter-orbital hopping between nearest neighbors along the y direction", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t3", "intra-orbital hopping between nearest neighbors", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model. The n-body interaction strength with largest n is set to unity");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
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
      cout << "see man page for option syntax or type FCIChern2SquareLatticeTwoOrbitalModel -h" << endl;
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

  char* FilePrefix = new char [512];
  int lenFilePrefix=0;
  if (Manager.GetBoolean("single-band") == false)
    {
      cout << "NotImplemented: only the single-band case is supported at the moment." <<endl;
      return 1;
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
          lenFilePrefix += sprintf (FilePrefix, "fermions_singleband_chern2_n_%d_x_%d_y_%d",  NbrParticles, NbrSiteX, NbrSiteY);
      else
      {
	  if (Manager.GetBoolean("three-body") == true)
              lenFilePrefix += sprintf (FilePrefix, "fermions_singleband_threebody_chern2_n_%d_x_%d_y_%d",  NbrParticles, NbrSiteX, NbrSiteY);
          else
              lenFilePrefix += sprintf (FilePrefix, "fermions_singleband_fourbody_chern2_n_%d_x_%d_y_%d",  NbrParticles, NbrSiteX, NbrSiteY);
      }
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetBoolean("flat-band") == false)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_w_%f", Manager.GetDouble("w-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetDouble("s-potential") != 0.0)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_s_%f", Manager.GetDouble("s-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) || Manager.GetBoolean("flat-band") == false)
          lenFilePrefix += sprintf (FilePrefix + lenFilePrefix, "_u_%f", Manager.GetDouble("u-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_v_%f", Manager.GetDouble("v-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_t1_%f_t2_%f_t3_%f_mus_%f_gx_%f_gy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
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
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
          ExportOneBody = true;
      TightBindingModelChern2TwoOrbitalSquareLattice TightBindingModel(NbrSiteX, NbrSiteY,
              Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), 
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

  TightBindingModelChern2TwoOrbitalSquareLattice TightBindingModel(NbrSiteX, NbrSiteY,
          Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), 
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
	      if ((NbrSiteX * NbrSiteY) <= 63)
		{
		  Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		}
	      else
		{
		  Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		}
 	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
 	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
 	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
		{
                  Hamiltonian = new ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, &TightBindingModel,
													Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
													Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
              else
		{
// 		  if (Manager.GetBoolean("three-body") == true)
//                   {
//                       Hamiltonian = new ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
//                               Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
//                               Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
//                               Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
//                   }
//                   else
//                   {
//                       Hamiltonian = new ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
//                               Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
//                               Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
//                               Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
//                   }
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

