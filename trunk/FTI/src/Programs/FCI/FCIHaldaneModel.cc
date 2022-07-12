#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"

#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHaldaneHoneycombLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"

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


// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the x direction
// nnHopping = nearest neighbor hoping amplitude
// nnnHopping =  next nearest neighbor hoping amplitude
// phi =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSiteX, int nbrSiteY, double nnHopping, double nnnHopping, double nnnnHopping, double phi, double mus);


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIHaldaneModel" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive two-body nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive two-body nearest next neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive three-body nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "s-potential", "repulsive three-body next-to-nearest neighbor potential strength", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "three-body", "use a three-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new BooleanOption  ('\n', "four-body", "use a four-body interaction in addition to a two-body interaction");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t3", "next to next nearest neighbor hoping amplitude", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "phi", "Haldane phase on nnn hopping (multiples of pi)", 1.0/3.0);
  (*SystemGroup) += new BooleanOption  ('\n', "phase-in-pi", "Haldane phase on nnn hopping given in multiples of pi");  
  (*SystemGroup) += new BooleanOption  ('\n', "WannierHilbertSpace", "Wannier Hilbert Space");  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model. The n-body interaction strength with largest n is set to unity");
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
      cout << "see man page for option syntax or type FCIHaldaneModel -h" << endl;
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

  double HaldanePhi;

  if (Manager.GetBoolean("phase-in-pi"))
    HaldanePhi = M_PI*Manager.GetDouble("phi");
  else
    HaldanePhi = Manager.GetDouble("phi");


  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* FilePrefix = new char [512];
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t1_%f_t2_%f_t3_%f_phi_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("phi"));
  int lenFilePrefix=0;
  if (Manager.GetBoolean("single-band") == false)
    {
      cout << "NotImplemented: only the single-band case is supported at the moment." <<endl;
      return 1;
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
          lenFilePrefix += sprintf (FilePrefix, "%s_singleband_haldane_n_%d_x_%d_y_%d",  StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY);
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_threebody_haldane_n_%d_x_%d_y_%d",  StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY);
          else
	    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_fourbody_haldane_n_%d_x_%d_y_%d",  StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY);
	}
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetBoolean("flat-band") == false)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_w_%f", Manager.GetDouble("w-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetDouble("s-potential") != 0.0)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_s_%f", Manager.GetDouble("s-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) || Manager.GetBoolean("flat-band") == false)
          lenFilePrefix += sprintf (FilePrefix + lenFilePrefix, "_u_%f", Manager.GetDouble("u-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_v_%f", Manager.GetDouble("v-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_%s_gx_%f_gy_%f", FileParameterString, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      if (Manager.GetDouble("mu-s") != 0.0)
          lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_mus_%f", Manager.GetDouble("mu-s"));
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
       TightBindingModelHaldaneHoneycombLattice TightBindingModel(NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
								  HaldanePhi, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), ExportOneBody);
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
  Abstract2DTightBindingModel* TightBindingModel;
  if (Manager.GetString("import-onebody") == 0)
    {
      TightBindingModel = new TightBindingModelHaldaneHoneycombLattice (NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), 
									HaldanePhi, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture());
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }
   
  if(Manager.GetBoolean("WannierHilbertSpace"))
    {
      MinKx=0;
      MaxKx=0;
    }
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
		  if(!Manager.GetBoolean("WannierHilbertSpace"))
		    Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		  else
		    Space =  new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, j);
		}
 	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
 	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
 	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
              {
                  Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
										       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
										       TightBindingModel, 
										       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
              }
              else
              {
		  if (Manager.GetBoolean("three-body") == true)
                  {
                      Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
												    Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
												    TightBindingModel, Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
                  }
                  else
                  {
                      Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
												   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
												   TightBindingModel, 
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

// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the x direction
// nnHopping = nearest neighbor hoping amplitude
// nnnHopping =  next nearest neighbor hoping amplitude
// phase =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSiteX, int nbrSiteY, double nnHopping, double nnnHopping, double nnnnHopping, double phase, double mus)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_-    E_-" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -1000.0;
  double MinEPlus = 1000.0;
  double MaxEPlus = 0.0;
  for (int kx = 0; kx < nbrSiteX; ++kx)
  {
      double x=2*M_PI*((double)kx)/nbrSiteX;
      for (int ky = 0; ky < nbrSiteY; ++ky)
      {
          double y=2*M_PI*((double)ky)/nbrSiteY;

	  Complex B1 = - nnHopping * Complex(1 + cos(x+y) + cos(y), + sin(x+y) + sin(y));
	  Complex B2 = - nnnnHopping * Complex(2* cos(x) + cos(x+2*y),  sin(x+2*y));
          double d0 = - 2.0 * nnnHopping * cos(phase) * (cos(x) + cos(y) + cos(x+y));
          double d3 = - 2.0 * nnnHopping * sin(phase) * (sin(x) + sin(y) - sin(x+y)) + mus;

	  // My Convention
          // Complex B1 = - nnHopping * Complex(1 + cos(x) + cos(y), - sin(x) - sin(y));
	  // Complex B2 = - nnnnHopping * Complex(cos(x+y)+2*cos(x-y),-sin(x+y) );
          // double d0 = - 2.0 * nnnHopping * cos(phase) * (cos(x) + cos(y) + cos(x-y));
          // double d3 = - 2.0 * nnnHopping * sin(phase) * (sin(x) - sin(y) - sin(x-y)) + mus;


	  HermitianMatrix TmpOneBobyHamiltonian(2, true);
	  TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d0 + d3);
	  TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1+B2);
	  TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d0 - d3);
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	  TmpOneBobyHamiltonian.Diagonalize(TmpDiag);
#endif   
	  if (MaxEMinus < TmpDiag(0, 0))
	    {
	      MaxEMinus = TmpDiag(0, 0);
	    }
	  if (MinEMinus > TmpDiag(0, 0))
	    {
	      MinEMinus = TmpDiag(0, 0);
	    }
	  if (MaxEPlus < TmpDiag(1, 1))
	    {
	      MaxEPlus = TmpDiag(1, 1);
	    }
	  if (MinEPlus > TmpDiag(1, 1))
	    {
	      MinEPlus = TmpDiag(1, 1);
	    }
	  File << (2.0 * M_PI * ((double) kx) / ((double) nbrSiteX)) << " " << (2.0 * M_PI * ((double) ky) / ((double) nbrSiteY)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
	}
      File << endl;
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MinEPlus - MaxEMinus) /(MaxEMinus - MinEMinus) ) << endl;
}
