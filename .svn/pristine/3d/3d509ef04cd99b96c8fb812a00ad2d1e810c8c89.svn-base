#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"
#include "HilbertSpace/BosonOnTorusShort.h"


#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "MainTask/GenericRealMainTask.h"

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
  OptionManager Manager ("FQHEHaldaneModelWannier" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed). This option has an effect only if a=0 was chosen, i.e. if we calculate in the exactly block diagonal limit", -1);
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "phi", "Haldane phase on nnn hopping", M_PI/3);
  (*SystemGroup) += new BooleanOption  ('\n', "phase-in-pi", "Haldane phase on nnn hopping given in multiples of pi");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "single-band", "project onto the lowest enregy band");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model. The n-body interaction strength with largest n is set to unity");
  (*SystemGroup) += new BooleanOption  ('\n', "old-aspect", "use aspect ratio Ly/Lx");
  (*SystemGroup) += new BooleanOption  ('\n', "gaugeB", "Takes into account the gauge on B sites");
  (*SystemGroup) += new BooleanOption  ('\n', "NoWannier", "No Wannier");
  (*SystemGroup) += new SingleDoubleOption  ('a', "aParam", "parameter multiplying the off-block diagonal elements. a=0 corresponds to the exactly block diagonal limit",1.0);
  (*SystemGroup) += new SingleDoubleOption  ('b', "bParam", "parameter doing the cross over between the block diagonal elements of FCI and those of FQHE. b=0 correspons to FQHE block diagonal elements", 1.0);
  (*SystemGroup) += new BooleanOption  ('\n', "twisted", "Connect to rectangular torus (default: rectangular torus)");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of the hamiltonian");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEHaldaneModelWannier -h" << endl;
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
  int lenFilePrefix=0;
  if(Manager.GetBoolean("flat-band") == false)
    {
      cout << "Not Implemented: only the flat-band case is supported at the moment." <<endl;
      return 1;
    }
  if (Manager.GetBoolean("single-band") == false)
    {
      cout << "Not Implemented: only the single-band case is supported at the moment." <<endl;
      return 1;
    }
  else
    {
      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
	lenFilePrefix += sprintf (FilePrefix, "%s_singleband_haldane_n_%d_x_%d_y_%d",  StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY);
      else
	{
	  if (Manager.GetBoolean("three-body") == true)
	    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_threebody_haldane_Wannier_a_%f_b_%f_n_%d_x_%d_y_%d",  StatisticPrefix, Manager.GetDouble("aParam"), Manager.GetDouble("bParam"), NbrParticles, NbrSiteX, NbrSiteY);
          else
	    lenFilePrefix += sprintf (FilePrefix, "%s_singleband_fourbody_haldane_n_%d_x_%d_y_%d",  StatisticPrefix, NbrParticles, NbrSiteX, NbrSiteY);
	}
      if (Manager.GetBoolean("gaugeB"))
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_gaugeB");
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetBoolean("flat-band") == false)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_w_%f", Manager.GetDouble("w-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) && Manager.GetDouble("s-potential") != 0.0)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_s_%f", Manager.GetDouble("s-potential"));
      if ((Manager.GetBoolean("three-body") == true || Manager.GetBoolean("four-body") == true) || Manager.GetBoolean("flat-band") == false)
	lenFilePrefix += sprintf (FilePrefix + lenFilePrefix, "_u_%f", Manager.GetDouble("u-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_v_%f", Manager.GetDouble("v-potential"));
      lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_t1_%f_t2_%f_t3_%f_phi_%f_gx_%f_gy_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), Manager.GetDouble("phi"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
      if (Manager.GetDouble("mu-s") != 0.0)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_mus_%f", Manager.GetDouble("mu-s"));
      if (Manager.GetBoolean("flat-band") == true)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_flatband");
      if (Manager.GetBoolean("old-aspect") == false)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_aspect");
      if (Manager.GetBoolean("twisted") == true)
	lenFilePrefix += sprintf(FilePrefix + lenFilePrefix, "_twist");
    }
  char* CommentLine = new char [256];
  if(Manager.GetDouble("aParam")!=0.0)
    {
      sprintf (CommentLine, "eigenvalues\n# ky ");
    }
  else
    {
      sprintf (CommentLine, "eigenvalues\n# x ky ");
    }
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      if (Manager.GetInteger("only-ky")>=0)
	{
	  sprintf (EigenvalueOutputFile, "%s_ky_%ld.dat", FilePrefix, Manager.GetInteger("only-ky"));
	}
      else
	sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSiteX, NbrSiteY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), HaldanePhi, Manager.GetDouble("mu-s"));
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
  bool FirstRunFlag = true;

  // If we keep off block diagonal terms
  if ((Manager.GetDouble("aParam")!=0.0) || ((Manager.GetDouble("aParam")==0.0)&&(Manager.GetDouble("bParam")==0.0)))
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "ky=" << j << endl;
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
		      //Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		      return 0;
		    }
		  else
		    {
		      //Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, i, j);
		      return 0;
		    }
		}
	      else
		{
		  Space = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
		Memory = Architecture.GetArchitecture()->GetLocalMemory();
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
		{
		  return 0;
		  // Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
		  // 									   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
		  // 									   Manager.GetDouble("t1"), Manager.GetDouble("t2"), HaldanePhi, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
		  // 									   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  if (Manager.GetBoolean("three-body") == true)
		    {
		      Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier(Space, NbrParticles, NbrSiteX, NbrSiteY, 
													   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
													   Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), HaldanePhi, Manager.GetDouble("mu-s"), 
													   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
													   Manager.GetBoolean("flat-band"), Manager.GetBoolean("gaugeB"), Manager.GetDouble("aParam"), Manager.GetDouble("bParam"), !Manager.GetBoolean("twisted"), !Manager.GetBoolean("old-aspect"), Manager.GetBoolean("NoWannier"), Architecture.GetArchitecture(), Memory);
		    }
		  else
		    {
		      return 0;
		      // Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
		      // 									       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
		      // 									       Manager.GetDouble("t1"), Manager.GetDouble("t2"), HaldanePhi, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
		      // 									       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		    }
		}

	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d", j);
	      char* EigenstateOutputFile = new char [512];
	      if (Manager.GetString("eigenstate-file")!=0)
		sprintf (EigenstateOutputFile, "%s_ky_%d", Manager.GetString("eigenstate-file"), j);
	      else
		sprintf (EigenstateOutputFile, "%s_ky_%d", FilePrefix, j);
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
  // else // if we are in the exactly block diagonal limit
  //   {
  //     for(int i=MinKx; i<= MaxKx; ++i)
  // 	{
  // 	  for (int j = MinKy; j <= MaxKy; ++j)
  // 	    {
  // 	      cout << "(kx=" << i << ",ky=" << j << ") : " << endl;
  // 	      if (Manager.GetBoolean("single-band") == false)
  // 		{
  // 		}
  // 	      else
  // 		{
  // 		  ParticleOnSphere* Space = 0;
  // 		  if (Manager.GetBoolean("boson") == false)
  // 		    {
  // 		      if ((NbrSiteX * NbrSiteY) <= 63)
  // 			{
  // 			  //Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
  // 			  return 0;
  // 			}
  // 		      else
  // 			{
  // 			  //Space = new FermionOnSquareLatticeMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, i, j);
  // 			  return 0;
  // 			}
  // 		    }
  // 		  else
  // 		    {
  // 		      Space = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, j, i);
  // 		      //Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, i, j);
  // 		      //Space = new BosonOnTorusShort(NbrParticles, NbrSiteX * NbrSiteY, i*NbrSiteY+j);	
  // 		    }
  // 		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
  // 		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
  // 		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  // 		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
  // 		  AbstractQHEHamiltonian* Hamiltonian = 0;
  // 		  if ((Manager.GetBoolean("three-body") == false) && (Manager.GetBoolean("four-body") == false))
  // 		    {
  // 		      return 0;
  // 		      // Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
  // 		      // 									   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
  // 		      // 									   Manager.GetDouble("t1"), Manager.GetDouble("t2"), HaldanePhi, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
  // 		      // 									   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
  // 		    }
  // 		  else
  // 		    {
  // 		      if (Manager.GetBoolean("three-body") == true)
  // 			{
  // 			  Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier(Space, NbrParticles, NbrSiteX, NbrSiteY, 
  // 													       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
  // 													       Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("t3"), HaldanePhi, Manager.GetDouble("mu-s"), 
  // 													       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
  // 													       Manager.GetBoolean("flat-band"), Manager.GetBoolean("gaugeB"), Manager.GetDouble("aParam"), Manager.GetDouble("bParam"), Manager.GetBoolean("NoWannier"), Architecture.GetArchitecture(), Memory);
  // 			}
  // 		      else
  // 			{
  // 			  return 0;
  // 			  // Hamiltonian = new ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian(Space, NbrParticles, NbrSiteX, NbrSiteY, 
  // 			  // 									       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("s-potential"),
  // 			  // 									       Manager.GetDouble("t1"), Manager.GetDouble("t2"), HaldanePhi, Manager.GetDouble("mu-s"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),
  // 			  // 									       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
  // 			}
  // 		    }

  // 		  char* ContentPrefix = new char[256];
  // 		  sprintf (ContentPrefix, "%d %d", i, j);
  // 		  char* EigenstateOutputFile = new char [512];
  // 		  if (Manager.GetString("eigenstate-file")!=0)
  // 		    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", Manager.GetString("eigenstate-file"), i, j);
  // 		  else
  // 		    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d", FilePrefix, i, j);
  // 		  // if(Manager.GetDouble("aParam")!=0)
  // 		  //   {
  // 		      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
  // 		      FirstRunFlag = false;
  // 		      MainTaskOperation TaskOperation (&Task);
  // 		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  // 		  //   }
  // 		  // else
  // 		  //   {
  // 		  //     GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile, true);
  // 		  //     FirstRunFlag = false;
  // 		  //     MainTaskOperation TaskOperation (&Task);
  // 		  //     TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  // 		  //   }
  // 		  cout << "------------------------------------" << endl;
  // 		  delete Hamiltonian;
  // 		  delete Space;
  // 		  delete[] EigenstateOutputFile;
  // 		  delete[] ContentPrefix;
  // 		}
  // 	    }
  // 	}
  //   }
  delete [] StatisticPrefix;
  delete [] FilePrefix;
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
