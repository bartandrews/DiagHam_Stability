#include "Options/Options.h"

#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU4SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian.h"
#include "Hamiltonian/ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian.h"
#include "Hamiltonian/ExplicitHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelFuKaneMeleLattice.h"

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


// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbit = amplitude of the spin orbit coupling
ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHopingDistortion111, double spinOrbit);


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FTI3DFuKaneMele" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kz", "only evalute a given y momentum sector (negative if all kz sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive on-site potential strength between opposite spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive nearest neighbor potential strength between opposite spins", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('d', "deltat-111", "distortion of the nearest neighbor coupling in the (111) direction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('l', "lambda-so", "spin orbit coupling (in t unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-z", "boundary condition twisting angle along z (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "four-bands", "perform the calculations within the full four band model");
  (*SystemGroup) += new BooleanOption ('\n', "project-fourbands", "project the hamiltonian from the four band model to the two band model");
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
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTI3DFuKaneMele -h" << endl;
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
  int TotalNbrSites = NbrSitesX * NbrSitesY * NbrSitesZ;
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
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky kz ");

  char* FilePrefix = new char [512];
  if (Manager.GetBoolean("four-bands") == true)
    {
      if (Manager.GetDouble("mu-s") == 0.0)
	sprintf (FilePrefix, "%s_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
      else
	sprintf (FilePrefix, "%s_quantumspinhall3d_fukanemele_fourbands_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f_mus_%f", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Manager.GetDouble("mu-s"));
    }
  else
    {
      if (Manager.GetDouble("mu-s") == 0.0)
	sprintf (FilePrefix, "%s_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"));
      else
	sprintf (FilePrefix, "%s_quantumspinhall3d_fukanemele_n_%d_x_%d_y_%d_z_%d_u_%f_v_%f_w_%f_dt111_%f_so_%f_gx_%f_gy_%f_gz_%f_mus_%f", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), Manager.GetDouble("mu-s"));
    }

  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
  else
    {
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);
    }

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {      
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModelFuKaneMeleLattice TightBindingModel(NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"),
							   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), 
							   Architecture.GetArchitecture(), ExportOneBody);
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

  TightBindingModelFuKaneMeleLattice TightBindingModel(NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"),
						       Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), 
						       Architecture.GetArchitecture());
						       
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  for (int k = MinKz; k <= MaxKz; ++k)
	    {
	      cout << "(kx=" << i << ",ky=" << j << ",kz=" << k << ") " << endl;
	      if (Manager.GetBoolean("four-bands") == false)
		{
		  ParticleOnSphereWithSpin* Space = 0;
		  if (Manager.GetBoolean("boson") == false)
		    {
		      Space = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		    }
		  else
		    {
		      Space = new BosonOnCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		    }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
		  AbstractQHEHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ,
										       Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
										       &TightBindingModel, 		     
										       Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, k);
		  char* EigenstateOutputFile = new char [512];
		  char* TmpExtention = new char [512];
		  sprintf (TmpExtention, "_kx_%d_ky_%d_kz_%d", i, j, k);
		  if (Manager.GetString("eigenstate-file")!=0)
		    sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d_kz_%d", Manager.GetString("eigenstate-file"), i, j, k);
		  else
		    {
		      char* TmpExtention = new char [512];
		      sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
		      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
		    }

		  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
		  FirstRunFlag = false;
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  cout << "------------------------------------" << endl;
		  delete Hamiltonian;
		  delete[] EigenstateOutputFile;
		  delete[] ContentPrefix;
		}
	      else
		{
		  ParticleOnSphere* Space = 0;
		  if (Manager.GetBoolean("boson") == false)
		    {
#ifdef __128_BIT_LONGLONG__
		      if (TotalNbrSites <= 15)
#else
			if (TotalNbrSites <= 7)
#endif
			  {
			    Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
			  }
			else
			  {
			    Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
			  }
		    }
		  else
		    {
		      Space = new BosonOnCubicLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
		    }
		  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
		  AbstractHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian((ParticleOnSphereWithSU4Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ,
											Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 1.0, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"),
											Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("gamma-z"), 		     
											Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);

		  if (Manager.GetBoolean("project-fourbands") == true)
		    {
		      ComplexMatrix* OneBodyBasis = ComputeSingleParticleTransformationMatrices(NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetDouble("deltat-111"), Manager.GetDouble("lambda-so"));
		      ComplexMatrix NBodyTransformationMatrix;
		      if (Manager.GetBoolean("boson") == false)
			{
			  NBodyTransformationMatrix = ((FermionOnCubicLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
			}
		      else
			{
			  NBodyTransformationMatrix = ((BosonOnCubicLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
			}
		      ComplexMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
		      Hamiltonian->GetHamiltonian(HRep);
		      ComplexMatrix TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
		      ParticleOnSphereWithSpin* TargetSpace;
		      if (Manager.GetBoolean("boson") == false)
			{
			  TargetSpace = new FermionOnCubicLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
			}
		      else
			{
			  TargetSpace = new BosonOnCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, NbrSitesZ, i, j, k);
			}
		      ComplexMatrix SU4SU2TransformationMatrix;
		      if (Manager.GetBoolean("boson") == false)
			{
			  SU4SU2TransformationMatrix = ((FermionOnCubicLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixSU4ToSU2(TargetSpace, 2, 3);
			}
		      else
			{
			  SU4SU2TransformationMatrix = ((BosonOnCubicLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixSU4ToSU2(TargetSpace, 0, 1);
			}
		      ComplexMatrix TransformedHRep2 = TransformedHRep.InvConjugate(SU4SU2TransformationMatrix);
		      if (Manager.GetDouble("u-potential") != 0.0)
			TransformedHRep2 /= Manager.GetDouble("u-potential");
		      
		      RealDiagonalMatrix TmpDiag;
		      HermitianMatrix HRep2(TransformedHRep2);
		      delete Hamiltonian;
		      delete Space;
		      delete[] OneBodyBasis;
		      Hamiltonian = new ExplicitHamiltonian(TargetSpace, &HRep2);
		      Space = TargetSpace;
		    }


		  char* ContentPrefix = new char[256];
		  sprintf (ContentPrefix, "%d %d %d", i, j, k);
		  char* EigenstateOutputFile = new char [512];
		  char* TmpExtention = new char [512];
		  sprintf (TmpExtention, "_kx_%d_ky_%d_kz_%d", i, j, k);
		  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);

		  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
		  FirstRunFlag = false;
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  cout << "------------------------------------" << endl;
		  delete Hamiltonian;
		  delete[] EigenstateOutputFile;
		  delete[] ContentPrefix;
		}
	    }
	}
    }
  return 0;
}

// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// nbrSitesZ = number of sites in the z direction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbit = amplitude of the spin orbit coupling

ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, int nbrSitesZ, double nnHopingDistortion111, double spinOrbit)
{
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[nbrSitesX * nbrSitesY * nbrSitesZ];
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  double KzFactor = 2.0 * M_PI / ((double) nbrSitesZ);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  for (int kz = 0; kz < nbrSitesZ; ++kz)
	    {
	      double TmpKx = ((double) kx) * KxFactor;
	      double TmpKy = ((double) ky) * KyFactor;
	      double TmpKz = ((double) kz) * KzFactor;
	      int Index = ((kx * nbrSitesY) + ky) * nbrSitesZ + kz;
	      HermitianMatrix TmpOneBodyHamiltonian(4, true);
// 	      Complex B1 = 1.0 + nnHopingDistortion111 + Phase(0.5 * (TmpKy + TmpKz))   + Phase(0.5 * (TmpKx + TmpKz))  + Phase(0.5 * (TmpKx + TmpKy)) ;
// 	      double d3 = spinOrbit * (sin (0.5 * (TmpKx + TmpKz))
// 				       - sin (0.5 * (TmpKx + TmpKy))
// 				       - sin (0.5 * (TmpKx - TmpKy))
// 				       + sin (0.5 * (TmpKx - TmpKz)));
// 	      double d4 = spinOrbit * (sin (0.5 * (TmpKx + TmpKz))
// 				       - sin (0.5 * (TmpKy + TmpKz))
// 				       - sin (0.5 * (TmpKy - TmpKz))
// 				       + sin (0.5 * (TmpKy - TmpKx)));
// 	      double d5 = spinOrbit * (sin (0.5 * (TmpKy + TmpKx))
// 				       - sin (0.5 * (TmpKx + TmpKz))
// 				       - sin (0.5 * (TmpKz - TmpKx))
// 				       + sin (0.5 * (TmpKz - TmpKy)));

	      Complex B1 = 1.0 + nnHopingDistortion111 + Phase(1.0 * TmpKx)   + Phase(1.0 * TmpKy)  + Phase(1.0 * TmpKz) ;
	      double d3 = spinOrbit * (sin (1.0 * TmpKy)
				       - sin (1.0 * TmpKz)
				       - sin (1.0 * (TmpKy - TmpKx))
				       + sin (1.0 * (TmpKz - TmpKx)));
	      double d4 = spinOrbit * (sin (1.0 * TmpKz)
				       - sin (1.0 * TmpKx)
				       - sin (1.0 * (TmpKz - TmpKy))
				       + sin (1.0 * (TmpKx - TmpKy)));
	      double d5 = spinOrbit * (sin (1.0 * TmpKx)
				       - sin (1.0 * TmpKy)
				       - sin (1.0 * (TmpKx - TmpKz))
				       + sin (1.0 * (TmpKy - TmpKz)));

	      Complex B2 = d3 + I() * d4;
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, -d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d5);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 2, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 3, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B2);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, -B2);
	      RealDiagonalMatrix TmpDiag;
	      OneBodyBasis[Index] = ComplexMatrix (4, 4);
	      OneBodyBasis[Index].SetToIdentity();
#ifdef __LAPACK__
	      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis[Index]);
#else
	      TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis[Index]);
#endif   
	    }
	}
    }
  return OneBodyBasis;
}

