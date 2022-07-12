#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian.h"
#include "Hamiltonian/ParticleOnSquareLatticeFourBandSimpleTIHamiltonian.h"
#include "Hamiltonian/ExplicitHamiltonian.h"

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


// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// mass = mass term of the simple TI model
void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double mass);

// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// mass = mass term of the simple TI model
ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double mass);


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHEQuantumSpinHall2DSimpleTI" , "0.01");
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site potential strength between different orbitals", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive on-site potential strength between opposite spins", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mass", "mass parameter of the simple TI model", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new BooleanOption ('\n', "four-bands", "perform the calculations within the full four band model");
  (*SystemGroup) += new BooleanOption ('\n', "project-fourbands", "project the hamiltonian from the four band model to the two band model");
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
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export hamiltonian in an ASCII column formatted file", 0);
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEQuantumSpinHall2DSimpleTI -h" << endl;
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
  int TotalNbrSites = NbrSitesX * NbrSitesY;
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
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetBoolean("four-bands") == true)
    {
      if (Manager.GetBoolean("project-fourbands") == true)
	{
	  sprintf (EigenvalueOutputFile, "%s_quantumspinhall3d_simpleti_projectedfourbands_n_%d_x_%d_y_%d_u_%f_v_%f_m_%f_gx_%f_gy_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("mass"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
      else
	{
	  sprintf (EigenvalueOutputFile, "%s_quantumspinhall3d_simpleti_fourbands_n_%d_x_%d_y_%d_u_%f_v_%f_m_%f_gx_%f_gy_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("mass"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
	}
    }
  else
    {
      sprintf (EigenvalueOutputFile, "%s_quantumspinhall_simpleti_n_%d_x_%d_y_%d_u_%f_v_%f_m_%f_gx_%f_gy_%f.dat", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("mass"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      ComputeSingleParticleSpectrum(EigenvalueOutputFile, NbrSitesX, NbrSitesY, Manager.GetDouble("mass"));
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
  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") " << endl;
	  if (Manager.GetBoolean("four-bands") == false)
	    {
	      ParticleOnSphereWithSpin* Space = 0;
	      if (Manager.GetBoolean("boson") == false)
		{
		  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      else
		{
		  Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      AbstractQHEHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,
										  Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("mass"),
										  Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), 		     
										  Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);
	      char* ContentPrefix = new char[256];
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      char* TmpExtention = new char [512];
	      sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
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
	  else
	    {
	      ParticleOnSphere* Space = 0;
	      cout << TotalNbrSites << endl;
	      if (Manager.GetBoolean("boson") == false)
		{
#ifdef __128_BIT_LONGLONG__
		  if (TotalNbrSites <= 15)
#else
		    if (TotalNbrSites <= 7)
#endif
		      {
			Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		      }
		    else
		      {
			Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		      }
		}
	      else
		{
		  Space = new BosonOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		}
	      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	      AbstractHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ParticleOnSquareLatticeFourBandSimpleTIHamiltonian((ParticleOnSphereWithSU4Spin*) Space, NbrParticles, NbrSitesX, NbrSitesY,
										   Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("mass"),
										   Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"),  		     
										   Manager.GetBoolean("flat-band"), Architecture.GetArchitecture(), Memory);

	      if (Manager.GetBoolean("project-fourbands") == true)
		{
		  ComplexMatrix* OneBodyBasis = ComputeSingleParticleTransformationMatrices(NbrSitesX, NbrSitesY, Manager.GetDouble("mass"));
		  ComplexMatrix NBodyTransformationMatrix;
		  if (Manager.GetBoolean("boson") == false)
		    {
		      NBodyTransformationMatrix = ((FermionOnSquareLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
		    }
		  else
		    {
		      NBodyTransformationMatrix = ((BosonOnSquareLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixOneBodyBasis(OneBodyBasis);
		    }
		  ComplexMatrix HRep (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
		  Hamiltonian->GetHamiltonian(HRep);
		  ComplexMatrix TransformedHRep = HRep.Conjugate(NBodyTransformationMatrix);
		  ParticleOnSphereWithSpin* TargetSpace;
		  if (Manager.GetBoolean("boson") == false)
		    {
		      TargetSpace = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  else
		    {
		      TargetSpace = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, i, j);
		    }
		  ComplexMatrix SU4SU2TransformationMatrix;
		  if (Manager.GetBoolean("boson") == false)
		    {
		      SU4SU2TransformationMatrix = ((FermionOnSquareLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixSU4ToSU2(TargetSpace, 2, 3);
		    }
		  else
		    {
		      SU4SU2TransformationMatrix = ((BosonOnSquareLatticeWithSU4SpinMomentumSpace*) Space)->TransformationMatrixSU4ToSU2(TargetSpace, 2, 3);
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
	      sprintf (ContentPrefix, "%d %d", i, j);
	      char* EigenstateOutputFile = new char [512];
	      char* TmpExtention = new char [512];
	      sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
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
  return 0;
}

// compute the single particle spectrum 
//
// outputFileName = name of the output file
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// mass = mass term of the simple TI model

void ComputeSingleParticleSpectrum(char* outputFileName, int nbrSitesX, int nbrSitesY, double mass)
{
  ofstream File;
  File.open(outputFileName);
  File << "# kx    ky     E_{-,1}    E_{-,2}    E_{+,1}    E_{+,2}" << endl;
  double MinEMinus = 0.0;
  double MaxEMinus = -10.0;
  double MinEPlus = 10.0;
  double MaxEPlus = 0.0;
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  HermitianMatrix TmpOneBodyHamiltonian(4, true);
	  Complex d2 (sin (((double) ky) * KyFactor), 0.0);
	  double d1 = sin (((double) kx) * KxFactor);
	  double d3 = (mass - cos (((double) kx) * KxFactor) - cos (((double) ky) * KyFactor));
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, d1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, d1);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, d2);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, d2);
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif   
	  if (MaxEMinus < TmpDiag(0, 0))
	    {
	      MaxEMinus = TmpDiag(0, 0);
	    }
	  if (MinEMinus > TmpDiag(0, 0))
	    {
	      MinEMinus = TmpDiag(0, 0);
	    }
	  if (MaxEPlus < TmpDiag(2, 2))
	    {
	      MaxEPlus = TmpDiag(2, 2);
	    }
	  if (MinEPlus > TmpDiag(2, 2))
	    {
	      MinEPlus = TmpDiag(2, 2);
	    }
	  File << (KxFactor * ((double) kx)) << " " << (KyFactor * ((double) ky)) << " " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) <<  " " << TmpDiag(2, 2) << " " << TmpDiag(3, 3) << endl;
	}
      File << endl;
    }
  cout << "Spread = " << (MaxEMinus - MinEMinus) << "  Gap = " <<  (MinEPlus - MaxEMinus) << "  Flatening = " << ((MaxEMinus - MinEMinus) / (MinEPlus - MaxEMinus)) << endl;
}

// compute the single particle tranformation matrices 
//
// nbrSitesX = number of sites in the x direction
// nbrSitesY = number of sites in the y direction
// mass = mass term of the simple TI model

ComplexMatrix* ComputeSingleParticleTransformationMatrices(int nbrSitesX, int nbrSitesY, double mass)
{
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[nbrSitesX * nbrSitesY];
  double KxFactor = 2.0 * M_PI / ((double) nbrSitesX);
  double KyFactor = 2.0 * M_PI / ((double) nbrSitesY);
  for (int kx = 0; kx < nbrSitesX; ++kx)
    {
      for (int ky = 0; ky < nbrSitesY; ++ky)
	{
	  int Index = (kx * nbrSitesY) + ky;
	  HermitianMatrix TmpOneBodyHamiltonian(4, true);
	  Complex d2 (sin (((double) ky) * KyFactor), 0.0);
	  double d1 = sin (((double) kx) * KxFactor);
	  double d3 = (mass - cos (((double) kx) * KxFactor) - cos (((double) ky) * KyFactor));
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, d1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, d1);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, d2);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, d2);
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
  return OneBodyBasis;
}

