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

#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceEtaPairingJ2Hamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian.h"


#include "Tools/FTITightBinding/TightBindingModelSimpleSquareLattice.h"
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
  OptionManager Manager ("HubbardSquareLatticeModelJ2S2" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz", "twice the spin Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1, 0 if both sectors have to be computed", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "discard-momentum", "do not use the momentum as a good quantum number");
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
      cout << "see man page for option syntax or type HubbardSquareLatticeModelJ2S2 -h" << endl;
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
  int TotalSz = Manager.GetInteger("sz");
  int NbrSites = NbrSitesX * NbrSitesY; 
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  int* SitesA = 0;
  int* SitesB = 0;
  int* BondTypes = 0;
  int NbrBonds = 0;
  
  
 
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [64];
  if (Manager.GetBoolean("boson") == false)
    {
      if (SzSymmetryFlag == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "fermions_hubbard");
	  else
	    sprintf (StatisticPrefix, "fermions_hubbard_gutzwiller");
	}
      else
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "fermions_hubbard_szsym");
	  else
	    sprintf (StatisticPrefix, "fermions_hubbard_gutzwiller_szsym");
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_hubbard");
	  else
	    sprintf (StatisticPrefix, "bosons_hubbard_gutzwiller");
	}
      else
	{
	  if (GutzwillerFlag == false)
	    sprintf (StatisticPrefix, "bosons_hubbard_szsym");
	  else
	    sprintf (StatisticPrefix, "bosons_hubbard_gutzwiller_szsym");
	}
    }
    
  

  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_square_j2s2_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
  
  char* CommentLine = new char [256];
  if (Manager.GetBoolean("discard-momentum"))
    {
      if (SzSymmetryFlag == false)
	{
	  sprintf (CommentLine, "sz");
	}
      else
	{
	  sprintf (CommentLine, "sz szp");
	}
    }
  else
    {
      if (SzSymmetryFlag == false)
	{
	  sprintf (CommentLine, "kx ky sz");
	}
      else
	{
	  sprintf (CommentLine, "kx ky sz szp");
	}
    }

  char* EigenvalueOutputFile = new char [512];
  sprintf(EigenvalueOutputFile, "%s_sz_%d.dat", FilePrefix, TotalSz);

  Abstract2DTightBindingModel* TightBindingModel;
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModel = new TightBindingModelSimpleSquareLattice (NbrSitesX, NbrSitesY, 1.0, 0.0, 0.0, 0.0,
								    Architecture.GetArchitecture(), ExportOneBody);
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
      TightBindingModel = new TightBindingModelSimpleSquareLattice (NbrSitesX, NbrSitesY, 1.0, 0.0, 0.0, 0.0,
								    Architecture.GetArchitecture(), true);
      char* BandStructureOutputFile = new char [1024];
      sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
      TightBindingModel->WriteBandStructure(BandStructureOutputFile);
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }

  bool FirstRunFlag = true;

  if (Manager.GetBoolean("discard-momentum"))
    {
      int SzParitySector = -1;
      int MaxSzParitySector = 1;
      if (SzSymmetryFlag == false)
	{
	  SzParitySector = 1;
	}
      else
	{
	  if ((Manager.GetInteger("sz-parity") != 0) && (TotalSz == 0))
	    {
	      SzParitySector = Manager.GetInteger("sz-parity");
	      MaxSzParitySector = SzParitySector;
	    }
	}
      ParticleOnSphereWithSpin* Space = 0;
      AbstractHamiltonian* Hamiltonian = 0;
      if (SzSymmetryFlag == false)
	{
	  cout << "Sz = " << TotalSz <<  endl;
	  if (GutzwillerFlag == false)
	    Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, TotalSz, NbrSites);
	  else
	    Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (NbrParticles, TotalSz, NbrSites);
	}
      else
	{
	  bool MinusParitySector = true;
	  if (SzParitySector == 1)
	    MinusParitySector = false;
	  cout << "Sz = " << TotalSz << "  SzParity = " << SzParitySector<< endl;
	  if (GutzwillerFlag == false)
	    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, TotalSz, NbrSites, MinusParitySector);
	  else
	    Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (NbrParticles, TotalSz, NbrSites, MinusParitySector);
	}
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      
      char* ContentPrefix = new char[256];
      if (SzSymmetryFlag == false)
	{
	  sprintf (ContentPrefix, "%d", TotalSz);
	}
      else
	{
	  sprintf (ContentPrefix, "%d %d", TotalSz, SzParitySector);
	}
      char* EigenstateOutputFile;
      char* TmpExtention = new char [512];
      if (SzSymmetryFlag == false)
	{
	  sprintf (TmpExtention, "_sz_%d", TotalSz);
	}
      else
	{
	  sprintf (TmpExtention, "_szp_%d_sz_%d", SzParitySector, TotalSz);
	}
      Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceEtaPairingJ2Hamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY, 1.0,
										  Architecture.GetArchitecture(), Memory);
      HermitianMatrix J2HRep (Hamiltonian->GetHilbertSpaceDimension(), true);
      ComplexMatrix J2Basis (Hamiltonian->GetHilbertSpaceDimension(), Hamiltonian->GetHilbertSpaceDimension());
      RealDiagonalMatrix J2Eigenvalues = RealDiagonalMatrix(Hamiltonian->GetHilbertSpaceDimension());
      Hamiltonian->GetHamiltonian(J2HRep);
      J2HRep.LapackDiagonalize(J2Eigenvalues, J2Basis);
      J2Eigenvalues.SortMatrixUpOrder(J2Basis);
      delete Hamiltonian;

      int* J2Values = new int[J2Eigenvalues.GetNbrRow()];
      double Error = 1e-10;
      int MaxJ2Value = NbrSites;
      int* NbrJ2Values = new int[MaxJ2Value + 1];
      int* J2ValueFirstIndex = new int[MaxJ2Value + 1];
      for (int i = 0; i <= MaxJ2Value; ++i)
	{
	  NbrJ2Values[i] = 0;
	  J2ValueFirstIndex[i] = -1;
	}
      for (int i = 0; i < J2Eigenvalues.GetNbrRow(); ++i)
	{
	  if (J2Eigenvalues[i] < -Error)
	    {
	      cout << "error with the J^2 eigenvalue " << i << " : " << J2Eigenvalues[i] << endl;
	      J2Values[i] = -1;
	    }
	  else
	    {
	      if (J2Eigenvalues[i] < Error)
		{
		  J2Values[i] = 0;
		  if (J2ValueFirstIndex[J2Values[i]] < 0)
		    {
		      J2ValueFirstIndex[J2Values[i]] = i;
		    }
		  NbrJ2Values[J2Values[i]]++;
		}	      
	      else
		{
		  double Tmp = - 1.0 + sqrt(4.0 * J2Eigenvalues[i] + 1.0);
		  if (fabs(Tmp - round(Tmp)) > (Error * Tmp))
		    {
		      cout << "error with the J^2 eigenvalue " << i << " : " << J2Eigenvalues[i] << " " << Tmp << " " << fabs(Tmp - round(Tmp)) << endl;
		      J2Values[i] = -1;		      
		    }
		  else
		    {
		      J2Values[i] = (int) round(Tmp);
		      if (J2ValueFirstIndex[J2Values[i]] < 0)
			{
			  J2ValueFirstIndex[J2Values[i]] = i;
			}
		      NbrJ2Values[J2Values[i]]++;
		    }
		}
	    }
	}
      Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian(Space, NbrParticles, NbrSites, 1.0, true,
									Architecture.GetArchitecture(), Memory);
//       for (int i = 0; i < Hamiltonian->GetHilbertSpaceDimension(); ++i)
// 	{
// 	  cout << i << " " << J2Values[i] << endl;
// 	}
      for (int TmpJ2 = 0; TmpJ2 <= MaxJ2Value; ++TmpJ2)
	{
	  if (NbrJ2Values[TmpJ2] > 0)
	    {
	      cout << "Nbr states with J=";
	      if ((TmpJ2 & 1) == 0)
		cout << (TmpJ2 >> 1);
	      else
		cout << TmpJ2 << "/2";
	      cout << " : " << NbrJ2Values[TmpJ2] << endl;
	      
	      HermitianMatrix S2HRep (NbrJ2Values[TmpJ2], true);
 	      RealDiagonalMatrix S2Eigenvalues = RealDiagonalMatrix(NbrJ2Values[TmpJ2]);
	      ComplexVector TmpVector (Hamiltonian->GetHilbertSpaceDimension());
	      for (int i = 0; i < NbrJ2Values[TmpJ2]; ++i)
		{
		  Hamiltonian->Multiply(J2Basis[i + J2ValueFirstIndex[TmpJ2]], TmpVector);
		  Complex Tmp = J2Basis[i + J2ValueFirstIndex[TmpJ2]] * TmpVector;
		  S2HRep.SetMatrixElement(i, i, Tmp);
		  for (int j = i + 1; j < NbrJ2Values[TmpJ2]; ++j)
		    {
		      Tmp = J2Basis[j + J2ValueFirstIndex[TmpJ2]] * TmpVector;
		      S2HRep.SetMatrixElement(i, j, Tmp);		      
		    }
		}

	      S2HRep.LapackDiagonalize(S2Eigenvalues);
	      S2Eigenvalues.SortMatrixUpOrder();

	      int* S2Values = new int[S2Eigenvalues.GetNbrRow()];
	      int MaxS2Value = NbrParticles;
	      int* NbrS2Values = new int[MaxS2Value + 1];
	      for (int i = 0; i <= MaxS2Value; ++i)
		{
		  NbrS2Values[i] = 0;
		}
	      for (int i = 0; i < S2Eigenvalues.GetNbrRow(); ++i)
		{
		  if (S2Eigenvalues[i] < -Error)
		    {
		      cout << "error with the S^2 eigenvalue " << i << " : " << S2Eigenvalues[i] << endl;
		      S2Values[i] = -1;
		    }
		  else
		    {
		      if (S2Eigenvalues[i] < Error)
			{
			  S2Values[i] = 0;
			  NbrS2Values[S2Values[i]]++;
			}	      
		      else
			{
			  double Tmp = - 1.0 + sqrt(4.0 * S2Eigenvalues[i] + 1.0);
			  if (fabs(Tmp - round(Tmp)) > (Error * Tmp))
			    {
			      cout << "error with the S^2 eigenvalue " << i << " : " << S2Eigenvalues[i] << " " << Tmp << " " << fabs(Tmp - round(Tmp)) << endl;
			      S2Values[i] = -1;		      
			    }
			  else
			    {
			      S2Values[i] = (int) round(Tmp);
			      NbrS2Values[S2Values[i]]++;
			    }
			}
		    }
		}
	      for (int TmpS2 = 0; TmpS2 <= MaxS2Value; ++TmpS2)
		{
		  if (NbrS2Values[TmpS2] > 0)
		    {
		      cout << "Nbr states with S=";
		      if ((TmpS2 & 1) == 0)
			cout << (TmpS2 >> 1);
		      else
			cout << TmpS2 << "/2";
		      cout << " : " << NbrS2Values[TmpS2] << endl;
		    }
		}		      
	    }
	}
      
      delete Hamiltonian;

      cout << "------------------------------------" << endl;
      delete Space;
      delete[] ContentPrefix;
    }
  else
    {
 //      int MinXMomentum = 0;
//       int MaxXMomentum = NbrSitesX - 1;
//       if (Manager.GetInteger("only-kx") >= 0)
// 	{
// 	  MaxXMomentum = Manager.GetInteger("only-kx");
// 	  MinXMomentum = MaxXMomentum;
// 	}
//       for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
// 	{
// 	  int MinYMomentum = 0;
// 	  int MaxYMomentum = NbrSitesY - 1;
// 	  if (Manager.GetInteger("only-ky") >= 0)
// 	    {
// 	      MaxYMomentum = Manager.GetInteger("only-ky");
// 	      MinYMomentum = MaxYMomentum;
// 	    }
// 	  for (int YMomentum = MinYMomentum; YMomentum <= MaxYMomentum; ++YMomentum)
// 	    {
// 	      int SzParitySector = -1;
// 	      int MaxSzParitySector = 1;
// 	      if (SzSymmetryFlag == false)
// 		{
// 		  SzParitySector = 1;
// 		}
// 	      else
// 		{
// 		  if ((Manager.GetInteger("sz-parity") != 0) && (TotalSz == 0))
// 		    {
// 		      SzParitySector = Manager.GetInteger("sz-parity");
// 		      MaxSzParitySector = SzParitySector;
// 		    }
// 		}
// 	      for (; SzParitySector <= MaxSzParitySector; SzParitySector += 2)
// 		{
// 		  ParticleOnSphereWithSpin* Space = 0;
// 		  AbstractHamiltonian* Hamiltonian = 0;
// 		  if (SzSymmetryFlag == false)
// 		    {
// 		      cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  Sz = " << TotalSz <<  endl;
// 		      if (GutzwillerFlag == false)
// 			Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, XMomentum, NbrSitesX,
// 										       YMomentum, NbrSitesY);
// 		      else
// 			Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, XMomentum, NbrSitesX,
// 													      YMomentum, NbrSitesY);
// 		    }
// 		  else
// 		    {
// 		      bool MinusParitySector = true;
// 		      if (SzParitySector == 1)
// 			MinusParitySector = false;
// 		      cout << "Kx = " << XMomentum << "  Ky = " << YMomentum << "  Sz = " << TotalSz << "  SzParity = " << SzParitySector<< endl;
// 		      if (GutzwillerFlag == false)
// 			Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, MinusParitySector, 
// 												 XMomentum, NbrSitesX,
// 												 YMomentum, NbrSitesY);
// 		      else
// 			Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, TotalSz, NbrSites, MinusParitySector, 
// 															XMomentum, NbrSitesX,
// 															YMomentum, NbrSitesY);
// 		    }
// 		  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
// 		    Memory = Architecture.GetArchitecture()->GetLocalMemory();
// 		  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		  
// 		  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
// 		  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites,XMomentum, NbrSitesX,
// 												  YMomentum, NbrSitesY, TightBindingMatrix, TightBindingMatrix,
// 												  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
// 												  DensityDensityInteractionupdown, 
// 												  Architecture.GetArchitecture(), Memory);
		  
// 		  char* ContentPrefix = new char[256];
// 		  if (SzSymmetryFlag == false)
// 		    {
// 		      sprintf (ContentPrefix, "%d %d %d", XMomentum, YMomentum, TotalSz);
// 		    }
// 		  else
// 		    {
// 		      sprintf (ContentPrefix, "%d %d %d %d", XMomentum, YMomentum, TotalSz, SzParitySector);
// 		    }
// 		  char* EigenstateOutputFile;
// 		  char* TmpExtention = new char [512];
// 		  if (SzSymmetryFlag == false)
// 		    {
// 		      sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", XMomentum, YMomentum, TotalSz);
// 		    }
// 		  else
// 		    {
// 		      sprintf (TmpExtention, "_szp_%d_kx_%d_ky_%d_sz_%d", SzParitySector, XMomentum, YMomentum, TotalSz);
// 		    }
// 		  char* TmpExtentionSpectrum = new char [32];
// 		  sprintf (TmpExtentionSpectrum, "_sz_%d.dat", TotalSz);
// 		  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, TmpExtentionSpectrum, TmpExtention);
		  
// 		  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
// 		  FirstRunFlag = false;
// 		  MainTaskOperation TaskOperation (&Task);
// 		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
// 		  cout << "------------------------------------" << endl;
// 		  delete Hamiltonian;
// 		  delete Space;
// 		  delete[] EigenstateOutputFile;
// 		  delete[] ContentPrefix;
// 		}
// 	    }
//	}
    }  
  return 0;
}
