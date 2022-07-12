#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Tools/FTITightBinding/TightBindingModelOFLGenericLatticeWithSymmetry.h"

#include "Hamiltonian/ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian.h"
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
  OptionManager Manager ("FCIOFLGenericLatticeWithSymmetry" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);
  Manager += SystemGroup;
  TightBindingModelOFLGenericLatticeWithSymmetry::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-points1", "number of unit cells along the 'x' (G1)-direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-points2", "number of unit cells along the 'y' (G2)-direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength (all/equal spin)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "spin-anisotropy", "relative strength of inter-spin interactions", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "depth", "depth of the optical lattice", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "cut-off-mode", "cut-off for calculation: 0:Periodic, 1:Square, 2:circular ", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cut-off-momentum", "maximum absolute momentum to be considered (circular case)", 20.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "precision-threshold", "precision required in equality of symmetry related bands", 1e-6);  
  (*SystemGroup) += new MultipleIntegerOption  ('\n', "nmax", "cut-off for square or periodic case (units of enlarged unit cells)", ',', ',', "10,5");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleIntegerOption  ('b', "nbr-bands", "project onto the lowest #n energy bands (1 or 2)", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "no-dispersion", "use a model without dispersion and a contant gap of 10 between the two lowest bands");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-curvature", "compute the Berry curvature of the lowest band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the Chern number of the fully filled band (only available in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "output-nbr-bands", "number of bands to include in output files", 10);

  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption ('\n',"show-hamiltonian", "show Hamiltonian matrix, and exit");
  (*SystemGroup) += new BooleanOption  ('\n', "test-hermitian", "Check hermitian symmetry, and exit");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 1000);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption ('\n',"show-hamiltonian", "show Hamiltonian matrix, and exit");
  (*MiscGroup) += new BooleanOption  ('\n', "test-hermitian", "Check hermitian symmetry, and exit");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-points1"); 
  int NbrSitesY = Manager.GetInteger("nbr-points2");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int BandIndex = 0;
      
  double LatticeDepth = Manager.GetDouble("depth");
  
  char* FilePrefix = new char [512];
  char* EigenvalueOutputFile = new char [512];

  // set cut-off mode
  TightBindingModelOFLGenericLatticeWithSymmetry::CutOffModes MyCutOffMode = static_cast<TightBindingModelOFLGenericLatticeWithSymmetry::CutOffModes>( Manager.GetInteger("cut-off-mode"));

  // get cut-off for reciprocal space hopping model
  int NMax1, NMax2, l; 
  int *NMax = Manager.GetIntegers("nmax",l);
  if (l<2)
    {
     NMax1 = NMax[0];
     NMax2 = NMax[0];
    }
  else
    {
      NMax1 = NMax[0];
      NMax2 = NMax[1];
    }
  delete NMax;

  bool ExportOneBody = false;
  if ((Manager.GetBoolean("singleparticle-spectrum") == false)||((Manager.GetBoolean("export-onebody")) || (Manager.GetBoolean("export-onebodytext")) || (Manager.GetBoolean("singleparticle-chernnumber")) || (Manager.GetBoolean("singleparticle-curvature")) ))
    ExportOneBody = true;

  cout << "Creating TB model"<<endl;
  TightBindingModelOFLGenericLatticeWithSymmetry TightBindingModel(Manager.GetInteger("nbr-points1"), Manager.GetInteger("nbr-points2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), MyCutOffMode, 
								       Manager.GetDouble("cut-off-momentum"), NMax1, NMax2, LatticeDepth, Manager.GetInteger("output-nbr-bands"), Manager.GetDouble("precision-threshold"), ExportOneBody);


  // assign filename
    // generate filename according to TightBindingModel
  if (Manager.GetString("eigenvalue-file")!=0)
    {
      char *TmpExtension = GetExtensionFromFileName(Manager.GetString("eigenvalue-file"),4);
      if (TmpExtension!=0)
	{
	  FilePrefix = RemoveExtensionFromFileName(Manager.GetString("eigenvalue-file"), TmpExtension);
	  delete [] TmpExtension;
	}
      else
	strcpy(FilePrefix,Manager.GetString("eigenvalue-file"));
    }
  else
    {
      int Offset=0;

      if (Manager.GetBoolean("boson") == false)
	Offset+=sprintf (FilePrefix, "fermions");
      else
	Offset+=sprintf (FilePrefix, "bosons");
      
      switch (Manager.GetInteger("nbr-bands"))
	{
	case 1:
	  Offset+=sprintf (FilePrefix+Offset, "_singleband");
	  break;
	case 2:
	  Offset+=sprintf (FilePrefix+Offset, "_twobands");
	  break;
	default:
	  Offset+=sprintf (FilePrefix+Offset, "_%ldbands",Manager.GetInteger("nbr-bands"));
	  break;
	}
      
      if (TightBindingModel.GetDescriptor()!=NULL)
	Offset+=sprintf (FilePrefix+Offset, "_%s",TightBindingModel.GetDescriptor());
      else
	Offset+=sprintf (FilePrefix+Offset, "_unnamedOFL");
      
      // size / character of unit cell
      Offset+=sprintf (FilePrefix+Offset, "_flav_%d_sl_%d_sym_%d_%d", TightBindingModel.GetNbrSubLatticeFlavours(), TightBindingModel.GetNbrSubLatticesPerFlavour(), TightBindingModel.GetSymmetryMultiplier1(), TightBindingModel.GetSymmetryMultiplier2());
      
      // number of particles / size of system
      Offset+=sprintf (FilePrefix+Offset, "_n_%d_x_%d_y_%d", NbrParticles, NbrSitesX, NbrSitesY);

      // interaction parameters
      if (Manager.GetBoolean("flat-band") == false)
	Offset+=sprintf (FilePrefix+Offset, "_u_%g", Manager.GetDouble("u-potential"));
      if (Manager.GetDouble("spin-anisotropy")!=1.0)
	Offset+=sprintf (FilePrefix+Offset, "_iso_%g", Manager.GetDouble("spin-anisotropy"));
      Offset+=sprintf (FilePrefix+Offset, "_las_%g_gx_%g_gy_%g", Manager.GetDouble("depth"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"));
    }

  

  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {            
      if (Manager.GetBoolean("singleparticle-chernnumber") == true)      
	{
	  cout << "Chern number = " << TightBindingModel.ComputeChernNumber(0) << endl;
	}
      sprintf(EigenvalueOutputFile,"%s_singleparticle.dat", FilePrefix);
      
      TightBindingModel.WriteAsciiSpectrum(EigenvalueOutputFile);
      double BandSpread = TightBindingModel.ComputeBandSpread(0);
      double DirectBandGap = TightBindingModel.ComputeDirectBandGap(0);
      cout << "Spread = " << BandSpread << "  Direct Gap = " << DirectBandGap  << "  Flattening = " << (BandSpread / DirectBandGap) << endl;

      if (Manager.GetBoolean("singleparticle-curvature"))
	{
	  char CurvatureFileName[512];
	  sprintf(CurvatureFileName,"%s_curvature.dat", FilePrefix);
	  TightBindingModel.ComputeBerryCurvature(/*band */ 0, CurvatureFileName);
	}
      
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
  int MaxKx = NbrSitesX*TightBindingModel.GetSymmetryMultiplier1() - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSitesY*TightBindingModel.GetSymmetryMultiplier2() - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }  
  
  sprintf(EigenvalueOutputFile, "%s.dat", FilePrefix);

  RealSymmetricMatrix UPotentialMatrix(TightBindingModel.GetNbrSubLatticeFlavours());
  UPotentialMatrix.SetAllEntries(Manager.GetDouble("u-potential"));
  if (Manager.GetDouble("spin-anisotropy")!=1.0)
    UPotentialMatrix.SetOffDiagonalEntries(Manager.GetDouble("u-potential")*Manager.GetDouble("spin-anisotropy"));

  bool FirstRunFlag = true;
  for (int i = MinKx; i <= MaxKx; ++i)
    {
      for (int j = MinKy; j <= MaxKy; ++j)
	{
	  cout << "(kx=" << i << ",ky=" << j << ") " << endl;
	  ParticleOnSphere* Space = 0;
	  AbstractQHEHamiltonian* Hamiltonian;
	  switch ( Manager.GetInteger("nbr-bands"))
	    {
	    case 1: 
	      {
		if (Manager.GetBoolean("boson") == false)
		  {
		    Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSitesX*TightBindingModel.GetSymmetryMultiplier1(), NbrSitesY*TightBindingModel.GetSymmetryMultiplier2(), i, j);
		  }
		else
		  {
		    Space = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSitesX*TightBindingModel.GetSymmetryMultiplier1(), NbrSitesY*TightBindingModel.GetSymmetryMultiplier2(), i, j);
		  }
		cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		Hamiltonian = new ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian((ParticleOnSphere*)Space, NbrParticles, &TightBindingModel, /* RealSymmetricMatrix &uPotentialMatrix */ Manager.GetDouble("u-potential"), Manager.GetBoolean("flat-band"), Manager.GetBoolean("no-dispersion"), Architecture.GetArchitecture(), Memory);
		
		break;
	      }
	    case 2:
	      {
		if (Manager.GetBoolean("boson") == false)
		  {
		    Space = new FermionOnSquareLatticeWithSpinMomentumSpace(NbrParticles, NbrSitesX*TightBindingModel.GetSymmetryMultiplier1(), NbrSitesY*TightBindingModel.GetSymmetryMultiplier2(), i, j);
		  }
		else
		  {
		    Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace(NbrParticles, NbrSitesX*TightBindingModel.GetSymmetryMultiplier1(), NbrSitesY*TightBindingModel.GetSymmetryMultiplier2(), i, j);
		  }
		cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
		Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
		Hamiltonian = new ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian((ParticleOnSphereWithSpin*)Space, NbrParticles, &TightBindingModel, /* RealSymmetricMatrix &uPotentialMatrix */ Manager.GetDouble("u-potential"), Manager.GetBoolean("flat-band"), Manager.GetBoolean("no-dispersion"), Architecture.GetArchitecture(), Memory);
		break;
	      }
	    default:
	      {
		cout << "Single band case not implemented at the moment"<<endl;
		exit(1);
		break;
	      }
	    }
	  
	  char* CommentLine = new char [256];
	  sprintf (CommentLine, "eigenvalues\n# kx ky ");
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d %d", i, j);
	  char* EigenstateOutputFile = new char [512];
	  //char* TmpExtention = new char [512];
	  //sprintf (TmpExtention, "_kx_%d_ky_%d", i, j);
	  //EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	  sprintf(EigenstateOutputFile,"%s_kx_%d_ky_%d", FilePrefix, i, j);
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	  delete Space;
	}
    }
  
  return 0;
}

