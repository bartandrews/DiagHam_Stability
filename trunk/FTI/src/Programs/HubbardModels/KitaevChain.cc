#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceFixedParity.h"

#include "Hamiltonian/ParticleOnLatticeRealSpacePairingHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelKitaevChain.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericRealMainTask.h"
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
  OptionManager Manager ("KitaevChain" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t-hopping", "nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "delta", "nearest neighbor superconducting coupling", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu", "on-site chemical potential", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
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
      cout << "see man page for option syntax or type KitaevChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSites = NbrSitesX ; 
  
  
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

 

  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "fermions_kitaevchain_ns_%d", NbrSitesX);
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t_%.6f_delta_%.6f_mu_%.6f", Manager.GetDouble("t-hopping"), Manager.GetDouble("delta"), 
	   Manager.GetDouble("mu"));

  char* CommentLine = new char [256];
  sprintf (CommentLine, "kx");

  char* EigenvalueOutputFile = new char [512];
//   if (Manager.GetDouble("u-potential") == 0.0)
    sprintf(EigenvalueOutputFile, "%s_%s.dat", FilePrefix, FileParameterString);
//   else
//     sprintf(EigenvalueOutputFile, "%s_%s_u_%f.dat", FilePrefix, FileParameterString, Manager.GetDouble("u-potential"));

  Abstract1DTightBindingModel* TightBindingModel;
  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true))
	ExportOneBody = true;
      TightBindingModel = new TightBindingModelKitaevChain(NbrSitesX, Manager.GetDouble("t-hopping"), Manager.GetDouble("delta"), Manager.GetDouble("mu"), 
							   Architecture.GetArchitecture(), ExportOneBody);
//       TightBindingModel->TestRealSpaceTightBindingHamiltonian();
//       cout << TightBindingModel->GetRealSpaceTightBindingNonHermitianHamiltonian() << endl;
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

  TightBindingModel = new  TightBindingModelKitaevChain(NbrSitesX, Manager.GetDouble("t-hopping"), Manager.GetDouble("delta"), 
							Manager.GetDouble("mu"), Architecture.GetArchitecture(), true);
  char* BandStructureOutputFile = new char [1024];
  sprintf (BandStructureOutputFile, "%s_%s_tightbinding.dat", FilePrefix, FileParameterString);
  TightBindingModel->WriteBandStructure(BandStructureOutputFile);


  RealSymmetricMatrix DensityDensityInteraction((TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand()) / 2, true);
//   if (Manager.GetDouble("u-potential") != 0.0)
//     {
//       double UPotential = Manager.GetDouble("u-potential");
//       for (int i = 0; i < NbrSites; ++i)
// 	{
// 	  DensityDensityInteraction.SetMatrixElement(i, i, UPotential);
// 	}
//     }
  bool FirstRunFlag = true;

  int MinXMomentum = 0;
  int MaxXMomentum = NbrSitesX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {
      MaxXMomentum = Manager.GetInteger("only-kx");
      MinXMomentum = MaxXMomentum;
    }

  if (Manager.GetBoolean("use-periodic"))
    {
      Lanczos.SetComplexAlgorithms();
      for (int XMomentum = MinXMomentum; XMomentum <= MaxXMomentum; ++XMomentum)
	{
// 	  int FermionParity = 0;
// 	  int MaxFermionParity = 1;
// 	  for (; FermionParity <= MaxFermionParity; ++FermionParity)
// 	    {
// 	      FermionOnLatticeRealSpaceFixedParity* Space = new FermionOnLatticeRealSpaceFixedParity(NbrSitesX, FermionParity);
// 	      cout << "Kx = " << XMomentum << "  Parity = " << FermionParity << endl;
// 	      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
// 		Memory = Architecture.GetArchitecture()->GetLocalMemory();
// 	      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	      
// 	      HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
// 	      ParticleOnLatticeRealSpacePairingHamiltonian* Hamiltonian = new ParticleOnLatticeRealSpacePairingHamiltonian(Space, NbrSites, TightBindingMatrix,
// 															   DensityDensityInteraction,
// 															   Architecture.GetArchitecture(), Memory);
	      
// 	      char* ContentPrefix = new char[256];
// 	      sprintf (ContentPrefix, "%d %d", XMomentum, FermionParity);
// 	      char* EigenstateOutputFile;
// 	      char* TmpExtention = new char [512];
// 	      sprintf (TmpExtention, "_kx_%d_i_%d", XMomentum, FermionParity);
// 	      EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	      
// 	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
// 	      FirstRunFlag = false;
// 	      MainTaskOperation TaskOperation (&Task);
// 	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
// 	      cout << "------------------------------------" << endl;
// 	      delete Hamiltonian;
// 	      delete Space;
// 	      delete[] EigenstateOutputFile;
// 	      delete[] ContentPrefix;
// 	    }
	}  
    }
  else
    {
      int FermionParity = 0;
      int MaxFermionParity = 1;
      for (; FermionParity <= MaxFermionParity; ++FermionParity)
	{
	  FermionOnLatticeRealSpaceFixedParity* Space = new FermionOnLatticeRealSpaceFixedParity(NbrSitesX, FermionParity);
	  cout << "Parity = " << FermionParity << endl;
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
	  ParticleOnLatticeRealSpacePairingHamiltonian* Hamiltonian = new ParticleOnLatticeRealSpacePairingHamiltonian(Space, NbrSites, TightBindingMatrix,
														       DensityDensityInteraction,
														       Architecture.GetArchitecture(), Memory);
	  
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d", FermionParity);
	  char* EigenstateOutputFile;
	  char* TmpExtention = new char [512];
	  sprintf (TmpExtention, "_par_%d", FermionParity);
	  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	  
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
