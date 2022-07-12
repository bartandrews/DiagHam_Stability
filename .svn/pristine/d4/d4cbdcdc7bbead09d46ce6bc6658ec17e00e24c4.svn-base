#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian.h"
 
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBandStructure.h"
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBlochHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "MainTask/GenericRealMainTask.h"

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
  OptionManager Manager ("FTIGenericInteractionFromFileTwoBands" , "0.02");
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
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "full-momentum", "compute the spectrum for all momentum sectors, disregarding symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-file", "name of the file containing the two-body interaction matrix elements");
  (*SystemGroup) += new BooleanOption  ('\n', "real-interaction", "assume that the two-body interaction matrix elements are real");
  (*SystemGroup) += new SingleStringOption  ('\n', "interaction-name", "name of the two-body interaction", "noname");
  (*SystemGroup) += new SingleDoubleOption  ('u', "interaction-rescaling", "global rescaling of the two-body intearction", 1.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "singleparticle-file", "optional name of the file containing the one-body matrix elements");
  (*SystemGroup) += new BooleanOption  ('\n', "full-singleparticle", "the one-body matrix element file contains off-diagonal inter-band contributions");
  (*SystemGroup) += new BooleanOption  ('\n', "complex-singlebody", "the one-body matrix element file contains complex entries (only valid when using --full-singleparticle)");
  (*SystemGroup) += new BooleanOption  ('\n', "add-valley", "add valley-like degree of freedom (i.e. U(1) symmetry) included in --interaction-file");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pz-value", "twice the valley Pz value", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ez-value", "twice the Ez =1/2(N_{1u}+N_{2d}-N_{1d}-N_{2u}) value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "add-spin", "add spin 1/2 degree of freedom while assuming an SU(2) invariant interaction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-value", "twice the spin Sz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-valleyspin", "use the spin per valley instead of --sz-value and --ez-value");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz1-value", "twice the Sz value in valley 1", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz2-value", "twice the Sz value in valley 2", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "conserve-bandoccuption", "assume that the interaction conserves the number of particles per band");
  (*SystemGroup) += new SingleStringOption ('\n', "selected-sectors", "provide an ascii file that indicates which symmetry sectors have to be computed");
  (*SystemGroup) += new BooleanOption  ('\n', "disable-pzsymmetry", "disable the valley Pz<->-Pz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "disable-szsymmetry", "disable the valley Sz<->-Sz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a file");
  (*SystemGroup) += new BooleanOption  ('\n', "flat-band", "use flat band model");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "apply a temporary energy shift during the diagonalization", 0.0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FTIGenericInteractionFromFileTwoBands -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "no interaction file defined" << endl;
      cout << "see man page for option syntax or type FTIGenericInteractionFromFileTwoBands -h" << endl;
      return -1;
    }
  if (!(IsFile(Manager.GetString("interaction-file"))))
    {
      cout << "interaction file " << Manager.GetString("interaction-file")<< " does not exist" << endl;
      return -1;
    }

  if ((Manager.GetBoolean("real-interaction") == true) && (Manager.GetBoolean("conserve-bandoccuption") == true))
    {
      cout << "warning, be sure that your interaction preserves band occupation when using --real-interaction" << endl;
    }
    
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSites = 2 * NbrSitesX * NbrSitesY;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  double EnergyShift = Manager.GetDouble("energy-shift");
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
  int NbrMomentumSectors = (MaxKx - MinKx + 1) * (MaxKy - MinKy + 1);
  int MinSz = NbrParticles & 1;
  int MaxSz = MinSz;
  if (Manager.GetBoolean("add-spin") == true)
    {
      MinSz = Manager.GetInteger("sz-value") | (NbrParticles & 1);
      if (Manager.GetBoolean("use-valleyspin") == true)
	{
	  MinSz = (Manager.GetInteger("sz1-value") + Manager.GetInteger("sz2-value")) | (NbrParticles & 1);
	}
      MaxSz = MinSz;
    }  
  int MinPz = NbrParticles & 1;
  int MaxPz = MinPz;
  if (Manager.GetBoolean("add-valley") == true)
    {
      MinPz = Manager.GetInteger("pz-value") | (NbrParticles & 1);
      MaxPz = MinPz;
    }  
  int MinEz = NbrParticles & 1;
  int MaxEz = MinEz;
  if ((Manager.GetBoolean("add-valley") == true) && (Manager.GetBoolean("add-spin") == true))
    {
      MinEz = Manager.GetInteger("ez-value") | (NbrParticles & 1);
      if (Manager.GetBoolean("use-valleyspin") == true)
	{
	  MinEz = (Manager.GetInteger("sz1-value") - Manager.GetInteger("sz2-value")) | (NbrParticles & 1);
	}
      MaxEz = MinEz;
    }  
  bool DisablePzMinusPzSymmetry = Manager.GetBoolean("disable-pzsymmetry");
  bool DisableSzMinusSzSymmetry = Manager.GetBoolean("disable-szsymmetry");
  bool UsePzMinusPzSymmetry = !DisablePzMinusPzSymmetry;
  bool UseSzMinusSzSymmetry = !DisableSzMinusSzSymmetry;
    
  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }


  char* FileSystemGeometry = new char [512];
  char* CommentLine = new char [256];
  if (Manager.GetBoolean("add-valley") == false)
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	      sprintf (CommentLine, "eigenvalues\n# kx ky");
	    }
	  else
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY);
	      sprintf (CommentLine, "eigenvalues\n# kx ky n1 n2");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	      if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		{
		  UseSzMinusSzSymmetry = true;
		  sprintf (CommentLine, "eigenvalues\n# Sz Szsym kx ky");
		}
	      else
		{
		  sprintf (CommentLine, "eigenvalues\n# Sz kx ky");
		}
	    }
	  else
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinSz);
	      sprintf (CommentLine, "eigenvalues\n# Sz kx ky n1up n1down n2up n2down");
	    }
	}
   }
  else
    {
      if (Manager.GetBoolean("add-spin") == false)
	{
	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz);
	      if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		{
		  UsePzMinusPzSymmetry = true;
		  sprintf (CommentLine, "eigenvalues\n# Pz Pzsym kx ky");
		}
	      else
		{
		  sprintf (CommentLine, "eigenvalues\n# Pz kx ky");
		}
	    }
	  else
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz);
	      sprintf (CommentLine, "eigenvalues\n# Pz kx ky n1plus n1minus n2plus n2minus");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("conserve-bandoccuption") == false)
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d_ez_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz, MinEz, MinSz);
	      if (((MinSz == 0) && (DisableSzMinusSzSymmetry == false)) || ((MinPz == 0) && (DisablePzMinusPzSymmetry == false)))
		{
		  UseSzMinusSzSymmetry = true;
		  UsePzMinusPzSymmetry = true;
		  sprintf (CommentLine, "eigenvalues\n# Pz Sz Ez Pzsym Szsym kx ky");
		}
	      else
		{
		  sprintf (CommentLine, "eigenvalues\n# Pz Sz Ez kx ky");
		}
	    }
	  else
	    {
	      sprintf (FileSystemGeometry, "n_%d_ns_%d_x_%d_y_%d_pz_%d_ez_%d_sz_%d", NbrParticles, NbrSites, NbrSitesX, NbrSitesY, MinPz, MinEz, MinSz);
	      sprintf (CommentLine, "eigenvalues\n# Pz Sz Ez kx ky n1upplus n2upplus n1upminus n2upminus n1downplus n2downplus n1downminus n2downminus");
	    }
	}
    }
  char* FilePrefix = new char [512 + strlen(FileSystemGeometry)];
  if (Manager.GetBoolean("flat-band"))
    {
      sprintf (FilePrefix, "%s_twoband_flatband_%s_%s", StatisticPrefix, Manager.GetString("interaction-name"), FileSystemGeometry);
    }
  else
    {
      sprintf (FilePrefix, "%s_twoband_u_%.3f_%s_%s", StatisticPrefix, Manager.GetDouble("interaction-rescaling"),
	       Manager.GetString("interaction-name"), FileSystemGeometry);
    }
  
  char* EigenvalueOutputFile = new char [512 + strlen(FilePrefix)];
  
  if (Manager.GetString("eigenvalue-file") != 0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
    }
  else
    {
      sprintf (EigenvalueOutputFile, "%s.dat", FilePrefix);
    }
  

  Abstract2DTightBindingModel* TightBindingModel;
  
  if (Manager.GetString("import-onebody") == 0)
    {
      if (Manager.GetString("singleparticle-file") == 0)
	{
	  double* DummyChemicalPotentials = new double[2];
	  DummyChemicalPotentials[0] = 0.0;
	  DummyChemicalPotentials[1] = 0.0;
	  TightBindingModel = new TightBindingModel2DAtomicLimitLattice (NbrSitesX, NbrSitesY, 2, DummyChemicalPotentials,
									 0.0, 0.0, Architecture.GetArchitecture(), true);
	}
      else
	{
	  double* DummyChemicalPotentials = new double[2];
	  DummyChemicalPotentials[0] = 0.0;
	  DummyChemicalPotentials[1] = 0.0;
	  int* TmpKxValues = 0;
	  int* TmpKyValues = 0;
	  int* TmpValleyIndices = 0;
	  int* TmpBandIndices1 = 0;	  
	  int* TmpBandIndices2 = 0;	  
	  double* TmpOneBodyEnergies = 0;
	  MultiColumnASCIIFile OneBodyEnergyFile;
	  if (OneBodyEnergyFile.Parse(Manager.GetString("singleparticle-file")) == false)
	    {
	      OneBodyEnergyFile.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if (OneBodyEnergyFile.GetNbrLines() == 0)
	    {
	      cout << Manager.GetString("singleparticle-file") << " is an empty file" << endl;
	      return 0;
	    }
	  int NbrEnergies = OneBodyEnergyFile.GetNbrLines();
	  if (Manager.GetBoolean("full-singleparticle") == false)
	    {
	      if (((Manager.GetBoolean("add-valley") == false) && (NbrEnergies != (2 * NbrSitesX * NbrSitesY)))
		  || ((Manager.GetBoolean("add-valley") == true) && (NbrEnergies != (4 * NbrSitesX * NbrSitesY))))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of lines (has "
		       << NbrEnergies << ", should be " << (2 * NbrSitesX * NbrSitesY)
		       << " without valley, " << (4 * NbrSitesX * NbrSitesY) << " with vallley)" << endl;
		  return 0;
		}
	      if (((Manager.GetBoolean("add-valley") == false) && (OneBodyEnergyFile.GetNbrColumns() < 4))
		  || ((Manager.GetBoolean("add-valley") == true) && (OneBodyEnergyFile.GetNbrColumns() < 5)))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of columns (has "
		       << OneBodyEnergyFile.GetNbrColumns() << ", should be at least 4 without valley, 5 with valley)" << endl;
		  return 0;
		}
	    }
	  else
	    {
	      if (((Manager.GetBoolean("add-valley") == false) && (NbrEnergies != (4 * NbrSitesX * NbrSitesY)))
		  || ((Manager.GetBoolean("add-valley") == true) && (NbrEnergies != (8 * NbrSitesX * NbrSitesY))))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of lines (has "
		       << NbrEnergies << ", should be " << (4 * NbrSitesX * NbrSitesY)
		       << " without valley and --full-singleparticle, " << (8 * NbrSitesX * NbrSitesY) << " with vallley and --full-singleparticle)" << endl;
		  return 0;
		}
	      if (((Manager.GetBoolean("add-valley") == false) && (OneBodyEnergyFile.GetNbrColumns() < 5))
		  || ((Manager.GetBoolean("add-valley") == true) && (OneBodyEnergyFile.GetNbrColumns() < 6)))
		{
		  cout << Manager.GetString("singleparticle-file") << " has a wrong number of columns (has "
		       << OneBodyEnergyFile.GetNbrColumns() << ", should be at least 5 without valley, 6 with valley when using --full-singleparticle)" << endl;
		  return 0;
		}
	    }

	  int TmpNbrBands = 2;
	  if (Manager.GetBoolean("add-valley") == true)
	    {
	      TmpNbrBands = 4;
	    }
	  if (Manager.GetBoolean("add-spin") == true)
	    {
	      TmpNbrBands *= 2;
	    }
	  int* KxValues = new int[NbrSitesX * NbrSitesY];
	  int* KyValues = new int[NbrSitesX * NbrSitesY];
	  int TmpIndex = 0;
	  int* IndexToKIndex = new int[NbrSitesX * NbrSitesY];
	  int* NbrBandsPerKSector = new int[NbrSitesX * NbrSitesY];
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
		{
		  KxValues[TmpIndex] = i;
		  KyValues[TmpIndex] = j;
		  IndexToKIndex[(i * NbrSitesY) + j] = TmpIndex;
		  ++TmpIndex;
		}
	    }
	  TmpKxValues = OneBodyEnergyFile.GetAsIntegerArray(0);
	  TmpKyValues = OneBodyEnergyFile.GetAsIntegerArray(1);

	  if (Manager.GetBoolean("full-singleparticle") == false)
	    {
	      double** Energies = new double*[NbrSitesX * NbrSitesY];
	      TmpIndex = 0;
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      Energies[TmpIndex] = new double[TmpNbrBands];
		      ++TmpIndex;
		    }
		}
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(3);
		  if (TmpOneBodyEnergies == 0)
		    {
		      OneBodyEnergyFile.DumpErrors(cout) << endl;
		      return 0;
		    }
		  
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i]] = TmpOneBodyEnergies[i];
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i]] = TmpOneBodyEnergies[i];
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][2 + TmpBandIndices1[i]] = TmpOneBodyEnergies[i];
			}
		    }			      
		}
	      else
		{
		  TmpValleyIndices = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(3);
		  TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(4);
		  if (TmpOneBodyEnergies == 0)
		    {
		      OneBodyEnergyFile.DumpErrors(cout) << endl;
		      return 0;
		    }
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i] + (TmpValleyIndices[i] + 1)] = TmpOneBodyEnergies[i];
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][TmpBandIndices1[i] + (TmpValleyIndices[i] + 1)] = TmpOneBodyEnergies[i];
			  Energies[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]][4 + TmpBandIndices1[i] + (TmpValleyIndices[i] + 1)] = TmpOneBodyEnergies[i];
			}
		    }			      
		}
	      TightBindingModel = new TightBindingModel2DExplicitBandStructure (NbrSitesX, NbrSitesY, TmpNbrBands, 0.0, 0.0,
										KxValues, KyValues, Energies,
										Architecture.GetArchitecture(), true);
	      TmpIndex = 0;
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      delete[] Energies[TmpIndex];
		      ++TmpIndex;
		    }
		}
	      delete[] Energies;
	    }
	  else
	    {
	      HermitianMatrix* BlochHamiltonian = new HermitianMatrix[NbrSitesX * NbrSitesY];
	      TmpIndex = 0;
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      BlochHamiltonian[TmpIndex] = HermitianMatrix(TmpNbrBands, true);
		      ++TmpIndex;
		    }
		}
	      Complex* TmpComplexOneBodyEnergies = 0;
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpBandIndices2 = OneBodyEnergyFile.GetAsIntegerArray(3);
		  if (Manager.GetBoolean("complex-singlebody") == false)
		    {
		      TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(4);
		      TmpComplexOneBodyEnergies = 0;
		    }
		  else
		    {
		      TmpComplexOneBodyEnergies = OneBodyEnergyFile.GetAsComplexArray(4);
		      TmpOneBodyEnergies = 0;
		    }
		  if ((TmpOneBodyEnergies == 0) && (TmpComplexOneBodyEnergies == 0))
		    {
		      OneBodyEnergyFile.DumpErrors(cout) << endl;
		      return 0;
		    }
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  if (TmpOneBodyEnergies != 0)
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i], TmpBandIndices2[i], TmpOneBodyEnergies[i]);
			    }
			  else
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i], TmpBandIndices2[i], TmpComplexOneBodyEnergies[i]);
			    }
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  if (TmpOneBodyEnergies != 0)
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i], TmpBandIndices2[i], TmpOneBodyEnergies[i]);
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(2 + TmpBandIndices1[i], 2 + TmpBandIndices2[i], TmpOneBodyEnergies[i]);
			    }
			  else
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i], TmpBandIndices2[i], TmpComplexOneBodyEnergies[i]);
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(2 + TmpBandIndices1[i], 2 + TmpBandIndices2[i], TmpComplexOneBodyEnergies[i]);
			    }
			}
		    }			      
		}
	      else
		{
		  TmpBandIndices1 = OneBodyEnergyFile.GetAsIntegerArray(2);
		  TmpBandIndices2 = OneBodyEnergyFile.GetAsIntegerArray(3);
		  TmpValleyIndices = OneBodyEnergyFile.GetAsIntegerArray(4);
		  if (Manager.GetBoolean("complex-singlebody") == false)
		    {
		      TmpOneBodyEnergies = OneBodyEnergyFile.GetAsDoubleArray(5);
		    }
		  else
		    {
		      TmpComplexOneBodyEnergies = OneBodyEnergyFile.GetAsComplexArray(5);
		    }
		  if ((TmpOneBodyEnergies == 0) && (TmpComplexOneBodyEnergies == 0))
		    {
		      OneBodyEnergyFile.DumpErrors(cout) << endl;
		      return 0;
		    }
		  if (Manager.GetBoolean("add-spin") == false)
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  if (TmpOneBodyEnergies != 0)
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpOneBodyEnergies[i]);
			    }
			  else
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpComplexOneBodyEnergies[i]);
			    }
			}
		    }
		  else
		    {
		      for (int i = 0; i < NbrEnergies; ++i)
			{
			  if (TmpOneBodyEnergies != 0)
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpOneBodyEnergies[i]);
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(4 + TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), 4 + TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpOneBodyEnergies[i]);
			    }
			  else
			    {
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpComplexOneBodyEnergies[i]);
			      BlochHamiltonian[IndexToKIndex[(TmpKxValues[i] * NbrSitesY) + TmpKyValues[i]]].SetMatrixElement(4 + TmpBandIndices1[i] + (TmpValleyIndices[i] + 1), 4 + TmpBandIndices2[i] + (TmpValleyIndices[i] + 1), TmpComplexOneBodyEnergies[i]);
			    }
			}
		    }
		}
	      TightBindingModel = new TightBindingModel2DExplicitBlochHamiltonian (NbrSitesX, NbrSitesY, TmpNbrBands, 0.0, 0.0,
										   KxValues, KyValues, BlochHamiltonian,
										   Architecture.GetArchitecture(), true);
	      delete[] BlochHamiltonian;
	    }
	  char* BandStructureOutputFile = new char [64 + strlen(FilePrefix)];
	  sprintf (BandStructureOutputFile, "%s_tightbinding.dat", FilePrefix);
	  TightBindingModel->WriteBandStructure(BandStructureOutputFile);

	  delete[] IndexToKIndex;
	  delete[] KxValues;
	  delete[] KyValues;
	}      
    }
  else
    {
      TightBindingModel = new Generic2DTightBindingModel(Manager.GetString("import-onebody")); 
    }

  bool FirstRunFlag = true;

  if (Manager.GetBoolean("real-interaction"))
    {
      Lanczos.SetRealAlgorithms();
    }

  int* NbrParticlesBand1UpPlus = 0;
  int* NbrParticlesBand2UpPlus = 0;
  int* NbrParticlesBand1UpMinus = 0;
  int* NbrParticlesBand2UpMinus = 0;
  int* NbrParticlesBand1DownPlus = 0;
  int* NbrParticlesBand2DownPlus = 0;
  int* NbrParticlesBand1DownMinus = 0;
  int* NbrParticlesBand2DownMinus = 0;
  int* KxMomenta = 0;
  int* KyMomenta = 0;
  int* SzValues = 0;
  int* PzValues = 0;
  int* EzValues = 0;
  int* SzParityValues1 = 0;
  int* SzParityValues2 = 0;
  int* PzParityValues1 = 0;
  int* PzParityValues2 = 0;
  int NbrSymmetrySectors = NbrMomentumSectors;
  if (Manager.GetString("selected-sectors") == 0)
    {
      if (Manager.GetBoolean("conserve-bandoccuption") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  NbrSymmetrySectors = 1;
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = NbrParticles;
		  NbrParticlesBand2UpPlus[0] = 0;
		  NbrParticlesBand1UpMinus[0] = 0;
		  NbrParticlesBand2UpMinus[0] = 0;
		  NbrParticlesBand1DownPlus[0] = 0;
		  NbrParticlesBand2DownPlus[0] = 0;
		  NbrParticlesBand1DownMinus[0] = 0;
		  NbrParticlesBand2DownMinus[0] = 0;	      
		  SzParityValues1[0] = 0;
		  SzParityValues2[0] = 0;
		  PzParityValues1[0] = 0;
		  PzParityValues2[0] = 0;
		}
	      else
		{
		  NbrSymmetrySectors = 1;
		  if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		    {
		      NbrSymmetrySectors = 2;
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = (NbrParticles + MinPz) / 2;
		  NbrParticlesBand2UpPlus[0] = 0;
		  NbrParticlesBand1UpMinus[0] = (NbrParticles - MinPz) / 2;
		  NbrParticlesBand2UpMinus[0] = 0;
		  NbrParticlesBand1DownPlus[0] = 0;
		  NbrParticlesBand2DownPlus[0] = 0;
		  NbrParticlesBand1DownMinus[0] = 0;
		  NbrParticlesBand2DownMinus[0] = 0;
		  SzParityValues1[0] = 0;
		  SzParityValues2[0] = 0;
		  PzParityValues2[0] = 0;
		  if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		    {
		      PzParityValues1[0] = 1;
		      NbrParticlesBand1UpPlus[NbrMomentumSectors] = (NbrParticles + MinPz) / 2;
		      NbrParticlesBand2UpPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1UpMinus[NbrMomentumSectors] = (NbrParticles - MinPz) / 2;
		      NbrParticlesBand2UpMinus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1DownPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand2DownPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1DownMinus[NbrMomentumSectors] = 0;
		      NbrParticlesBand2DownMinus[NbrMomentumSectors] = 0;
		      SzParityValues1[NbrMomentumSectors] = 0;
		      PzParityValues1[NbrMomentumSectors] = -1;
		      SzParityValues2[NbrMomentumSectors] = 0;
		      PzParityValues2[NbrMomentumSectors] = 0;
		    }
		  else
		    {
		      PzParityValues1[0] = 0;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  NbrSymmetrySectors = 1;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      NbrSymmetrySectors = 2;
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz) / 2;
		  NbrParticlesBand2UpPlus[0] = 0;
		  NbrParticlesBand1UpMinus[0] = 0;
		  NbrParticlesBand2UpMinus[0] = 0;
		  NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz) / 2;
		  NbrParticlesBand2DownPlus[0] = 0;
		  NbrParticlesBand1DownMinus[0] = 0;
		  NbrParticlesBand2DownMinus[0] = 0;	      
		  PzParityValues1[0] = 0;
		  PzParityValues2[0] = 0;
		  SzParityValues2[0] = 0;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      SzParityValues1[0] = 1;
		      NbrParticlesBand1UpPlus[NbrMomentumSectors] = (NbrParticles + MinPz) / 2;
		      NbrParticlesBand2UpPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1UpMinus[NbrMomentumSectors] = (NbrParticles - MinPz) / 2;
		      NbrParticlesBand2UpMinus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1DownPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand2DownPlus[NbrMomentumSectors] = 0;
		      NbrParticlesBand1DownMinus[NbrMomentumSectors] = 0;
		      NbrParticlesBand2DownMinus[NbrMomentumSectors] = 0;
		      PzParityValues1[NbrMomentumSectors] = 0;
		      SzParityValues1[NbrMomentumSectors] = -1;
		      PzParityValues2[NbrMomentumSectors] = 0;
		      SzParityValues2[NbrMomentumSectors] = 0;
		    }
		  else
		    {
		      SzParityValues1[0] = 0;
		    }
		}
	      else
		{
		  NbrSymmetrySectors = 1;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      NbrSymmetrySectors *= 2;
		    }
		  if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
		    {
		      NbrSymmetrySectors *= 2;
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus[0] = (NbrParticles + MinSz + MinPz + MinEz);
		  NbrParticlesBand1UpMinus[0] = (NbrParticles + MinSz - MinPz - MinEz);
		  NbrParticlesBand1DownPlus[0] = (NbrParticles - MinSz + MinPz - MinEz);
		  NbrParticlesBand1DownMinus[0] = (NbrParticles - MinSz - MinPz + MinEz);			  
		  NbrParticlesBand2UpPlus[0] = 0;
		  NbrParticlesBand2UpMinus[0] = 0;
		  NbrParticlesBand2DownPlus[0] = 0;
		  NbrParticlesBand2DownMinus[0] = 0;			  
		  if ((NbrParticlesBand1UpPlus[0] < 0) || (NbrParticlesBand1UpMinus[0] < 0) || (NbrParticlesBand1DownPlus[0] < 0) || (NbrParticlesBand1DownMinus[0] < 0)
		      || ((NbrParticlesBand1UpPlus[0] & 3) != 0) ||  ((NbrParticlesBand1UpMinus[0] & 3) != 0)
		      || ((NbrParticlesBand1DownPlus[0] & 3) != 0) ||  ((NbrParticlesBand1DownMinus[0] & 3) != 0))
		    {
		      cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to 4N_{up,+}=" << NbrParticlesBand1UpPlus[0]
			   << " 4N_{up,-}=" << NbrParticlesBand1UpMinus[0] << " 4N_{down,+}=" << NbrParticlesBand1DownPlus[0]
			   << " 4N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		      return 0;
		    }
		  NbrParticlesBand1UpPlus[0] /= 4;
		  NbrParticlesBand1UpMinus[0] /= 4;
		  NbrParticlesBand1DownPlus[0] /= 4;
		  NbrParticlesBand1DownMinus[0] /= 4;
		  PzParityValues1[0] = 0;
		  SzParityValues1[0] = 0;
		  PzParityValues2[0] = 0;
		  SzParityValues2[0] = 0;
		  if ((NbrParticlesBand1UpPlus[0] > NbrParticles) || (NbrParticlesBand1UpMinus[0] > NbrParticles)
		      || (NbrParticlesBand1DownPlus[0] > NbrParticles) || (NbrParticlesBand1DownMinus[0] > NbrParticles))
		    {
		      cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to N_{up,+}=" << NbrParticlesBand1UpPlus[0]
			   << " N_{up,-}=" << NbrParticlesBand1UpMinus[0] << " N_{down,+}=" << NbrParticlesBand1DownPlus[0]
			   << " N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		      return 0;
		    }
		  cout << "N_{up,+}=" << NbrParticlesBand1UpPlus[0] << " N_{up,-}=" << NbrParticlesBand1UpMinus[0]
		       << " N_{down,+}=" << NbrParticlesBand1DownPlus[0] << " N_{down,-}=" << NbrParticlesBand1DownMinus[0] << endl;
		  if ((MinSz == 0) && (DisableSzMinusSzSymmetry == false))
		    {
		      if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
			{
			  PzParityValues1[0] = 1;
			  SzParityValues1[0] = 1;
			  PzParityValues1[NbrMomentumSectors] = -1;
			  SzParityValues1[NbrMomentumSectors] = 1;
			  PzParityValues1[2 * NbrMomentumSectors] = 1;
			  SzParityValues1[2 * NbrMomentumSectors] = -1;
			  PzParityValues1[3 * NbrMomentumSectors] = -1;
			  SzParityValues1[3 * NbrMomentumSectors] = -1;
			  for (int i = NbrMomentumSectors; i < NbrSymmetrySectors; i += NbrMomentumSectors)
			    {
			      NbrParticlesBand1UpPlus[i] = NbrParticlesBand1UpPlus[0];
			      NbrParticlesBand2UpPlus[i] = NbrParticlesBand2UpPlus[0];
			      NbrParticlesBand1UpMinus[i] = NbrParticlesBand1UpMinus[0];
			      NbrParticlesBand2UpMinus[i] = NbrParticlesBand2UpMinus[0];
			      NbrParticlesBand1DownPlus[i] = NbrParticlesBand1DownPlus[0];
			      NbrParticlesBand2DownPlus[i] = NbrParticlesBand2DownPlus[0];
			      NbrParticlesBand1DownMinus[i] = NbrParticlesBand1DownMinus[0];
			      NbrParticlesBand2DownMinus[i] = NbrParticlesBand2DownMinus[0];
			      PzParityValues2[i] = 0;
			      SzParityValues2[i] = 0;
			    }
			}
		      else
			{
			  SzParityValues1[0] = 1;
			  NbrParticlesBand1UpPlus[NbrMomentumSectors] = NbrParticlesBand1UpPlus[0];
			  NbrParticlesBand2UpPlus[NbrMomentumSectors] = NbrParticlesBand2UpPlus[0];
			  NbrParticlesBand1UpMinus[NbrMomentumSectors] = NbrParticlesBand1UpMinus[0];
			  NbrParticlesBand2UpMinus[NbrMomentumSectors] = NbrParticlesBand2UpMinus[0];
			  NbrParticlesBand1DownPlus[NbrMomentumSectors] = NbrParticlesBand1DownPlus[0];
			  NbrParticlesBand2DownPlus[NbrMomentumSectors] = NbrParticlesBand2DownPlus[0];
			  NbrParticlesBand1DownMinus[NbrMomentumSectors] = NbrParticlesBand1DownMinus[0];
			  NbrParticlesBand2DownMinus[NbrMomentumSectors] = NbrParticlesBand2DownMinus[0];
			  PzParityValues1[NbrMomentumSectors] = 0;
			  SzParityValues1[NbrMomentumSectors] = -1;
			  PzParityValues2[NbrMomentumSectors] = 0;
			  SzParityValues2[NbrMomentumSectors] = 0;
			}
		    }
		  else
		    {
		      if ((MinPz == 0) && (DisablePzMinusPzSymmetry == false))
			{
			  PzParityValues1[0] = 1;
			  NbrParticlesBand1UpPlus[NbrMomentumSectors] = NbrParticlesBand1UpPlus[0];
			  NbrParticlesBand2UpPlus[NbrMomentumSectors] = NbrParticlesBand2UpPlus[0];
			  NbrParticlesBand1UpMinus[NbrMomentumSectors] = NbrParticlesBand1UpMinus[0];
			  NbrParticlesBand2UpMinus[NbrMomentumSectors] = NbrParticlesBand2UpMinus[0];
			  NbrParticlesBand1DownPlus[NbrMomentumSectors] = NbrParticlesBand1DownPlus[0];
			  NbrParticlesBand2DownPlus[NbrMomentumSectors] = NbrParticlesBand2DownPlus[0];
			  NbrParticlesBand1DownMinus[NbrMomentumSectors] = NbrParticlesBand1DownMinus[0];
			  NbrParticlesBand2DownMinus[NbrMomentumSectors] = NbrParticlesBand2DownMinus[0];
			  PzParityValues1[NbrMomentumSectors] = -1;
			  SzParityValues1[NbrMomentumSectors] = 0;
			  PzParityValues2[NbrMomentumSectors] = 0;
			  SzParityValues2[NbrMomentumSectors] = 0;
			}
		    }
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  NbrSymmetrySectors = NbrMomentumSectors * (NbrParticles + 1);
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; i += NbrMomentumSectors)
		    {
		      NbrParticlesBand1UpPlus[i] = i;
		      NbrParticlesBand2UpPlus[i] = NbrParticles - i;
		      NbrParticlesBand1UpMinus[i] = 0;
		      NbrParticlesBand2UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = 0;
		      NbrParticlesBand2DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = 0;
		      NbrParticlesBand2DownMinus[i] = 0;
		    }
		}
	      else
		{
		  NbrSymmetrySectors = 0;
		  int MaxNbrParticlesPerSector = NbrParticles;
		  if (MaxNbrParticlesPerSector > (NbrSitesX * NbrSitesY))
		    {
		      MaxNbrParticlesPerSector = NbrSitesX * NbrSitesY;
		    }
		  for (int TmpPlus1 = 0; TmpPlus1 <= MaxNbrParticlesPerSector; ++TmpPlus1)
		    {
		      int TmpPlus2 = ((MinPz + NbrParticles) / 2) - TmpPlus1;
		      if ((TmpPlus2 >= 0) && (TmpPlus2 <= MaxNbrParticlesPerSector))
			{
			  for (int TmpMinus1 = 0; TmpMinus1 <= MaxNbrParticlesPerSector; ++TmpMinus1)
			    {
			      int TmpMinus2 = ((NbrParticles - MinPz) / 2) - TmpMinus1;
			      if ((TmpMinus2 >= 0) && (TmpMinus2 <= MaxNbrParticlesPerSector))
				{
				  ++NbrSymmetrySectors;
				}
			    }
			}
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  NbrSymmetrySectors = 0;
		  for (int TmpPlus1 = 0; TmpPlus1 <= MaxNbrParticlesPerSector; ++TmpPlus1)
		    {
		      int TmpPlus2 = ((MinPz + NbrParticles) / 2) - TmpPlus1;
		      if ((TmpPlus2 >= 0) && (TmpPlus2 <= MaxNbrParticlesPerSector))
			{
			  for (int TmpMinus1 = 0; TmpMinus1 <= MaxNbrParticlesPerSector; ++TmpMinus1)
			    {
			      int TmpMinus2 = ((NbrParticles - MinPz) / 2) - TmpMinus1;
			      if ((TmpMinus2 >= 0) && (TmpMinus2 <= MaxNbrParticlesPerSector))
				{
				  NbrParticlesBand1UpPlus[NbrSymmetrySectors] = TmpPlus1;
				  NbrParticlesBand2UpPlus[NbrSymmetrySectors] = TmpPlus2;
				  NbrParticlesBand1UpMinus[NbrSymmetrySectors] = TmpMinus1;
				  NbrParticlesBand2UpMinus[NbrSymmetrySectors] = TmpMinus2;
				  NbrParticlesBand1DownPlus[NbrSymmetrySectors] = 0;
				  NbrParticlesBand2DownPlus[NbrSymmetrySectors] = 0;
				  NbrParticlesBand1DownMinus[NbrSymmetrySectors] = 0;
				  NbrParticlesBand2DownMinus[NbrSymmetrySectors] = 0;
				  NbrSymmetrySectors += NbrMomentumSectors;
				}
			    }
			}
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  NbrSymmetrySectors = 0;
		  int MaxNbrParticlesPerSector = NbrParticles;
		  if (MaxNbrParticlesPerSector > (NbrSitesX * NbrSitesY))
		    {
		      MaxNbrParticlesPerSector = NbrSitesX * NbrSitesY;
		    }
		  for (int TmpUp1 = 0; TmpUp1 <= MaxNbrParticlesPerSector; ++TmpUp1)
		    {
		      int TmpUp2 = ((MinSz + NbrParticles) / 2) - TmpUp1;
		      if ((TmpUp2 >= 0) && (TmpUp2 <= MaxNbrParticlesPerSector))
			{
			  for (int TmpDown1 = 0; TmpDown1 <= MaxNbrParticlesPerSector; ++TmpDown1)
			    {
			      int TmpDown2 = ((NbrParticles - MinSz) / 2) - TmpDown1;
			      if ((TmpDown2 >= 0) && (TmpDown2 <= MaxNbrParticlesPerSector))
				{
				  ++NbrSymmetrySectors;
				}
			    }
			}
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  PzParityValues1 = new int[NbrSymmetrySectors];
		  SzParityValues1 = new int[NbrSymmetrySectors];
		  PzParityValues2 = new int[NbrSymmetrySectors];
		  SzParityValues2 = new int[NbrSymmetrySectors];
		  NbrSymmetrySectors = 0;
		  for (int TmpUp1 = 0; TmpUp1 <= MaxNbrParticlesPerSector; ++TmpUp1)
		    {
		      int TmpUp2 = ((MinSz + NbrParticles) / 2) - TmpUp1;
		      if ((TmpUp2 >= 0) && (TmpUp2 <= MaxNbrParticlesPerSector))
			{
			  for (int TmpDown1 = 0; TmpDown1 <= MaxNbrParticlesPerSector; ++TmpDown1)
			    {
			      int TmpDown2 = ((NbrParticles - MinSz) / 2) - TmpDown1;
			      if ((TmpDown2 >= 0) && (TmpDown2 <= MaxNbrParticlesPerSector))
				{
				  NbrParticlesBand1UpPlus[NbrSymmetrySectors] = TmpUp1;
				  NbrParticlesBand2UpPlus[NbrSymmetrySectors] = TmpUp2;
				  NbrParticlesBand1UpMinus[NbrSymmetrySectors] = 0;
				  NbrParticlesBand2UpMinus[NbrSymmetrySectors] = 0;
				  NbrParticlesBand1DownPlus[NbrSymmetrySectors] = TmpDown1;
				  NbrParticlesBand2DownPlus[NbrSymmetrySectors] = TmpDown2;
				  NbrParticlesBand1DownMinus[NbrSymmetrySectors] = 0;
				  NbrParticlesBand2DownMinus[NbrSymmetrySectors] = 0;
				  PzParityValues1[NbrSymmetrySectors] = 0;
				  PzParityValues2[NbrSymmetrySectors] = 0;
				  SzParityValues1[NbrSymmetrySectors] = 0;
				  SzParityValues2[NbrSymmetrySectors] = 0;
				  NbrSymmetrySectors += NbrMomentumSectors;
				}
			    }
			}
		    }
		}
	      else
		{
		  int NbrParticlesUpPlus = (NbrParticles + MinSz + MinPz + MinEz);
		  int NbrParticlesUpMinus = (NbrParticles + MinSz - MinPz - MinEz);
		  int NbrParticlesDownPlus = (NbrParticles - MinSz + MinPz - MinEz);
		  int NbrParticlesDownMinus = (NbrParticles - MinSz - MinPz + MinEz);			  
		  if ((NbrParticlesUpPlus < 0) || (NbrParticlesUpMinus < 0) || (NbrParticlesDownPlus < 0) || (NbrParticlesDownMinus < 0)
		      || ((NbrParticlesUpPlus & 3) != 0) ||  ((NbrParticlesUpMinus & 3) != 0)
		      || ((NbrParticlesDownPlus & 3) != 0) ||  ((NbrParticlesDownMinus & 3) != 0))
		    {
		      cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to 4N_{up,+}=" << NbrParticlesUpPlus
			   << " 4N_{up,-}=" << NbrParticlesUpMinus << " 4N_{down,+}=" << NbrParticlesDownPlus
			   << " 4N_{down,-}=" << NbrParticlesDownMinus << endl;
		      return 0;
		    }
		  NbrParticlesUpPlus /= 4;
		  NbrParticlesUpMinus /= 4;
		  NbrParticlesDownPlus /= 4;
		  NbrParticlesDownMinus /= 4;
		  if ((NbrParticlesUpPlus > NbrParticles) || (NbrParticlesUpMinus > NbrParticles)
		      || (NbrParticlesDownPlus > NbrParticles) || (NbrParticlesDownMinus > NbrParticles))
		    {
		      cout << "Incompatible values of N, 2Sz, 2Pz and 2Ez, lead to N_{up,+}=" << NbrParticlesUpPlus
			   << " N_{up,-}=" << NbrParticlesUpMinus << " N_{down,+}=" << NbrParticlesDownPlus
			   << " N_{down,-}=" << NbrParticlesDownMinus << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = 0;
		  int MaxNbrParticlesPerSector = NbrParticles;
		  if (MaxNbrParticlesPerSector > (NbrSitesX * NbrSitesY))
		    {
		      MaxNbrParticlesPerSector = NbrSitesX * NbrSitesY;
		    }
		  for (int TmpUpPlus1 = 0; TmpUpPlus1 <= MaxNbrParticlesPerSector; ++TmpUpPlus1)
		    {
		      int TmpUpPlus2 = NbrParticlesUpPlus - TmpUpPlus1;
		      if ((TmpUpPlus2 >= 0) && (TmpUpPlus2 <= MaxNbrParticlesPerSector))
			{
			  for (int TmpUpMinus1 = 0; TmpUpMinus1 <= MaxNbrParticlesPerSector; ++TmpUpMinus1)
			    {
			      int TmpUpMinus2 = NbrParticlesUpMinus - TmpUpMinus1;
			      if ((TmpUpPlus2 >= 0) && (TmpUpPlus2 <= MaxNbrParticlesPerSector))
				{
				  for (int TmpDownPlus1 = 0; TmpDownPlus1 <= MaxNbrParticlesPerSector; ++TmpDownPlus1)
				    {
				      int TmpDownPlus2 = NbrParticlesDownPlus - TmpDownPlus1;
				      if ((TmpDownPlus2 >= 0) && (TmpDownPlus2 <= MaxNbrParticlesPerSector))
					{
					  for (int TmpDownMinus1 = 0; TmpDownMinus1 <= MaxNbrParticlesPerSector; ++TmpDownMinus1)
					    {
					      int TmpDownMinus2 = NbrParticlesDownMinus - TmpDownMinus1;
					      if ((TmpDownMinus2 >= 0) && (TmpDownMinus2 <= MaxNbrParticlesPerSector))
						{
						  ++NbrSymmetrySectors;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		  NbrSymmetrySectors *= NbrMomentumSectors;
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  NbrSymmetrySectors = 0;
		  for (int TmpUpPlus1 = 0; TmpUpPlus1 <= MaxNbrParticlesPerSector; ++TmpUpPlus1)
		    {
		      int TmpUpPlus2 = NbrParticlesUpPlus - TmpUpPlus1;
		      if ((TmpUpPlus2 >= 0) && (TmpUpPlus2 <= MaxNbrParticlesPerSector))
			{
			  for (int TmpUpMinus1 = 0; TmpUpMinus1 <= MaxNbrParticlesPerSector; ++TmpUpMinus1)
			    {
			      int TmpUpMinus2 = NbrParticlesUpMinus - TmpUpMinus1;
			      if ((TmpUpPlus2 >= 0) && (TmpUpPlus2 <= MaxNbrParticlesPerSector))
				{
				  for (int TmpDownPlus1 = 0; TmpDownPlus1 <= MaxNbrParticlesPerSector; ++TmpDownPlus1)
				    {
				      int TmpDownPlus2 = NbrParticlesDownPlus - TmpDownPlus1;
				      if ((TmpDownPlus2 >= 0) && (TmpDownPlus2 <= MaxNbrParticlesPerSector))
					{
					  for (int TmpDownMinus1 = 0; TmpDownMinus1 <= MaxNbrParticlesPerSector; ++TmpDownMinus1)
					    {
					      int TmpDownMinus2 = NbrParticlesDownMinus - TmpDownMinus1;
					      if ((TmpDownMinus2 >= 0) && (TmpDownMinus2 <= MaxNbrParticlesPerSector))
						{
						  NbrParticlesBand1UpPlus[NbrSymmetrySectors] = TmpUpPlus1;
						  NbrParticlesBand2UpPlus[NbrSymmetrySectors] = TmpUpPlus2;
						  NbrParticlesBand1UpMinus[NbrSymmetrySectors] = TmpUpMinus1;
						  NbrParticlesBand2UpMinus[NbrSymmetrySectors] = TmpUpMinus2;
						  NbrParticlesBand1DownPlus[NbrSymmetrySectors] = TmpDownPlus1;
						  NbrParticlesBand2DownPlus[NbrSymmetrySectors] = TmpDownPlus2;
						  NbrParticlesBand1DownMinus[NbrSymmetrySectors] = TmpDownMinus1;
						  NbrParticlesBand2DownMinus[NbrSymmetrySectors] = TmpDownMinus2;
						  NbrSymmetrySectors += NbrMomentumSectors;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
  
      KxMomenta = new int[NbrSymmetrySectors];
      KyMomenta = new int[NbrSymmetrySectors];
      SzValues = new int[NbrSymmetrySectors];
      PzValues = new int[NbrSymmetrySectors];
      EzValues = new int[NbrSymmetrySectors];
      int TmpIndex = 0;
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  for (int j = MinKy; j <= MaxKy; ++j)
	    {
	      for (int k = 0; k < NbrSymmetrySectors; k += NbrMomentumSectors)
		{
		  NbrParticlesBand1UpPlus[k + TmpIndex] = NbrParticlesBand1UpPlus[k];
		  NbrParticlesBand2UpPlus[k + TmpIndex] = NbrParticlesBand2UpPlus[k];
		  NbrParticlesBand1UpMinus[k + TmpIndex] = NbrParticlesBand1UpMinus[k];
		  NbrParticlesBand2UpMinus[k + TmpIndex] = NbrParticlesBand2UpMinus[k];
		  NbrParticlesBand1DownPlus[k + TmpIndex] = NbrParticlesBand1DownPlus[k];
		  NbrParticlesBand2DownPlus[k + TmpIndex] = NbrParticlesBand2DownPlus[k];
		  NbrParticlesBand1DownMinus[k + TmpIndex] = NbrParticlesBand1DownMinus[k];
		  NbrParticlesBand2DownMinus[k + TmpIndex] = NbrParticlesBand2DownMinus[k];
		  KxMomenta[k + TmpIndex] = i;
		  KyMomenta[k + TmpIndex] = j;
		  SzValues[k + TmpIndex] = MinSz;
		  PzValues[k + TmpIndex] = MinPz;
		  EzValues[k + TmpIndex] = MinEz;
		  SzParityValues1[k + TmpIndex] = SzParityValues1[k];
		  SzParityValues2[k + TmpIndex] = SzParityValues2[k];
		  PzParityValues1[k + TmpIndex] = PzParityValues1[k];
		  PzParityValues2[k + TmpIndex] = PzParityValues2[k];
		}
	      ++TmpIndex;
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile SymmetrySectorsFile;
      if (SymmetrySectorsFile.Parse(Manager.GetString("selected-sectors")) == false)
	{
	  SymmetrySectorsFile.DumpErrors(cout) << endl;
	  return 0;
	}
      if (SymmetrySectorsFile.GetNbrLines() == 0)
	{
	  cout << Manager.GetString("selected-sectors") << " is an empty file" << endl;
	  return 0;
	}
      if (Manager.GetBoolean("conserve-bandoccuption") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if (SymmetrySectorsFile.GetNbrColumns() < 2)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least two)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		  SzValues = new int [NbrSymmetrySectors];
		  PzValues = new int [NbrSymmetrySectors];
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      SzValues[i] = NbrParticles;
		      PzValues[i] = NbrParticles;
		      EzValues[i] = NbrParticles;
		      NbrParticlesBand1UpPlus[i] = NbrParticles;
		      NbrParticlesBand2UpPlus[i] = 0;
		      NbrParticlesBand1UpMinus[i] = 0;
		      NbrParticlesBand2UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = 0;
		      NbrParticlesBand2DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = 0;
		      NbrParticlesBand2DownMinus[i] = 0;
		    }
		}
	      else
		{
		  if (SymmetrySectorsFile.GetNbrColumns() < 3)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least three when using --add-valley)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		  SzValues = new int [NbrSymmetrySectors];
		  PzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      SzValues[i] = NbrParticles;
		      EzValues[i] = NbrParticles;
		      NbrParticlesBand1UpPlus[i] = (NbrParticles + PzValues[i]) / 2;
		      NbrParticlesBand2UpPlus[i] = 0;
		      NbrParticlesBand1UpMinus[i] = (NbrParticles - PzValues[i]) / 2;
		      NbrParticlesBand2UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = 0;
		      NbrParticlesBand2DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = 0;
		      NbrParticlesBand2DownMinus[i] = 0;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if (SymmetrySectorsFile.GetNbrColumns() < 3)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least three when using --add-spin)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		  SzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		  PzValues = new int [NbrSymmetrySectors];
		  EzValues = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      PzValues[i] = NbrParticles;
		      EzValues[i] = NbrParticles;
		      NbrParticlesBand1UpPlus[i] = (NbrParticles + SzValues[i]) / 2;
		      NbrParticlesBand2UpPlus[i] = 0;
		      NbrParticlesBand1UpMinus[i] = 0;
		      NbrParticlesBand2UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = (NbrParticles - PzValues[i]) / 2;
		      NbrParticlesBand2DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = 0;
		      NbrParticlesBand2DownMinus[i] = 0;
		    }
		}
	      else
		{
		  if (SymmetrySectorsFile.GetNbrColumns() < 5)
		    {
		      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least five when using --add-spin and --add-valley)" << endl;
		      return 0;
		    }
		  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
		  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
		  KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
		  SzValues = SymmetrySectorsFile.GetAsIntegerArray(2);
		  PzValues = SymmetrySectorsFile.GetAsIntegerArray(3);
		  EzValues = SymmetrySectorsFile.GetAsIntegerArray(4);
		  NbrParticlesBand1UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2UpMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownPlus = new int [NbrSymmetrySectors];
		  NbrParticlesBand1DownMinus = new int [NbrSymmetrySectors];
		  NbrParticlesBand2DownMinus = new int [NbrSymmetrySectors];
		  for (int i = 0; i < NbrSymmetrySectors; ++i)
		    {
		      NbrParticlesBand1UpPlus[i] = (NbrParticles + SzValues[i] + PzValues[i] + EzValues[i]) / 4;
		      NbrParticlesBand2UpPlus[i] = 0;
		      NbrParticlesBand1UpMinus[i] = (NbrParticles + SzValues[i] - PzValues[i] - EzValues[i]) / 4;
		      NbrParticlesBand2UpMinus[i] = 0;
		      NbrParticlesBand1DownPlus[i] = (NbrParticles - SzValues[i] + PzValues[i] - EzValues[i]) / 4;
		      NbrParticlesBand2DownPlus[i] = 0;
		      NbrParticlesBand1DownMinus[i] = (NbrParticles - SzValues[i] - PzValues[i] + EzValues[i]) / 4;
		      NbrParticlesBand2DownMinus[i] = 0;
		    }
		}
	    }
	}
      else
	{
	  if (SymmetrySectorsFile.GetNbrColumns() < 10)
	    {
	      cout << Manager.GetString("selected-sectors") << " has a wrong number of columns (should be at least ten when using --conserve-bandoccuption)" << endl;
	      return 0;
	    }
	  NbrSymmetrySectors = SymmetrySectorsFile.GetNbrLines();
	  KxMomenta = SymmetrySectorsFile.GetAsIntegerArray(0);
	  KyMomenta = SymmetrySectorsFile.GetAsIntegerArray(1);
	  NbrParticlesBand1UpPlus = SymmetrySectorsFile.GetAsIntegerArray(2);
	  NbrParticlesBand1UpMinus = SymmetrySectorsFile.GetAsIntegerArray(3);
	  NbrParticlesBand1DownPlus = SymmetrySectorsFile.GetAsIntegerArray(4);
	  NbrParticlesBand1DownMinus = SymmetrySectorsFile.GetAsIntegerArray(5);
	  NbrParticlesBand2UpPlus = SymmetrySectorsFile.GetAsIntegerArray(6);
	  NbrParticlesBand2UpMinus = SymmetrySectorsFile.GetAsIntegerArray(7);
	  NbrParticlesBand2DownPlus = SymmetrySectorsFile.GetAsIntegerArray(8);
	  NbrParticlesBand2DownMinus = SymmetrySectorsFile.GetAsIntegerArray(9);
	  SzValues = new int [NbrSymmetrySectors];
	  PzValues = new int [NbrSymmetrySectors];
	  EzValues = new int [NbrSymmetrySectors];
	  for (int i = 0; i < NbrSymmetrySectors; ++i)
	    {
	      SzValues[i] = (NbrParticlesBand1UpPlus[i] + NbrParticlesBand1UpMinus[i]
			     - NbrParticlesBand1DownPlus[i] - NbrParticlesBand1DownMinus[i]
			     + NbrParticlesBand2UpPlus[i] + NbrParticlesBand2UpMinus[i]
			     - NbrParticlesBand2DownPlus[i] - NbrParticlesBand2DownMinus[i]);
	      PzValues[i] = (NbrParticlesBand1UpPlus[i] - NbrParticlesBand1UpMinus[i]
			     + NbrParticlesBand1DownPlus[i] - NbrParticlesBand1DownMinus[i]
			     + NbrParticlesBand2UpPlus[i] - NbrParticlesBand2UpMinus[i]
			     + NbrParticlesBand2DownPlus[i] - NbrParticlesBand2DownMinus[i]);
	      EzValues[i] = (NbrParticlesBand1UpPlus[i] - NbrParticlesBand1UpMinus[i]
			     - NbrParticlesBand1DownPlus[i] + NbrParticlesBand1DownMinus[i]
			     + NbrParticlesBand2UpPlus[i] - NbrParticlesBand2UpMinus[i]
			     - NbrParticlesBand2DownPlus[i] + NbrParticlesBand2DownMinus[i]);
	    }
	}
    }
  
  int TotalDim = 0;
  for (int SymmetrySectorIndex = 0; SymmetrySectorIndex < NbrSymmetrySectors; ++SymmetrySectorIndex)
    {
      if (Manager.GetBoolean("add-valley") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("conserve-bandoccuption") == false)
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ") : " << endl;
		}
	      else
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" <<  KyMomenta[SymmetrySectorIndex]
		       << ",n1=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
		       << ",n2=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex] << ") : " << endl;
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("conserve-bandoccuption") == false)
		{
		  if (UseSzMinusSzSymmetry == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",Sz<->-Sz=" << SzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	      else
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex]
		       << ",n1up=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
		       << ",n2up=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex]
		       << ",n1down=" << NbrParticlesBand1DownPlus[SymmetrySectorIndex]
		       << ",n2down=" << NbrParticlesBand2DownPlus[SymmetrySectorIndex] << ") : " << endl;
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("conserve-bandoccuption") == false)
		{
		  if (UsePzMinusPzSymmetry == false)
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ") : " << endl;
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",Pz<->-Pz=" << PzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	      else
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex]
		       << ",n1plus=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
		       << ",n2plus=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex]
		       << ",n1minus=" << NbrParticlesBand1UpMinus[SymmetrySectorIndex]
		       << ",n2minus=" << NbrParticlesBand2UpMinus[SymmetrySectorIndex] << ") : " << endl;
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("conserve-bandoccuption") == false)
		{
		  if ((UsePzMinusPzSymmetry == false) && (UseSzMinusSzSymmetry == false))
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",2ez=" << EzValues[SymmetrySectorIndex]<< ") : " << endl;		      
		    }
		  else
		    {
		      cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex] << ",2sz=" << SzValues[SymmetrySectorIndex] << ",2ez=" << EzValues[SymmetrySectorIndex]<< ",Sz<->-Sz=" << SzParityValues1[SymmetrySectorIndex] << ",Pz<->-Pz=" << PzParityValues1[SymmetrySectorIndex] << ") : " << endl;
		    }
		}
	      else
		{
		  cout << "(kx=" << KxMomenta[SymmetrySectorIndex] << ",ky=" << KyMomenta[SymmetrySectorIndex] << ",2pz=" << PzValues[SymmetrySectorIndex]
		       << ",n1upplus=" << NbrParticlesBand1UpPlus[SymmetrySectorIndex]
		       << ",n2upplus=" << NbrParticlesBand2UpPlus[SymmetrySectorIndex]
		       << ",n1upminus=" << NbrParticlesBand1UpMinus[SymmetrySectorIndex]
		       << ",n2upminus=" << NbrParticlesBand2UpMinus[SymmetrySectorIndex]
		       << ",n1downplus=" << NbrParticlesBand1DownPlus[SymmetrySectorIndex]
		       << ",n2downplus=" << NbrParticlesBand2DownPlus[SymmetrySectorIndex]
		       << ",n1downminus=" << NbrParticlesBand1DownMinus[SymmetrySectorIndex]
		       << ",n2downminus=" << NbrParticlesBand2DownMinus[SymmetrySectorIndex] << ") : " << endl;
		}
	    }
	}
      ParticleOnSphereWithSpin* Space = 0;
      AbstractQHEHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (Manager.GetBoolean("boson") == false)
	{
	  if (Manager.GetBoolean("add-spin") == false)
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 32)
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			}
		      else
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      if ((NbrSitesX * NbrSitesY) <= 32)
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpace (NbrParticles, NbrParticlesBand1UpPlus[SymmetrySectorIndex], NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			}
		      else
			{
			  Space = new FermionOnSquareLatticeWithSpinMomentumSpaceLong (NbrParticles, NbrParticlesBand1UpPlus[SymmetrySectorIndex], NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			}
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {			  
		      if ((NbrSitesX * NbrSitesY) <= 16)
			{
			  if (PzParityValues1[SymmetrySectorIndex] == 0)
			    {
			      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
											  PzValues[SymmetrySectorIndex], 10000000ul);
			    }
			  else
			    {
			      // Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (NbrParticles, NbrSitesX, NbrSitesY,
			      // 									    KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
			      // 									    PzValues[SymmetrySectorIndex],
			      // 									    (PzParityValues1[SymmetrySectorIndex] == -1), 10000000ul);
			      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry (NbrParticles, NbrSitesX, NbrSitesY,
														  KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
														  PzValues[SymmetrySectorIndex],
														  (PzParityValues1[SymmetrySectorIndex] == -1), 10000000ul);
			    }			    
			}
		      else
			{
			  cout << "SU(4) not supported with more than 16 momenta" << endl;
			  Space = 0;
			  //		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      int FakePz = (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]
				    + NbrParticlesBand1UpMinus[SymmetrySectorIndex] - NbrParticlesBand2UpMinus[SymmetrySectorIndex]);
		      int FakeEz = (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]
				    - NbrParticlesBand1UpMinus[SymmetrySectorIndex] + NbrParticlesBand2UpMinus[SymmetrySectorIndex]);
		      if ((NbrSitesX * NbrSitesY) <= 16)
			{
			  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
										      PzValues[SymmetrySectorIndex], FakePz, FakeEz, 10000000ul);
			}
		      else
			{
			  cout << "SU(4) not supported with more than 16 momenta" << endl;
			  Space = 0;
			  //		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 16)
			{
			  if (SzParityValues1[SymmetrySectorIndex] == 0)
			    {
			      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
											  SzValues[SymmetrySectorIndex], 10000000ul);
			    }
			  else
			    {
			      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (NbrParticles, NbrSitesX, NbrSitesY,
												    KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
												    SzValues[SymmetrySectorIndex],
												    (SzParityValues1[SymmetrySectorIndex] == -1), 10000000ul);
			    }
			}
		      else
			{
			  cout << "SU(4) not supported with more than 16 momenta" << endl;
			  Space = 0;
			  //		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      int FakePz = (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]
				    + NbrParticlesBand1DownPlus[SymmetrySectorIndex] - NbrParticlesBand2DownPlus[SymmetrySectorIndex]);
		      int FakeEz = (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]
				    - NbrParticlesBand1DownPlus[SymmetrySectorIndex] + NbrParticlesBand2DownPlus[SymmetrySectorIndex]);
		      if ((NbrSitesX * NbrSitesY) <= 16)
			{
			  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
										      SzValues[SymmetrySectorIndex], FakePz, FakeEz, 10000000ul);
			}
		      else
			{
			  cout << "SU(4) not supported with more than 16 momenta" << endl;
			  Space = 0;
			  //		      Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      if ((NbrSitesX * NbrSitesY) <= 8)
			{
			  Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
										      NbrParticlesBand1DownMinus[SymmetrySectorIndex], NbrParticlesBand1DownPlus[SymmetrySectorIndex],
										      NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex], 10000000ul);
			}
		      else
			{
			  if ((NbrSitesX * NbrSitesY) <= 16)
			    {			  
			      Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong(NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
											     NbrParticlesBand1DownMinus[SymmetrySectorIndex], NbrParticlesBand1DownPlus[SymmetrySectorIndex],
											     NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex], 10000000ul);
			    }
			  else
			    {
			      cout << "SU(8) not supported with more than 16 momenta" << endl;
			      Space = 0;			      
			    }
			}
		    }
		  else
		    {
		      int TmpNbrParticles[8];
		      TmpNbrParticles[0] = NbrParticlesBand1DownMinus[SymmetrySectorIndex];
		      TmpNbrParticles[1] = NbrParticlesBand2DownMinus[SymmetrySectorIndex];
		      TmpNbrParticles[2] = NbrParticlesBand1DownPlus[SymmetrySectorIndex];
		      TmpNbrParticles[3] = NbrParticlesBand2DownPlus[SymmetrySectorIndex];
		      TmpNbrParticles[4] = NbrParticlesBand1UpMinus[SymmetrySectorIndex];
		      TmpNbrParticles[5] = NbrParticlesBand2UpMinus[SymmetrySectorIndex];
		      TmpNbrParticles[6] = NbrParticlesBand1UpPlus[SymmetrySectorIndex];
		      TmpNbrParticles[7] = NbrParticlesBand2UpPlus[SymmetrySectorIndex];
		      if ((NbrSitesX * NbrSitesY) <= 8)
			{
			  Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], TmpNbrParticles, 10000000ul);
			}
		      else
			{
			  if ((NbrSitesX * NbrSitesY) <= 16)
			    {			  
			      Space = new FermionOnSquareLatticeWithSU8SpinMomentumSpaceLong(NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], TmpNbrParticles, 10000000ul);
			    }
			  else
			    {
			      cout << "SU(8) not supported with more than 16 momenta" << endl;
			      Space = 0;			      
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY, KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
	}
      cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
      TotalDim += Space->GetHilbertSpaceDimension();
      if (Space->GetHilbertSpaceDimension() > 0)
	{
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
	  
	  if (Manager.GetBoolean("real-interaction"))
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
												Manager.GetString("interaction-file"),
												TightBindingModel, Manager.GetBoolean("flat-band"), 
												Manager.GetDouble("interaction-rescaling"),
												Manager.GetBoolean("add-spin"),
												Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
													Manager.GetString("interaction-file"),
													TightBindingModel, Manager.GetBoolean("flat-band"), 
													Manager.GetDouble("interaction-rescaling"),
													Manager.GetBoolean("add-spin"),
													Architecture.GetArchitecture(), Memory);
		}		
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-valley") == false)
		{
		  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
											    Manager.GetString("interaction-file"),
											    TightBindingModel, Manager.GetBoolean("flat-band"), 
											    Manager.GetDouble("interaction-rescaling"),
											    Manager.GetBoolean("add-spin"),
											    Architecture.GetArchitecture(), Memory);
		}
	      else
		{
		  Hamiltonian = new ParticleOnLatticeFromFileInteractionTwoBandWithSpinHamiltonian (Space, NbrParticles, NbrSitesX, NbrSitesY,
												    Manager.GetString("interaction-file"),
												    TightBindingModel, Manager.GetBoolean("flat-band"), 
												    Manager.GetDouble("interaction-rescaling"),
												    Manager.GetBoolean("add-spin"),
												    Architecture.GetArchitecture(), Memory);
		}				
	    }
	  
	  
	  char* ContentPrefix = new char[256];
	  char* EigenstateOutputFile = new char [512];
	  char* TmpExtention = new char[256];
	  if (Manager.GetBoolean("add-valley") == false)
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      sprintf (ContentPrefix, "%d %d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d_ky_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d %d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex],
			       NbrParticlesBand2UpPlus[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d_ky_%d_bz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], (NbrParticlesBand1UpPlus[SymmetrySectorIndex] - NbrParticlesBand2UpPlus[SymmetrySectorIndex]));
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      if (UseSzMinusSzSymmetry == false)
			{
			  sprintf (ContentPrefix, "%d %d %d", SzValues[SymmetrySectorIndex], KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d %d", SzValues[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d_szsym_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   SzValues[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d %d %d %d %d", SzValues[SymmetrySectorIndex], KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex],
			       NbrParticlesBand1DownPlus[SymmetrySectorIndex], NbrParticlesBand2UpPlus[SymmetrySectorIndex], 
			       NbrParticlesBand2DownPlus[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d_ky_%d_sz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("add-spin") == false)
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      if (UsePzMinusPzSymmetry == false)
			{
			  sprintf (ContentPrefix, "%d %d %d", PzValues[SymmetrySectorIndex], KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], PzValues[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d %d", PzValues[SymmetrySectorIndex], PzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_pzsym_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], PzParityValues1[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d %d %d %d %d", PzValues[SymmetrySectorIndex], KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], NbrParticlesBand1UpPlus[SymmetrySectorIndex],
			       NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand2UpPlus[SymmetrySectorIndex], 
			       NbrParticlesBand2UpMinus[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex], PzValues[SymmetrySectorIndex]);
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("conserve-bandoccuption") == false)
		    {
		      if ((UsePzMinusPzSymmetry == false) && (UseSzMinusSzSymmetry == false))
			{
			  sprintf (ContentPrefix, "%d %d %d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex], EzValues[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_ez_%d_sz_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], EzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex]);
			}
		      else
			{
			  sprintf (ContentPrefix, "%d %d %d %d %d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex], EzValues[SymmetrySectorIndex],
				   PzParityValues1[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex],
				   KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex]);
			  sprintf (TmpExtention, "_kx_%d_ky_%d_pz_%d_ez_%d_sz_%d_pzsym_%d_szsym_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
				   PzValues[SymmetrySectorIndex], EzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex],
				   PzParityValues1[SymmetrySectorIndex], SzParityValues1[SymmetrySectorIndex]);
			}
		    }
		  else
		    {
		      sprintf (ContentPrefix, "%d %d %d %d %d %d %d %d %d %d %d %d %d", PzValues[SymmetrySectorIndex], SzValues[SymmetrySectorIndex], EzValues[SymmetrySectorIndex],
			       KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
			       NbrParticlesBand1UpPlus[SymmetrySectorIndex], NbrParticlesBand2UpPlus[SymmetrySectorIndex],
			       NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand2UpMinus[SymmetrySectorIndex],
			       NbrParticlesBand1DownPlus[SymmetrySectorIndex], NbrParticlesBand2DownPlus[SymmetrySectorIndex],
			       NbrParticlesBand1DownMinus[SymmetrySectorIndex], NbrParticlesBand2DownMinus[SymmetrySectorIndex]);
		      sprintf (TmpExtention, "_kx_%d_ky_%d_s1_%d_s2_%d_s3_%d_s4_%d_s5_%d_s6_%d_s7_%d_s8_%d", KxMomenta[SymmetrySectorIndex], KyMomenta[SymmetrySectorIndex],
			       NbrParticlesBand1UpPlus[SymmetrySectorIndex], NbrParticlesBand2UpPlus[SymmetrySectorIndex],
			       NbrParticlesBand1UpMinus[SymmetrySectorIndex], NbrParticlesBand2UpMinus[SymmetrySectorIndex],
			       NbrParticlesBand1DownPlus[SymmetrySectorIndex], NbrParticlesBand2DownPlus[SymmetrySectorIndex],
			       NbrParticlesBand1DownMinus[SymmetrySectorIndex], NbrParticlesBand2DownMinus[SymmetrySectorIndex]);
			      
		    }
		}
	    }
	  EigenstateOutputFile = ReplaceExtensionToFileName(EigenvalueOutputFile, ".dat", TmpExtention);
	  delete[] TmpExtention;
	  Hamiltonian->ShiftHamiltonian(EnergyShift);
	  if (Manager.GetBoolean("real-interaction"))
	    {
	      GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix,
				       CommentLine, EnergyShift,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	  else
	    {
	      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix,
					  CommentLine, EnergyShift,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	    }
	  FirstRunFlag = false;
	  
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
	}
      delete Space;
    }
  cout << "Total dim=" << TotalDim << endl;
  return 0;
}

