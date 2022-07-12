#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "Hamiltonian/ExplicitHamiltonian.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereWithSU2SpinEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  
  (*SystemGroup) += new BooleanOption ('\n', "single-sza", "focus on a single Sza sector");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-sza", "if --single-sza is used, provides twice the Sza value", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption ('\n', "no-sz", "indicates that the input states are not Sz eigenstates");
  (*SystemGroup) += new BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
  (*SystemGroup) += new SingleStringOption  ('\n', "realspace-cut", "use real space partition instead of particle partition, providing the orbital weights");
  (*SystemGroup) += new BooleanOption  ('\n', "symbreak-patch", "when using real space partition, assume a patch that breaks all spatial symmetries");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
  (*SystemGroup) += new SingleStringOption ('\n', "selected-blocks", "provide a column formatted ascii file that indicates which block of the reduced density matrix should be computed");
  (*SystemGroup) += new BooleanOption ('\n', "partial-es", "only compute the first few entanglement energues per quantum number sector");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "sza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Sz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*OutputGroup) += new BooleanOption ('\n', "export-densitymatrix", "write a single block of the reduced density matrix from a file and exit");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereWithSU2SpinEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }


  int NbrParticles = 0; 
  int LzMax = 0; 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool SVDFlag = Manager.GetBoolean("use-svd");
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int FilterSza = Manager.GetInteger("sza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  bool NoSzFlag = Manager.GetBoolean("no-sz");
  bool ComplexFlag = Manager.GetBoolean("complex");
  bool RealSpaceCut = false;
  bool SymmetryBreakingPatch = false;
  if (Manager.GetString("realspace-cut") != 0)
    {
      RealSpaceCut = true;
      if (Manager.GetBoolean("symbreak-patch"))
	{
	  SymmetryBreakingPatch = true;
	}
    }
  else
    {
      if (Manager.GetBoolean("symbreak-patch"))
	{
	  cout << "error, --symbreak-patch should be used with --realspace-cut" << endl;
	  return -1;
	}
    }
  bool PartialDiagonalization = Manager.GetBoolean("partial-es");
  int* TotalLz = 0;
  int* TotalSz = 0;
  int* LzSymmetry = 0;
  int* SzSymmetry = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphereWithSpin** Spaces = 0;
  RealVector* GroundStates = 0;
  ComplexVector* ComplexGroundStates = 0;
  char** GroundStateFiles = 0;

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      TotalSz = new int[1];
      LzSymmetry = new int[1];
      SzSymmetry = new int[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       NbrSpaces = DegeneratedFile.GetNbrLines();
       GroundStateFiles = new char* [NbrSpaces];
       TotalLz = new int[NbrSpaces];
       TotalSz = new int[NbrSpaces];
       LzSymmetry = new int[NbrSpaces];
       SzSymmetry = new int[NbrSpaces];
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
    }

  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      TotalSz[i] = 0;
      LzSymmetry[i] = 0;
      SzSymmetry[i] = 0;
      if (NoSzFlag == false)
	{
	  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(GroundStateFiles[i],
								   NbrParticles, LzMax, TotalLz[i], TotalSz[i], SzSymmetry[i], LzSymmetry[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
      else
	{
	  if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
							   NbrParticles, LzMax, TotalLz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }

  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("max-na") > 0)
    {
      MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
    }
  else
    {
      if (RealSpaceCut == true)
        MaxSubsystemNbrParticles = NbrParticles;
    }
  int MinSubsystemNbrParticles = Manager.GetInteger("min-na");
  if ((RealSpaceCut == true) && (MinSubsystemNbrParticles == 1))
    {
      MinSubsystemNbrParticles = 0;
    }
  int* NbrNUp = new int[NbrSpaces];
  int* NbrNDown = new int[NbrSpaces];
  int MaxNbrNUp = 0;
  int MaxNbrNDown = 0;
  for (int i = 0; i < NbrSpaces; ++i)
    { 
      NbrNUp[i] = (NbrParticles + TotalSz[i]);
      NbrNDown[i] = (NbrParticles - TotalSz[i]);
      NbrNUp[i] >>= 1;
      NbrNDown[i] >>= 1;
      if (NbrNUp[i] > MaxNbrNUp)
	MaxNbrNUp = NbrNUp[i];
      if (NbrNDown[i] > MaxNbrNDown)
	MaxNbrNDown = NbrNDown[i];
    }
  if (NoSzFlag == true)
    {
      MaxNbrNUp = 0;
      MaxNbrNDown = 0;
   }


  int TotalNbrReducedDensityMatrixBlocks = 0;
  int* SubsystemNbrParticleSectors = 0;
  int* SubsystemTotalSzSectors = 0;
  int* SubsystemTotalLzSectors = 0; 
  int* SubsystemSzSymmetrySectors = 0;
  if (Manager.GetString("selected-blocks") == 0)
    {
      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int LocalMinNbrNUp = 0;
	  int LocalMaxNbrNUp = MaxNbrNUp;
	  if (Manager.GetBoolean("single-sza") == true)
	    {
	      LocalMinNbrNUp = (SubsystemNbrParticles + Manager.GetInteger("only-sza")) /2;
	      LocalMaxNbrNUp = LocalMinNbrNUp;
	      if (LocalMaxNbrNUp > MaxNbrNUp)
		{
		  LocalMaxNbrNUp = LocalMinNbrNUp - 1;
		}
	    }
	  for (int SubsystemNbrNUp = LocalMinNbrNUp; SubsystemNbrNUp <= LocalMaxNbrNUp; ++SubsystemNbrNUp)
	    {
	      int SubsystemNbrNDown = SubsystemNbrParticles - SubsystemNbrNUp;
	      if (((SubsystemNbrNDown >= 0) && (SubsystemNbrNDown <= MaxNbrNDown)) || (NoSzFlag == true))
		{
		  int SubsystemTotalSz = SubsystemNbrNUp - SubsystemNbrNDown;
		  
		  int SubsystemMaxTotalLz = 0;
		  if (Statistics == true)
		    {
		      SubsystemMaxTotalLz = (((SubsystemNbrNUp * LzMax) - (SubsystemNbrNUp * (SubsystemNbrNUp - 1)))
					     + ((SubsystemNbrNDown * LzMax) - (SubsystemNbrNDown * (SubsystemNbrNDown - 1))));
		    }
		  else
		    {
		      SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
		    }
		  int SubsystemTotalLz = -SubsystemMaxTotalLz; 
		  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
		    {
		      ++TotalNbrReducedDensityMatrixBlocks;
		    }
		}
	    }
	}
      SubsystemNbrParticleSectors = new int[TotalNbrReducedDensityMatrixBlocks];
      SubsystemTotalSzSectors = new int[TotalNbrReducedDensityMatrixBlocks];
      SubsystemTotalLzSectors = new int[TotalNbrReducedDensityMatrixBlocks];
      SubsystemSzSymmetrySectors = new int[TotalNbrReducedDensityMatrixBlocks];
      TotalNbrReducedDensityMatrixBlocks = 0;
      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int LocalMinNbrNUp = 0;
	  int LocalMaxNbrNUp = MaxNbrNUp;
	  if (Manager.GetBoolean("single-sza") == true)
	    {
	      LocalMinNbrNUp = (SubsystemNbrParticles + Manager.GetInteger("only-sza")) /2;
	      LocalMaxNbrNUp = LocalMinNbrNUp;
	      if (LocalMaxNbrNUp > MaxNbrNUp)
		{
		  LocalMaxNbrNUp = LocalMinNbrNUp - 1;
		}
	    }
	  for (int SubsystemNbrNUp = LocalMinNbrNUp; SubsystemNbrNUp <= LocalMaxNbrNUp; ++SubsystemNbrNUp)
	    {
	      int SubsystemNbrNDown = SubsystemNbrParticles - SubsystemNbrNUp;
	      if (((SubsystemNbrNDown >= 0) && (SubsystemNbrNDown <= MaxNbrNDown)) || (NoSzFlag == true))
		{
		  int SubsystemTotalSz = SubsystemNbrNUp - SubsystemNbrNDown;
		  
		  int SubsystemMaxTotalLz = 0;
		  if (Statistics == true)
		    {
		      SubsystemMaxTotalLz = (((SubsystemNbrNUp * LzMax) - (SubsystemNbrNUp * (SubsystemNbrNUp - 1)))
					     + ((SubsystemNbrNDown * LzMax) - (SubsystemNbrNDown * (SubsystemNbrNDown - 1))));
		    }
		  else
		    {
		      SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
		    }
		  int SubsystemTotalLz = -SubsystemMaxTotalLz; 
		  int LocalMinSzSymmetrySector = -1;
		  int LocalMaxSzSymmetrySector = 1;
		  if ((SubsystemTotalSz != 0) || (SzSymmetry[0] == 0) || (Manager.GetBoolean("disable-szsymmetry") == true))
		    {
		      LocalMinSzSymmetrySector = 0;
		      LocalMaxSzSymmetrySector = 0;
		    }
		  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
		    {
		      for (int SubsystemSzSymmetrySector = LocalMinSzSymmetrySector; SubsystemSzSymmetrySector <= LocalMinSzSymmetrySector; SubsystemSzSymmetrySector += 2)
			{
			  SubsystemNbrParticleSectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemNbrParticles;
			  SubsystemTotalSzSectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemTotalSz;
			  SubsystemTotalLzSectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemTotalLz;
			  SubsystemSzSymmetrySectors[TotalNbrReducedDensityMatrixBlocks] = SubsystemSzSymmetrySector;
			  ++TotalNbrReducedDensityMatrixBlocks;
			}
		    }
		}
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile BlockFile;
      if (BlockFile.Parse(Manager.GetString("selected-blocks")) == false)
	{
	  BlockFile.DumpErrors(cout);
	  return -1;
	}
      TotalNbrReducedDensityMatrixBlocks = BlockFile.GetNbrLines();
      SubsystemNbrParticleSectors = BlockFile.GetAsIntegerArray(0);
      MinSubsystemNbrParticles = 100000;
      MaxSubsystemNbrParticles = 0;
      for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
	{
	  if (SubsystemNbrParticleSectors[i] < MinSubsystemNbrParticles)
	    MinSubsystemNbrParticles = SubsystemNbrParticleSectors[i];
	  if (SubsystemNbrParticleSectors[i] > MaxSubsystemNbrParticles)
	    MaxSubsystemNbrParticles = SubsystemNbrParticleSectors[i];
	}
      SubsystemTotalLzSectors = BlockFile.GetAsIntegerArray(1);
      SubsystemTotalSzSectors = BlockFile.GetAsIntegerArray(2);
      if (BlockFile.GetNbrColumns() >= 4)
	{
	  SubsystemSzSymmetrySectors = BlockFile.GetAsIntegerArray(3);
	}
      else
	{
	  SubsystemSzSymmetrySectors = new int[TotalNbrReducedDensityMatrixBlocks];
	  for (int i = 0; i < TotalNbrReducedDensityMatrixBlocks; ++i)
	    {
	      SubsystemSzSymmetrySectors[i] = 0;
	    }
	}
    }

  if (ComplexFlag == false)
    {
      GroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	  {
	    cout << "can't open vector file " << GroundStateFiles[i] << endl;
	    return -1;      
	  }
    }
  else
    {
      ComplexGroundStates = new ComplexVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	  {
	    cout << "can't open vector file " << GroundStateFiles[i] << endl;
	    return -1;      
	  }
    }
  
  Spaces = new ParticleOnSphereWithSpin* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Statistics == true)
	{
	  if (NoSzFlag == false)
	    {
	      Spaces[i] = new FermionOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax, TotalSz[i]);
	    }
	  else
	    {
	      cout << "error : bosons without Sz are not yet supported" << endl;
	      return 0;
//	      Spaces[i] = new FermionOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax);
	    }
	  
	}
      else
	{
	  if (NoSzFlag == false)
	    {
	      if (Manager.GetBoolean("use-alt") == false)
		{
		  Spaces[i] = new BosonOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax, TotalSz[i]);
		}
	      else
		{
		  if (LzSymmetry[i] == 0)
		    {
		      if (SzSymmetry[i] == 0)
			{		      
			  if (Manager.GetString("load-hilbert") == 0)
			    {
			      Spaces[i] = new BosonOnSphereWithSU2Spin (NbrParticles, TotalLz[i], LzMax, TotalSz[i]);
			    }
			  else
			    {
			      Spaces[i] = new BosonOnSphereWithSU2Spin (Manager.GetString("load-hilbert"));
			    }
			}
		      else
			{		      
			  if (Manager.GetString("load-hilbert") == 0)
			    {
			      Spaces[i] = new BosonOnSphereWithSU2SpinSzSymmetry (NbrParticles, TotalLz[i], LzMax, TotalSz[i], (SzSymmetry[i] == -1));
			    }
			  else
			    {
			      Spaces[i] = new BosonOnSphereWithSU2SpinSzSymmetry (Manager.GetString("load-hilbert"));
			    }
			}

		    }
		  else
		    {
		      if (SzSymmetry[i] == 0)
			{		      
			  if (Manager.GetString("load-hilbert") == 0)
			    {
			      Spaces[i] = new BosonOnSphereWithSU2SpinLzSymmetry (NbrParticles, LzMax, TotalSz[i], (LzSymmetry[i] == -1));
			    }
			  else
			    {
			      Spaces[i] = new BosonOnSphereWithSU2SpinLzSymmetry (Manager.GetString("load-hilbert"));
			    }
			}
		      else
			{		      
			  if (Manager.GetString("load-hilbert") == 0)
			    {
			      Spaces[i] = new BosonOnSphereWithSU2SpinLzSzSymmetry (NbrParticles, LzMax, TotalSz[i], (SzSymmetry[i] == -1), (LzSymmetry[i] == -1));
			    }
			  else
			    {
			      Spaces[i] = new BosonOnSphereWithSU2SpinLzSzSymmetry (Manager.GetString("load-hilbert"));
			    }
			}
		    }
		}
	    }
	  else
	    {
	      cout << "error : bosons without Sz are not yet supported" << endl;
	      return 0;
//	      Spaces[i] = new BosonOnSphereWithSpin (NbrParticles, TotalLz[i], LzMax);
	    }
	}
      if(ComplexFlag == false)
	{
	  if (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
      else
	{
	  if (Spaces[i]->GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and ground state" << endl;
	      return 0;
	    }
	}
    }
  
  if ((DensityMatrixFileName != 0) && (Architecture.GetArchitecture()->CanWriteOnDisk()))
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      if (SymmetryBreakingPatch == false)
	{
	  if (NoSzFlag == false)
	    {
	      DensityMatrixFile << "#  N    Sz    Sz<->-Sz    Nup    Ndown    Lz    lambda";
	    }
	  else
	    {
	      DensityMatrixFile << "#  N    Lz    lambda";
	    }
	}
      else
	{
	  if (NoSzFlag == false)
	    {
	      DensityMatrixFile << "#  N    Sz    Nup    Ndown    lambda";
	    }
	  else
	    {
	      DensityMatrixFile << "#  N    lambda";
	    }
	}
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  char* OutputFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy(OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      if (RealSpaceCut == false)
	{
	 OutputFileName  = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
	}
      else
	{
	 OutputFileName  = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "realent");	  
	}
      if (OutputFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk())  
    {
      ofstream File;
      File.open(OutputFileName, ios::binary | ios::out);
      File.precision(14);
      File.close();
    }
  cout.precision(14);

  double TotalTrace = 0.0;
  double TotalEntanglementEntropy = 0.0;

  double* WeightAOrbitalsUp = 0;
  double* WeightBOrbitalsUp = 0;
  double* WeightAOrbitalsDown = 0;
  double* WeightBOrbitalsDown = 0;
  int NbrAOrbitals = LzMax + 1;
  int NbrBOrbitals = LzMax + 1;
  if ((Manager.GetString("realspace-cut") != 0) && (SymmetryBreakingPatch == false))
    {
      ConfigurationParser RealSpaceWeights;
      if (RealSpaceWeights.Parse(Manager.GetString("realspace-cut")) == false)
	{
	  RealSpaceWeights.DumpErrors(cout) << endl;
	  return -1;
	}
      double* TmpSquareWeights = 0;
      int TmpNbrOrbitals = 0;
      if (RealSpaceWeights.GetAsDoubleArray("OrbitalSquareWeights", ' ', TmpSquareWeights, TmpNbrOrbitals) == false)
	{
	  cout << "OrbitalSquareWeights is not defined or as a wrong value" << endl;
	  return -1;
	}
      if (TmpNbrOrbitals > (LzMax + 1))
	{
	  cout << "error, the number of weights (" << TmpNbrOrbitals << ") cannot exceed the number of orbitals (" << (LzMax + 1) << ")" << endl;
	  return -1;
	}
      NbrAOrbitals = (LzMax + 1 + TmpNbrOrbitals) / 2;
      WeightAOrbitalsUp = new double [NbrAOrbitals];
      WeightAOrbitalsDown = new double [NbrAOrbitals];
      for (int i = 0; i < (NbrAOrbitals - TmpNbrOrbitals); ++i)
	{
	  WeightAOrbitalsUp[i] = 1.0;
	  WeightAOrbitalsDown[i] = 1.0;
	}
      for (int i = NbrAOrbitals - TmpNbrOrbitals; i < NbrAOrbitals; ++i)
	{
	  WeightAOrbitalsUp[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	  WeightAOrbitalsDown[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	}
      NbrBOrbitals = (LzMax + 1 + TmpNbrOrbitals) / 2;
      WeightBOrbitalsUp = new double [NbrBOrbitals];
      WeightBOrbitalsDown = new double [NbrBOrbitals];
      for (int i = 0; i < TmpNbrOrbitals; ++i)
	{
	  WeightBOrbitalsUp[i] = sqrt(1.0 - TmpSquareWeights[i]);
	  WeightBOrbitalsDown[i] = sqrt(1.0 - TmpSquareWeights[i]);
	}
      for (int i = TmpNbrOrbitals; i < NbrBOrbitals; ++i)
	{
	  WeightBOrbitalsUp[i] = 1.0;
	  WeightBOrbitalsDown[i] = 1.0;
	}      
    }  


  double* EntanglementEntropies = new double[MaxSubsystemNbrParticles - MinSubsystemNbrParticles + 1];
  double* DensitySums = new double[MaxSubsystemNbrParticles - MinSubsystemNbrParticles + 1];
  for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] = 0.0;
      DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] = 0.0;
    }


  if (SymmetryBreakingPatch == true)
    {
      int* NbrConnectedOrbitalAUp = new int[NbrAOrbitals];
      int* NbrConnectedOrbitalADown = new int[NbrAOrbitals]; 
      int** ConnectedOrbitalAUp = new int*[NbrAOrbitals];
      int** ConnectedOrbitalADown = new int*[NbrAOrbitals];
      double** FullWeightAOrbitalsUp = new double*[NbrAOrbitals]; 
      double** FullWeightAOrbitalsDown = new double*[NbrAOrbitals];
      int* NbrConnectedOrbitalBUp = new int[NbrBOrbitals];
      int* NbrConnectedOrbitalBDown = new int[NbrBOrbitals]; 
      int** ConnectedOrbitalBUp = new int*[NbrBOrbitals];
      int** ConnectedOrbitalBDown = new int*[NbrBOrbitals];
      double** FullWeightBOrbitalsUp = new double*[NbrBOrbitals]; 
      double** FullWeightBOrbitalsDown = new double*[NbrBOrbitals];

      MultiColumnASCIIFile RealSpaceWeightFile;
      if (RealSpaceWeightFile.Parse(Manager.GetString("realspace-cut")) == false)
        {
          RealSpaceWeightFile.DumpErrors(cout);
          return -1;
        }
      if (RealSpaceWeightFile.GetNbrColumns() != 4)
        {
          cout << "error, wrong number of columns in " << Manager.GetString("realspace-cut") << endl;
          return -1;
        }
      int TmpNbrWeigths = RealSpaceWeightFile.GetNbrLines();
      int* TmpIndices1 = RealSpaceWeightFile.GetAsIntegerArray(0);
      if (TmpIndices1 == 0)
        {
          RealSpaceWeightFile.DumpErrors(cout);
          return -1;      
        }
      int* TmpIndices2 = RealSpaceWeightFile.GetAsIntegerArray(1);
      if (TmpIndices1 == 0)
        {
          RealSpaceWeightFile.DumpErrors(cout);
          return -1;      
        }
      double* TmpAWeights = RealSpaceWeightFile.GetAsDoubleArray(2);
      if (TmpAWeights == 0)
        {
          RealSpaceWeightFile.DumpErrors(cout);
          return -1;      
        }
      double* TmpBWeights = RealSpaceWeightFile.GetAsDoubleArray(3);
      if (TmpBWeights == 0)
        {
          RealSpaceWeightFile.DumpErrors(cout);
          return -1;      
        }
      for (int i = 0; i < NbrAOrbitals; ++i)
	{
  	  NbrConnectedOrbitalAUp[i] = 0;
	}
      
      for (int i = 0; i < TmpNbrWeigths; ++i)
	{
  	  NbrConnectedOrbitalAUp[TmpIndices1[i]]++;
	}
      for (int i = 0; i < NbrAOrbitals; ++i)
	{
	  NbrConnectedOrbitalADown[i] = NbrConnectedOrbitalAUp[i];
	  if (NbrConnectedOrbitalAUp[i] > 0)
	    {
	      ConnectedOrbitalAUp[i] = new int [NbrConnectedOrbitalAUp[i]];
	      ConnectedOrbitalADown[i] = new int [NbrConnectedOrbitalADown[i]];
	      FullWeightAOrbitalsUp[i] = new double [NbrConnectedOrbitalAUp[i]];
	      FullWeightAOrbitalsDown[i] = new double [NbrConnectedOrbitalADown[i]];
	      NbrConnectedOrbitalAUp[i] = 0;
	      NbrConnectedOrbitalADown[i] = 0;
	    }
	  else
	    {
	      ConnectedOrbitalAUp[i] = 0;
	      ConnectedOrbitalADown[i] = 0;
	      FullWeightAOrbitalsUp[i] = 0;
	      FullWeightAOrbitalsDown[i] = 0;
	    }
	}
      for (int i = 0; i < TmpNbrWeigths; ++i)
	{
	  ConnectedOrbitalAUp[TmpIndices1[i]][NbrConnectedOrbitalAUp[TmpIndices1[i]]] = TmpIndices2[i];
	  FullWeightAOrbitalsUp[TmpIndices1[i]][NbrConnectedOrbitalAUp[TmpIndices1[i]]] = TmpAWeights[i];
	  NbrConnectedOrbitalAUp[TmpIndices1[i]]++;
	  ConnectedOrbitalADown[TmpIndices1[i]][NbrConnectedOrbitalADown[TmpIndices1[i]]] = TmpIndices2[i];
	  FullWeightAOrbitalsDown[TmpIndices1[i]][NbrConnectedOrbitalADown[TmpIndices1[i]]] = TmpAWeights[i];
	  NbrConnectedOrbitalADown[TmpIndices1[i]]++;
	}
      for (int i = 0; i < NbrBOrbitals; ++i)
	{
	  NbrConnectedOrbitalBUp[i] = 0;
	}
      
      for (int i = 0; i < TmpNbrWeigths; ++i)
	{
	  NbrConnectedOrbitalBUp[TmpIndices1[i]]++;
	}
      for (int i = 0; i < NbrBOrbitals; ++i)
	{
	  NbrConnectedOrbitalBDown[i] = NbrConnectedOrbitalBUp[i];
	  if (NbrConnectedOrbitalBUp[i] > 0)
	    {
	      ConnectedOrbitalBUp[i] = new int [NbrConnectedOrbitalBUp[i]];
	      ConnectedOrbitalBDown[i] = new int [NbrConnectedOrbitalBDown[i]];
	      FullWeightBOrbitalsUp[i] = new double [NbrConnectedOrbitalBUp[i]];
	      FullWeightBOrbitalsDown[i] = new double [NbrConnectedOrbitalBDown[i]];
	      NbrConnectedOrbitalBUp[i] = 0;
	      NbrConnectedOrbitalBDown[i] = 0;
	    }
	  else
	    {
	      ConnectedOrbitalBUp[i] = 0;
	      ConnectedOrbitalBDown[i] = 0;
	      FullWeightBOrbitalsUp[i] = 0;
	      FullWeightBOrbitalsDown[i] = 0;
	    }
	}
      for (int i = 0; i < TmpNbrWeigths; ++i)
	{
	  ConnectedOrbitalBUp[TmpIndices1[i]][NbrConnectedOrbitalBUp[TmpIndices1[i]]] = TmpIndices2[i];
	  FullWeightBOrbitalsUp[TmpIndices1[i]][NbrConnectedOrbitalBUp[TmpIndices1[i]]] = TmpBWeights[i];
  	  NbrConnectedOrbitalBUp[TmpIndices1[i]]++;
	  ConnectedOrbitalBDown[TmpIndices1[i]][NbrConnectedOrbitalBDown[TmpIndices1[i]]] = TmpIndices2[i];
	  FullWeightBOrbitalsDown[TmpIndices1[i]][NbrConnectedOrbitalBDown[TmpIndices1[i]]] = TmpBWeights[i];
  	  NbrConnectedOrbitalBDown[TmpIndices1[i]]++;
	}


//       for (int i = 0; i < NbrAOrbitals; ++i)
// 	{
// 	  NbrConnectedOrbitalAUp[i] = 1;
// 	  NbrConnectedOrbitalADown[i] = 1;
// 	  ConnectedOrbitalAUp[i] = new int [NbrConnectedOrbitalAUp[i]];
// 	  ConnectedOrbitalADown[i] = new int [NbrConnectedOrbitalADown[i]];
// 	  FullWeightAOrbitalsUp[i] = new double [NbrConnectedOrbitalAUp[i]];
// 	  FullWeightAOrbitalsDown[i] = new double [NbrConnectedOrbitalADown[i]];
// 	  for (int j = 0; j < NbrConnectedOrbitalAUp[i]; ++j)
// 	    {
// 	      ConnectedOrbitalAUp[i][j] = i;
// 	      ConnectedOrbitalADown[i][j] = i;
// 	      FullWeightAOrbitalsUp[i][j] = WeightAOrbitalsUp[i];
// 	      FullWeightAOrbitalsDown[i][j] = WeightAOrbitalsDown[i];
// 	    }
// 	}
//       for (int i = 0; i < NbrBOrbitals; ++i)
// 	{
// 	  NbrConnectedOrbitalBUp[i] = 1;
// 	  NbrConnectedOrbitalBDown[i] = 1;
// 	  ConnectedOrbitalBUp[i] = new int [NbrConnectedOrbitalBUp[i]];
// 	  ConnectedOrbitalBDown[i] = new int [NbrConnectedOrbitalBDown[i]];
// 	  FullWeightBOrbitalsUp[i] = new double [NbrConnectedOrbitalBUp[i]];
// 	  FullWeightBOrbitalsDown[i] = new double [NbrConnectedOrbitalBDown[i]];
// 	  for (int j = 0; j < NbrConnectedOrbitalBUp[i]; ++j)
// 	    {
// 	      ConnectedOrbitalBUp[i][j] = i;
// 	      ConnectedOrbitalBDown[i][j] = i;
// 	      FullWeightBOrbitalsUp[i][j] = WeightBOrbitalsUp[i];
// 	      FullWeightBOrbitalsDown[i][j] = WeightBOrbitalsDown[i];
// 	    }
// 	}


      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int LocalMinNbrNUp = 0;
	  int LocalMaxNbrNUp = MaxNbrNUp;
	  if (Manager.GetBoolean("single-sza") == true)
	    {
	      LocalMinNbrNUp = (SubsystemNbrParticles + Manager.GetInteger("only-sza")) /2;
	      LocalMaxNbrNUp = LocalMinNbrNUp;
	      if (LocalMaxNbrNUp > MaxNbrNUp)
		{
		  LocalMaxNbrNUp = LocalMinNbrNUp - 1;
		}
	    }
	  for (int SubsystemNbrNUp = LocalMinNbrNUp; SubsystemNbrNUp <= LocalMaxNbrNUp; ++SubsystemNbrNUp)
	    {
	      int SubsystemNbrNDown = SubsystemNbrParticles - SubsystemNbrNUp;
	      if (((SubsystemNbrNDown >= 0) && (SubsystemNbrNDown <= MaxNbrNDown)) || (NoSzFlag == true))
		{
		  int SubsystemTotalSz = SubsystemNbrNUp - SubsystemNbrNDown;
		  
		  int SubsystemMaxTotalLz = 0;
		  if (Statistics == true)
		    {
		      SubsystemMaxTotalLz = (((SubsystemNbrNUp * LzMax) - (SubsystemNbrNUp * (SubsystemNbrNUp - 1)))
					     + ((SubsystemNbrNDown * LzMax) - (SubsystemNbrNDown * (SubsystemNbrNDown - 1))));
		    }
		  else
		    {
		      SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax;
		    }
		  int SubsystemTotalLz = -SubsystemMaxTotalLz; 
		  RealMatrix* TmpEntanglementMatrices = new RealMatrix[SubsystemMaxTotalLz + 1];
		  int TmpIndex = 0;
		  int TmpNbrNonZeroEntanglementMatrices = 0;
		  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
		    {
		      cout << "computing PES entanglement matrix for NA=" << SubsystemNbrParticles << " NAup=" << SubsystemNbrNUp << " NAdown=" << SubsystemNbrNDown << " 2LzA=" << SubsystemTotalLz  << endl;
		      timeval SVDTotalStartingTime;
		      timeval SVDTotalEndingTime;
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(SVDTotalStartingTime), 0);
			}
		      TmpEntanglementMatrices[TmpIndex] = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
															SubsystemNbrNUp - SubsystemNbrNDown, 
															GroundStates[0], true, Architecture.GetArchitecture());
		      if (ShowTimeFlag == true)
			{
			  gettimeofday (&(SVDTotalEndingTime), 0);
			  double Dt = (double) ((SVDTotalEndingTime.tv_sec - SVDTotalStartingTime.tv_sec) + 
						((SVDTotalEndingTime.tv_usec - SVDTotalStartingTime.tv_usec) / 1000000.0));		      
			  cout << "particle entanglement matrix evaluated in " << Dt << "s" << endl;
			}
		      if (TmpEntanglementMatrices[TmpIndex].GetNbrRow() > 0)
			{
			  ++TmpNbrNonZeroEntanglementMatrices;
			}
		      ++TmpIndex;
		    }
		  RealMatrix* TmpEntanglementMatrices2 = new RealMatrix[TmpNbrNonZeroEntanglementMatrices];
		  int* TmpEntanglementMatrixLzSectors = new int [TmpNbrNonZeroEntanglementMatrices];
		  TmpNbrNonZeroEntanglementMatrices = 0;
		  TmpIndex = 0;
		  SubsystemTotalLz = -SubsystemMaxTotalLz; 
		  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
		    {
		      if (TmpEntanglementMatrices[TmpIndex].GetNbrRow() > 0)
			{
			  TmpEntanglementMatrices2[TmpNbrNonZeroEntanglementMatrices] = TmpEntanglementMatrices[TmpIndex];
			  TmpEntanglementMatrixLzSectors[TmpNbrNonZeroEntanglementMatrices] = SubsystemTotalLz;
			  ++TmpNbrNonZeroEntanglementMatrices;
			}
		      ++TmpIndex;		      
		    }
		  delete[] TmpEntanglementMatrices;
		  TmpEntanglementMatrices = TmpEntanglementMatrices2;
		  cout << "computing RSES entanglement matrix for NA=" << SubsystemNbrParticles << " NAup=" << SubsystemNbrNUp << " NAdown=" << SubsystemNbrNDown << endl;
		  RealMatrix PartialEntanglementMatrix = Spaces[0]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemNbrNUp - SubsystemNbrNDown, 
																		      NbrAOrbitals, NbrConnectedOrbitalAUp, NbrConnectedOrbitalADown, 
																		      ConnectedOrbitalAUp, ConnectedOrbitalADown, 
																		      FullWeightAOrbitalsUp, FullWeightAOrbitalsDown,
																		      NbrBOrbitals, NbrConnectedOrbitalBUp, NbrConnectedOrbitalBDown, 
																		      ConnectedOrbitalBUp, ConnectedOrbitalBDown,
																		      FullWeightBOrbitalsUp, FullWeightBOrbitalsDown,
																		      TmpNbrNonZeroEntanglementMatrices, TmpEntanglementMatrixLzSectors, TmpEntanglementMatrices);
		  timeval TotalStartingTime;
		  timeval TotalEndingTime;
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalStartingTime), 0);
		    }
//		  cout << PartialEntanglementMatrix << endl;
		  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "singular value decomposition done in " << Dt << "s" << endl;
		    }
		  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		    {
		      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		    }
		  for (int i = 0; i < TmpDimension; ++i)
		    {
		      TmpValues[i] *= TmpValues[i];
		    }
		  RealDiagonalMatrix TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);  
		  
		  TmpDiag.SortMatrixDownOrder();
		  if ((DensityMatrixFileName != 0) && (Architecture.GetArchitecture()->CanWriteOnDisk()))
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      if (NoSzFlag == false)
			{
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemNbrNUp << " " 
					      <<  SubsystemNbrNDown << " " << TmpDiag[i] << endl;
			}
		      else
			{
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << TmpDiag[i] << endl;
			}
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] +=TmpDiag[i];
			}
		    }
		}
	    }
	}
    }
  else
    {
      for (int BlockIndex = 0; BlockIndex < TotalNbrReducedDensityMatrixBlocks; ++BlockIndex)
	{
	  int SubsystemNbrParticles = SubsystemNbrParticleSectors[BlockIndex];
	  int SubsystemTotalSz = SubsystemTotalSzSectors[BlockIndex];
	  int SubsystemTotalLz = SubsystemTotalLzSectors[BlockIndex];
	  int SubsystemSzSymmetrySector = SubsystemSzSymmetrySectors[BlockIndex];
	  int SubsystemNbrNUp = (SubsystemNbrParticles + SubsystemTotalSz) / 2;
	  int SubsystemNbrNDown = (SubsystemNbrParticles - SubsystemTotalSz) / 2;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  if (NoSzFlag == false)
	    {
	      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Nup=" << SubsystemNbrNUp 
		   << " subsystem total Ndown=" << SubsystemNbrNDown << " subsystem total Lz=" << SubsystemTotalLz;
	      if (SubsystemSzSymmetrySector != 0)
		{
		  cout << " subsystem Sz<->Sz =" << SubsystemSzSymmetrySector << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
	    }
	  RealSymmetricMatrix PartialDensityMatrix;
	  HermitianMatrix ComplexPartialDensityMatrix;
	  RealMatrix PartialEntanglementMatrix; 
	  ComplexMatrix ComplexPartialEntanglementMatrix; 
	  
	  for (int i = 0; i < NbrSpaces; ++i)
	    {
	      int ComplementaryNbrNUp = NbrNUp[i] - SubsystemNbrNUp;
	      int ComplementaryNbrNDown = NbrNDown[i] - SubsystemNbrNDown;
	      int ComplementaryMaxTotalLz = 0;
	      if (Statistics == true)
		{
		  ComplementaryMaxTotalLz = (((ComplementaryNbrNUp * LzMax) - (ComplementaryNbrNUp * (ComplementaryNbrNUp - 1)))
					     + ((ComplementaryNbrNDown * LzMax) - ( ComplementaryNbrNDown * ( ComplementaryNbrNDown - 1))));
		}
	      else
		{
		  ComplementaryMaxTotalLz = (ComplementaryNbrNUp + ComplementaryNbrNDown) * LzMax;
		}
	      
	      if ((SubsystemNbrNUp <= NbrNUp[i]) && (SubsystemNbrNDown <= NbrNDown[i]) && (abs(TotalLz[i] - SubsystemTotalLz) <= ComplementaryMaxTotalLz ))
		{
		  RealSymmetricMatrix TmpMatrix;
		  HermitianMatrix ComplexTmpMatrix;
		  
		  if (Statistics == true)
		    {
		      if (ComplexFlag == false)
			{
			  if ((SVDFlag == true) || (RealSpaceCut == true) || (PartialDiagonalization == true))
			    {
			      if (RealSpaceCut == true)
				{
				  PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,  SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i] , true);
				  if(PartialEntanglementMatrix.GetNbrRow() != 0)
				    {
				      Spaces[i]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, SubsystemNbrNUp - SubsystemNbrNDown ,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), PartialEntanglementMatrix);
				      if ((SVDFlag == false) && (PartialDiagonalization == false))
					{
					  if (PartialEntanglementMatrix.GetNbrRow() >= PartialEntanglementMatrix.GetNbrColumn())
					    {
					      PartialDensityMatrix = RealSymmetricMatrix(PartialEntanglementMatrix);
					    }
					  else
					    {
					      PartialDensityMatrix = RealSymmetricMatrix(PartialEntanglementMatrix, true);
					    }
					  PartialEntanglementMatrix = RealMatrix();
					}
				    }
				}
			      else
				{
				  PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz , 
															    SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i] , false);
				}
			    }
			  else
			    {
			      TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
												   SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i]);
			    }
			}
		    }
		  else
		    {
		      if (ComplexFlag == false)
			{
			  if ((SVDFlag == true) || (RealSpaceCut == true) || (PartialDiagonalization == true))
			    {
			      if (RealSpaceCut == true)
				{
				  timeval SVDTotalStartingTime;
				  timeval SVDTotalEndingTime;
				  if (ShowTimeFlag == true)
				    {
				      gettimeofday (&(SVDTotalStartingTime), 0);
				    }
				  if (SubsystemSzSymmetrySector == 0)
				    {
				      PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
																SubsystemNbrNUp - SubsystemNbrNDown, 
																GroundStates[i] , true, Architecture.GetArchitecture());
				    }
				  else
				    {
				      PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
																SubsystemNbrNUp - SubsystemNbrNDown, SubsystemSzSymmetrySector,
																GroundStates[i] , true, Architecture.GetArchitecture());
				    }
				  if (ShowTimeFlag == true)
				    {
				      gettimeofday (&(SVDTotalEndingTime), 0);
				      double Dt = (double) ((SVDTotalEndingTime.tv_sec - SVDTotalStartingTime.tv_sec) + 
							    ((SVDTotalEndingTime.tv_usec - SVDTotalStartingTime.tv_usec) / 1000000.0));		      
				      cout << "particle entanglement matrix evaluated in " << Dt << "s" << endl;
				    }
				  if(PartialEntanglementMatrix.GetNbrRow() != 0)
				    {
				      if (ShowTimeFlag == true)
					{
					  gettimeofday (&(SVDTotalStartingTime), 0);
					}
				      if (SubsystemSzSymmetrySector == 0)
					{
					  Spaces[i]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, 
																       SubsystemNbrNUp - SubsystemNbrNDown, 
																       NbrAOrbitals, WeightAOrbitalsUp, WeightAOrbitalsDown,
																       NbrBOrbitals, WeightBOrbitalsUp, WeightBOrbitalsDown, 
																       PartialEntanglementMatrix);
					}
				      else
					{
					  Spaces[i]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, 
																       SubsystemNbrNUp - SubsystemNbrNDown, 
																       SubsystemSzSymmetrySector,
																       NbrAOrbitals, WeightAOrbitalsUp, WeightAOrbitalsDown,
																       NbrBOrbitals, WeightBOrbitalsUp, WeightBOrbitalsDown, 
																       PartialEntanglementMatrix);
					}
				      if (ShowTimeFlag == true)
					{
					  gettimeofday (&(SVDTotalEndingTime), 0);
				      double Dt = (double) ((SVDTotalEndingTime.tv_sec - SVDTotalStartingTime.tv_sec) + 
							    ((SVDTotalEndingTime.tv_usec - SVDTotalStartingTime.tv_usec) / 1000000.0));		      
				      cout << "conversion from particle entanglement matrix to real space entanglement matrix done in " << Dt << "s" << endl;
					}					      
				      
				      if ((SVDFlag == false) || (PartialDiagonalization == true))
					{
					  if (ShowTimeFlag == true)
					    {
					      gettimeofday (&(SVDTotalStartingTime), 0);
					    }
					  if (PartialEntanglementMatrix.GetNbrRow() >= PartialEntanglementMatrix.GetNbrColumn())
					    {
					      PartialDensityMatrix = RealSymmetricMatrix(PartialEntanglementMatrix, Architecture.GetArchitecture());
					    }
					  else
					    {
					      PartialDensityMatrix = RealSymmetricMatrix(PartialEntanglementMatrix, true);
					    }
					  PartialEntanglementMatrix = RealMatrix();
					  if (ShowTimeFlag == true)
					    {
					      gettimeofday (&(SVDTotalEndingTime), 0);
					      double Dt = (double) ((SVDTotalEndingTime.tv_sec - SVDTotalStartingTime.tv_sec) + 
								    ((SVDTotalEndingTime.tv_usec - SVDTotalStartingTime.tv_usec) / 1000000.0));		      
					      cout << "conversion from entanglement matrix to reduced density matrix done in " << Dt << "s" << endl;
					    }					      
					}
				      
				    }
				}
			      else
				{
				  if (SubsystemSzSymmetrySector == 0)
				    {
				      PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz , 
																SubsystemNbrNUp - SubsystemNbrNDown, GroundStates[i],
																false, Architecture.GetArchitecture());
				    }
				  else
				    {
				      PartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz , 
																SubsystemNbrNUp - SubsystemNbrNDown, SubsystemSzSymmetrySector,
																GroundStates[i],
																false, Architecture.GetArchitecture());
				    }
				}
			    }
			  else
			    {
			      if (NoSzFlag == false)
				{
				  if (ComplexFlag == false)
				    TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
													 SubsystemNbrNUp, SubsystemNbrNDown, GroundStates[i], Architecture.GetArchitecture());
				  //else
				  //ComplexTmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,   SubsystemNbrNUp, SubsystemNbrNDown, ComplexGroundStates[i], Architecture.GetArchitecture());
				}
			      else
				{
				  if (ComplexFlag == false)
				    TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 
													 GroundStates[i], Architecture.GetArchitecture());
				  else
				    ComplexTmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,  ComplexGroundStates[i], Architecture.GetArchitecture());
				  
				}
			    }
			}
		    }
		  if (SVDFlag == false)
		    {
		      if (ComplexFlag == false)
			{
			  if (PartialDensityMatrix.GetNbrRow() > 0)
			    PartialDensityMatrix += TmpMatrix;
			  else
			    PartialDensityMatrix = TmpMatrix;
			}
		      else
			{
			  if (ComplexPartialDensityMatrix.GetNbrRow() > 0)
			    ComplexPartialDensityMatrix += ComplexTmpMatrix;
			  else
			    ComplexPartialDensityMatrix = ComplexTmpMatrix;
			}
		    }
		}
	    }
	  if (ComplexFlag == false)
	    {
	      if (NbrSpaces > 1)
		{ 
		  if (SVDFlag == true)
		    {
		      cout <<"SVD is not possible for more than one state"<<endl;
		      return -1;
		    }
		  else
		    PartialDensityMatrix /= ((double) NbrSpaces);
		}
	    }
	  else
	    {
	      if (NbrSpaces > 1)
		ComplexPartialDensityMatrix /= ((double) NbrSpaces);
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
	    }
	  if (ComplexFlag == false)
	    {
	      if ((PartialDensityMatrix.GetNbrRow() > 1)||(PartialEntanglementMatrix.GetNbrRow() >= 1))
		{
		  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
		  if (SVDFlag == false)
		    {
		      if (Manager.GetBoolean("export-densitymatrix") == true)
			{
			  if (Architecture.GetArchitecture()->CanWriteOnDisk())
			    {
			      char* TmpBlockFileName = 0;
			      if (DensityMatrixFileName == 0)
				{
				  TmpBlockFileName = new char[256];
				  sprintf (TmpBlockFileName, "densitymatrix_na_%d_lza_%d_sza_%d_szasym_%d.mat", SubsystemNbrParticles, SubsystemTotalLz, 
					   SubsystemTotalSz, SubsystemSzSymmetrySector);		  
				}
			      else
				{
				  TmpBlockFileName = new char[strlen(DensityMatrixFileName) + 256];
				  sprintf (TmpBlockFileName, "%s_na_%d_lza_%d_sza_%d_szasym_%d.mat", DensityMatrixFileName, SubsystemNbrParticles, SubsystemTotalLz, 
					   SubsystemTotalSz, SubsystemSzSymmetrySector);
				}
			      if (PartialDensityMatrix.WriteMatrix(TmpBlockFileName) == false)
				{
				  cout << "error, can't write the reduced density matrix block " << TmpBlockFileName << endl;
				  return -1;
				}
			      delete[] TmpBlockFileName;
			    }
			  return 0;
			}
		      char* TmpEigenstateNamePrefix = new char[512];
		      if (NoSzFlag == false)
			{
			  if (Statistics == true)
			    {
			      sprintf (TmpEigenstateNamePrefix,
				       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d",
				       NbrParticles, LzMax, TotalLz[0], 
				       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz);
			    }
			  else
			    {
			      sprintf (TmpEigenstateNamePrefix,
				       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d",
				       NbrParticles, LzMax, TotalLz[0], 
				       SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz);
			    }
			}
		      else
			{
			  if (Statistics == true)
			    {
			      sprintf (TmpEigenstateNamePrefix,
				       "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d",
				       NbrParticles, LzMax, TotalLz[0], 
				       SubsystemNbrParticles, SubsystemTotalLz);
			    }
			  else
			    {
			      sprintf (TmpEigenstateNamePrefix,
				       "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d",
				       NbrParticles, LzMax, TotalLz[0], 
				       SubsystemNbrParticles, SubsystemTotalLz);
			    }
			}
		      if (PartialDiagonalization == true)
			{
			  UndescribedHilbertSpace* DummyHilbertSpace = new UndescribedHilbertSpace(PartialDensityMatrix.GetNbrRow());
			  PartialDensityMatrix *= -1.0;
			  ExplicitHamiltonian* Hamiltonian = new ExplicitHamiltonian(DummyHilbertSpace, &PartialDensityMatrix);
			  GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, "", "", 0.0, 
						   0, true, TmpEigenstateNamePrefix);
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());		      
			  TmpDiag = Task.GetEigenvalues();
			  TmpDiag *= -1.0;
			  delete Hamiltonian;
			  delete DummyHilbertSpace;
			}
		      else
			{
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    {
			      if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  char* TmpEigenstateName = new char[64 + strlen(TmpEigenstateNamePrefix)];
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf(TmpEigenstateName, "%s.%d.vec", TmpEigenstateNamePrefix, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
				}
			    }
			  else
			    {
			      if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[64 + strlen(TmpEigenstateNamePrefix)];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf(TmpEigenstateName, "%s.%d.vec", TmpEigenstateNamePrefix, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.Diagonalize(TmpDiag);
				}
			    }
#else
			  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif 
			}
		    }
		  else
		    {
		      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			{
			  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			}
		      for (int i = 0; i < TmpDimension; ++i)
			{
			  TmpValues[i] *= TmpValues[i];
			}
		      TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);  
		    }
		  
		  TmpDiag.SortMatrixDownOrder();
		  if ((DensityMatrixFileName != 0) &&(Architecture.GetArchitecture()->CanWriteOnDisk()))
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      if (NoSzFlag == false)
			{
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemSzSymmetrySector << " " << SubsystemNbrNUp << " " 
					      <<  SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			}
		      else
			{
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			}
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] +=TmpDiag[i];
			}
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "diagonalization done in " << Dt << "s" << endl;
		    }
		}
	      else
		{
		  if (PartialDensityMatrix.GetNbrRow() == 1)
		    {
		      double TmpValue = PartialDensityMatrix(0,0);
		      if ((DensityMatrixFileName != 0) &&(Architecture.GetArchitecture()->CanWriteOnDisk()))
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  if (NoSzFlag == false)
			    {
			      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemSzSymmetrySector << " " << SubsystemNbrNUp << " " 
						<< SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpValue << endl;
			    }
			  else
			    {
			      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			    }
			  DensityMatrixFile.close();
			}		  
		      if (TmpValue > 1e-14)
			{
			  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpValue * log(TmpValue);
			  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpValue;
			}
		    }
		}
	    }
	  else
	    {
	      if (ComplexPartialDensityMatrix.GetNbrRow() > 1)
		{
		  RealDiagonalMatrix TmpDiag (ComplexPartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    {
		      if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
			{
			  ComplexMatrix TmpEigenstates(ComplexPartialDensityMatrix.GetNbrRow(),
						       ComplexPartialDensityMatrix.GetNbrRow(), true);
			  for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
			    TmpEigenstates[i][i] = 1.0;
			  ComplexPartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
			  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
			  char* TmpEigenstateName = new char[512];
			  int MaxNbrEigenstates = NbrEigenstates;
			  if (NbrEigenstates == 0)
			    MaxNbrEigenstates = ComplexPartialDensityMatrix.GetNbrRow();
			  for (int i = 0; i < MaxNbrEigenstates; ++i)
			    {
			      if (TmpDiag[i] > 1e-14)
				{
				  if (NoSzFlag == false)
				    {
				      if (Statistics == true)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
					       NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					}
				      else
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					}
				    }
				  else
				    {
				      if (Statistics == true)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					}
				      else
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					}
				    }
				  TmpEigenstates[i].WriteVector(TmpEigenstateName);
				}
			    }
			  delete[] TmpEigenstateName;
			}
		      else
			{
			  ComplexPartialDensityMatrix.LapackDiagonalize(TmpDiag);
			}
		    }
		  else
		    {
		      if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && (FilterSza == SubsystemTotalSz))
			{
			  ComplexMatrix TmpEigenstates(ComplexPartialDensityMatrix.GetNbrRow(),
						       ComplexPartialDensityMatrix.GetNbrRow(), true);
			  for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
			    TmpEigenstates[i][i] = 1.0;
			  ComplexPartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
			  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
			  char* TmpEigenstateName = new char[512];
			  int MaxNbrEigenstates = NbrEigenstates;
			  if (NbrEigenstates == 0)
			    MaxNbrEigenstates = ComplexPartialDensityMatrix.GetNbrRow();
			  for (int i = 0; i < MaxNbrEigenstates; ++i)
			    {
			      if (TmpDiag[i] > 1e-14)
				{
				  if (NoSzFlag == false)
				    {
				      if (Statistics == true)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					}
				      else
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_sza_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalSz, SubsystemTotalLz, i);
					}
				    }
				  else
				    {
				      if (Statistics == true)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					}
				      else
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_su2_density_n_%d_2s_%d_lz_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], 
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					}
				    }
				  TmpEigenstates[i].WriteVector(TmpEigenstateName);
				}
			    }
			  delete[] TmpEigenstateName;
			}
		      else
			{
			  ComplexPartialDensityMatrix.Diagonalize(TmpDiag);
			}
		    }
#else
		  ComplexPartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		  TmpDiag.SortMatrixDownOrder();
		  if ((DensityMatrixFileName != 0) && (Architecture.GetArchitecture()->CanWriteOnDisk()))
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      if (NoSzFlag == false)
			{
			  for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemSzSymmetrySector  << " " << SubsystemNbrNUp << " " 
					      <<  SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			}
		      else
			{
			  for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			}
		      DensityMatrixFile.close();
		    }
		  for (int i = 0; i < ComplexPartialDensityMatrix.GetNbrRow(); ++i)
		    {
		      if (TmpDiag[i] > 1e-14)
			{
			  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpDiag[i] * log(TmpDiag[i]);
			  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] +=TmpDiag[i];
			}
		    }
		  if (ShowTimeFlag == true)
		    {
		      gettimeofday (&(TotalEndingTime), 0);
		      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		      cout << "diagonalization done in " << Dt << "s" << endl;
		    }
		}
	      else
		{
		  if (ComplexPartialDensityMatrix.GetNbrRow() == 1)
		    {
		      double TmpValue = ComplexPartialDensityMatrix(0,0);
		      if ((DensityMatrixFileName != 0) && (Architecture.GetArchitecture()->CanWriteOnDisk()))
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  if (NoSzFlag == false)
			    {
			      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalSz << " " << SubsystemNbrNUp << " " 
						<< SubsystemNbrNDown << " " << SubsystemTotalLz << " " << TmpValue << endl;
			    }
			  else
			    {
			      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			    }
			  DensityMatrixFile.close();
			}		  
		      if (TmpValue > 1e-14)
			{
			  EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpValue * log(TmpValue);
			  DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] += TmpValue;
			}
		    }
		} 
	    }		  
	}
    }
  if (Architecture.GetArchitecture()->CanWriteOnDisk())  
    {
      ofstream File;
      File.open(OutputFileName, ios::binary | ios::out | ios::app);
      File.precision(14);
      for (int SubsystemNbrParticles = MinSubsystemNbrParticles; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  File << SubsystemNbrParticles << " " << (-EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles]) 
	       << " " << DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles] << " " << (1.0 - DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles]) << endl;
	  if (RealSpaceCut == true)
	    {
	      TotalTrace += DensitySums[SubsystemNbrParticles - MinSubsystemNbrParticles];
	      TotalEntanglementEntropy -= EntanglementEntropies[SubsystemNbrParticles - MinSubsystemNbrParticles];
	    }
	}
      File.close();
    }

  if ((RealSpaceCut == true) && (Architecture.GetArchitecture()->CanWriteOnDisk()))
    {
      ofstream File;
      File.open(OutputFileName, ios::binary | ios::out | ios::app);
      File.precision(14);
      cout <<" Total density matrix trace is equal to " << TotalTrace << endl;
      File << "# Total density matrix trace = " << TotalTrace << endl;
      File << "# Total entanglement entropy = " << TotalEntanglementEntropy << endl;
      File.close();
    }
}

