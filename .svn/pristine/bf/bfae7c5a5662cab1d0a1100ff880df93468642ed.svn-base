#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum values of Lz whose sectors has to be evaluated", 0);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-lza", "minimum values of Lz whose sectors has to be evaluated", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-lza", "maximum values of Lz whose sectors has to be evaluated (0 if equal to half the total system size)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleStringOption  ('\n', "states", "single column file describing states whose entanglement spectrum has to be individually computed. All states must belong to the same Hilbert space");
  (*SystemGroup) += new BooleanOption  ('\n', "multiple-density", "store the density matrix eigenvalues for multiple states (names are infered from those given in 'states', unless specified in 'density-prefix')");
  (*SystemGroup) += new SingleStringOption  ('\n', "density-prefix", "prefix for the name of the file where the density matrix eigenvalues have to be stored");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-lvalue", "compute the L value of each reduced density matrix eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "largest-lz", "only compute the largest block of the reduced density matrix (Lz=0 or 1/2)");
  (*SystemGroup) += new BooleanOption  ('\n', "positive-lz", "only compute the positive Lz sectors");

  (*SystemGroup) += new BooleanOption  ('\n', "realspace-cut", "use real space partition instead of particle partition");
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-theta-top", "inclination angle that defines the top of the real space parition (in degrees)", 0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-theta-bot", "inclination angle that defines the bottom of the real space parition (in degrees)", 90);
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-phi-range", "angle between the 2 longitudes that defines the real space parition (in degrees)", 360);
  (*SystemGroup) += new BooleanOption  ('\n', "realspace-cylinder", "use real space partition instead of particle partition in the cylinder geometry");
  (*SystemGroup) += new SingleDoubleOption ('\n', "realspace-cylindercut", "x coordinate of the cut on the cylinder", 0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "realspace-generic", "use a generic real space partition instead of particle partition (geometrical weight has to be provided through this external file)");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");

  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension)");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*OutputGroup) += new SingleStringOption ('\n', "save-matrix", "save the entanglement matrix (in SVD mode) or the reduced density matrix in a ascii file");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*OutputGroup) += new BooleanOption ('c', "complex", "compute the entanglement spectrum of a complex state");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0) && (Manager.GetString("states") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereFermionEntanglementEntropyParticlePartition -h" << endl;
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


  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
  unsigned long MemorySpace = Manager.GetInteger("fast-search") << 20;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char** DensityMatrixFileName = 0;
  if (Manager.GetString("density-matrix") != 0)
  {
    DensityMatrixFileName = new char*[1];
    DensityMatrixFileName[0] = Manager.GetString("density-matrix");
  }
  bool MultipleDensityFlag = Manager.GetBoolean("multiple-density");
  bool ComputeLValueFlag = Manager.GetBoolean("compute-lvalue");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  bool LargestLSector = Manager.GetBoolean("largest-lz");
  bool PositiveLzSectors = Manager.GetBoolean("positive-lz");
  bool RealSpaceCut = Manager.GetBoolean("realspace-cut");
  bool RealSpaceCutCylinder = Manager.GetBoolean("realspace-cylinder");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  int* TotalLz = 0;
  bool Statistics = true;
  bool SVDFlag = Manager.GetBoolean("use-svd");
  int NbrSpaces = 1;
  int NbrStates = 1;
  ParticleOnSphere** Spaces = 0;
  RealVector* GroundStates = 0;
  ComplexVector* ComplexGroundStates = 0;
  char** GroundStateFiles = 0;
  int MaxLzA = Manager.GetInteger("max-lza");
  int MinLzA = Manager.GetInteger("min-lza");
  timeval AllSectorsStartingTime;
  timeval AllSectorsEndingTime;
  if (ShowTimeFlag == true)
   {
      gettimeofday (&(AllSectorsStartingTime), 0);
    }
 
  if ((ComputeLValueFlag == true) && (Manager.GetString("density-matrix") == 0) && (MultipleDensityFlag == false))
    {
      cout << "compute-lvalue only valid when density-matrix is activated" << endl;
      return - 1;
    }
  if ((SVDFlag == false) && (Manager.GetString("states") != 0))
    {
      cout << "--states option only supported for SVD method" << endl;
      return -1;
    }
  if ((Manager.GetString("realspace-generic") == 0) && (Manager.GetString("states") != 0))
    {
      cout << "--states option only supported for real space generic cut" << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") == 0) && (Manager.GetString("states") == 0))
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if ((Manager.GetString("degenerated-groundstate") != 0) and (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false))
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
	if ((Manager.GetString("states") != 0) and (DegeneratedFile.Parse(Manager.GetString("states")) == false))
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
       
      if (Manager.GetString("degenerated-groundstate") != 0)
      {
	NbrSpaces = DegeneratedFile.GetNbrLines();
	NbrStates = NbrSpaces;
	GroundStateFiles = new char* [NbrSpaces];
	TotalLz = new int[NbrSpaces];
	for (int i = 0; i < NbrSpaces; ++i)
	  {
	    GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	    strcpy (GroundStateFiles[i], DegeneratedFile(0, i));     
	  }
      }
      else
      {
	NbrStates = DegeneratedFile.GetNbrLines();
	GroundStateFiles = new char* [NbrStates];
	TotalLz = new int[1];
	for (int i = 0; i < NbrStates; ++i)
	  {
	    GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	    strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	  }
      }
    }

  int TmpNbrEntanglementMatrices = NbrStates;
  cout << "nbr ent matrices = " << TmpNbrEntanglementMatrices << endl;
  if (NbrSpaces > 1)
    TmpNbrEntanglementMatrices = 1;
  for (int i = 0; i < NbrStates; ++i)
    {
      int TmpTotalLz = 0;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticles, LzMax, TmpTotalLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}      
      if (Statistics == false)
	{
	  cout << GroundStateFiles[i] << " is not a fermionic state" << endl;
	  return -1;
	}
      if (i < NbrSpaces)
      {
	TotalLz[i] = TmpTotalLz;
    
	if (((NbrParticles * LzMax) & 1) != (TotalLz[i] & 1))
	{
	  cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFiles[i] << endl;
	  return -1;
	}
      }
    }

  if (Manager.GetBoolean("complex") == false)
  {
    GroundStates = new RealVector [NbrStates]; 
    for (int i = 0; i < NbrStates; ++i)
      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
  }
  else
  {
    ComplexGroundStates = new ComplexVector [NbrStates];
    for (int i = 0; i < NbrStates; ++i)
    {
      if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	{
	  cout << "can't open vector file " << GroundStateFiles[i] << endl;
	  return -1;      
	}
// 	cout << ComplexGroundStates[i].Norm() << endl;
    }
  }


  Spaces = new ParticleOnSphere* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Manager.GetBoolean("haldane") == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 63)
#else
	    if (LzMax <= 31)
#endif
	      {
		Spaces[i] = new FermionOnSphere(NbrParticles, TotalLz[i], LzMax, MemorySpace);
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    Spaces[i] = new FermionOnSphereLong(NbrParticles, TotalLz[i], LzMax, MemorySpace);
		  }
		else
		  Spaces[i] = new FermionOnSphereUnlimited(NbrParticles, TotalLz[i], LzMax, MemorySpace);
	}
      else
	{
	  int* ReferenceState = 0;
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	      return -1;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return -1;     
	    }
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  Spaces[i] = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  Spaces[i] = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz[i], LzMax, ReferenceState, MemorySpace);
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (Manager.GetString("load-hilbert") != 0)
		      Spaces[i] = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		    else
		      Spaces[i] = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz[i], LzMax, ReferenceState, MemorySpace);
		  } 
	}
      if (Manager.GetBoolean("complex"))
	{
	  if (NbrSpaces == NbrStates)
	    {
	      if (Spaces[i]->GetHilbertSpaceDimension() != ComplexGroundStates[i].GetVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space and ground state" << endl;
		  return 0;
		}
	    }
	  else
	    {
	      for (int j = 0; j < NbrStates; ++j)
		{
		  if (Spaces[i]->GetHilbertSpaceDimension() != ComplexGroundStates[j].GetVectorDimension())
		    {
		      cout << "dimension mismatch between Hilbert space and ground state" << endl;
		      return 0;
		    }
		}
	    }
	}
      else
	{
	  if (NbrSpaces == NbrStates)
	    {
	      if (Spaces[i]->GetHilbertSpaceDimension() != GroundStates[i].GetVectorDimension())
		{
		  cout << "dimension mismatch between Hilbert space and ground state" << endl;
		  return 0;
		}
	    }
	  else
	    {
	      for (int j = 0; j < NbrStates; ++j)
		{
		  if (Spaces[i]->GetHilbertSpaceDimension() != GroundStates[j].GetVectorDimension())
		    {
		      cout << "dimension mismatch between Hilbert space and ground state" << endl;
		      return 0;
		    }
		}
	    }
	}
    }
  
  double* WeightAOrbitals = 0;
  double* WeightBOrbitals = 0;
  int NbrAOrbitals = LzMax + 1;
  int NbrBOrbitals = LzMax + 1;
  char* CutName = 0;
  if (Manager.GetString("realspace-generic") != 0)
    {
      ConfigurationParser RealSpaceWeights;
      if (RealSpaceWeights.Parse(Manager.GetString("realspace-generic")) == false)
	{
	  RealSpaceWeights.DumpErrors(cout) << endl;
	  return -1;
	}
      if (RealSpaceWeights["Name"] != 0)
	{
	  CutName = new char [strlen(RealSpaceWeights["Name"]) + 1];
	  strcpy (CutName, RealSpaceWeights["Name"]);
	}
      else
	{
	  CutName = new char [32];
	  sprintf (CutName, "genericcut");
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
      WeightAOrbitals = new double [NbrAOrbitals];
      for (int i = 0; i < (NbrAOrbitals - TmpNbrOrbitals); ++i)
	WeightAOrbitals[i] = 1.0;
      for (int i = NbrAOrbitals - TmpNbrOrbitals; i < NbrAOrbitals; ++i)
	{
	  WeightAOrbitals[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	}
      NbrBOrbitals = (LzMax + 1 + TmpNbrOrbitals) / 2;
      WeightBOrbitals = new double [NbrBOrbitals];
      for (int i = 0; i < TmpNbrOrbitals; ++i)
	{
	  WeightBOrbitals[i] = sqrt(1.0 - TmpSquareWeights[i]);
	}
      for (int i = TmpNbrOrbitals; i < NbrBOrbitals; ++i)
	{
	  WeightBOrbitals[i] = 1.0;
	}      
    }

  cout << "number of orbitals in A = " << NbrAOrbitals << endl;
  cout << "number of orbitals in B = " << NbrBOrbitals << endl;

  ofstream File;
  char** TmpFileName;
  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
      TmpFileName = new char* [1];
      TmpFileName[0] = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (TmpFileName[0], Manager.GetString("output-file"));
    }
  else
    {
      TmpFileName = new char* [TmpNbrEntanglementMatrices];
      for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
	{
	  if (RealSpaceCut == false)
	    {
	      TmpFileName[l] = ReplaceExtensionToFileName(GroundStateFiles[l], "vec", "partent");
	      if (TmpFileName == 0)
		{
		  cout << "no vec extension was find in " << GroundStateFiles[l] << " file name" << endl;
		  return 0;
		}
	    }
	  else
	    {
	      char* TmpExtension = new char [512];
	      if (CutName != 0)
		{
		  sprintf(TmpExtension, "%s.realent", CutName);
		}
	      else
		{
		  if (RealSpaceCutCylinder == false)
		    {
		      sprintf(TmpExtension, "_thetabot_%.6f_thetabot_%.6f.realent", Manager.GetDouble("realspace-theta-top"), 
			      Manager.GetDouble("realspace-theta-bot"));
		    }
		  else
		    {
		      sprintf(TmpExtension, "_x_%.6f.realent", Manager.GetDouble("realspace-cylindercut"));
		    }
		}
	      TmpFileName[l] = ReplaceExtensionToFileName(GroundStateFiles[l], "vec", TmpExtension);
	      if (TmpFileName == 0)
		{
		  cout << "no vec extension was find in " << GroundStateFiles[l] << " file name" << endl;
		  return 0;
		}
	    }
	}
    }
  cout.precision(14);
  
  if (MultipleDensityFlag)
  {
    DensityMatrixFileName = new char*[TmpNbrEntanglementMatrices];
    for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
    {
	DensityMatrixFileName[l] = new char[512];
	if (RealSpaceCut == false)
	{
	  cout << "Not implemented" << endl;
	}
	else
	{
	  char* TmpExtension = new char[512];
	   if (CutName != 0)
	    {
	      sprintf(TmpExtension, "%s.real.full.ent", CutName);
	    }
	  else
	    {
	      if (RealSpaceCutCylinder == false)
		{
		  sprintf(TmpExtension, "_thetabot_%.6f_thetabot_%.6f.real.full.ent", Manager.GetDouble("realspace-theta-top"), 
			  Manager.GetDouble("realspace-theta-bot"));
		}
	      else
		{
		  sprintf(TmpExtension, "_x_%.6f.real.full.ent", Manager.GetDouble("realspace-cylindercut"));
		}
	    }
	  if (Manager.GetString("density-prefix") == 0)
	  {
	    DensityMatrixFileName[l] = ReplaceExtensionToFileName(GroundStateFiles[l], "vec", TmpExtension);
	    if (DensityMatrixFileName[l] == 0)
	      {
		cout << "no vec extension was find in " << GroundStateFiles[l] << " file name" << endl;
		return 0;
	      }
	  }
	  else
	  {
	     sprintf(DensityMatrixFileName[l], "%s.%d.%s", Manager.GetString("density-prefix"), l, TmpExtension);
	  }
	}
    }
  }
  
  if ((MultipleDensityFlag) || (Manager.GetString("density-matrix") != 0))
    {
      ofstream DensityMatrixFile;
      for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
      {
	DensityMatrixFile.open(DensityMatrixFileName[l], ios::binary | ios::out); 
	DensityMatrixFile << "#  N    Lz    lambda";
	if (ComputeLValueFlag == true)
	  DensityMatrixFile << " L^2 L";
	DensityMatrixFile << endl;
	DensityMatrixFile.close();
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
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  if ((RealSpaceCut == true) && (SubsystemNbrParticles == 1))
    {
      SubsystemNbrParticles = 0;
    }

  double Ratio, Perimeter, Height;
  if (RealSpaceCutCylinder)
    {
      Ratio = Manager.GetDouble("ratio");
      cout << "Cylinder geometry"<<endl;
      Perimeter = sqrt(2.0 * M_PI * Ratio * (LzMax + 1));
      Height = Perimeter/Ratio;
      cout<<"L= "<<Perimeter<<" H= "<<Height<<endl;
    }

  double* TotalTrace = new double[TmpNbrEntanglementMatrices];
  double* TotalEntanglementEntropy = new double[TmpNbrEntanglementMatrices];
  double* EntanglementEntropy = new double[TmpNbrEntanglementMatrices];
  double* DensitySum = new double[TmpNbrEntanglementMatrices];  
  for (int i = 0; i < TmpNbrEntanglementMatrices; ++i)
  {
    TotalTrace[i] = 0.0;
    TotalEntanglementEntropy[i] = 0.0;
  }
  
  bool FirstRun = true;
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
	for (int i = 0; i < TmpNbrEntanglementMatrices; ++i)
	{
	  EntanglementEntropy[i] = 0.0;
	  DensitySum[i] = 0.0;	 
	}
      
      int ComplementarySubsystemNbrParticles = NbrParticles - SubsystemNbrParticles;
      int SubsystemMaxTotalLz = SubsystemNbrParticles * (NbrAOrbitals - 1) - (SubsystemNbrParticles * (SubsystemNbrParticles - 1));
      int ComplementaryMaxTotalLz = ComplementarySubsystemNbrParticles * (NbrBOrbitals - 1) - (ComplementarySubsystemNbrParticles * (ComplementarySubsystemNbrParticles - 1));
      cout << "SubsystemMaxTotalLz = " << SubsystemMaxTotalLz << "    ComplementaryMaxTotalLz = " << ComplementaryMaxTotalLz << endl;
      int SubsystemTotalLz = -SubsystemMaxTotalLz;
      while ((SubsystemMaxTotalLz - ComplementaryMaxTotalLz) > TotalLz[0])
	SubsystemMaxTotalLz -= 2;
      while ((SubsystemTotalLz + ComplementaryMaxTotalLz) < TotalLz[0])
	SubsystemTotalLz += 2;
	
      if ((MinLzA != -1) && (MinLzA > SubsystemTotalLz))
	SubsystemTotalLz = MinLzA;
      if ((MaxLzA != -1) && (MaxLzA < SubsystemMaxTotalLz))
	SubsystemMaxTotalLz = MaxLzA;
      
      if ((LargestLSector == true) && (RealSpaceCut == false))
	{
	  if (((LzMax * NbrParticles) & 1) == 0)
	    {
	      SubsystemTotalLz = 0;
	      SubsystemMaxTotalLz = 0;
	    }
	  else
	    {
	      SubsystemTotalLz = 1;
	      SubsystemMaxTotalLz = 1;
	    }
	}
      if (PositiveLzSectors == true)
	{
	  if (((LzMax * NbrParticles) & 1) == 0)
	    {
	      SubsystemTotalLz = 0;
	    }
	  else
	    {
	      SubsystemTotalLz = 1;
	    }
	}
      cout << "SubsystemMaxTotalLz = " << SubsystemMaxTotalLz << "    ComplementaryMaxTotalLz = " << ComplementaryMaxTotalLz << endl;
      for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
	{
	  cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  
	  RealSymmetricMatrix PartialDensityMatrix;
	  RealMatrix PartialEntanglementMatrix;
	  ComplexMatrix ComplexPartialEntanglementMatrix;
	  RealMatrix* MultiplePartialEntanglementMatrix = 0;
	  ComplexMatrix* MultipleComplexPartialEntanglementMatrix = 0;
	  
	  if (RealSpaceCut == false)
	    {
	      if (SVDFlag == false)
		{
		  PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0]);
		}
	      else
		{
		  if (Manager.GetBoolean("complex") == false)
		    {
		      PartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0],false);
		    }
		  else
		    {
		      ComplexPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[0],false);
		    }
		}
	    }
	  else //real space
	    {
	      if (SVDFlag == false)
		{
                  if (RealSpaceCutCylinder == false)
		    {
		      if (WeightAOrbitals == 0)
			{
			  PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, 
													   Manager.GetDouble("realspace-theta-top"), 
													   Manager.GetDouble("realspace-theta-bot"), 
													   Manager.GetDouble("realspace-phi-range"), GroundStates[0]);
			}
		      else
			{
			  PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixGenericRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, 
														  NbrAOrbitals, WeightAOrbitals, NbrBOrbitals, WeightBOrbitals, 
														  GroundStates[0]);
			}
		    }
                  else //cylinder
		    PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrixRealSpacePartitionCylinder(SubsystemNbrParticles, SubsystemTotalLz, Perimeter, Height, 
													     Manager.GetDouble("realspace-cylindercut"), GroundStates[0]);
		}
	      else
		{
                  if (RealSpaceCutCylinder == false)
                    {
		      if (WeightAOrbitals == 0)
			{
			  if (Manager.GetBoolean("complex") == false)
			    {
			      PartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0],true);
			      if(PartialEntanglementMatrix.GetNbrRow() != 0)
				{
				  Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), PartialEntanglementMatrix);
				}
			    }
			  else
			    {
			      ComplexPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[0],true);
			      if(ComplexPartialEntanglementMatrix.GetNbrRow() != 0)
				{
				  Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), ComplexPartialEntanglementMatrix);
				}
			    }			    
			}
		      else
			{
			  if (Manager.GetBoolean("complex") == false)
			    {
			      if (TmpNbrEntanglementMatrices == 1)
				{
				  PartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,     NbrAOrbitals, NbrBOrbitals, GroundStates[0], true);
				  // 			    cout << PartialEntanglementMatrix << endl;
				  if(PartialEntanglementMatrix.GetNbrRow() != 0)
				    {
				      Spaces[0]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, 
																   NbrAOrbitals, WeightAOrbitals, NbrBOrbitals, WeightBOrbitals, 
																   PartialEntanglementMatrix);
				      // 				cout << PartialEntanglementMatrix << endl;
				    }
				}
			      else
				{
				  MultiplePartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,     NbrAOrbitals, NbrBOrbitals, GroundStates, TmpNbrEntanglementMatrices, true);
				  // 			      cout << MultiplePartialEntanglementMatrix[0] << endl;
				  if(MultiplePartialEntanglementMatrix[0].GetNbrRow() != 0)
				    {
				      MultiplePartialEntanglementMatrix = Spaces[0]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, 
																				       NbrAOrbitals, WeightAOrbitals, NbrBOrbitals, WeightBOrbitals, 
																				       MultiplePartialEntanglementMatrix, TmpNbrEntanglementMatrices);
				      // 				cout << MultiplePartialEntanglementMatrix[0] << endl;
				    }
				}
			    }
			  else
			    {
			      if (TmpNbrEntanglementMatrices == 1)
				{
				  ComplexPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, 	    NbrAOrbitals, NbrBOrbitals, ComplexGroundStates[0], true);
				  // 			    cout << ComplexPartialEntanglementMatrix << endl;
				  if(ComplexPartialEntanglementMatrix.GetNbrRow() != 0)
				    {
				      Spaces[0]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, NbrAOrbitals, WeightAOrbitals, NbrBOrbitals, WeightBOrbitals, ComplexPartialEntanglementMatrix);
				      // 				cout << ComplexPartialEntanglementMatrix << endl;
				    }
				}
			      else
				{
				  MultipleComplexPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz,     NbrAOrbitals, NbrBOrbitals, ComplexGroundStates, TmpNbrEntanglementMatrices, true);
				  // 			      cout << MultiplePartialEntanglementMatrix[0] << endl;
				  if(MultipleComplexPartialEntanglementMatrix[0].GetNbrRow() != 0)
				    {
				      MultipleComplexPartialEntanglementMatrix = Spaces[0]->EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, 
																					      NbrAOrbitals, WeightAOrbitals, NbrBOrbitals, WeightBOrbitals, 
																					      MultipleComplexPartialEntanglementMatrix, TmpNbrEntanglementMatrices);
				      // 				cout << MultiplePartialEntanglementMatrix[0] << endl;
				    }
				}
			    }
			}
		    }
                  else //cylinder
                    {
		      if (Manager.GetBoolean("complex") == false)
			{
			PartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0],true);
			if(PartialEntanglementMatrix.GetNbrRow() != 0)
			  {
			    Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrixCylinder(SubsystemNbrParticles, SubsystemTotalLz, Perimeter, Height, Manager.GetDouble("realspace-cylindercut"), PartialEntanglementMatrix);
			  }
		      }
		      else
		      {
			cout << "Warning, cylinder not implemented for complex vectors" << endl;
// 			ComplexPartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[0],true);
// 			if(ComplexPartialEntanglementMatrix.GetNbrRow() != 0)
// 			  {
// 			    Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrixCylinder(SubsystemNbrParticles, SubsystemTotalLz, Perimeter, Height, Manager.GetDouble("realspace-cylindercut"), ComplexPartialEntanglementMatrix);
// 		      }
                     }
		    }
		}
	    }
	  
	  for (int i = 1; i < NbrSpaces; ++i)
	    {
	      RealSymmetricMatrix TmpMatrix;
	      RealMatrix TmpEntanglementMatrix;
	      ComplexMatrix ComplexTmpEntanglementMatrix;
	      if (RealSpaceCut == false)
		{
		  if (SVDFlag == false)
		    {
		      TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
		    }
		  else
		    {
		      if (Manager.GetBoolean("complex") == false)
			TmpEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
		      else
			ComplexTmpEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[i]);
		    }
		}
	      else //real space cut
		{
		  if (SVDFlag == false)
		    {
                       if (RealSpaceCutCylinder == false)		      
                           TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), GroundStates[i]);
                       else //cylinder
                           TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrixRealSpacePartitionCylinder(SubsystemNbrParticles, SubsystemTotalLz, Perimeter, Height, Manager.GetDouble("realspace-cylindercut"), GroundStates[i]);
		    }
		  else
		    {
                      if (RealSpaceCutCylinder == false)
                        {
			  if (Manager.GetBoolean("complex") == false)
			  {
			    TmpEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i],true);
			    if(PartialEntanglementMatrix.GetNbrRow() != 0)
			      {
				Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), TmpEntanglementMatrix);
			      }
			  }
			  else
			  {
			    ComplexTmpEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[i],true);
			    if(ComplexPartialEntanglementMatrix.GetNbrRow() != 0)
			      {
				Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz,Manager.GetDouble("realspace-theta-top"), Manager.GetDouble("realspace-theta-bot"), Manager.GetDouble("realspace-phi-range"), ComplexTmpEntanglementMatrix);
			      }
			  }
			}
                      else //cylinder
                       {
			   if (Manager.GetBoolean("complex") == false) 
			   {
			    TmpEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i],true);
			    if(PartialEntanglementMatrix.GetNbrRow() != 0)
			      {
				Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrixCylinder(SubsystemNbrParticles, SubsystemTotalLz, Perimeter, Height, Manager.GetDouble("realspace-cylindercut"), TmpEntanglementMatrix);
			      }
			   }
			   else
			   {
			     cout << "Waning, cylinder not implemented for complex vectors" << endl;
// 			     ComplexTmpEntanglementMatrix = Spaces[0]-	>EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[i],true);
// 			    if(ComplexPartialEntanglementMatrix.GetNbrRow() != 0)
// 			      {
// 				Spaces[0]->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrixCylinder(SubsystemNbrParticles, SubsystemTotalLz, Perimeter, Height, Manager.GetDouble("realspace-cylindercut"), ComplexTmpEntanglementMatrix);
// 			      }
			   }
                       }  
		    }
		}
	      
	      if (SVDFlag == false)
		{
		  PartialDensityMatrix += TmpMatrix;
		}
	      else
		{
		  if (Manager.GetBoolean("complex") == false)
		    PartialEntanglementMatrix += TmpEntanglementMatrix;
		  else
		    ComplexPartialEntanglementMatrix += ComplexTmpEntanglementMatrix;
		}
	    }
	  
	  if (NbrSpaces > 1)
	    {
	      if (SVDFlag == false)
		{
		  PartialDensityMatrix /= ((double) NbrSpaces);
		}
	      else
		{
		  if (Manager.GetBoolean("complex") == false)
		    PartialEntanglementMatrix /= sqrt((double) NbrSpaces);
		  else
		    ComplexPartialEntanglementMatrix /= sqrt((double) NbrSpaces);
		}
	    }
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "reduced density matrix evaluated in " << Dt << "s" << endl;
	    }
	  if ((PartialDensityMatrix.GetNbrRow() > 1) || 
	      ((SVDFlag == true) && (PartialEntanglementMatrix.GetNbrColumn() >= 1) && (PartialEntanglementMatrix.GetNbrRow() >= 1)) || 
	      ((SVDFlag == true) && (Manager.GetBoolean("complex")) && (ComplexPartialEntanglementMatrix.GetNbrColumn() >= 1) && (ComplexPartialEntanglementMatrix.GetNbrRow() >= 1)) || (SVDFlag == true && (MultiplePartialEntanglementMatrix != 0) && (MultiplePartialEntanglementMatrix[0].GetNbrColumn() >= 1) && (MultiplePartialEntanglementMatrix[0].GetNbrRow() >= 1)) || (SVDFlag == true && (MultipleComplexPartialEntanglementMatrix != 0) && (MultipleComplexPartialEntanglementMatrix[0].GetNbrColumn() >= 1) && (MultipleComplexPartialEntanglementMatrix[0].GetNbrRow() >= 1)))
	    {
	      if (Manager.GetString("save-matrix") != 0)
		{	
		  ofstream OutputDensityMatrixFile;
		  OutputDensityMatrixFile.open(Manager.GetString("save-matrix"), ios::binary | ios::out); 
		  OutputDensityMatrixFile.precision(14);
		  if (SVDFlag == false)
		    {
		      OutputDensityMatrixFile << PartialDensityMatrix;
		    }
		  else
		    {
		      if (Manager.GetBoolean("complex"))
			OutputDensityMatrixFile << ComplexPartialEntanglementMatrix;
		      else
			OutputDensityMatrixFile << PartialEntanglementMatrix;
		    }
		  OutputDensityMatrixFile.close();
		}
	      if (ShowTimeFlag == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	     RealDiagonalMatrix* TmpDiag;
	     TmpDiag = new RealDiagonalMatrix[TmpNbrEntanglementMatrices];
	     for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
		  TmpDiag[l] = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
	      
	      if (ComputeLValueFlag == false)
		{
		  if (SVDFlag == false)
		    {
#ifdef __LAPACK__
		      if (LapackFlag == true)
			PartialDensityMatrix.LapackDiagonalize(TmpDiag[0]);
		      else
			PartialDensityMatrix.Diagonalize(TmpDiag[0]);
#else
		      PartialDensityMatrix.Diagonalize(TmpDiag[0]);
#endif		  
		    }
		  else
		    {
// 		      RealMatrix TmpMat;
// 		      RealMatrix TmpMat2;
// 		      TmpMat2.Copy(PartialEntanglementMatrix);
// 		      TmpMat = PartialEntanglementMatrix.DuplicateAndTranspose();
// 		      TmpMat2.Multiply(TmpMat);
// 		      RealSymmetricMatrix TmpMat3 ((Matrix&) TmpMat2);
// 		      RealDiagonalMatrix TmpDiag2(TmpMat3.GetNbrRow());
// 		      TmpMat3.LapackDiagonalize(TmpDiag2);
// 		      TmpDiag2.SortMatrixDownOrder();

		      for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
		      {
			double* TmpValues = 0;
			int TmpDimension;
			if (Manager.GetBoolean("complex") == false)
			{
			  if (TmpNbrEntanglementMatrices == 1)
			  {
			    TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			    TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			    if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
			    {
			    TmpDimension = PartialEntanglementMatrix.GetNbrRow();
			    }
			  }
			  else
			  {
			    TmpValues = MultiplePartialEntanglementMatrix[l].SingularValueDecomposition();
			    TmpDimension = MultiplePartialEntanglementMatrix[0].GetNbrColumn();
			    if (TmpDimension > MultiplePartialEntanglementMatrix[0].GetNbrRow())
			    {
			      TmpDimension = MultiplePartialEntanglementMatrix[0].GetNbrRow();
			    }
			  }
			}
			else
			{
			  if (TmpNbrEntanglementMatrices == 1)
			  {
			    TmpValues = ComplexPartialEntanglementMatrix.SingularValueDecomposition();
			    TmpDimension = ComplexPartialEntanglementMatrix.GetNbrColumn();
			    if (TmpDimension > ComplexPartialEntanglementMatrix.GetNbrRow())
			    {
			      TmpDimension = ComplexPartialEntanglementMatrix.GetNbrRow();
			    }
			  }
			  else
			  {
			    TmpValues = MultipleComplexPartialEntanglementMatrix[l].SingularValueDecomposition();
			    TmpDimension = MultipleComplexPartialEntanglementMatrix[0].GetNbrColumn();
			    if (TmpDimension > MultipleComplexPartialEntanglementMatrix[0].GetNbrRow())
			    {
			      TmpDimension = MultipleComplexPartialEntanglementMatrix[0].GetNbrRow();
			    }
			      
			  }
			}
			for (int i = 0; i < TmpDimension; ++i)
			{
			  TmpValues[i] *= TmpValues[i];
			}
			TmpDiag[l] = RealDiagonalMatrix(TmpValues, TmpDimension);
		      }
		    }
		    for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
		      {
			TmpDiag[l].SortMatrixDownOrder();
		    
			if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName[l], ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			
			  for (int i = 0; i < TmpDiag[l].GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[l][i] << endl;
			  DensityMatrixFile.close();
			}
		      }
		  }
	      else
		{
		  if (Manager.GetBoolean("complex"))
		    cout << "Warning, not implemented for complex vectors" << endl;
		  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
					    PartialDensityMatrix.GetNbrRow(), true);
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    TmpEigenstates[i][i] = 1.0;
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    PartialDensityMatrix.LapackDiagonalize(TmpDiag[0], TmpEigenstates);
		  else
		    PartialDensityMatrix.Diagonalize(TmpDiag[0], TmpEigenstates);
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag[0], TmpEigenstates);
#endif
		  TmpDiag[0].SortMatrixDownOrder(TmpEigenstates);
		  FermionOnSphere TmpDestinationHilbertSpace(SubsystemNbrParticles, SubsystemTotalLz, LzMax);
		  ParticleOnSphereSquareTotalMomentumOperator OperMomentum (&TmpDestinationHilbertSpace, LzMax);
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName[0], ios::binary | ios::out | ios::app); 
		  DensityMatrixFile.precision(14);
		  char* TmpEigenstateName = new char[512];
		  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		    {
		      double TmpSqrMomentum = (OperMomentum.MatrixElement(TmpEigenstates[i], TmpEigenstates[i])).Re;
		      double TmpMomentum = 0.0;
		      if (TmpSqrMomentum > 0.0)
			TmpMomentum = (0.5 * (sqrt ((4.0 * TmpSqrMomentum) + 1.0) - 1.0));
		      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[0][i] << " " << TmpSqrMomentum << " " << TmpMomentum << endl;
		      if ((EigenstateFlag == true) && (FilterLza == SubsystemTotalLz) && ((NbrEigenstates == 0) || (NbrEigenstates > i)))
			{
			  sprintf (TmpEigenstateName,
				   "fermions_particlereduceddensity_na_%d_n_%d_2s_%d_lz_%d.%d.vec",
				   SubsystemNbrParticles, NbrParticles, LzMax, SubsystemTotalLz, i);
			  TmpEigenstates[i].WriteVector(TmpEigenstateName);
			}
		    }
		  delete[] TmpEigenstateName;
		  DensityMatrixFile.close();
		}
		
		for (int l = 0; l < TmpNbrEntanglementMatrices; ++l)
		{
		  for (int i = 0; i < TmpDiag[l].GetNbrRow(); ++i)
		    {
		      if (TmpDiag[l][i] > 1e-14)
			{
			  EntanglementEntropy[l] += TmpDiag[l][i] * log(TmpDiag[l][i]);
			  DensitySum[l] +=TmpDiag[l][i];
			}
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
	      if ((SVDFlag == false) && (PartialDensityMatrix.GetNbrRow() == 1))
		{
		  double TmpValue = PartialDensityMatrix(0,0);
		  if (DensityMatrixFileName != 0)
		    {
		      ofstream DensityMatrixFile;
		      DensityMatrixFile.open(DensityMatrixFileName[0], ios::binary | ios::out | ios::app); 
		      DensityMatrixFile.precision(14);
		      if (ComputeLValueFlag == false)
			{
			  DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			}
		      else		      
			{
			  if (SubsystemNbrParticles == 1)
			    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << " " << ((LzMax * (LzMax + 2)) / 4.0) << " " << (LzMax / 2.0) << endl;
			  else
			    {
			      FermionOnSphere TmpDestinationHilbertSpace(SubsystemNbrParticles, SubsystemTotalLz, LzMax);
			      ParticleOnSphereSquareTotalMomentumOperator OperMomentum (&TmpDestinationHilbertSpace, LzMax);
			      RealVector TmpEigenstate(1);
			      TmpEigenstate[0] = 1.0;
			      double TmpSqrMomentum = (OperMomentum.MatrixElement(TmpEigenstate, TmpEigenstate)).Re;
			      double TmpMomentum = 0.0;
			      if (TmpSqrMomentum > 0.0)
				TmpMomentum = (0.5 * (sqrt ((4.0 * TmpSqrMomentum) + 1.0) - 1.0));
			      DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << " " << TmpSqrMomentum << " " << TmpMomentum << endl;
			    }
			}
		      DensityMatrixFile.close();
		    }		  
		  if (TmpValue > 1e-14)
		    {
		      EntanglementEntropy[0] += TmpValue * log(TmpValue);
		      DensitySum[0] += TmpValue;
		    }
		}
	    }
	  if (MultiplePartialEntanglementMatrix != 0)
	    delete[] MultiplePartialEntanglementMatrix;
	  if (MultipleComplexPartialEntanglementMatrix != 0)
	    delete[] MultipleComplexPartialEntanglementMatrix;
	  if (TmpNbrEntanglementMatrices == 1)
	    cout << "DensitySum = " << DensitySum[0] << endl;
	}
      for (int i = 0; i < TmpNbrEntanglementMatrices; ++i)
      {
	if (Manager.GetString("output-file") == 0)
	  if (FirstRun)
	    File.open(TmpFileName[i], ios::binary | ios::out);
	  else
	    File.open(TmpFileName[i], ios::binary | ios::out | ios::app);
	  
	File.precision(14);
	File << SubsystemNbrParticles << " " << (-EntanglementEntropy[i]) << " " << DensitySum[i] << " " << (1.0 - DensitySum[i]) << endl;
	File.close();
// 	cout << "trace = " << DensitySum[i] << endl;
	TotalEntanglementEntropy[i] += (-EntanglementEntropy[i]);
	TotalTrace[i] += DensitySum[i];
      }
      FirstRun = false;
    }
  if (RealSpaceCut == true)
    {
      for (int i = 0; i < TmpNbrEntanglementMatrices; ++i)
      {
	cout << "Entanglement entropy = " << TotalEntanglementEntropy[i] << endl;
	cout << "Total trace = " << TotalTrace[i] << endl;
	if (Manager.GetString("output-file") == 0)
	  File.open(TmpFileName[i], ios::binary | ios::out | ios::app);	  
	File.precision(14);
	File << "# Entanglement entropy = " << TotalEntanglementEntropy[i] << endl;
	File << "# Total trace = " << TotalTrace[i] << endl;
	File.close();
      }
   }
  if (GroundStates != 0)
    delete[] GroundStates;
  if (ComplexGroundStates != 0)
    delete[] ComplexGroundStates;
  delete[] Spaces;
  delete[] TmpFileName;
  if (DensityMatrixFileName != 0)
    delete[] DensityMatrixFileName;
  if (ShowTimeFlag == true)
   {
      gettimeofday (&(AllSectorsEndingTime), 0);
      double Dt = (double) ((AllSectorsEndingTime.tv_sec - AllSectorsStartingTime.tv_sec) + 
				    ((AllSectorsEndingTime.tv_usec - AllSectorsStartingTime.tv_usec) / 1000000.0));		      
      cout << "total calculation time = " << Dt << "s" << endl;
    }
  return 0;
}

