#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
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
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the reduced density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "na-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed number of particles", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new SingleDoubleOption  ('\n', "diag-precision", "convergence precision in non LAPACK mode", 1e-7);
#endif

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereEntanglementEntropy -h" << endl;
      return -1;
    }

  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterNa = Manager.GetInteger("na-eigenstate");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  int* TotalLz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphere** Spaces = 0;
  RealVector* GroundStates = 0;
  ComplexVector* ComplexGroundStates = 0;
  char** GroundStateFiles = 0;
  double* Weights =0;
  bool WeightFlag = false;
  bool ComplexFlag = false;
  bool SVDFlag = Manager.GetBoolean("use-svd");

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      Weights = new double[1];
      Weights[0] = 1.0;
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
       for (int i = 0; i < NbrSpaces; ++i)
	 {
	   GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	   strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	 }
       if (DegeneratedFile.GetNbrColumns() > 1)
	 {
	   Weights = DegeneratedFile.GetAsDoubleArray(1);
	   WeightFlag = true;
	 }
       else
	 {
	   Weights = new double[NbrSpaces];
	   for (int i = 0; i < NbrSpaces; ++i)
	     Weights[i] = 1.0;
	 }
    }
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticles, LzMax, TotalLz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == false)
	{
	  cout << GroundStateFiles[i] << " is not a fermionic state" << endl;
	  return -1;
	}
      if (((NbrParticles * LzMax) & 1) != (TotalLz[i] & 1))
	{
	  cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }


  GroundStates = new RealVector [NbrSpaces];  
  ComplexGroundStates = new ComplexVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (GroundStates[i].ReadVectorTest (GroundStateFiles[i]) == true)
	{
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	}
      else
	{
	  ComplexFlag = true;
	  if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }

  Spaces = new ParticleOnSphere* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      if (Manager.GetBoolean("huge-basis") == true)
	{
	  if (Manager.GetString("load-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  Spaces[i] = new  FermionOnSphereHaldaneHugeBasis (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    {
#ifdef __64_BITS__
	      if (LzMax <= 63)
#else
		if (LzMax <= 31)
#endif
		  {
		    Spaces[i] = new FermionOnSphere(NbrParticles, TotalLz[i], LzMax, MemorySpace);
		    if ((SymmetrizedBasis == true) && (TotalLz[i] == 0))
		      {
			FermionOnSphereSymmetricBasis TmpSpace(NbrParticles, LzMax, MemorySpace);
			RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundStates[i], *((FermionOnSphere*) Spaces[i]));
			GroundStates[i] = OutputState;
		      }
		  }
		else
#ifdef __128_BIT_LONGLONG__
		  if (LzMax <= 126)
#else
		    if (LzMax <= 62)
#endif
		      {
			Spaces[i] = new FermionOnSphereLong(NbrParticles, TotalLz[i], LzMax, MemorySpace);
			if ((SymmetrizedBasis == true) && (TotalLz[i] == 0))
			  {
			    FermionOnSphereSymmetricBasisLong TmpSpace(NbrParticles, LzMax, MemorySpace);
			    RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundStates[i], *((FermionOnSphereLong*) Spaces[i]));
			    GroundStates[i] = OutputState;
			  }
		      }
		    else
		      Spaces[i] = new FermionOnSphereUnlimited(NbrParticles, TotalLz[i], LzMax, MemorySpace);
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (Manager.GetString("reference-file") == 0)
		{
		  ReferenceState = new int[LzMax + 1];
		  for (int i = 0; i <= LzMax; ++i)
		    ReferenceState[i] = 0;
		  if (strcasecmp(Manager.GetString("reference-state"), "laughlin") == 0)
		    for (int i = 0; i <= LzMax; i += 3)
		      ReferenceState[i] = 1;
		  else
		    if (strcasecmp(Manager.GetString("reference-state"), "pfaffian") == 0)
		      for (int i = 0; i <= LzMax; i += 4)
			{
			  ReferenceState[i] = 1;
			  ReferenceState[i + 1] = 1;
			}
		    else
		      if (strcasecmp(Manager.GetString("reference-state"), "readrezayi3") == 0)
			for (int i = 0; i <= LzMax; i += 5)
			  {
			    ReferenceState[i] = 1;
			    ReferenceState[i + 1] = 1;
			    ReferenceState[i + 2] = 1;
			  }
		      else
			{
			  cout << "unknown reference state " << Manager.GetString("reference-state") << endl;
			  return -1;
			}
		}
	      else
		{
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
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereHaldaneBasis*) Spaces[i])->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		    if ((SymmetrizedBasis == true) && (TotalLz == 0))
		      {
			FermionOnSphereHaldaneSymmetricBasis TmpSpace(NbrParticles, LzMax, ReferenceState, MemorySpace);
			RealVector OutputState = TmpSpace.ConvertToHaldaneNbodyBasis(GroundStates[i], * ((FermionOnSphereHaldaneBasis*) Spaces[i]));
			GroundStates[i] = OutputState;
		      }
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
			if (Manager.GetString("save-hilbert") != 0)
			  {
			    ((FermionOnSphereHaldaneBasisLong*) Spaces[i])->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			    return 0;
			  }
			if ((SymmetrizedBasis == true) && (TotalLz == 0))
			  {
			    FermionOnSphereHaldaneSymmetricBasisLong TmpSpace(NbrParticles, LzMax, ReferenceState, MemorySpace);
			    RealVector OutputState = TmpSpace.ConvertToHaldaneNbodyBasis(GroundStates[i], * ((FermionOnSphereHaldaneBasisLong*) Spaces[i]));
			    GroundStates[i] = OutputState;
			  }
		      } 
	    }
	}
      if (((ComplexFlag == false) && (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension()))
	  || ((ComplexFlag == true) && (Spaces[i]->GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension())))
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# l_a    N    Lz    lambda" << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName;
      TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  int MeanSubsystemSize = LzMax >> 1;
  if ((LzMax & 1) != 0)
    ++MeanSubsystemSize;
  if (Manager.GetInteger("max-la") > 0)
    {
      MeanSubsystemSize = Manager.GetInteger("max-la");
      if (MeanSubsystemSize > LzMax)
	MeanSubsystemSize = LzMax;
    }
  int SubsystemSize = Manager.GetInteger("min-la");
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  
  BinomialCoefficients Coefs(MeanSubsystemSize);
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int MaxSubsystemNbrParticles = NbrParticles;
      if (MaxSubsystemNbrParticles > SubsystemSize)
	MaxSubsystemNbrParticles = SubsystemSize;
      int SubsystemNbrParticles = NbrParticles - (LzMax + 1 - SubsystemSize);
      if (SubsystemNbrParticles < 0)
	SubsystemNbrParticles = 0;
      long MaximumSize = 0;
      for (int i = SubsystemNbrParticles; i <= MaxSubsystemNbrParticles; ++i)
	MaximumSize += Coefs(SubsystemSize, i);
      
      for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
	{
	  int SubsystemTotalLz = 0;
	  int SubsystemLzMax = SubsystemSize - 1;
	  int SubsystemMaxTotalLz = (SubsystemNbrParticles * (SubsystemLzMax - SubsystemNbrParticles + 1));
	  SubsystemTotalLz = -SubsystemMaxTotalLz; 
	  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
	    {
	      if (ComplexFlag == false)
		{
		  RealSymmetricMatrix PartialDensityMatrix;
		  RealMatrix PartialEntanglementMatrix;
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      if ((TotalLz[i] - SubsystemTotalLz) <= (((LzMax + 1) * (NbrParticles - 2 *SubsystemNbrParticles)) + (SubsystemSize * SubsystemNbrParticles) - 
							      ((NbrParticles - SubsystemNbrParticles) * (NbrParticles - SubsystemNbrParticles))) &&
			  ((EigenstateFlag == false) || ((FilterNa == SubsystemNbrParticles) && (FilterLza == SubsystemTotalLz))))
			{
			  cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
			  if (SVDFlag == false)
			    {
			      RealSymmetricMatrix TmpPartialDensityMatrix = Spaces[i]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
			      if (WeightFlag == true)
				TmpPartialDensityMatrix *= Weights[i];
			      if (PartialDensityMatrix.GetNbrRow() == 0)
				PartialDensityMatrix = TmpPartialDensityMatrix;
			      else
				PartialDensityMatrix += TmpPartialDensityMatrix;
			    }
			  else
			    {
			      RealMatrix TmpPartialEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
			      if (WeightFlag == true)
				TmpPartialEntanglementMatrix*= sqrt(Weights[i]);
			      if (PartialEntanglementMatrix.GetNbrRow() == 0)
				PartialEntanglementMatrix = TmpPartialEntanglementMatrix;
			      else
				PartialEntanglementMatrix += TmpPartialEntanglementMatrix;				
			    }
			}
		    }
		  
		  if(SVDFlag == false)
		    {
		      if ((NbrSpaces > 1) && (WeightFlag == false))
			PartialDensityMatrix /= ((double) NbrSpaces);
		    }
		  else
		    {
		      if ((NbrSpaces > 1) && (WeightFlag == false))
			PartialEntanglementMatrix /= sqrt(((double) NbrSpaces));
		    }
		  
		  if ((PartialDensityMatrix.GetNbrRow() > 1) || ((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))
		    {
		      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
		      if (SVDFlag == false)
			{
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    {
			      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				  && (FilterLza == SubsystemTotalLz ))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
				  TmpDiag.SortMatrixDownOrder();
				}
			    }
			  else
			    {
			      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				  && (FilterLza == SubsystemTotalLz ))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder();
				}
			    }
#else
			  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
			      && (FilterLza == SubsystemTotalLz ))
			    {
			      if (PartialDensityMatrix.GetNbrRow() == 1)
				{
				  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder();
				}
			      else
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   NbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			    }
			  else
			    {
			      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
			      TmpDiag.SortMatrixDownOrder();
			    }
#endif		  		
			}
		      else
			{
			  if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
			    {
			      cout << "PartialEntanglementMatrix = " << PartialEntanglementMatrix.GetNbrRow() << " x " << PartialEntanglementMatrix.GetNbrColumn() << endl;
			      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
			      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
				{
				  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
				}
			      for (int i = 0; i < TmpDimension; ++i)
				TmpValues[i] *= TmpValues[i];
			      TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
			      TmpDiag.SortMatrixDownOrder();
			    }
			  else
			    {
			      double TmpValue = 0.0;
			      if (PartialEntanglementMatrix.GetNbrRow() == 1)
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
				    TmpValue += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
				}
			      else
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
				    TmpValue += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
				}
			      TmpDiag = RealDiagonalMatrix(1, 1);
			      TmpDiag[0] = TmpValue;
			    }
			}
		      
		      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			{
			  if (TmpDiag[i] > 1e-14)
			    {
			      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			      DensitySum += TmpDiag[i];
			    }
			}
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  for (int i = 0; i <TmpDiag.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			  DensityMatrixFile.close();
			}
		    }
		  else
		    if (PartialDensityMatrix.GetNbrRow() == 1)
		      {
			double TmpValue = PartialDensityMatrix(0,0);
			if (TmpValue > 1e-14)
			  {
			    EntanglementEntropy += TmpValue * log(TmpValue);
			    DensitySum += TmpValue;
			  }
			if (DensityMatrixFileName != 0)
			  {
			    ofstream DensityMatrixFile;
			    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			    DensityMatrixFile.precision(14);
			    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			    DensityMatrixFile.close();
			  }		  
		      }
		} //complex
	      else
		{
		  HermitianMatrix PartialDensityMatrix;
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      if ((TotalLz[i] - SubsystemTotalLz) <= (((LzMax + 1) * (NbrParticles - 2 *SubsystemNbrParticles)) + (SubsystemSize * SubsystemNbrParticles) - 
							      ((NbrParticles - SubsystemNbrParticles) * (NbrParticles - SubsystemNbrParticles))) &&
			  ((EigenstateFlag == false) || ((FilterNa == SubsystemNbrParticles) && (FilterLza == SubsystemTotalLz))))
			{
			  cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
			  HermitianMatrix TmpPartialDensityMatrix = Spaces[i]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[i]);
			  if (WeightFlag == true)
			    TmpPartialDensityMatrix *= Weights[i];
			  if (PartialDensityMatrix.GetNbrRow() == 0)
			    PartialDensityMatrix = TmpPartialDensityMatrix;
			  else
			    PartialDensityMatrix += TmpPartialDensityMatrix;
			}
		    }
		  if ((NbrSpaces > 1) && (WeightFlag == false))
		    PartialDensityMatrix /= ((double) NbrSpaces);
		  if (PartialDensityMatrix.GetNbrRow() > 1)
		    {
		      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		      if (LapackFlag == true)
			{
			  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
			      && (FilterLza == SubsystemTotalLz ))
			    {
			      ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							   PartialDensityMatrix.GetNbrRow(), true);
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				TmpEigenstates[i][i] = 1.0;
			      PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
			      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
			      char* TmpEigenstateName = new char[512];
			      int MaxNbrEigenstates = NbrEigenstates;
			      if (NbrEigenstates == 0)
				MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
			      for (int i = 0; i < MaxNbrEigenstates; ++i)
				{
				  if (TmpDiag[i] > 1e-14)
				    {
				      sprintf (TmpEigenstateName,
					       "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
					       NbrParticles, LzMax, TotalLz[0], SubsystemSize,
					       SubsystemNbrParticles, SubsystemTotalLz, i);
				      TmpEigenstates[i].WriteVector(TmpEigenstateName);
				    }
				}
			      delete[] TmpEigenstateName;
			    }
			  else
			    {
			      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
			      TmpDiag.SortMatrixDownOrder();
			    }
			}
		      else
			{
			  if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
			      && (FilterLza == SubsystemTotalLz ))
			    {
			      ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							   PartialDensityMatrix.GetNbrRow(), true);
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				TmpEigenstates[i][i] = 1.0;
			      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
			      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
			      char* TmpEigenstateName = new char[512];
			      int MaxNbrEigenstates = NbrEigenstates;
			      if (NbrEigenstates == 0)
				MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
			      for (int i = 0; i < MaxNbrEigenstates; ++i)
				{
				  if (TmpDiag[i] > 1e-14)
				    {
				      sprintf (TmpEigenstateName,
						   "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
					       NbrParticles, LzMax, TotalLz[0], SubsystemSize,
					       SubsystemNbrParticles, SubsystemTotalLz, i);
				      TmpEigenstates[i].WriteVector(TmpEigenstateName);
				    }
				}
			      delete[] TmpEigenstateName;
			    }
			  else
			    {
			      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
			      TmpDiag.SortMatrixDownOrder();
			    }
			}
#else
		      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
			  && (FilterLza == SubsystemTotalLz ))
			{
			  if (PartialDensityMatrix.GetNbrRow() == 1)
			    {
			      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
			      TmpDiag.SortMatrixDownOrder();
			    }
			  else
			    {
			      ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							   PartialDensityMatrix.GetNbrRow(), true);
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				TmpEigenstates[i][i] = 1.0;
			      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
			      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
			      char* TmpEigenstateName = new char[512];
			      int MaxNbrEigenstates = NbrEigenstates;
			      if (NbrEigenstates == 0)
				MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
			      for (int i = 0; i < MaxNbrEigenstates; ++i)
				{
				  if (TmpDiag[i] > 1e-14)
				    {
				      sprintf (TmpEigenstateName,
						   "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
					       NbrParticles, LzMax, TotalLz[0], SubsystemSize,
					       SubsystemNbrParticles, SubsystemTotalLz, i);
				      TmpEigenstates[i].WriteVector(TmpEigenstateName);
				    }
				}
			      delete[] TmpEigenstateName;
			    }
			}
		      else
			{
			  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
			  TmpDiag.SortMatrixDownOrder();
			}
#endif		  			
		      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			{
			  if (TmpDiag[i] > 1e-14)
			    {
			      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
			      DensitySum += TmpDiag[i];
			    }
			}
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			  DensityMatrixFile.close();
			}
		    }
		  else
		    if (PartialDensityMatrix.GetNbrRow() == 1)
		      {
			double TmpValue = PartialDensityMatrix(0,0);
			if (TmpValue > 1e-14)
			  {
			    EntanglementEntropy += TmpValue * log(TmpValue);
			    DensitySum += TmpValue;
			  }
			if (DensityMatrixFileName != 0)
			  {
				ofstream DensityMatrixFile;
				DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
				DensityMatrixFile.precision(14);
				DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
				DensityMatrixFile.close();
			  }		  
		      }
		}
	    }
	}
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
    }
  File.close();
}

