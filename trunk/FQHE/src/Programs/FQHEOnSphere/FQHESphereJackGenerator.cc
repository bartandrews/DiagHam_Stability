#include "Vector/RealVector.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneLargeBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// export the Jack polynomial in a text file
//
// outputTxtFileName = output file name
// outputState = Jack polynomial decomposition
// initialSpace = pointer to the Hilbert space 
// manager = reference to the option manager
// nbrParticles = number of particles
// alphaNumerator= numerator of the alpha factor
// alphaDenominator= denominator of the alpha factor
// statistics = particle statistics (true for fermions)
void FQHESphereJackTxtExportPolynomial(char* outputTxtFileName, Vector& outputState, ParticleOnSphere& initialSpace, OptionManager& manager, int nbrParticles, long alphaNumerator, long alphaDenominator, bool statistics);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereJackGenerator" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleDoubleOption  ('a', "alpha", "alpha coefficient of the Jack polynomial", -2.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "numerator-alpha", "numerator of the alpha coefficient of the Jack polynomial (enable in rational mode)", -2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "denominator-alpha", "denominator of the alpha coefficient of the Jack polynomial (enable in rational mode)", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0) to speed up calculations");
  (*SystemGroup) += new SingleStringOption  ('\n', "initial-state", "use an optional state where some of the components have already been computed, improving computation time");
  (*SystemGroup) += new BooleanOption  ('\n', "resume", "resume Jack calculation (only available in huge mode)");
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-index", "compute the Jack polynomial from the min-index-th component (require an initial state)", 0l);
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-index", "compute the Jack polynomial from the max-index-th component (require an initial state, 0 if it has computed up to the end)", 0l);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "compute the slater decomposition of the Jack polynomial times Vandermonde");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new BooleanOption  ('\n', "large-basis", "use large Hilbert space support (i.e. handle non-squeezed Hilbert space larger than 2^31 without hard-drive storage)");
  (*SystemGroup) += new SingleIntegerOption ('\n', "huge-fulldim", "indicate the full Hilbert space dimension (i.e. without squeezing) when using huge Hilbert space (0 if it has to be computed)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "file-size", "maximum file size (in MBytes) when using huge mode", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n' ,"huge-vector", "maximum memory (in MBytes) that can allocated for buffering vector when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n' ,"huge-blocks", "maximum memory (in MBytes) that can allocated for buffering indices when using huge mode (useful to improve parallelization speed-up)", 100);
  (*SystemGroup) += new BooleanOption ('\n', "disk-storage", "use disk storage in huge mode both for the Hilbert space and vectors");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "large-memory", "maximum memory (in kBytes) that can allocated for precalculations when using huge mode", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "check-singularity", "display configurations which may produce singularities");
  (*SystemGroup) += new BooleanOption  ('\n', "check-connected", "display lowest configuration connected to each squeezed paritition");
  (*SystemGroup) += new BooleanOption  ('\n', "use-symbolic", "use symbolic calculation to solve singular coefficient (only available in rational mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "full-symbolic", "use symbolic calculation for every case (only available in rational mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "display-symbolic", "display symbolic calculation");
#ifdef __GMP__
  (*SystemGroup) += new BooleanOption  ('\n', "use-gmp", "use arbitrary precision integers instead of fixed precision integers in rational mode");
#else
  (*SystemGroup) += new BooleanOption  ('\n', "use-longlong", "use 128bit(64bits) integers instead of 64bits(32bits) integers in rational mode");
#endif
  (*SystemGroup) += new BooleanOption  ('\n', "p-truncated", "use a p-truncated basis instead of the full squeezed basis");
  (*SystemGroup) += new SingleIntegerOption ('\n', "p-truncation", "p-truncation for the p-truncated basis (if --p-truncated is used)", 0);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  (*OutputGroup) += new BooleanOption ('\n', "fortran-header", "add a header to the text file to make the file fortan compatible");
  (*OutputGroup) += new BooleanOption ('\n', "column-partition", "use column formatted output for the partitions");
  (*OutputGroup) += new BooleanOption ('\n', "integer-occupation", "display the occupation numbers as an integer");
  (*OutputGroup) += new BooleanOption ('n', "normalize", "express the Jack polynomial in the normalized basis");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereJackGenerator -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0; 
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  double Alpha = ((SingleDoubleOption*) Manager["alpha"])->GetDouble();
  long AlphaDenominator = Manager.GetInteger("denominator-alpha");
  long AlphaNumerator = Manager.GetInteger("numerator-alpha");
  int TotalLz = 0;
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  int SymbolicDepth = 0;
  if (Manager.GetBoolean("use-symbolic") == true)
    {
      SymbolicDepth = 1;
    }
  if (Manager.GetBoolean("full-symbolic") == true)
    {
      SymbolicDepth = 2;
    }
  if (Manager.GetBoolean("display-symbolic") == true)
    {
      SymbolicDepth = 3;
    }

  long MinIndex = Manager.GetInteger("min-index");
  long MaxIndex = Manager.GetInteger("max-index");
  if (((MinIndex != 0l) || (MaxIndex != 0)) && (Manager.GetString("initial-state") == 0))
    {
      cout << "error, min-index/max-index options require an inital state" << endl;
      return 0;
    }

  if ((OutputTxtFileName == 0) && (OutputFileName == 0) && (Manager.GetString("save-hilbert")==0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  int* ReferenceState = 0;
  if (Manager.GetString("reference-file") == 0)
    {
      cout << "error, a reference file is needed" << endl;
      return 0;
    }
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return 0;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return 0;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return 0;
    }
  int MaxNbrLz;
  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
    {
      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
      return 0;     
    }
  if (MaxNbrLz != (NbrFluxQuanta + 1))
    {
      cout << "wrong LzMax value in ReferenceState" << endl;
      return 0;     
    }

  if (Manager.GetBoolean("fermion") == false)
    {
      if (Manager.GetBoolean("huge-basis") == true)
	{
	  BosonOnSphereHaldaneHugeBasisShort* InitialSpace = 0;
	  if (Manager.GetString("save-hilbert") != 0)
	    {
	      InitialSpace = new BosonOnSphereHaldaneHugeBasisShort (NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("file-size"), ReferenceState, ((unsigned long) Manager.GetInteger("memory")) << 20, false);
	      InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
	      return 0;
	    }
	  if (Manager.GetString("load-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  InitialSpace = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	  cout << "dimension = " << InitialSpace->GetLargeHilbertSpaceDimension() << endl;
	  bool DiskStorageFlag = InitialSpace->CheckDiskStorage();	  
	  if (DiskStorageFlag == true)
	    {
	      cout << "using disk storage for the Hilbert space and output state, some options might be disabled" << endl;
	    }
	  if (Manager.GetBoolean("disk-storage") == true)
	    DiskStorageFlag = true;
	  if (DiskStorageFlag == false)
	    {
	      if  (Manager.GetBoolean("rational") == false)
		{
		  RealVector OutputState;
		  if (Manager.GetString("initial-state") == 0) 
		    {
		      OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		    }
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false) 
		    {
		      InitialSpace->GenerateJackPolynomial(OutputState, Alpha, MinIndex, MaxIndex, OutputFileName);
		    }
		  else
		    InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha, MinIndex, MaxIndex, OutputFileName);
		  
		  if (Manager.GetBoolean("normalize"))
		    InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		}
	      else
		{
#ifdef __GMP__
		  if (Manager.GetBoolean("use-gmp") == false)
#else
		    if (Manager.GetBoolean("use-longlong") == false)
#endif
		      {
			cout << "only use-gmp/use-longlong rational mode are available for bosons in huge mode" << endl;
			return 0;
		      }
		  LongRationalVector OutputState;
		  if (Manager.GetString("initial-state") == 0) 
		    {
		      OutputState = LongRationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		    }
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false) 
		    {
		      InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, 
							   MinIndex, MaxIndex, OutputFileName);
		    }
		  else
		    {
		      cout << "symmetrized basis in huge and rational modes is not available" << endl;
		      return -1;
		    }
		  
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		  return 0;
		}
	    }
	  else
	    {
	      if (OutputFileName == 0)
		{
		  cout << "a binary output file name has to be provided when using disk storage mode" << endl;
		  return -1;
		}
	      if (SymmetrizedBasis == false)    
		{
		  InitialSpace->GenerateJackPolynomialSparse(Alpha, Architecture.GetArchitecture(), OutputFileName, MinIndex, MaxIndex, Manager.GetInteger("huge-vector") << 20, Manager.GetInteger("huge-blocks") << 20, Manager.GetBoolean("resume"));
		}
	      else
		InitialSpace->GenerateSymmetrizedJackPolynomialSparse(Alpha, Architecture.GetArchitecture(), OutputFileName, MinIndex, MaxIndex, Manager.GetInteger("huge-vector") << 20, Manager.GetInteger("huge-blocks") << 20, Manager.GetBoolean("resume"));
	    }
	  return 0;
	}
#ifdef __64_BITS__
      if ((NbrFluxQuanta + NbrParticles - 1) < 63)
#else
      if ((NbrFluxQuanta + NbrParticles - 1) < 31)
#endif
	{
	  BosonOnSphereHaldaneBasisShort* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    InitialSpace = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
	  else
	    {
	      InitialSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	      if (Manager.GetString("save-hilbert") != 0)
		{
		  InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		  return 0;
		}
	    }
	  if (Manager.GetBoolean("rational") == false)
	    {
	      RealVector OutputState;
	      if (Manager.GetBoolean("check-singularity") == true)
		{
		  OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  InitialSpace->CheckPossibleSingularCoefficientsInJackPolynomial(OutputState, Alpha, 1e-14);
		  cout << "partitions that may lead to singular coefficients : " << endl;
		  for (long i = 1l; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		    if (OutputState[i] != 0.0)
		      {
			InitialSpace->PrintStateMonomial(cout, i) << " = ";
			InitialSpace->PrintState(cout, i) << " (pos = " << i<< ")" << endl;
		      }
		  return 0;
		}
	      if (Manager.GetBoolean("check-connected") == true)
		{
		  ((BosonOnSphereHaldaneBasisShort*) InitialSpace)->CheckMaximumConnectedStateInJackPolynomial();
		  return 0;
		}
	      if (Manager.GetString("initial-state") == 0)
		OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
	      else
		if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		  {
		    cout << "can't open " << Manager.GetString("initial-state") << endl;
		    return -1;
		  }
	      if (SymmetrizedBasis == false)    
		InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
	      else
		InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
	      if (Manager.GetBoolean("normalize"))
		InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);      
	      if (OutputTxtFileName != 0)
		{
		  FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		}
	      if (OutputFileName != 0)
		{
		  OutputState.WriteVector(OutputFileName);
		}
	    }
	  else
	    {
#ifdef __GMP__
	      if (Manager.GetBoolean("use-gmp") == false)
#else
	      if (Manager.GetBoolean("use-longlong") == false)
#endif
		{
		  RationalVector OutputState;
		  if (Manager.GetBoolean("check-singularity") == true)
		    {
		      OutputState = RationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);		  
		      InitialSpace->CheckPossibleSingularCoefficientsInJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, true);
		      cout << "partitions that may lead to singular coefficients : " << endl;
		      Rational Zero = 0l;
		      for (long i = 1l; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
			if (OutputState[i] != Zero)
			  {
			    InitialSpace->PrintStateMonomial(cout, i) << " = ";
			    InitialSpace->PrintState(cout, i) << " (pos = " << i << ", connectivity = " << OutputState[i] << ") " << endl;
			  }
		      return 0;
		    }
		  if (Manager.GetString("initial-state") == 0)
		    OutputState = RationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false)    
		    InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, SymbolicDepth);
		  if (Manager.GetBoolean("normalize"))
		    {
		      cout << "calculations have been done with rational numbers, normalization will not be done" << endl;
		    }
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		}
	      else
		{
		  LongRationalVector OutputState;
		  if (Manager.GetString("initial-state") == 0)
		    OutputState = LongRationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false)    
		    InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		  else
		    InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		  if (Manager.GetBoolean("normalize"))
		    {
		      cout << "calculations have been done with rational numbers, normalization will not be done" << endl;
		    }
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		}
	    }
	}
      else
	{
	  BosonOnSphereHaldaneBasisLong* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    InitialSpace = new BosonOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert"));
	  else
	    {
	      InitialSpace = new BosonOnSphereHaldaneBasisLong(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	      if (Manager.GetString("save-hilbert") != 0)
		{
		  InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		  return 0;
		}
	    }
	  if (Manager.GetBoolean("rational") == false)
	    {
	      RealVector OutputState;
	      if (Manager.GetBoolean("check-singularity") == true)
		{
		  OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  InitialSpace->CheckPossibleSingularCoefficientsInJackPolynomial(OutputState, Alpha, 1e-14);
		  cout << "partitions that may lead to singular coefficients : " << endl;
		  for (long i = 1l; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		    if (OutputState[i] != 0.0)
		      {
			InitialSpace->PrintStateMonomial(cout, i) << " = ";
			InitialSpace->PrintState(cout, i) << endl;
		      }
		  return 0;
		}
	      if (Manager.GetString("initial-state") == 0)
		OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
	      else
		if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		  {
		    cout << "can't open " << Manager.GetString("initial-state") << endl;
		    return -1;
		  }
	      if (Manager.GetBoolean("normalize"))
		InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
	      if (SymmetrizedBasis == false)    
		InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
	      else
		InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
	      if (Manager.GetBoolean("normalize"))
		InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);      
	      if (OutputTxtFileName != 0)
		{
		  FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		}
	      if (OutputFileName != 0)
		{
		  OutputState.WriteVector(OutputFileName);
		}
	    }
	  else
	    {
#ifdef __GMP__
	      if (Manager.GetBoolean("use-gmp") == false)
#else
	      if (Manager.GetBoolean("use-longlong") == false)
#endif
		{
		  cout << "rational mode without --use-gmp or --use-longlong is not available for bosons with ";
#ifdef __64_BITS__
		  cout << "(LzMax + NbrParticles - 1) < 63" << endl;
#else
		  cout << "( LzMax+ NbrParticles - 1) < 31" << endl;
#endif
		  return 0;
		}
	      else
		{
		  LongRationalVector OutputState;
		  if (Manager.GetString("initial-state") == 0)
		    OutputState = LongRationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false)    
		    InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		  else
		    InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		  if (Manager.GetBoolean("normalize"))
		    {
		      cout << "calculations have been done with rational numbers, normalization will not be done" << endl;
		    }
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		}
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("huge-basis") == true)
	{
	  FermionOnSphereHaldaneHugeBasis* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    InitialSpace = new FermionOnSphereHaldaneHugeBasis(Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	  else
	    {
	      InitialSpace = new FermionOnSphereHaldaneHugeBasis (NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("file-size"), ReferenceState, ((unsigned long) Manager.GetInteger("memory")) << 20, false, Manager.GetInteger("huge-fulldim"));
	      if (Manager.GetString("save-hilbert") != 0)
		{		      
		  InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		  return 0;
		}
	    }
	  bool DiskStorageFlag = InitialSpace->CheckDiskStorage();	  
	  if (DiskStorageFlag == true)
	    {
	      cout << "using disk storage for the Hilbert space and output state, some options might be disabled" << endl;
	    }
	  if (Manager.GetBoolean("disk-storage") == true)
	    DiskStorageFlag = true;
	  if (DiskStorageFlag == false)
	    {
	      if  (Manager.GetBoolean("rational") == false)
		{
		  RealVector OutputState;
		  if (Manager.GetString("initial-state") == 0)
		    OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false)    
		    InitialSpace->GenerateJackPolynomial(OutputState, Alpha, MinIndex, MaxIndex, OutputFileName);
		  else
		    InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha, MinIndex, MaxIndex, OutputFileName);
		  if (Manager.GetBoolean("normalize"))
		    InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		}
	      else
		{
#ifdef __GMP__
		  if (Manager.GetBoolean("use-gmp") == false)
#else
		    if (Manager.GetBoolean("use-longlong") == false)
#endif
		      {
			cout << "only use-gmp/use-longlong rational mode are available for fermions" << endl;
			return 0;
		      }
		  LongRationalVector OutputState;
		  if (Manager.GetString("initial-state") == 0)
		    OutputState = LongRationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		  else
		    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		      {
			cout << "can't open " << Manager.GetString("initial-state") << endl;
			return -1;
		      }
		  if (SymmetrizedBasis == false)    
		    {
		      InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, MinIndex, MaxIndex, OutputFileName);
		    }
		  else
		    {
		      cout << "symmetrized basis in huge and rational modes is not available" << endl;
		    }
		  if (OutputTxtFileName != 0)
		    {
		      FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		    }
		  if (OutputFileName != 0)
		    {
		      OutputState.WriteVector(OutputFileName);
		    }
		  return 0;
		}	      
	    }
	  else
	    {
	      if (SymmetrizedBasis == false)    
		InitialSpace->GenerateJackPolynomialSparse(Alpha, Architecture.GetArchitecture(), OutputFileName, MinIndex, MaxIndex, Manager.GetInteger("huge-vector") << 20, Manager.GetInteger("huge-blocks") << 20, Manager.GetBoolean("resume"));
	      else
		InitialSpace->GenerateSymmetrizedJackPolynomialSparse(Alpha, Architecture.GetArchitecture(), OutputFileName, MinIndex, MaxIndex, Manager.GetInteger("huge-vector") << 20, Manager.GetInteger("huge-blocks") << 20, Manager.GetBoolean("resume"));
	    }
	  return 0;
	}
#ifdef __64_BITS__
      if (NbrFluxQuanta < 63)
#else
	if (NbrFluxQuanta < 31)
#endif
	  {
	    if (Manager.GetBoolean("p-truncated") == true)
	      {
		FermionOnSpherePTruncated* InitialSpace = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
		RealVector OutputState;
		OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
		if (Manager.GetBoolean("normalize"))
		  InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
		if (OutputTxtFileName != 0)
		  {
		    FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		  }
		if (OutputFileName != 0)
		  {
		    OutputState.WriteVector(OutputFileName);
		  }
		return 0;
	      }
	    
	    FermionOnSphereHaldaneBasis* InitialSpace;
	    if (Manager.GetBoolean("large-basis") == true)
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  InitialSpace = new FermionOnSphereHaldaneLargeBasis(Manager.GetString("load-hilbert"), Manager.GetInteger("large-memory") << 10);
		else
		  {
		    InitialSpace = new FermionOnSphereHaldaneLargeBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState, Manager.GetInteger("large-memory") << 10);	  
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
	      }
	    else
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  InitialSpace = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"));
		else
		  {
		    InitialSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
	      }
	    if (Manager.GetBoolean("rational") == false)
	      {
		RealVector OutputState;
		if (Manager.GetBoolean("check-singularity") == true)
		  {
		    OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		    InitialSpace->CheckPossibleSingularCoefficientsInJackPolynomial(OutputState, Alpha, 1e-14);
		    cout << "partitions that may lead to singular coefficients : " << endl;
		    for (long i = 1l; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		      if (OutputState[i] != 0.0)
			{
			  InitialSpace->PrintStateMonomial(cout, i) << " = ";
			  InitialSpace->PrintState(cout, i) << endl;
			}
		    return 0;
		  }
		if (Manager.GetString("initial-state") == 0)
		  OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		else
		  if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		    {
		      cout << "can't open " << Manager.GetString("initial-state") << endl;
		      return -1;
		    }
		if (SymmetrizedBasis == false)    
		  InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
		else
		  InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
		if (Manager.GetBoolean("normalize"))
		  InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
		if (OutputTxtFileName != 0)
		  {
		    FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		  }
		if (OutputFileName != 0)
		  {
			OutputState.WriteVector(OutputFileName);
		  }
	      }
	    else
	      {
#ifdef __GMP__
		if (Manager.GetBoolean("use-gmp") == false)
#else
		  if (Manager.GetBoolean("use-longlong") == false)
#endif
		    {
		      RationalVector OutputState;
		      if (Manager.GetBoolean("check-singularity") == true)
			{
			  OutputState = RationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);		  
			  InitialSpace->CheckPossibleSingularCoefficientsInJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, true);
			  cout << "partitions that may lead to singular coefficients : " << endl;
			  Rational Zero = 0l;
			  for (long i = 1l; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
			    if (OutputState[i] != Zero)
			      {
				InitialSpace->PrintStateMonomial(cout, i) << " = ";
				InitialSpace->PrintState(cout, i) << " (pos = " << i << ", connectivity = " << OutputState[i] << ") " << endl;
			      }
			  return 0;
			}
		      cout << "only use-gmp/use-longlong rational mode are available for fermions" << endl;
		      return 0;
		    }
		  else
		    {
		      LongRationalVector OutputState;
		      if (Manager.GetString("initial-state") == 0)
			OutputState = LongRationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		      else
			if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
			  {
			    cout << "can't open " << Manager.GetString("initial-state") << endl;
			    return -1;
			  }
		      if (SymmetrizedBasis == false)    
			InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		      else
			InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		      if (Manager.GetBoolean("normalize"))
			{
			  cout << "calculations have been done with rational numbers, normalization will not be done" << endl;
			    }
		      if (OutputTxtFileName != 0)
			{
			  FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
			}
		      if (OutputFileName != 0)
			{
			  OutputState.WriteVector(OutputFileName);
			}
		    }
	      }
	  }
	else
	  {
	    FermionOnSphereHaldaneBasisLong* InitialSpace;
	    if (Manager.GetString("load-hilbert") != 0)
	      InitialSpace = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert"));
	    else
	      {
		InitialSpace = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
		if (Manager.GetString("save-hilbert") != 0)
		  {
		    InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		    return 0;
		  }
	      }
	    if (Manager.GetBoolean("rational") == false)
	      {
		RealVector OutputState;
		if (Manager.GetBoolean("check-singularity") == true)
		  {
		    OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		    InitialSpace->CheckPossibleSingularCoefficientsInJackPolynomial(OutputState, Alpha, 1e-14);
		    cout << "partitions that may lead to singular coefficients : " << endl;
		    for (long i = 1l; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		      if (OutputState[i] != 0.0)
			{
			  InitialSpace->PrintStateMonomial(cout, i) << " = ";
			  InitialSpace->PrintState(cout, i) << endl;
			}
		    return 0;
		  }
		if (Manager.GetString("initial-state") == 0)
		  OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		else
		  if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
		    {
		      cout << "can't open " << Manager.GetString("initial-state") << endl;
		      return -1;
		    }
		if (SymmetrizedBasis == false)    
		  InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
		else
		  InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
		if (Manager.GetBoolean("normalize"))
		  InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
		if (OutputTxtFileName != 0)
		  {
		    FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
		  }
		if (OutputFileName != 0)
		  {
		    OutputState.WriteVector(OutputFileName);
		  }
	      }
	    else
	      {
#ifdef __GMP__
		if (Manager.GetBoolean("use-gmp") == false)
#else
		  if (Manager.GetBoolean("use-longlong") == false)
#endif
		    {
		      cout << "only use-gmp/use-longlong rational mode are available for fermions" << endl;
		      return 0;
		    }
		  else
		    {
		      LongRationalVector OutputState;
		      if (Manager.GetString("initial-state") == 0)
			OutputState = LongRationalVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
		      else
			if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
			  {
			    cout << "can't open " << Manager.GetString("initial-state") << endl;
			    return -1;
			  }
		      if (SymmetrizedBasis == false)    
			InitialSpace->GenerateJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		      else
			InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, AlphaNumerator, AlphaDenominator, Architecture.GetArchitecture(), SymbolicDepth, MinIndex, MaxIndex, OutputFileName);
		      if (Manager.GetBoolean("normalize"))
			{
			  cout << "calculations have been done with rational numbers, normalization will not be done" << endl;
			    }
		      if (OutputTxtFileName != 0)
			{
			  FQHESphereJackTxtExportPolynomial(OutputTxtFileName, OutputState, *InitialSpace, Manager, NbrParticles, AlphaNumerator, AlphaDenominator, Manager.GetBoolean("fermion"));
			}
		      if (OutputFileName != 0)
			{
			  OutputState.WriteVector(OutputFileName);
			}
		    }
	      }
	  }
    }
  return 0;
}

// export the Jack polynomial in a text file
//
// outputTxtFileName = output file name
// outputState = Jack polynomial decomposition
// initialSpace = pointer to the Hilbert space 
// manager = reference to the option manager
// nbrParticles = number of particles
// alphaNumerator= numerator of the alpha factor
// alphaDenominator= denominator of the alpha factor
// statistics = particle statistics (true for fermions)

void FQHESphereJackTxtExportPolynomial(char* outputTxtFileName, Vector& outputState, ParticleOnSphere& initialSpace, OptionManager& manager, int nbrParticles, long alphaNumerator, long alphaDenominator, bool statistics)
{
  ofstream File;
  File.open(outputTxtFileName, ios::binary | ios::out);
  File.precision(14);
  if (manager.GetBoolean("fortran-header"))
    {
      File << initialSpace.GetLargeHilbertSpaceDimension() << " " << alphaNumerator << " " << alphaDenominator;
      if (statistics == true)
	File << " 1" << endl;
      else
	File << " 0" << endl;
    }
   if (manager.GetBoolean("integer-occupation"))
     {			  
      for (long i = 0; i < initialSpace.GetLargeHilbertSpaceDimension(); ++i)
	{
	  initialSpace.PrintCompactState(File, i) << " ";
	  outputState.PrintComponent(File, i) << endl;
	}
    }
  else
    {
      if (manager.GetBoolean("column-partition"))
	{
	  for (long i = 0; i < initialSpace.GetLargeHilbertSpaceDimension(); ++i)
	    {
	      outputState.PrintComponent(File, i) << " ";
	      initialSpace.PrintColumnFormattedStateMonomial(File, i) << endl;
	    }	  
	}
      else
	{			  
	  for (long i = 0; i < initialSpace.GetLargeHilbertSpaceDimension(); ++i)
	    {
	      outputState.PrintComponent(File, i) << " ";
	      initialSpace.PrintStateMonomial(File, i) << endl;
	    }
	}
    }
  File.close();
}
