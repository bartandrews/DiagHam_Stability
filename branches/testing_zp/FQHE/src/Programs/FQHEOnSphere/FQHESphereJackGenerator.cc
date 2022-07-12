#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneLargeBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereJackGenerator" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleDoubleOption  ('a', "alpha", "alpha coefficient of the Jack polynomial", -2.0);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0) to speed up calculations");
  (*SystemGroup) += new BooleanOption  ('\n', "sym-storage", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0), both for speed and storage");
  (*SystemGroup) += new SingleStringOption  ('\n', "initial-state", "use an optional state where some of the components have already been computed, improving computation time");
  (*SystemGroup) += new SingleIntegerOption ('\n', "min-index", "compute the Jack polynomial from the min-index-th component (require an initial state)", 0l);
  (*SystemGroup) += new SingleIntegerOption ('\n', "max-index", "compute the Jack polynomial from the max-index-th component (require an initial state, 0 if it has computed up to the end)", 0l);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "compute the slater decomposition of the Jack polynomial times Vandermonde");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new BooleanOption  ('\n', "large-basis", "use large Hilbert space support (i.e. handle non-squeezed Hilbert space larger than 2^31 without hard-drive storage)");
  (*SystemGroup) += new SingleIntegerOption ('\n', "huge-fulldim", "indicate the full Hilbert space dimension (i.e. without squeezing) when using huge Hilbert space (0 if it has to be computed)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "file-size", "maximum file size (in MBytes) when using huge mode", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "large-memory", "maximum memory (in kBytes) that can allocated for precalculations when using huge mode", 1);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  (*OutputGroup) += new BooleanOption ('n', "normalize", "express the Jack polynomial in the normalized basis");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereJackGenerator -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0; 
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  double Alpha = ((SingleDoubleOption*) Manager["alpha"])->GetDouble();
  int TotalLz = 0;
  char* OutputFileName = ((SingleStringOption*) Manager["bin-output"])->GetString();
  char* OutputTxtFileName = ((SingleStringOption*) Manager["txt-output"])->GetString();

  long MinIndex = Manager.GetInteger("min-index");
  long MaxIndex = Manager.GetInteger("max-index");
  if (((MinIndex != 0l) || (MaxIndex != 0)) && (Manager.GetString("initial-state") == 0))
    {
      cout << "error, min-index/max-index options require an inital state" << endl;
      return 0;
    }

  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  int* ReferenceState = 0;
  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
    {
      cout << "error, a reference file is needed" << endl;
      return 0;
    }
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
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
      cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
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
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);
	      for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		{
		  File << OutputState[i] << " ";
		  InitialSpace->PrintStateMonomial(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      OutputState.WriteVector(OutputFileName);
	    }
	  return 0;
	}
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
      RealVector OutputState;
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
	  ofstream File;
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputState[i] << " ";
	      InitialSpace->PrintStateMonomial(File, i) << endl;
	    }
	  File.close();
	}
      if (OutputFileName != 0)
	{
	  OutputState.WriteVector(OutputFileName);
	}
    }
  else
    {
      if (Manager.GetBoolean("sym-storage") == false)
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
		  ofstream File;
		  File.open(OutputTxtFileName, ios::binary | ios::out);
		  File.precision(14);
		  for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		    {
		      File << OutputState[i] << " ";
		      InitialSpace->PrintStateMonomial(File, i) << endl;
		    }
		  File.close();
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
	    InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
	  else
	    InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
	  if (Manager.GetBoolean("normalize"))
            InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);
	      for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		{
		  File << OutputState[i] << " ";
		  InitialSpace->PrintStateMonomial(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      OutputState.WriteVector(OutputFileName);
	    }
	}
      else
	{
	  FermionOnSphereHaldaneSymmetricBasis* InitialSpace;
	  if (Manager.GetString("load-hilbert") != 0)
	    InitialSpace = new FermionOnSphereHaldaneSymmetricBasis(Manager.GetString("load-hilbert"));
	  else
	    {
	      InitialSpace = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, NbrFluxQuanta, ReferenceState);	  
	      if (Manager.GetString("save-hilbert") != 0)
		{
		  InitialSpace->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		  return 0;
		}
	    }
	  RealVector OutputState;
	  if (Manager.GetString("initial-state") == 0)
	    OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
	  else
	    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
	      {
		cout << "can't open " << Manager.GetString("initial-state") << endl;
		return -1;
	      }
	  InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
	  if (Manager.GetBoolean("normalize"))
	    InitialSpace->ConvertFromUnnormalizedMonomial(OutputState);
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);
	      for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
		{
		  File << OutputState[i] << " ";
		  InitialSpace->PrintStateMonomial(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      OutputState.WriteVector(OutputFileName);
	    }
	}
    }    
  return 0;
}

