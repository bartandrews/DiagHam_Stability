#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

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
  OptionManager Manager ("FQHESphereWithSpinConvertHaldaneBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('r', "reverse", "convert a state from the n-body basis to the Haldane basis");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of haldane_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinConvertHaldaneBasis -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-file"])->GetString() == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereWithSpinConvertHaldaneBasis -h" << endl;
      return -1;
    }
  if (IsFile(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  bool ReverseFlag = ((BooleanOption*) Manager["reverse"])->GetBoolean();
  bool Statistics = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(),
							   NbrParticles, NbrFluxQuanta, TotalLz, TotalSz, SzSymmetrizedBasis, SzMinusParity, 
							   LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
  if (((SingleIntegerOption*) Manager["total-lz"])->GetInteger() >= 0)
    TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger(); 
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
    }

  char* OutputFileName = 0;
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    {
      OutputFileName = new char [strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
      strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
    }
  else
    {
      char* InputFileName = ((SingleStringOption*) Manager["input-file"])->GetString(); 
      char* TagPosition = strcasestr(InputFileName, "haldane_");
      if (TagPosition == 0)
	{
	  cout << "no default output name can be built from " << InputFileName << endl;
	  return -1;
	}
      OutputFileName = new char [strlen(InputFileName) - 7];
      strncpy (OutputFileName, InputFileName, TagPosition - InputFileName);
      strcpy (OutputFileName + (TagPosition - InputFileName), TagPosition + 8);
    }


  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;      
    }

  int** ReferenceStates = 0;
  int NbrReferenceStates;
  if (Manager.GetString("reference-file") == 0)
    {
      cout << "error, a reference file is needed for fermions in Haldane basis" << endl;
      return 0;
    }
  if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceStates, NbrReferenceStates) == false)
    {
      cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
      return 0;
    }      
  FermionOnSphereWithSpinHaldaneBasis InitialSpace(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, ReferenceStates, NbrReferenceStates);
  if ((ReverseFlag == false) && (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension()))
    {
      cout << "dimension mismatch between Hilbert space and input state" << endl;
    }
  FermionOnSphereWithSpin TargetSpace(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
  if ((ReverseFlag == true) && (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension()))
    {
      cout << "dimension mismatch between Hilbert space and input state" << endl;
    }	    
  RealVector OutputState;
  if (ReverseFlag == false)
    OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
  else
    OutputState = InitialSpace.ConvertFromNbodyBasis(State, TargetSpace);
  if (OutputState.WriteVector(OutputFileName) == false)
    {
      cout << "error while writing output state " << OutputFileName << endl;
      return -1;
    }	  
  return 0;
}
