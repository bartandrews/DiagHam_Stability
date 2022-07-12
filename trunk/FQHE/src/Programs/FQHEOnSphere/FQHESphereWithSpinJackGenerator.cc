#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLargeBasis.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

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
  OptionManager Manager ("FQHESphereWithSpinJackGenerator" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "initial-state", "use an optional state where some of the components have already been computed, improving computation time");
  (*SystemGroup) += new BooleanOption  ('\n', "large-basis", "use large Hilbert space support (i.e. handle non-squeezed Hilbert space larger than 2^31 without hard-drive storage)");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  (*OutputGroup) += new BooleanOption ('\n', "txt-separatespin", "for the text output, use the sign convention which separates spins");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinJackGenerator -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0; 
  int TotalSz = 0;
  double Alpha = ((SingleDoubleOption*) Manager["alpha"])->GetDouble();
  int TotalLz = 0;
  char* OutputFileName = ((SingleStringOption*) Manager["bin-output"])->GetString();
  char* OutputTxtFileName = ((SingleStringOption*) Manager["txt-output"])->GetString();

  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  int** ReferenceStates = 0;
  int NbrReferenceStates;
  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
    {
      cout << "error, a reference file is needed" << endl;
      return 0;
    }
  if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceStates, NbrReferenceStates) == false)
    {
      cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
      return 0;
    }

  FermionOnSphereWithSpinHaldaneBasis* InitialSpace;
  if (Manager.GetBoolean("large-basis") == false)
    {
      //   if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
      //     InitialSpace = new FermionOnSphereWithSpinHaldaneBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString());
      //   else
      //    {
      InitialSpace = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, ReferenceStates, NbrReferenceStates); 
      if (Manager.GetString("save-hilbert") != 0)
	{
	  InitialSpace->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
	  return 0;
	}
      //}
    }
  else
    {
      if (Manager.GetString("load-hilbert") != 0)
	{
	  InitialSpace = new FermionOnSphereWithSpinHaldaneLargeBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString());
	}
      else
	{
	  InitialSpace = new FermionOnSphereWithSpinHaldaneLargeBasis(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, ReferenceStates, NbrReferenceStates); 
	  if (Manager.GetString("save-hilbert") != 0)
	    {
	      ((FermionOnSphereWithSpinHaldaneLargeBasis*) InitialSpace)->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
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
  InitialSpace->GenerateJackPolynomial(OutputState, Alpha);

  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      if (Manager.GetBoolean("txt-separatespin") == false)
	{
	  for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputState[i] << " ";
	      InitialSpace->PrintStateMonomial(File, i) << endl;
	    }
	}
      else
	{
	  for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << (OutputState[i] * InitialSpace->GetSpinSeparationSignFromIndex(i)) << " ";
	      InitialSpace->PrintStateMonomialSeparatedSpin(File, i) << endl;
	    }
	}
      File.close();
    }
  if (OutputFileName != 0)
    {
      OutputState.WriteVector(OutputFileName);
    }

  return 0;
}

