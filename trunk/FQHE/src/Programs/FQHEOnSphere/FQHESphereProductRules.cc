#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

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
  OptionManager Manager ("FQHESphereProductRules" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that describes the smaller system");
  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a file as the definition of the reference state of the input state");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new SingleStringOption  ('\n', "common-pattern", "pattern that has to be shared between the two n-body states");
  (*SystemGroup) += new SingleStringOption  ('\n', "added-pattern", "additional pattern that has inserted to go from the smaller system to the larger one");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "factor", "multiplicative factor to use when going from the smaller system to the larger one", 1.0);
  (*SystemGroup) += new SingleStringOption  ('\n', "initial-state", "use an optional state where some of the components have already been computed, improving computation time");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "inputhuge-basis", "use huge Hilbert space support for the input state");
  (*SystemGroup) += new BooleanOption  ('\n', "outputhuge-basis", "use huge Hilbert space support for the output state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the fused vector that will be generated");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector into a text file");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "inputload-hilbert", "load Hilbert space description from the input file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "outputload-hilbert", "load Hilbert space description from the output file",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereProductRules -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  MultiColumnASCIIFile InputVectors;
  
  int InputNbrParticles = 0;
  int InputLzMax = 0;
  int InputTotalLz = 0;
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						   InputNbrParticles, InputLzMax, InputTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from input state name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  int CommonNbrParticles = 0;
  int CommonLzMax = 0;
  int* CommonPattern = 0;
  if (FQHEGetRootPartition(Manager.GetString("common-pattern"), CommonNbrParticles, CommonLzMax, CommonPattern) == false)
    return -1;
  int AddedNbrParticles = 0;
  int AddedLzMax = 0;
  int* AddedPattern = 0;
  if (FQHEGetRootPartition(Manager.GetString("added-pattern"), AddedNbrParticles, AddedLzMax, AddedPattern) == false)
    return -1;

  RealVector InputVector;
  if (InputVector.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }

  ParticleOnSphere* InputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("inputhuge-basis") == true)
	{
	  if (Manager.GetString("inputload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  InputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("inputload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("input-reference"), InputNbrParticles, InputLzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("inputload-hilbert") != 0)
	    InputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("inputload-hilbert"));	  
	  else
	    InputBasis = new BosonOnSphereHaldaneBasisShort(InputNbrParticles, InputTotalLz, InputLzMax, ReferenceState);	  
	}
    }
  else
    {
      int* ReferenceState = 0;
      if (FQHEGetRootPartition(Manager.GetString("input-reference"), InputNbrParticles, InputLzMax, ReferenceState) == false)
	return -1;
      if (Manager.GetString("inputload-hilbert") != 0)
	InputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("inputload-hilbert"));	  
      else
	InputBasis = new FermionOnSphereHaldaneBasis(InputNbrParticles, InputTotalLz, InputLzMax, ReferenceState);
    }

  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("outputhuge-basis") == true)
	{
	  if (Manager.GetString("outputload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  OutputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("outputload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("output-reference"), NbrParticles, LzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("outputload-hilbert") != 0)
	    OutputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("outputload-hilbert"));	  
	  else
	    OutputBasis = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);	  
	}
    }
  else
    {
      int* ReferenceState = 0;
      if (FQHEGetRootPartition(Manager.GetString("output-reference"), NbrParticles, LzMax, ReferenceState) == false)
	return -1;
      if (Manager.GetString("outputload-hilbert") != 0)
	OutputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("outputload-hilbert"));	  
      else
	OutputBasis = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState);
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [256];
      if (Statistics == false)
	sprintf (OutputFileName, "bosons_productrules_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
      else
	sprintf (OutputFileName, "fermions_productrules_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
    }

  RealVector OutputState;
  if (Manager.GetString("initial-state") == 0)
    OutputState = RealVector(OutputBasis->GetLargeHilbertSpaceDimension(), true);
  else
    if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
      {
	cout << "can't open " << Manager.GetString("initial-state") << endl;
	return -1;
      }

  OutputBasis->ProductRules(OutputState, InputVector, InputBasis, CommonPattern, CommonLzMax + 1, AddedPattern, AddedLzMax + 1, 
			    Manager.GetDouble("factor"), SymmetrizedBasis);

  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
	{
	  File << OutputState[i] << " ";
	  OutputBasis->PrintStateMonomial(File, i) << endl;
	}
      File.close();
    }
  if (OutputFileName != 0)
    {
      if (OutputState.WriteVector(OutputFileName) == false)
	{
	  cout << "error while writing output state " << OutputFileName << endl;
	  return -1;
	}	  
    }

  return 0;
}


