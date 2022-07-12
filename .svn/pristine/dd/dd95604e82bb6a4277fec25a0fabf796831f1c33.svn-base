#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"

#include "Vector/ComplexVector.h"

#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

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
  OptionManager Manager ("HubbardASCIIToState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "number of flux quanta", 20);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new SingleStringOption ('i', "ascii-state", "name of the input ASCII description of the state (should use the same convention than FQHESphereShowBasis output)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of trying to replace .txt extension with .vec or appending .vec extension");
  (*SystemGroup) += new BooleanOption  ('\n', "no-normalization", "do not normalize the final state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardASCIIToState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if (Manager.GetString("ascii-state") == 0)
    {
      cout << "An input file has to be provided" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSites = Manager.GetInteger("nbr-sites"); 
    
  MultiColumnASCIIFile InputFile(':');
  if (InputFile.Parse(Manager.GetString("ascii-state")) == false)
    {
      InputFile.DumpErrors(cout) << endl;
      return -1;
    }
  Complex* Coefficients = InputFile.GetAsComplexArray(1);
  if (Coefficients == 0)
    {
      InputFile.DumpErrors(cout) << endl;
      return -1;
    }

  ParticleOnSphere* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      cout << "bosonic Hubbard model not implemented" << endl;
      return -1;
    }
  else
    {
      if (Manager.GetBoolean("gutzwiller") == false)
	{
	  Space = new FermionOnLatticeWithSpinRealSpace(NbrParticles, NbrSites);
	}
      else
	{
	  Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(NbrParticles, NbrSites);
	}
    }
  
  ComplexVector State (Space->GetHilbertSpaceDimension(), true);

  for (int i = 0; i < InputFile.GetNbrLines(); ++i)
    {
      char* TmpString =  InputFile(0, i);
      CleanLine(TmpString);
      int TmpIndex = Space->FindStateIndex(TmpString);
      if (TmpIndex != -1)
	{
	  State[TmpIndex] = Coefficients[i];
	}
      else
	{
	  cout << "warning , invalid state |" <<  TmpString << ">" << endl;
	}
    }
  if (Manager.GetBoolean("no-normalization") == false)
    State /= State.Norm();

  char* OutputFile = Manager.GetString("output-file");
  if (OutputFile == 0)
    {
      OutputFile = ReplaceExtensionToFileName(Manager.GetString("ascii-state"), "txt", "vec");
      if (OutputFile == 0)
	OutputFile = AddExtensionToFileName(Manager.GetString("ascii-state"), "vec");
    }
  State.WriteVector(OutputFile);

  return 0;
}

