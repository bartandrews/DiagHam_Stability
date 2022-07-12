#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"


#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("FQHESphereDipolarPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle)", 8);
  (*SystemGroup) += new SingleDoubleOption ('v', "set-v0", "fix v0 to a given value", 0.0);  
  (*SystemGroup) += new SingleDoubleOption ('r', "set-ratio", "fix the ratio v2/v0 to a given value", 0.0);  
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_dipolar_2s_x.dat)");
  (*SystemGroup) += new BooleanOption ('\n', "std-output", "use standard output instead of an output file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereDipolarPseudopotentials -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFlux = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger();
  int NbrPseudopotentials = NbrFlux + 1;
  double* Pseudopotentials = EvaluateDipolarPseudopotentials(NbrFlux, !Manager.GetBoolean("std-output"));

  if (Manager.GetDouble("set-ratio") != 0.0)
    {
      Pseudopotentials[0] = Pseudopotentials[2] / Manager.GetDouble("set-ratio");
    }
  else
    {
      Pseudopotentials[0] = Manager.GetDouble("set-v0");
    }
  char* OutputFile;
  if (Manager.GetString("output") == 0l)
    OutputFile = Manager.GetFormattedString("pseudopotential_dipolar_2s_%nbr-flux%.dat");
  else
    {
      OutputFile = new char [strlen(((SingleStringOption*) Manager["output"])->GetString()) + 1];
      strcpy (OutputFile, ((SingleStringOption*) Manager["output"])->GetString());
    }

  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  File << "Pseudopotentials=";        
  for (int i = 0; i < NbrPseudopotentials; ++i)
    File << " " << Pseudopotentials[i];
  File << endl;      
  File.close();
  
  cout.precision(14);
  
  delete[] OutputFile;
  delete[] Pseudopotentials;

 return 0;
}
