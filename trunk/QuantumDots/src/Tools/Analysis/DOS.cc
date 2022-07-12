#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/DOSSpectra.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("DOS" , "0.01");
  OptionGroup* DOSGroup = new OptionGroup ("DOS");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += DOSGroup;
  Manager += MiscGroup;

  (*DOSGroup) += new SingleStringOption('\n', "input", "name of the input file", "eigenvalues");
  (*DOSGroup) += new SingleIntegerOption('n', "nbr-state", "number of states", 10);
  (*DOSGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in eV unit)", 0.0);
  (*DOSGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in eV unit)", 0.0);
  (*DOSGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in eV unit)", 0.01);
  (*DOSGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in eV unit) in the spectrum", 2e-4);
  (*DOSGroup) += new SingleStringOption('\n', "output", "name of the output file", "DOS.txt");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  char* InputFile = ((SingleStringOption*) Manager["input"])->GetString();
  int NbrState = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  int Number = 1;
  char** Files = new char* [Number]; int* State = new int[Number];
  for (int i = 0; i < Number; ++i)
    {
      State[i] = NbrState;
      Files[i] = new char[100];
      Files[0] = InputFile;
    }
  DOSSpectra DOS(Number, Files, State, Gamma, Min, Max, Step);
  DOS.WriteSpectra(OutputFile);
  
  return 1;
}
