#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/Spectra.h"

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
  OptionManager Manager ("Spectrum" , "0.01");
  OptionGroup* SpectrumGroup = new OptionGroup ("Spectrum");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += SpectrumGroup;
  Manager += MiscGroup;

  (*SpectrumGroup) += new SingleStringOption('\n', "input", "name of the input file", "PolarizationZ.txt");
  (*SpectrumGroup) += new SingleIntegerOption('\n', "nbr-state", "number of states", 10);
  (*SpectrumGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in eV unit)", 0.0);
  (*SpectrumGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in eV unit)", 0.0);
  (*SpectrumGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in eV unit)", 0.01);
  (*SpectrumGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in eV unit) in the spectrum", 2e-4);
  (*SpectrumGroup) += new SingleStringOption('\n', "output", "name of the output file", "Spectrum.txt");

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
  Spectra Spectrum(Number, Files, State, Gamma, Min, Max, Step);
  Spectrum.WriteSpectra(OutputFile);
  
  return 1;
}
