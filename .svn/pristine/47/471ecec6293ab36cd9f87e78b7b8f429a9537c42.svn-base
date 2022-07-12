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
#include <strstream>
#include <string>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("SumSpectrum" , "0.01");
  OptionGroup* SumSpectrumGroup = new OptionGroup ("SumSpectrum");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += SumSpectrumGroup;
  Manager += MiscGroup;

  (*SumSpectrumGroup) += new SingleStringOption('\n', "input", "name of the spectrum input file", "PolarizationZ.txt");
  (*SumSpectrumGroup) += new SingleStringOption('\n', "prefix", "prefix of the directories containing spectrum", "");
  (*SumSpectrumGroup) += new SingleIntegerOption('\n', "begin", "number of the first directory", 0);
  (*SumSpectrumGroup) += new SingleIntegerOption('\n', "end", "number of the last directory", 0);
  (*SumSpectrumGroup) += new SingleIntegerOption('n', "nbr-state", "number of states in each spectrum", 10);
  (*SumSpectrumGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in eV unit)", 0.0);
  (*SumSpectrumGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in eV unit)", 0.0);
  (*SumSpectrumGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in eV unit)", 0.01);
  (*SumSpectrumGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in eV unit) in the spectrum", 2e-4);
  (*SumSpectrumGroup) += new SingleStringOption('\n', "output", "name of the output file", "SumSpectrum.txt");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new BooleanOption ('v', "verbose", "verbose mode", false);

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
  char* Prefix = ((SingleStringOption*) Manager["prefix"])->GetString();
  int Begin = ((SingleIntegerOption*) Manager["begin"])->GetInteger();
  int End = ((SingleIntegerOption*) Manager["end"])->GetInteger();
  int NbrState = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  bool VerboseFlag = ((BooleanOption*) Manager["verbose"])->GetBoolean();

  int Number = End - Begin + 1;
  char* Prefixbis = new char [100]; char* InputFilebis = new char [100];
  AddString (Prefixbis, Prefix, 0, "");  strcpy(InputFilebis, "/"); strcat(InputFilebis, InputFile);
  char** Files = new char* [Number]; int* State = new int[Number];
  for (int i = Begin; i <= End; ++i)
    {
      State[i - Begin] = NbrState;
      Files[i - Begin] = new char[200];
      if (i < 10)
	AddString(Files[i - Begin], Prefixbis, i, InputFilebis);
      else
	AddString(Files[i - Begin], Prefix, i, InputFilebis);
      if (VerboseFlag)
	cout << Files[i - Begin] << endl;
    }
  Spectra SumSpectrum(Number, Files, State, Gamma, Min, Max, Step);
  SumSpectrum.WriteSpectra(OutputFile);
  
  return 1;
}
