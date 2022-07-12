#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/DOSSpectra.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <fstream>
#ifdef __SSTREAM_STYLE__
#include <sstream>
#else
#include <strstream>
#endif
#include <string>
#include <unistd.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("SumDOS2" , "0.01");
  OptionGroup* InputOptionGroup = new OptionGroup ("input options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += InputOptionGroup;
  Manager += PlotOptionGroup;
  Manager += MiscGroup;

  (*InputOptionGroup) += new SingleStringOption('d', "directory", "pattern to identify directories that contain states description (if pattern is mydir/run, any directory in mydir such that runx where x is an integer will be considered)");
  (*InputOptionGroup) += new SingleStringOption('\n', "input", "beginning of the name of file describing states (without relative/absolute path, files has to be of the form NAMEraw where x is an integer", "eigenvalues");
  (*InputOptionGroup) += new SingleIntegerOption('\n', "begin", "number of the first directory", 0);
  (*InputOptionGroup) += new SingleIntegerOption('\n', "end", "number of the last directory (0 if it has to be detected automatically)", 0);
  (*InputOptionGroup) += new SingleIntegerOption('n', "nbr-state", "number of initial states (0 if it has to be detected automatically)", 0);

  (*PlotOptionGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in meV unit)", 0.0);
  (*PlotOptionGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in meV unit)", 0.0);
  (*PlotOptionGroup) += new SingleDoubleOption('t', "temperature", "temperature (in Kelvin unit)", 10.0);
  (*PlotOptionGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in meV unit)", 0.01);
  (*PlotOptionGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in meV unit) in the spectrum", 2e-4);
  (*PlotOptionGroup) += new SingleStringOption('\n', "output", "name of the output file", "absorption.dat");

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

  char* StateDirectory = ((SingleStringOption*) Manager["directory"])->GetString();
  char* StateInputFile = ((SingleStringOption*) Manager["input"])->GetString();
  int Begin = ((SingleIntegerOption*) Manager["begin"])->GetInteger();
  int End = ((SingleIntegerOption*) Manager["end"])->GetInteger();

  int NbrState = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  bool VerboseFlag = ((BooleanOption*) Manager["verbose"])->GetBoolean();

  char** StateMatchedDirectories;
  int StateNbrMatchedDirectories = GetAllFilesDirectories(StateDirectory, StateMatchedDirectories);
  if (StateNbrMatchedDirectories == 0)
    {
      cout << "found no matching directories "  << StateDirectory << endl;
      return -1;
    }
  bool ErrorFlag = false;
  char** StateMatchedSpectra = new char* [StateNbrMatchedDirectories];
  char* StateSpectrumName = AddExtensionToFileName(StateInputFile, "raw");  
  for (int i = 0; i < StateNbrMatchedDirectories; ++i)
    {
      StateMatchedSpectra[i] = ConcatenatePathAndFileName(StateMatchedDirectories[i], StateSpectrumName);
      if (IsFile(StateMatchedSpectra[i]) == false)
	{
	  ErrorFlag = true;
	  cout << "can't open spectrum " << StateMatchedSpectra[i] << endl;
	}
      delete[] StateMatchedDirectories[i];
    }

  if (End == 0)
    {
      End = StateNbrMatchedDirectories - 1;
    }
  else
    {
      if (End > StateNbrMatchedDirectories)
	{
	  End = StateNbrMatchedDirectories;
	  cout << "warning : maximum requested sample has been reduced to " << End << endl;
	}      
    }
  if (Begin > End)
    {
      cout << "error : minimum requested sample (" << Begin << ") is higer than the maximum avalaible one (" << End << ")" << endl;
      ErrorFlag = true;
    }

  if (ErrorFlag == false)
    {  
      DOSSpectra SumDOS(End - Begin + 1, StateMatchedSpectra + Begin, NbrState, Gamma, Min, Max, Step);
      SumDOS.WriteSpectra(OutputFile);      
    }

   for (int i = 0; i < StateNbrMatchedDirectories; ++i)
   {
      delete[] StateMatchedSpectra[i];
    }
  delete[] StateMatchedDirectories;

  if (ErrorFlag == true)
    {
      return -1;
    }
  return 1;
}
