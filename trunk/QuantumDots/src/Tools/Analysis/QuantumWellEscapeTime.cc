#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/QuantumWellBFieldEscapeProbability.h"

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
  OptionManager Manager ("QuantumWellEscapeTime" , "0.01");
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
  (*InputOptionGroup) += new SingleIntegerOption('\n', "initial-index", "index of the initial state (-1 if probality has to evaluated for all possible states)", -1);
  (*InputOptionGroup) += new BooleanOption('\n', "subband", "compute the probability to stay in the same subband that the initial state instead of the probability to stay in the initial state");

  (*PlotOptionGroup) += new SingleDoubleOption('t', "time-step", "time step (in hbar/E unit)", 0.5);
  (*PlotOptionGroup) += new SingleIntegerOption('\n', "nbr-step", "number of time steps", 10000);
  (*PlotOptionGroup) += new SingleStringOption('\n', "output", "name of the output file", "absorption.dat");
  (*PlotOptionGroup) += new BooleanOption('\n', "log", "plot -log of the probability");

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

  int NbrStates = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double TimeStep = ((SingleDoubleOption*) Manager["time-step"])->GetDouble();
  int NbrTimeSteps  = ((SingleIntegerOption*) Manager["nbr-step"])->GetInteger();
  int InitialStateIndex  = ((SingleIntegerOption*) Manager["initial-index"])->GetInteger();
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
  char*** StateMatchedEigenvectors = new char** [StateNbrMatchedDirectories];
  int* StateNbrMatchedEigenvectors = new int [StateNbrMatchedDirectories];
  int StateNbrEigenvectors = 0;
  char* StateSpectrumName = AddExtensionToFileName(StateInputFile, "raw");  
  for (int i = 0; i < StateNbrMatchedDirectories; ++i)
    {
      char* TmpPattern = ConcatenatePathAndFileName(StateMatchedDirectories[i], StateInputFile);
      StateNbrMatchedEigenvectors[i] = GetAllFilesDirectories(TmpPattern, StateMatchedEigenvectors[i], ".vec"); 
      if (StateNbrEigenvectors == 0)
	{
	  StateNbrEigenvectors = StateNbrMatchedEigenvectors[i];
	}
      else
	if (StateNbrEigenvectors != StateNbrMatchedEigenvectors[i])
	  {
	    ErrorFlag = true;
	    cout << "wrong number of eigenvectors in " << StateMatchedDirectories[i] 
		 << "(found " << StateNbrMatchedEigenvectors[i] << ", need " << StateNbrEigenvectors << ")" << endl;
	  }
      StateMatchedSpectra[i] = ConcatenatePathAndFileName(StateMatchedDirectories[i], StateSpectrumName);
      if (IsFile(StateMatchedSpectra[i]) == false)
	{
	  ErrorFlag = true;
	  cout << "can't open spectrum " << StateMatchedSpectra[i] << endl;
	}
      delete[] StateMatchedDirectories[i];
      delete[] TmpPattern;
    }

  if (End == 0)
    End = StateNbrMatchedDirectories - 1;

  if (Begin > End)
    {
      cout << "error : minimum requested sample (" << Begin << ") is higer than the maximum avalaible one (" << End << ")" << endl;
      ErrorFlag = true;
    }

  if (NbrStates == 0)
    {
      NbrStates = StateNbrEigenvectors;      
    }
  else
    {
      if (StateNbrEigenvectors < NbrStates)
	{
	  cout << "not enough initial states" << endl;
	  ErrorFlag = true;
	}
    }

  cout << "found " << (End - Begin + 1) << " samples with " << NbrStates << " states each" << endl;
  if (ErrorFlag == false)
    {  
      QuantumWellBFieldEscapeProbability EscapeProbability(End - Begin + 1, NbrStates >> 1, StateMatchedSpectra, StateMatchedEigenvectors,
							   TimeStep, NbrTimeSteps, InitialStateIndex, ((BooleanOption*) Manager["subband"])->GetBoolean(),
							   ((BooleanOption*) Manager["log"])->GetBoolean());
      EscapeProbability.WriteSpectra(OutputFile);    
    }

  for (int i = 0; i < StateNbrMatchedDirectories; ++i)
    {
      for (int j = 0; j < StateNbrMatchedEigenvectors[i]; ++j)
	delete[] StateMatchedEigenvectors[i][j];
      delete[] StateMatchedEigenvectors[i];
      delete[] StateMatchedSpectra[i];
    }
  delete[] StateMatchedDirectories;
  delete[] StateNbrMatchedEigenvectors;

  if (ErrorFlag == true)
    {
      return -1;
    }
  return 1;
}
