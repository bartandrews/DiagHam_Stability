#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/Spectra/QuantumWellBFieldAbsorptionSpectra.h"

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
#include <dirent.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;


#define ME_KB 11.60451


int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("AbsorptionInQuantumWellBField" , "0.01");
  OptionGroup* InputOptionGroup = new OptionGroup ("input options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += InputOptionGroup;
  Manager += PlotOptionGroup;
  Manager += MiscGroup;

  (*InputOptionGroup) += new SingleStringOption('\n', "initial-dir", "pattern to identify directories that contain initial states description (if pattern is mydir/run, any directory in mydir such that runx where x is an integer will be considered)");
  (*InputOptionGroup) += new SingleStringOption('\n', "final-dir", "pattern to identify directories that contain final states description (using same conventions than initial-directory)");
  (*InputOptionGroup) += new SingleStringOption('\n', "initial-input", "beginning of the name of file describing initial states (without relative/absolute path, files has to be of the form NAMEraw for eigenvalues and NAMEx.vec for eigenvectors where x is an integer");
  (*InputOptionGroup) += new SingleStringOption('\n', "final-input", "beginning of the name of file describing final states (using same conventions than initial-input, use same name as initial-input if not defined)");
  (*InputOptionGroup) += new SingleIntegerOption('\n', "begin", "number of the first directory", 0);
  (*InputOptionGroup) += new SingleIntegerOption('\n', "end", "number of the last directory (0 if it has to be detected automatically)", 0);
  (*InputOptionGroup) += new SingleIntegerOption('n', "nbr-state", "number of initial states (0 if it has to be detected automatically)", 0);
  (*InputOptionGroup) += new SingleDoubleOption('z', "z-size", "sample size ine the z direction (in Angstrom)", 58.7);

  (*PlotOptionGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in meV unit)", 0.0);
  (*PlotOptionGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in meV unit)", 0.0);
  (*PlotOptionGroup) += new SingleDoubleOption('t', "temperature", "temperature (in Kelvin unit, any negative value for infinite temperature)", -1.0);
  (*PlotOptionGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in meV unit)", 0.01);
  (*PlotOptionGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in meV unit) in the spectrum", 2e-4);
  (*PlotOptionGroup) += new SingleStringOption('\n', "output", "name of the output file", "absorption.dat");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new BooleanOption ('v', "verbose", "verbose mode", false);

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type AsorptionInQuantumWellBFiled -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  char* InitialStateDirectory = ((SingleStringOption*) Manager["initial-dir"])->GetString();
  char* FinalStateDirectory = ((SingleStringOption*) Manager["final-dir"])->GetString();
  char* InitialStateInputFile = ((SingleStringOption*) Manager["initial-input"])->GetString();
  char* FinalStateInputFile = ((SingleStringOption*) Manager["final-input"])->GetString();
  if (FinalStateInputFile == 0)
    {
      FinalStateInputFile = InitialStateInputFile;
    }
  int Begin = ((SingleIntegerOption*) Manager["begin"])->GetInteger();
  int End = ((SingleIntegerOption*) Manager["end"])->GetInteger();


  int NbrStates = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Beta =  ME_KB / ((SingleDoubleOption*) Manager["temperature"])->GetDouble();
  if (Beta < 0.0)
    Beta = 0.0;
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  double ZSize = ((SingleDoubleOption*) Manager["z-size"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  char** InitialStateMatchedDirectories;
  int InitialStateNbrMatchedDirectories = GetAllFilesDirectories(InitialStateDirectory, InitialStateMatchedDirectories);
  if (InitialStateNbrMatchedDirectories == 0)
    {
      cout << "found no matching directories "  << InitialStateDirectory << endl;
      return -1;
    }
  char** FinalStateMatchedDirectories;
  int FinalStateNbrMatchedDirectories = GetAllFilesDirectories(FinalStateDirectory, FinalStateMatchedDirectories);
  if (FinalStateNbrMatchedDirectories == 0)
    {
      
      cout << "found no matching directories "  << FinalStateDirectory << endl;
      return -1;
    }

  bool ErrorFlag = false;
  char** InitialStateMatchedSpectra = new char* [InitialStateNbrMatchedDirectories];
  char*** InitialStateMatchedEigenvectors = new char** [InitialStateNbrMatchedDirectories];
  int* InitialStateNbrMatchedEigenvectors = new int [InitialStateNbrMatchedDirectories];
  int InitialStateNbrEigenvectors = 0;
  char* InitialStateSpectrumName = AddExtensionToFileName(InitialStateInputFile, "raw");  
  for (int i = 0; i < InitialStateNbrMatchedDirectories; ++i)
    {
      char* TmpPattern = ConcatenatePathAndFileName(InitialStateMatchedDirectories[i], InitialStateInputFile);
      InitialStateNbrMatchedEigenvectors[i] = GetAllFilesDirectories(TmpPattern, InitialStateMatchedEigenvectors[i], ".vec"); 
      if (InitialStateNbrEigenvectors == 0)
	{
	  InitialStateNbrEigenvectors = InitialStateNbrMatchedEigenvectors[i];
	}
      else
	if (InitialStateNbrEigenvectors != InitialStateNbrMatchedEigenvectors[i])
	  {
	    ErrorFlag = true;
	    cout << "wrong number of eigenvectors in " << InitialStateMatchedDirectories[i] << "(found " << InitialStateNbrMatchedEigenvectors[i] << ", need " << InitialStateNbrEigenvectors << ")" << endl;
	  }
      InitialStateMatchedSpectra[i] = ConcatenatePathAndFileName(InitialStateMatchedDirectories[i], InitialStateSpectrumName);
      if (IsFile(InitialStateMatchedSpectra[i]) == false)
	{
	  ErrorFlag = true;
	  cout << "can't open spectrum " << InitialStateMatchedSpectra[i] << endl;
	}
      delete[] InitialStateMatchedDirectories[i];
      delete[] TmpPattern;
    }
  char** FinalStateMatchedSpectra = new char* [FinalStateNbrMatchedDirectories];
  char*** FinalStateMatchedEigenvectors = new char** [FinalStateNbrMatchedDirectories];
  int* FinalStateNbrMatchedEigenvectors = new int [FinalStateNbrMatchedDirectories];
  int FinalStateNbrEigenvectors = 0;
  char* FinalStateSpectrumName = AddExtensionToFileName(FinalStateInputFile, "raw");  
  for (int i = 0; i < FinalStateNbrMatchedDirectories; ++i)
    {
      char* TmpPattern = ConcatenatePathAndFileName(FinalStateMatchedDirectories[i], FinalStateInputFile);
      FinalStateNbrMatchedEigenvectors[i] = GetAllFilesDirectories(TmpPattern, FinalStateMatchedEigenvectors[i], ".vec"); 
      if (FinalStateNbrEigenvectors == 0)
	{
	  FinalStateNbrEigenvectors = FinalStateNbrMatchedEigenvectors[i];
	}
      else
	if (FinalStateNbrEigenvectors != FinalStateNbrMatchedEigenvectors[i])
	  {
	    ErrorFlag = true;
	    cout << "wrong number of eigenvectors in " << FinalStateMatchedDirectories[i] << "(found " << FinalStateNbrMatchedEigenvectors[i] << ", need " << FinalStateNbrEigenvectors << ")" << endl;
	  }
      FinalStateMatchedSpectra[i] = ConcatenatePathAndFileName(FinalStateMatchedDirectories[i], FinalStateSpectrumName);
      if (IsFile(FinalStateMatchedSpectra[i]) == false)
	{
	  ErrorFlag = true;
	  cout << "can't open spectrum " << FinalStateMatchedSpectra[i] << endl;
	}
      delete[] FinalStateMatchedDirectories[i];
      delete[] TmpPattern;
    }

  if (End == 0)
    {
      if (FinalStateNbrMatchedDirectories >= InitialStateNbrMatchedDirectories)
	End = InitialStateNbrMatchedDirectories;
      else
	End = FinalStateNbrMatchedDirectories;
      --End;
    }
  else
    {
      if (End > FinalStateNbrMatchedDirectories)
	{
	  End = FinalStateNbrMatchedDirectories;
	  cout << "warning : maximum requested sample has been reduced to " << End << endl;
	}
      if (End > InitialStateNbrMatchedDirectories)
	{
	  End = InitialStateNbrMatchedDirectories;
	  cout << "warning : maximum requested sample has been reduced to " << End << endl;
	}      
    }
  if (Begin > End)
    {
      cout << "error : minimum requested sample (" << Begin << ") is higer than the maximum avalaible one (" << End << ")" << endl;
      ErrorFlag = true;
    }

  if (NbrStates == 0)
    {
      NbrStates = InitialStateNbrEigenvectors;      
      if (FinalStateNbrEigenvectors != (2 * InitialStateNbrEigenvectors))
	{
	  cout << "final states must have twice the number of initial states" << endl;
	  ErrorFlag = true;
	}
    }
  else
    {
      if (InitialStateNbrEigenvectors < NbrStates)
	{
	  cout << "not enough initial states" << endl;
	  ErrorFlag = true;
	}
      if (FinalStateNbrEigenvectors < (2 * NbrStates))
	{
	  cout << "not enough final states" << endl;
	  ErrorFlag = true;
	}
    }

  if (ErrorFlag == false)
    {  
      QuantumWellBFieldAbsorptionSpectra AbsorptionSpectrum(End - Begin + 1, NbrStates, InitialStateMatchedSpectra, InitialStateMatchedEigenvectors,
							    2 * NbrStates, FinalStateMatchedSpectra, FinalStateMatchedEigenvectors,
							    0.0, 0.0, ZSize, Gamma, Beta, Min, Max, Step);
      AbsorptionSpectrum.WriteSpectra(OutputFile);    
    }

  for (int i = 0; i < InitialStateNbrMatchedDirectories; ++i)
    {
      for (int j = 0; j < InitialStateNbrMatchedEigenvectors[i]; ++j)
	delete[] InitialStateMatchedEigenvectors[i][j];
      delete[] InitialStateMatchedEigenvectors[i];
      delete[] InitialStateMatchedSpectra[i];
    }
  delete[] InitialStateMatchedDirectories;
  delete[] InitialStateNbrMatchedEigenvectors;
  for (int i = 0; i < FinalStateNbrMatchedDirectories; ++i)
    {
      for (int j = 0; j < FinalStateNbrMatchedEigenvectors[i]; ++j)
	delete[] FinalStateMatchedEigenvectors[i][j];
      delete[] FinalStateMatchedEigenvectors[i];
      delete[] FinalStateMatchedSpectra[i];
    }
  delete[] FinalStateMatchedDirectories;
  delete[] FinalStateNbrMatchedEigenvectors;

  if (ErrorFlag == true)
    {
      return -1;
    }
  return 1;
}


