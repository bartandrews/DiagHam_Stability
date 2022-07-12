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
#ifdef __SSTREAM_STYLE__
#include <sstream>
#else
#include <strstream>
#endif
#include <cstring>
#include <unistd.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("SumDOS" , "0.01");
  OptionGroup* SumDOSGroup = new OptionGroup ("SumDOS");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += SumDOSGroup;
  Manager += MiscGroup;

  (*SumDOSGroup) += new SingleStringOption('\n', "input", "name of the input file (directory mode only)", "eigenvalues");
  (*SumDOSGroup) += new SingleStringOption('\n', "prefix", "prefix of the directories containing spectrum (directory mode only)", "");
  (*SumDOSGroup) += new SingleIntegerOption('\n', "begin", "number of the first directory (directory mode only)", 0);
  (*SumDOSGroup) += new SingleIntegerOption('\n', "end", "number of the last directory (directory mode only)", 0);
  (*SumDOSGroup) += new SingleStringOption('\n', "file-prefix", "prefix of the files containing spectrum (file mode only)", "");
  (*SumDOSGroup) += new SingleStringOption('\n', "file-suffix", "suffix of the files containing spectrum (file mode only)", "");
  (*SumDOSGroup) += new SingleIntegerOption('\n', "nbr-files", "maximum number of files to look at  (file mode only, files must be named as [file-prefix]xxxx[file-suffix], where xxxx belongs non-contiguous list of integer ranging from 0 to nbr-files)", 2000);
  (*SumDOSGroup) += new SingleIntegerOption('n', "nbr-state", "number of states", 10);
  (*SumDOSGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in eV unit)", 0.0);
  (*SumDOSGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in eV unit)", 0.0);
  (*SumDOSGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in eV unit)", 0.01);
  (*SumDOSGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in eV unit) in the spectrum", 2e-4);
  (*SumDOSGroup) += new SingleStringOption('\n', "output", "name of the output file", "SumDOS.txt");
  (*SumDOSGroup) += new BooleanOption ('f', "file-mode", "use file mode for spectrum storage (i.e. all spectra in the same directory)  instead of directory mode (i.e. one spectrum per directory)");

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

  int NbrState = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  bool VerboseFlag = ((BooleanOption*) Manager["verbose"])->GetBoolean();


  if (((BooleanOption*) Manager["file-mode"])->GetBoolean() == false)
    {
      char* InputFile = ((SingleStringOption*) Manager["input"])->GetString();
      char* Prefix = ((SingleStringOption*) Manager["prefix"])->GetString();
      int Begin = ((SingleIntegerOption*) Manager["begin"])->GetInteger();
      int End = ((SingleIntegerOption*) Manager["end"])->GetInteger();
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
      DOSSpectra SumDOS(Number, Files, State, Gamma, Min, Max, Step);
      SumDOS.WriteSpectra(OutputFile);
    }
  else
    {
      int MaxNbrFiles = ((SingleIntegerOption*) Manager["nbr-files"])->GetInteger();
      char* FilePrefix = ((SingleStringOption*) Manager["file-prefix"])->GetString();
      char* FileSuffix = ((SingleStringOption*) Manager["file-suffix"])->GetString();
      int MaxStringLength = strlen (((SingleStringOption*) Manager["file-prefix"])->GetString()) + 10 + strlen (((SingleStringOption*) Manager["file-suffix"])->GetString());
      char** Files = new char* [MaxNbrFiles];
      int* State = new int [MaxNbrFiles];
      int RealNbrFiles = 0;
      char* TmpFileName = new char[MaxStringLength];
      while (MaxNbrFiles >= 0)
	{
	  sprintf (TmpFileName, "%s%d%s", FilePrefix, MaxNbrFiles, FileSuffix);
	  if (access(TmpFileName, R_OK) == 0)
	    {
	      Files[RealNbrFiles] = TmpFileName;
	      State[RealNbrFiles] = NbrState;
	      if (VerboseFlag)
		cout << TmpFileName << endl;
	      TmpFileName = new char[MaxStringLength];
	      ++RealNbrFiles;
	    }
	  --MaxNbrFiles;
	}      
      if (RealNbrFiles > 0)
	{
	  DOSSpectra SumDOS(RealNbrFiles, Files, State, Gamma, Min, Max, Step);
	  SumDOS.WriteSpectra(OutputFile);      
	}      
    }
  return 1;
}
