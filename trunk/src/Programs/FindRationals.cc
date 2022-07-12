#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/Options.h"

#include "MathTools/LongRational.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("FindRationals" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-file", "name of the file containing the double numbers");
  (*SystemGroup) += new SingleStringOption  ('o', "output-file", "name of the file where the closest rational numbers will be stored");
  (*SystemGroup) += new SingleIntegerOption  ('c', "column-index", "index of the column that contains the double numbers (0 being the first column)", 0);
  (*SystemGroup) += new BooleanOption  ('s', "std-output", "use standard output instead of an output file");
  (*SystemGroup) += new SingleIntegerOption  ('n', "max-denominator", "maximum denominator allowed for a rational number", 1000l);
  (*SystemGroup) += new SingleDoubleOption ('f', "filter-error", "if non zero, only shows lines of the input file for which the error between the double number and the closest rational is lower than a give threshold", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "shift-numbers", "shift all the double numbers by a given amount", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "rescale", "rescale all the double numbers", 1.0);
  (*SystemGroup) += new BooleanOption  ('\n', "add-index", "add an index to each output line");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FindRationals -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, an input file should be provided" << endl;
      return 0;
    }
  if ((Manager.GetString("output-file") == 0) && (Manager.GetBoolean("std-output") == false))
    {
      cout << "error, an output file should be provided or the strandard output should be used" << endl;
      return 0;
    }

  bool AddIndexFlag = Manager.GetBoolean("add-index");
  MultiColumnASCIIFile InputFile;
  if (InputFile.Parse(Manager.GetString("input-file")) == false)
    {
      InputFile.DumpErrors(cout);
      return false;
    }
  int NbrLines = InputFile.GetNbrLines();
  int NbrColumns = InputFile.GetNbrColumns();
  double* TmpNumbers = InputFile.GetAsDoubleArray(Manager.GetInteger("column-index"));
  LongRational TmpRational;
  double MaxError = Manager.GetDouble("filter-error");
  long MaxDenominator = Manager.GetInteger("max-denominator");
  double Shift = Manager.GetDouble("shift-numbers");
  double Rescale = Manager.GetDouble("rescale");
  if (Manager.GetBoolean("std-output") == false)
    {
      ofstream File;
      File.open(Manager.GetString("output-file"), ios::binary | ios::out); 
      File.precision(14);
      for (int i = 0; i < NbrLines; ++i)
	{
	  double Error = TmpRational.GetClosestRational((TmpNumbers[i] + Shift) * Rescale, MaxDenominator);
	  if ((MaxError == 0.0) || (fabs(Error) < MaxError))
	    {
	      if (AddIndexFlag == true)
		{
		  cout << i << " ";
		}
	      for (int j = 0; j < NbrColumns; ++j)
		{
		  File << InputFile(j, i) << " ";
		}
	      File << TmpRational << " " << Error << endl;
	    }
	}
      File.close();
    }
  else
    {
      for (int i = 0; i < NbrLines; ++i)
	{
	  double Error = TmpRational.GetClosestRational((TmpNumbers[i] + Shift) * Rescale, MaxDenominator);
	  if ((MaxError == 0.0) || (fabs(Error) < MaxError))
	    {
	      if (AddIndexFlag == true)
		{
		  cout << i << " : ";
		}
	      for (int j = 0; j < NbrColumns; ++j)
		{
		  cout << InputFile(j, i) << " ";
		}
	      cout << TmpRational << " " << Error << endl;
	    }
	}
    }
  return 0;
}

