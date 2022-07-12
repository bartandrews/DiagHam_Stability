#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <fstream>
#ifdef __SSTREAM_STYLE__
#include <sstream>
#else
#include <strstream>
#endif
#include <string>
#include <unistd.h>
#include <math.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;


// read a two column formatted array of numbers from a file
// 
// filename = name of  the file that conatins the spectrum (with optional relative/absolute path)
// xValues = reference ont the array where x values will be stored
// yValues = reference ont the array where y values will be stored
// nbrValues = reference on the number of values per column to retrieve from the file (0 it has to be automatically determined, and so is both array allocation)
// return value = true if no error occured
bool ReadDatas(char* filename, double*& xValues, double*& yValues, int& nbrValues);

// search a value in an array in a given range of indices where the array valuesare monotonous
//
// array = array in which to search in
// value = value to found
// min = lowest index to start search
// max = highest index to start search
// return value = smallest index such that the value is between array[index]  and array[index+1]
int SearchValueInArray (double* array, double value, int min, int max);


int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("EvaluateBroadening" , "0.01");

  OptionGroup* InputOptionGroup = new OptionGroup ("input options");
  OptionGroup* PrecisionOptionGroup = new OptionGroup ("precision options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += InputOptionGroup;
  Manager += PrecisionOptionGroup;
  Manager += MiscGroup;

  (*InputOptionGroup) += new SingleStringOption('\0', "input", "data input file (two column formatted text)");
  (*PrecisionOptionGroup) += new SingleDoubleOption('\n', "maximum-precision", "relative precision allowed for maximum detection", 0.01);
  (*PrecisionOptionGroup) += new SingleDoubleOption('\n', "half-height", "multiplicative factor to use with maximum value to find broadening", 0.5);

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if ((Manager.ProceedOptions(argv, argc, cout) == false) || ((SingleStringOption*) Manager["input"])->GetString() == 0)
    {
      cout << "see man page for option syntax or type EvaluateBroadening -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  char* InputFile = ((SingleStringOption*) Manager["input"])->GetString();
  double MaxPrecision = ((SingleDoubleOption*) Manager["maximum-precision"])->GetDouble();
  double HalfHeight = ((SingleDoubleOption*) Manager["half-height"])->GetDouble();
  

  double* XValues = 0;
  double* YValues = 0;
  int NbrValues = 0;

  if (ReadDatas(InputFile, XValues, YValues, NbrValues) == false)
    {
      return 1;
    } 

  double TrueMaximum = YValues[0];
  for (int i = 0; i < NbrValues; ++i)
    if (YValues[i] > TrueMaximum)
      TrueMaximum = YValues[i];
  double MaximumLowerBound = (1.0 - MaxPrecision) * TrueMaximum;
  double Maximum = 0.0;
  double LocalSlope = fabs(TrueMaximum / (XValues[1] - XValues[0]));
  double TmpSlope = 0.0;
  int MaximumPosition = 0;
  int ReducedNbrValues = NbrValues - 1;
  for (int i = 1; i < ReducedNbrValues; ++i)
    {
      if (YValues[i] > MaximumLowerBound)
	{
	  TmpSlope = fabs((YValues[i + 1] - YValues[i - 1])  / (XValues[i + 1] - XValues[i - 1]));
	  if (TmpSlope < LocalSlope)
	    {
	      Maximum = YValues[i];
	      LocalSlope = TmpSlope;
	      MaximumPosition = i;
	    }
	}
    }
  cout << "true maximum = " << TrueMaximum << endl << "guessed maximum = " << Maximum << " at " << XValues[MaximumPosition] << "(pos = " << MaximumPosition << ") with slope = " << LocalSlope << endl;

  double HalfMaximum = HalfHeight * Maximum;
  int LeftmostPosition = SearchValueInArray(YValues, HalfMaximum, 0, MaximumPosition);
  int RightmostPosition = SearchValueInArray(YValues, HalfMaximum, MaximumPosition, ReducedNbrValues);
  
  double Broadening = 0.5 * (XValues[RightmostPosition + 1] + XValues[RightmostPosition] - XValues[LeftmostPosition + 1] - XValues[LeftmostPosition]);
  double BroadeningError = 0.5 * (XValues[RightmostPosition + 1] - XValues[RightmostPosition] + XValues[LeftmostPosition + 1] - XValues[LeftmostPosition]);
  cout << "broadening = " << Broadening << " +/- " << BroadeningError << endl;

  delete[] XValues;
  delete[] YValues;
}

// read a two column formatted array of numbers from a file
// 
// filename = name of  the file that conatins the spectrum (with optional relative/absolute path)
// xValues = reference ont the array where x values will be stored
// yValues = reference ont the array where y values will be stored
// nbrValues = reference on the number of values per column to retrieve from the file (0 it has to be automatically determined, and so is both array allocation)
// return value = true if no error occured

bool ReadDatas(char* filename, double*& xValues, double*& yValues, int& nbrValues)
{
  if (nbrValues == 0)
    {
      ifstream File2;
      File2.open(filename, ios::out);
      if (!File2.is_open())
	{
	  cout << "error while opening file : " << filename << endl;
	  return false;
	}
      double Dummy1;
      double Dummy2;
      while (File2.tellg() >= 0)
	{
	  File2 >> Dummy1 >> Dummy2;
	  ++nbrValues;
	}      
      --nbrValues;
      File2.close();
      xValues = new double[nbrValues];
      yValues = new double[nbrValues];
    }
  ifstream File;
  File.open(filename, ios::out);
  if (!File.is_open())
    {
      cout << "error while opening file : " << filename << endl;
      return false;
    }
  for (int j = 0; j < nbrValues; ++j)
    {
      if (File.tellg() < 0)
	{
	  cout << filename <<  " has to few eigenvalues" << endl;
	  File.close();
	  return false;
	}
      else
	File >> xValues[j] >> yValues[j];
    }
  File.close();
  return true;  
}


// search a value in an array in a given range of indices where the array valuesare monotonous
//
// array = array in which to search in
// value = value to found
// min = lowest index to start search
// max = highest index to start search
// return value = smallest index such that the value is between array[index]  and array[index+1]

int SearchValueInArray (double* array, double value, int min, int max)
{
  int Tmp;
  if (array[min] < array[max])
    {      
      while ((max - min) > 1)
	{
	  Tmp = (max + min) >> 1;
	  if (array[Tmp] < value)
	    min = Tmp;
	  else
	    max = Tmp;
	}
    }
  else
    {
      while ((max - min) > 1)
	{
	  Tmp = (max + min) >> 1;
	  if (array[Tmp] > value)
	    min = Tmp;
	  else
	    max = Tmp;
	}
    }
  return min;
}
