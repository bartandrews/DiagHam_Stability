#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("VectorBinary2Ascii" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the binary vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the vector will be stored in ASCII mode");
  (*SystemGroup) += new BooleanOption  ('\n', "add-index", "write ascii vector in a two-column formatted output (first column for the component index, second for the component value)");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorBinary2Ascii -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["input-vector"])->GetString() == 0)
    {
      cout << "VectorBinary2Ascii requires an input file" << endl << "see man page for option syntax or type VectorBinary2Ascii -h" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["input-vector"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-vector"])->GetString() << endl;
      return -1;      
    }
  if (((SingleStringOption*) Manager["output-vector"])->GetString() == 0)
    {
      cout << "VectorBinary2Ascii requires an output file" << endl << "see man page for option syntax or type VectorBinary2Ascii -h" << endl;
      return -1;
    }

  ofstream File;
  File.open(((SingleStringOption*) Manager["output-vector"])->GetString(), ios::out);
  File.precision(14);
  if (((BooleanOption*) Manager["add-index"])->GetBoolean() == true)
    for (int i = 0; i < State.GetVectorDimension(); ++i)
      File << i << " " << State[i] << endl;
  else
    for (int i = 0; i < State.GetVectorDimension(); ++i)
      File << State[i] << endl;
  File.close();
  
}
