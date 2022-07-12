#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/SingleStringOption.h"
#include "Options/BooleanOption.h"

#include "Vector/ComplexVector.h"

#include "MathTools/Complex.h"

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
  OptionManager Manager ("OverlapComplex" , "0.01");
  OptionGroup* FileGroup = new OptionGroup ("File");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += FileGroup;
  Manager += MiscGroup;

  (*FileGroup) += new SingleStringOption('\n', "file1", "name of the first input file");
  (*FileGroup) += new SingleStringOption('\n', "file2", "name of the second input file");

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

  char* File1 = ((SingleStringOption*) Manager["file1"])->GetString();
  char* File2 = ((SingleStringOption*) Manager["file2"])->GetString();
  
  ComplexVector vector1, vector2;
  vector1.ReadVector (File1); vector2.ReadVector (File2); 

  Complex c = vector1 * vector2;

  cout << Norm (c);

  return 1;
}
