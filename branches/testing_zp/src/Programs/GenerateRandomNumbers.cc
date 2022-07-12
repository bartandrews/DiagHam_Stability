#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"

#include "Options/Options.h"


#include <iostream>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("GenerateRandomNumbers" , "0.01");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += MonteCarloGroup;
  Manager += MiscGroup;

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-random", "number of random numbers to generate", 1000000);
  (*MonteCarloGroup) += new SingleStringOption ('o', "random-file", "name of the file where the random numbers will be stored", "random.dat");
  (*MonteCarloGroup) += new SingleStringOption  ('\n', "random-generator", "name of the random number to use (StdlibRandomNumberGenerator or NumRecRandomGenerator)", "StdlibRandomNumberGenerator");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenerateRandomNumbers -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  AbstractRandomNumberGenerator* RandomNumber = 0;
  if (strcasecmp(Manager.GetString("random-generator"), "NumRecRandomGenerator") == 0)
    {
      cout << "Using NumRecRandomGenerator"<<endl;
      RandomNumber = new NumRecRandomGenerator (29457);      
    }
  else
    {
      cout << "Using StdlibRandomNumberGenerator"<<endl;
      RandomNumber = new StdlibRandomNumberGenerator (29457);
    }
  FileRandomNumberGenerator RandomNumber2 (RandomNumber, Manager.GetInteger("nbr-random"), ((SingleStringOption*) Manager["random-file"])->GetString());
  return 0;  
}
