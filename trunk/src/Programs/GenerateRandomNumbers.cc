#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"

#include "Options/Options.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/time.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


// include code to test read-speed of vectors
#define TEST_SPEED

//#define TEST_COMPLEX


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

#ifdef TEST_SPEED
  OptionGroup* TestGroup = new OptionGroup ("test options");
  Manager += TestGroup;

  (*TestGroup) += new BooleanOption  ('\n', "test", "test I/O speed");
  (*TestGroup) += new SingleIntegerOption  ('\n', "size", "size for trial (in kb)", 100);
  (*TestGroup) += new SingleIntegerOption  ('\n', "trials", "number of times trial data is read", 10);
#endif

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


#ifdef TEST_SPEED
  if (Manager.GetBoolean("test"))
    {
#ifdef TEST_COMPLEX      
      long Size = (((long)Manager.GetInteger("size")) << 10)/sizeof(Complex);
      int Trials = Manager.GetInteger("trials");
      ComplexVector TestVector (Size);
      for (int i=0; i<Size; ++i)
	{
	  TestVector[i].Re=RandomNumber->GetRealRandomNumber();
	  TestVector[i].Im=RandomNumber->GetRealRandomNumber();
	}
#else
      long Size = (((long)Manager.GetInteger("size")) << 10)/sizeof(double);
      int Trials = Manager.GetInteger("trials");
      RealVector TestVector (Size);
      for (int i=0; i<Size; ++i)
	TestVector[i]=RandomNumber->GetRealRandomNumber();
#endif
      TestVector.Normalize();
      TestVector.WriteVector("test-vector.vec");
      TestVector.ReadVector("test-vector.vec");
      timeval TotalStartingTime2;
      timeval TotalEndingTime2;
      double Dt2, Dtb;
      cout << "------------------------------------------------------------------" << endl;
      cout << "start blocked" << endl;
      gettimeofday (&(TotalStartingTime2), 0);
      for (int t=0; t<Trials; ++t)
	{
	  TestVector.ReadVector("test-vector.vec");
	}
      gettimeofday (&(TotalEndingTime2), 0);
      Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
	((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
      cout << "time = " << Dt2 << endl;
      cout << "------------------------------------------------------------------" << endl;


      TestVector.ByteWriteVector("test-vector-b.vec");

      cout << "start bit-sized" << endl;
      gettimeofday (&(TotalStartingTime2), 0);
      for (int t=0; t<Trials; ++t)
	{
	  TestVector.ByteReadVector("test-vector.vec");
	}
      gettimeofday (&(TotalEndingTime2), 0);
      Dtb = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
	((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
      cout << "time = " << Dtb << endl;
      cout << "------------------------------------------------------------------" << endl;

      cout << "start blocked" << endl;
      gettimeofday (&(TotalStartingTime2), 0);
      for (int t=0; t<Trials; ++t)
	{
	  TestVector.ReadVector("test-vector.vec");
	}
      gettimeofday (&(TotalEndingTime2), 0);
      Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
	((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
      cout << "time = " << Dt2 << endl;
      cout << "------------------------------------------------------------------" << endl;


      TestVector.ByteWriteVector("test-vector-b.vec");

      gettimeofday (&(TotalStartingTime2), 0);
      cout << "start bit-sized" << endl;
      for (int t=0; t<Trials; ++t)
	{
	  TestVector.ByteReadVector("test-vector.vec");
	}
      gettimeofday (&(TotalEndingTime2), 0);
      Dtb = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
	((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
      cout << "time = " << Dtb << endl;
      cout << "------------------------------------------------------------------" << endl;

      
    }
#endif
  
  return 0;  
}
