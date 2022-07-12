#include "Tools/FQHEFiles/FQHEOnDiskFileTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/Options.h"

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);  

  OptionManager Manager ("FQHEDiskExcessCharge" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "density", "name of the file containing the sampled density");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nu-numerator", "numerator of the filling factor", 1);
  (*SystemGroup) += new SingleIntegerOption  ('d', "nu-denominator", "denominator of the filling factor", 1);
  (*SystemGroup) += new SingleIntegerOption  ('s', "shift", "shift bewteen nu^{-1} N and N_phi", 1);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the rho extension with rhorho extension)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskExcessCharge -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = 0;
  bool Statistics = true;

  if (FQHEOnDiskFindSystemInfoFromFileName(Manager.GetString("density"), NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("density") << endl;
      return -1;
    }

  cout << NbrParticles << " " << LzMax << endl;
  int FillingFactorNumerator =  Manager.GetInteger("nu-numerator");
  int FillingFactorDenominator =  Manager.GetInteger("nu-denominator");
  double FillingFactor = ((double) FillingFactorNumerator) / ((double) FillingFactorDenominator);
  int Shift = Manager.GetInteger("shift");


  MultiColumnASCIIFile Data;
  
  if (Data.Parse(Manager.GetString("density")) == false)
    {
      Data.DumpErrors(cout) << endl;
      return -1;
    } 
  double* XValues = Data.GetAsDoubleArray(0);
  double* YValues = Data.GetAsDoubleArray(1);
  int NbrValues = Data.GetNbrLines();

  double Factor = FillingFactor * M_PI * ((double) (LzMax)) / ((double) NbrParticles);
  double MeanValue = ((double) NbrParticles)  / (M_PI * ((double) (2 * LzMax))) * ((double) NbrParticles) / (FillingFactor * ((double) (LzMax)));
  double DropletRadius = sqrt((double) (2 * LzMax)); 
  for (int i = 0; i < NbrValues; ++i)
    {
      if (fabs(XValues[i]) < DropletRadius)
	YValues[i] = MeanValue - YValues[i];
      else
	YValues[i] *= -1.0;
    }
  double Integral = 0.0;
  cout << 0.0 << " " << 0.0 << endl;
  int i = 0;
  while ((i < NbrValues) && (XValues[i] < 0.0))
    ++i;
  for (; i < NbrValues; ++i)
    {
      // Integral += ((XValues[i - 1] * YValues[i - 1] * exp(-0.5 * XValues[i - 1] * XValues[i - 1])) + (XValues[i] * YValues[i] * exp(-0.5 * XValues[i] * XValues[i]))) * (XValues[i] - XValues[i - 1]);
      Integral += ((XValues[i - 1] * YValues[i - 1]) + (XValues[i] * YValues[i])) * (XValues[i] - XValues[i - 1]);
      cout << XValues[i] << " " << (Integral * Factor) << endl;
    }
  delete[] XValues;
  delete[] YValues;
  return 0;
}

