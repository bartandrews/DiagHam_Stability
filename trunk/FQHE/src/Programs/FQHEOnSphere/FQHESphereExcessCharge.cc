#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

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

  OptionManager Manager ("FQHESphereExcessCharge" , "0.01");
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
      cout << "see man page for option syntax or type FQHESphereExcessCharge -h" << endl;
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

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("density"), NbrParticles, LzMax, TotalLz, Statistics) == false)
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
  double MeanValue = ((double) (NbrParticles * (LzMax + Shift))) / (4.0 * M_PI * ((double) (LzMax)));
  for (int i = 0; i < NbrValues; ++i)
    YValues[i] = MeanValue - YValues[i];

  double Integral = 0.0;
  cout << 0.0 << " " << 0.0 << endl;
  for (int i = 1; i < NbrValues; ++i)
    {
      Integral += ((sin(XValues[i - 1]) * YValues[i - 1]) + (sin(XValues[i]) * YValues[i])) * (XValues[i] - XValues[i - 1]);
      cout << XValues[i] << " " << (Integral * Factor) << endl;
    }
  Integral *= M_PI;
  delete[] XValues;
  delete[] YValues;
  return 0;
}

