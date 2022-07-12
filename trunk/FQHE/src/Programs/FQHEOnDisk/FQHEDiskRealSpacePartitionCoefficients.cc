#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#ifdef __GSL__
#include <gsl/gsl_sf_gamma.h>
#endif

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital
// radius = radius of the region A centered at the origin
// return value = square of the orbital weight
double FQHEDiskComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double radius);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEDiskRealSpacePartitionCoefficients" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta", 10);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "radius", "disk radius of the region A centered around the origin ", 5.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "optional output file name (default is realspace_disk_radius_*_2s_*.dat)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskRealSpacePartitionCoefficients -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");

  double Radius = Manager.GetDouble("radius");
  double CutPosition = 0.0;
  char* OutputFile = 0;
  if (Manager.GetString("output-file") == 0)
    {
      OutputFile = new char[512];
      sprintf (OutputFile, "realspace_disk_radius_%.6f_2s_%d.dat", Radius, NbrFluxQuanta);
    }
  else
    {
      OutputFile = new char[strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFile, Manager.GetString("output-file"));
    }
  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  File << "# real space coefficients for disk cut of radius " << Radius << " on a disk with N_phi=" << NbrFluxQuanta << endl
       << "OrbitalSquareWeights =";
  int NbrCoefficients = 0;
  double* Coefficients = 0;
  Coefficients = new double [NbrFluxQuanta + 1];
  for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
    {
      Coefficients[NbrCoefficients] = FQHEDiskComputeSharpRealSpaceCutCoefficient(NbrCoefficients, Radius);
    }
  for (int i = 0; i < NbrCoefficients; ++i)
    File << " " << Coefficients[i];
  File << endl;
  File.close();
  delete[] Coefficients;
  return 0;
}


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital
// radius = radius of the region A centered at the origin
// return value = square of the orbital weight

double FQHEDiskComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double radius)
{  
  radius *= radius;
  radius *= 0.5;
#ifdef __GSL__
  return gsl_sf_gamma_inc_P((double) (orbitalIndex + 1), radius);
#else
  cout << "gsl is required" << endl;
  return 0.0;
#endif
}
