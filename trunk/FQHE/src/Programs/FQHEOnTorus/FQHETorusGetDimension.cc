#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusGetDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('q', "nbr-flux", "number of flux quanta", 8);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistics (default value))");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_torus_n_nbrparticles_q_nbrfluxquanta.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusGetDimension -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int MomentumModulo = FindGCD(NbrParticles, NbrFluxQuanta);
  long** Dimensions = new long* [MomentumModulo];
  for (int x = 0; x < MomentumModulo; ++x)
    Dimensions[x] = new long[MomentumModulo];
  long TotalDimension = 0l;

  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
      for (int x = 0; x < MomentumModulo; ++x)
	for (int y = 0; y < MomentumModulo; ++y)
	  {
	    BosonOnTorusWithMagneticTranslationsShort Space (NbrParticles, NbrFluxQuanta, x, y);
	    Dimensions[x][y] = Space.GetHilbertSpaceDimension();
	    TotalDimension += Space.GetHilbertSpaceDimension();
	  }
    }
  else
    {
      for (int x = 0; x < MomentumModulo; ++x)
	for (int y = 0; y < MomentumModulo; ++y)
	  {
	    FermionOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, x, y);
	    Dimensions[x][y] = Space.GetHilbertSpaceDimension();
	    TotalDimension += Space.GetHilbertSpaceDimension();
	  }
    }


  if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
    {
      char* OutputFileName = 0;
      if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
	{
	  OutputFileName = new char[256];
	  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	    sprintf (OutputFileName, "bosons_torus_n_%d_q_%d.dim", NbrParticles, NbrFluxQuanta);
	  else
	    sprintf (OutputFileName, "fermions_torus_n_%d_q_%d.dim", NbrParticles, NbrFluxQuanta);
	}
      else
	{
	  OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
	  strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
	}
      ofstream File;
      File.open(OutputFileName, ios::binary | ios::out);
      File << "# Hilbert space dimension in each momemtum sector for " << NbrParticles << " ";
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	File << "bosons";
      else
	File << "femions";
      File << " on the torus geometry with " << NbrFluxQuanta << " flux quanta" << endl;
      File << "# total Hilbert space dimension = " << TotalDimension << endl << "#" << endl << "#  kx  ky  dimension" << endl;
      for (int x = 0; x < MomentumModulo; ++x)
	for (int y = 0; y < MomentumModulo; ++y)
	  File << x << " "  << y << " " << Dimensions[x][y] << endl;
      File.close();
      delete[] OutputFileName;
    }
  else
    for (int x = 0; x < MomentumModulo; ++x)
      for (int y = 0; y < MomentumModulo; ++y)
	cout << " (k_x = " << x << ", k_y = " << y << ") : " << Dimensions[x][y] << endl;

  for (int x = 0; x < MomentumModulo; ++x)
    delete[] Dimensions[x];
  delete[] Dimensions;

  return 0;
}

