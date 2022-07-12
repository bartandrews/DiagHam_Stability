#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
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
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('q', "nbr-flux", "number of flux quanta", 8);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusShowBasis -h" << endl;
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

  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
      for (int x = 0; x < MomentumModulo; ++x)
	for (int y = 0; y < MomentumModulo; ++y)
	  {
	    BosonOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, x, y);
	    cout << " (k_x = " << x << ", k_y = " << y << ") : " << endl;
	    for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
	      Space.PrintState(cout, i) << endl;
	    cout << endl;
	  }
    }
  else
    {
      for (int x = 0; x < MomentumModulo; ++x)
	for (int y = 0; y < MomentumModulo; ++y)
	  {
	    FermionOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, x, y);
	    cout << " (k_x = " << x << ", k_y = " << y << ") : " << endl;
	    for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
	      Space.PrintState(cout, i) << endl;;
	    cout << endl;
	  }
    }
}

