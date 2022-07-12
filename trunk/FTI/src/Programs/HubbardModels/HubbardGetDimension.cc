// #include "HilbertSpace/FermionOnLatticeWithSpinMomentumSpace.h"
// #include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"

#include "MathTools/BinomialCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// nbrSites = number of sites
//gutzwillerFlag = false if double occupations are allowed
// return value = Hilbert space dimension
long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int nbrSites, bool gutzwillerFlag, int nbrSpinUp = -1);

// save dimensions in a given file
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrSites = number of sites
// statistics = true for bosons, false for fermions
// totalDimension = total Hilbert space dimension
bool WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrSites, bool statistics, long totalDimension);



int main(int argc, char** argv)
{
  OptionManager Manager ("HubbardGetDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "number of flux quanta", 20);
  (*SystemGroup) += new BooleanOption  ('\n', "conserve-sz", "Sz is a good quantum number");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-spin", "twice the total spin value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_n_nbrparticles_q_nbrfluxquanta.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereGetDimension -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int TotalSz = Manager.GetInteger("total-spin");
  if ((Manager.GetBoolean("conserve-sz")) && ((TotalSz % 2) != (NbrParticles % 2)))
  {
   cout << "Number of particles and total spin should have the same parity" << endl;
   return 0;
  }   
  int NbrSpinUp = (TotalSz + NbrParticles) >> 1;
  if ((Manager.GetBoolean("conserve-sz")) && (NbrSpinUp < 0))
  {
   cout << "Number of up spins is negative" << endl;
   return 0;
  }
  
  int NbrSites = Manager.GetInteger("nbr-sites"); 
  bool GutzwillerFlag = Manager.GetBoolean("gutzwiller");
  
  long TotalDimension;
  if (Manager.GetBoolean("conserve-sz") == false)
    TotalDimension = FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrSites, GutzwillerFlag);
  else
    TotalDimension = FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrSites, GutzwillerFlag, NbrSpinUp);
  
  if (Manager.GetBoolean("save-disk") == true)
    {
      char* OutputFileName = 0;
      if (Manager.GetString("output-file") == 0)
	{
	  OutputFileName = new char[256];
	  if (Manager.GetBoolean("boson") == true)
	    sprintf (OutputFileName, "bosons_hubbard_n_%d_2s_%d.dim", NbrParticles, NbrSites);
	  else
	    if (GutzwillerFlag == false)
	      if (Manager.GetBoolean("conserve-sz") == false)
		sprintf (OutputFileName, "fermions_hubbard_n_%d_2s_%d.dim", NbrParticles, NbrSites);
	      else
		sprintf (OutputFileName, "fermions_hubbard_n_%d_2s_%d_sz_%d.dim", NbrParticles, NbrSites, TotalSz);
	    else
	      sprintf (OutputFileName, "fermions_hubbard_gutzwiller_n_%d_2s_%d.dim", NbrParticles, NbrSites);
	}
      else
	{
	  OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFileName, Manager.GetString("output-file"));
	}
      WriteDimensionToDisk (OutputFileName, NbrParticles, NbrSites, Manager.GetBoolean("boson"),
				    TotalDimension);
      delete[] OutputFileName;
    }
  else
    {
      cout << "Hilbert space dimension = " << TotalDimension << endl;
    }
	  
  return 0;
}
      

// evaluate Hilbert space dimension for bosons
//
// nbrFermions = number of fermions
//nbrSites = number of sites
// return value = Hilbert space dimension

long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int nbrSites, bool gutzwillerFlag, int nbrSpinUp)
{
  long dimension;
  if (gutzwillerFlag == false)
  {
    if (nbrSpinUp == -1)
    {
      BinomialCoefficients binomials(2*nbrSites);
      dimension = binomials(2*nbrSites, nbrFermions); 
    }
    else
    {
     BinomialCoefficients binomials(nbrSites);
     dimension = binomials(nbrSites, nbrSpinUp); 
     dimension *= binomials(nbrSites, (nbrFermions - nbrSpinUp));
    }
  }
  else
  {
     if (nbrSpinUp == -1)
    {
      BinomialCoefficients binomials(nbrSites);
      int NbrHoles = nbrSites - nbrFermions;
      dimension = binomials(nbrSites, NbrHoles);
      for (int i = 0; i < nbrFermions; ++i)
	dimension *= 2l;
    }
    else
    {
      BinomialCoefficients binomials(nbrSites);
      int NbrHoles = nbrSites - nbrFermions;
      dimension = binomials(nbrSites, NbrHoles);
      BinomialCoefficients binomials1(nbrFermions);
      dimension *= binomials1(nbrFermions, nbrSpinUp);
    }
  }
  
  return dimension;
}


// save dimensions in a given file
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrsites = number of sites
// statistics = true for bosons, false for fermions
// totalDimension = total Hilbert space dimension

bool WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrSites, bool statistics,
			  long totalDimension)
{
  ofstream File;
  File.open(outputFileName, ios::binary | ios::out);
  File << "# Hilbert space dimension for " << nbrParticles << " ";
  if (statistics == true)
    File << "bosons";
  else
    File << "femions";
  File << " for the Hubbard model with " << nbrSites << " sites" << endl;
  File << "# total Hilbert space dimension = " << totalDimension << endl << endl 
       << "N = " << nbrParticles << endl
       << "2S = " << nbrSites << endl;
  File.close();
  return true;
}