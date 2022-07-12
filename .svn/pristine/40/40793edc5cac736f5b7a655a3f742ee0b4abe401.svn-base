#include "HilbertSpace/BosonOnSphere.h"

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

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension
long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for fermions
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension
long FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

// evaluate Hilbert space dimension for fermions with SU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension
long FermionSU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

// evaluate Hilbert space dimension for fermions with SU(2)xSU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// return value = Hilbert space dimension
long FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin);

// evaluate Hilbert space dimension for fermions with SU(4) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// totalEntanglement = number of particles with entanglement plus
// return value = Hilbert space dimension
long FermionSU4ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement);

// evaluate Hilbert space dimension for fermions with SU(3) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension
long FermionSU3ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3);

// save dimensions in a given file
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array taht contains dimension of each Lz sector
// lDimensions = array that contains dimension of each L sector (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value
// lzMax = twice the maximum Lz value
// totalDimension = total Hilbert space dimension
bool WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			  long* lzDimensions, long* lDimensions, int lzMin, int lzMax, long totalDimension);

// save dimensions in a given output stream for fermions with SU(2) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& FermionSU2WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for fermions with SU(3) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& FermionSU3WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for fermions with SU(4) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& FermionSU4WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereGetDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta", 20);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(3) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su2su2-spin", "consider particles with SU(2)xSU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "ground-only", "get the dimension only for the largest subspace");
  (*SystemGroup) += new BooleanOption ('\n', "use-files", "use dimension files that have been previously generated to increase speed. Files must be in current directory and obey the statistics_sphere_n_nbrparticles_q_nbrfluxquanta.dim naming convention");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_n_nbrparticles_q_nbrfluxquanta.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereGetDimension -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int  LzMin = 0;
  if (((NbrParticles * NbrFluxQuanta) & 1) != 0)
    LzMin = 1;
  if ((((BooleanOption*) Manager["su4-spin"])->GetBoolean() == false) && (((BooleanOption*) Manager["su2-spin"])->GetBoolean() == false) && 
      (((BooleanOption*) Manager["su3-spin"])->GetBoolean() == false) && (((BooleanOption*) Manager["su2su2-spin"])->GetBoolean() == false))
    if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
      {
	if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	  cout << BosonEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
	else
	  cout << FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
      }
    else
      {
	int LzMax = 0;
	if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	  LzMax = (NbrParticles * NbrFluxQuanta);
	else
	  LzMax = ((NbrFluxQuanta - NbrParticles + 1) * NbrParticles);
	long* LzDimensions = new long [1 + ((LzMax - LzMin) >> 1)];
	long* LDimensions = new long [1 + ((LzMax - LzMin) >> 1)];
	if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	  for (int x = LzMin; x <= LzMax; x += 2)
	    LzDimensions[(x - LzMin) >> 1] = BosonEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x);
	else
	  for (int x = LzMin; x <= LzMax; x += 2)
	    LzDimensions[(x - LzMin) >> 1] =  FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, x);
	LDimensions[(LzMax - LzMin) >> 1] = LzDimensions[(LzMax - LzMin) >> 1];
	long TotalDimension = LzDimensions[(LzMax - LzMin) >> 1];
	for (int x = LzMax - 2; x >= LzMin; x -= 2)
	  {
	    LDimensions[(x - LzMin) >> 1] =  LzDimensions[(x - LzMin) >> 1] - LzDimensions[((x - LzMin) >> 1) + 1];
	    TotalDimension += LzDimensions[(x - LzMin) >> 1];
	  }
	
	if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
	  {
	    char* OutputFileName = 0;
	    if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
	      {
		OutputFileName = new char[256];
		if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
		  sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		else
		  sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
	      }
	    else
	      {
		OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
		strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
	      }
	    WriteDimensionToDisk (OutputFileName, NbrParticles, NbrParticles, ((BooleanOption*) Manager["boson"])->GetBoolean(),
				  LzDimensions, LDimensions, LzMin, LzMax, TotalDimension);
	    delete[] OutputFileName;
	  }
	else
	  {
	    cout << "Lz =";
	    for (int x = LzMin; x <= LzMax; x += 2)
	      cout << " " << LzDimensions[(x - LzMin) >> 1];
	    cout << endl << "L =";
	    for (int x = LzMin; x <= LzMax; x += 2)
	      cout << " " << LDimensions[(x - LzMin) >> 1];	  
	    cout << endl;
	  }
	delete[] LzDimensions;
	delete[] LDimensions;
      }
  else
    if (((BooleanOption*) Manager["su3-spin"])->GetBoolean() == false)
      {
	int Sz = 0;
	if (NbrParticles & 1)
	  Sz = 1;
	if (((BooleanOption*) Manager["su4-spin"])->GetBoolean() == true)
	  {
	    if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	      {
		cout << "SU(4) mode for bosons not yet available" << endl;	
		return -1;
	      }
	    else
	      {
		if (NbrParticles > (((NbrFluxQuanta + 1) << 2)))
		  {
		    cout << "error : number of flux quanta is too low" << endl;
		    return -1;
		  }
		if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
		  cout << FermionSU4ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
									 (Sz + NbrParticles) >> 1, (Sz + NbrParticles) >> 1, (Sz + NbrParticles) >> 1) << endl;
		else
		  {
		    if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
		      {
			char* OutputFileName = 0;
			if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
			  {
			    OutputFileName = new char[256];
			    sprintf (OutputFileName, "fermions_sphere_su4_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
			  }
			else
			  {
			    OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
			    strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
			  }		  
			ofstream File;
			File.open(OutputFileName, ios::binary | ios::out);
			FermionSU4WriteDimension(File, NbrParticles, NbrFluxQuanta);
			File.close();
			delete[] OutputFileName;
		      }
		    else
		      FermionSU4WriteDimension (cout, NbrParticles, NbrFluxQuanta);
		  }
	      }
	  }
	else
	  {
	    if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	      {
		cout << "SU(2) mode not yet available" << endl;	
		return -1;
	      }
	    else
	      {
		if (NbrParticles > (((NbrFluxQuanta + 1) << 1)))
		  {
		    cout << "error : number of flux quanta is too low" << endl;
		    return -1;
		  }
		if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
		  cout << FermionSU2ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
									 (Sz + NbrParticles) >> 1) << endl;
		else
		  {
		    if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
		      {
			char* OutputFileName = 0;
			if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
			  {
			    OutputFileName = new char[256];
			    sprintf (OutputFileName, "fermions_sphere_su2_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
			  }
			else
			  {
			    OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
			    strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
			  }		  
			ofstream File;
			File.open(OutputFileName, ios::binary | ios::out);
			FermionSU2WriteDimension(File, NbrParticles, NbrFluxQuanta);
			File.close();
			delete[] OutputFileName;
		      }
		    else
		      FermionSU2WriteDimension (cout, NbrParticles, NbrFluxQuanta);
		  }	      
	      }
	  }
      }
    else
      {
	if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	  {
	    cout << "SU(3) mode not yet available" << endl;	
	    return -1;
	  }
	else
	  {
	    if (NbrParticles > (((NbrFluxQuanta + 1) * 3)))
	      {
		cout << "error : number of flux quanta is too low" << endl;
		return -1;
	      }
	    if (((BooleanOption*) Manager["ground-only"])->GetBoolean() == true)
	      {
		int MeanNbrParticles = NbrParticles / 3;
		cout << FermionSU3ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								       MeanNbrParticles, MeanNbrParticles, (NbrParticles - (2 * MeanNbrParticles))) << endl;
	      }
	    else
	      {
		if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
		  {
		    char* OutputFileName = 0;
		    if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
		      {
			OutputFileName = new char[256];
			sprintf (OutputFileName, "fermions_sphere_su3_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		      }
		    else
		      {
			OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
			strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
		      }		  
		    ofstream File;
		    File.open(OutputFileName, ios::binary | ios::out);
		    FermionSU3WriteDimension(File, NbrParticles, NbrFluxQuanta);
		    File.close();
		    delete[] OutputFileName;
		  }
		else
		  FermionSU3WriteDimension (cout, NbrParticles, NbrFluxQuanta);
	      }	      
	  }
      }
}

// evaluate Hilbert space dimension for bosons
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return BosonShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + lzMax * nbrBosons) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for bosons
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long BosonShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0l;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    {
      return 1l;
    }
  long TmpDim = 0;
  while ((totalLz >= 0) && (nbrBosons > 0))
    {
      TmpDim += BosonShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}

// evaluate Hilbert space dimension for fermions
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz for fermions
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return (long) 0;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return (long) 0;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return (long) 1;
  if (LzTotalMax == totalLz)
    return (long) 1;
  return  (FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   +  FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}

// evaluate Hilbert space dimension for fermions with SU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension

long FermionSU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  unsigned long Tmp = 0l;  
  if (nbrFermions > 2)    
    Tmp += FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1);
  else
    if ((totalLz == (2 * lzMax)) && (totalSpin == 1))
      ++Tmp;
  return  (Tmp + FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1)
	   + FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin)
	   + FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin));
}

// evaluate Hilbert space dimension for fermions with SU(2)xSU(2) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrFermions * lzMax)
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// return value = Hilbert space dimension

long FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) ||  
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  unsigned long Tmp = 0l;
  if (nbrFermions >= 3)    
    {
      Tmp += (FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, totalIsospin - 1)
	      + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 2)
	      + (2l * FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1))
	      + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin)
	      + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, totalIsospin - 1));

      if (nbrFermions > 3)
	{
	  Tmp += (FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2)
		  + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1)
		  + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2)
		  + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1));
	  if (nbrFermions == 4)
	    {
	      if ((totalLz == (4 * lzMax)) && (totalSpin == 2) && (totalIsospin == 2))
		++Tmp;      
	    }
	  else
	    Tmp += FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, totalIsospin - 2);
	}
      else
	if ((totalLz == (3 * lzMax)) && (((totalSpin == 2) || (totalSpin == 1)) && ((totalIsospin == 2) || (totalIsospin == 1))))
	  ++Tmp;
    }
  else
    if (totalLz == (2 * lzMax))
      {
 	switch (totalSpin)
 	  {
 	  case 2:
	    if (totalIsospin == 1)
	      ++Tmp;
 	    break;
 	  case 1:
	    switch (totalIsospin)
	      {
	      case 2:
		++Tmp;
		break;
	      case 1:
		Tmp += 2l;
		break;
	      case 0:
		++Tmp;
		break;
	      }
	    break;
 	  case 0:
	    if (totalIsospin == 1) 
	      ++Tmp;
	    break; 
 	  }
      }

  return  (Tmp + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin)
	   + FermionSU2SU2ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin));
}

// evaluate Hilbert space dimension for fermions with SU(4) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// totalEntanglement = number of particles with entanglement plus
// return value = Hilbert space dimension

long FermionSU4ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) || (totalEntanglement < 0) ||
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions) || (totalEntanglement > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < totalEntanglement) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalEntanglement)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) 
      || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz)
      || ((((2 * lzMax + nbrFermions + 1 - totalEntanglement) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if ((nbrFermions == 0) && (totalLz == 0) && (totalSpin == 0) && (totalIsospin == 0) && (totalEntanglement == 0))
    return 1l;
  if (nbrFermions == 1) 
    if ((lzMax >= totalLz) && (totalEntanglement != (totalSpin ^ totalIsospin)))
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  unsigned long Tmp = 0l;
  if (nbrFermions >= 3)    
    {
      Tmp += (FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 1)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 1)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 2)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin, totalEntanglement - 1)
	      + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, totalIsospin - 1, totalEntanglement - 1));
      
      if (nbrFermions > 3)
	{
 	  Tmp += (FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 1)
 		  + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 1)
 		  + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 2)
 		  + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 2));
 	  if (nbrFermions == 4)
 	    {
 	      if ((totalLz == (4 * lzMax)) && (totalSpin == 2) && (totalIsospin == 2) && (totalEntanglement == 2))
 		++Tmp;      
 	    }
 	  else
 	    Tmp += FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 2);
 	}
      else
	if ((totalLz == (3 * lzMax)) && (((totalSpin * totalIsospin * totalEntanglement) == 4) || ((totalSpin * totalIsospin * totalEntanglement) == 1)))
 	  ++Tmp;
    }
  else
    if (totalLz == (2 * lzMax))
      {
 	switch (totalSpin)
 	  {
 	  case 2:
	    if ((totalIsospin == 1) && (totalEntanglement == 1))
	      ++Tmp;
 	    break;
 	  case 1:
	    switch (totalIsospin)
	      {
	      case 2:
		if (totalEntanglement == 1)
		  ++Tmp;
		break;
	      case 1:
		if (totalEntanglement != 1)
		  ++Tmp;
		break;
	      case 0:
		if (totalEntanglement == 1)
		  ++Tmp;
		break;
	      }
	    break;
 	  case 0:
	    if ((totalIsospin == 1)  && (totalEntanglement == 1))
	      ++Tmp;
	    break; 
 	  }
      }

  return  (Tmp + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin, totalEntanglement)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1, totalEntanglement)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1, totalEntanglement - 1)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin, totalEntanglement - 1)
	   + FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin, totalEntanglement));

}

// evaluate Hilbert space dimension for fermions with SU(3) spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long FermionSU3ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrFermions < 0) || (totalLz < 0) || (lzMax < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0) || (lzMax < 0) || 
      ((nbrN1 - 1)> lzMax) || ((nbrN2 - 1)> lzMax) || ((nbrN3 - 1)> lzMax) ||
      ((nbrFermions * lzMax - (((nbrN1 * nbrN1) + (nbrN2 * nbrN2) + (nbrN3 * nbrN3) - nbrFermions) >> 1)) < totalLz))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;
  unsigned long Tmp = 0l;
  if (nbrFermions >= 3)
    {
      Tmp += (FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrN1 - 1, nbrN2 - 1, nbrN3)
	      + FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrN1 - 1, nbrN2, nbrN3 - 1)
	      + FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrN1, nbrN2 - 1, nbrN3 - 1));
      if (nbrFermions == 3)
	{
	  if ((totalLz == (3 * lzMax)) && (nbrN1 == 1) && (nbrN2 == 1) && (nbrN3 == 1))
	    ++Tmp;      
	}
      else
	Tmp += FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (lzMax * 3), nbrN1 - 1, nbrN2 - 1, nbrN3 -1);
    }
  else
    {
      if ((totalLz == (2 * lzMax)) && (nbrN1 <= 1) && (nbrN2 <= 1) && (nbrN3 <= 1))
	++Tmp;
    }
  return  (Tmp + FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1 - 1, nbrN2, nbrN3)
	   + FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1, nbrN2 - 1, nbrN3)
	   + FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1, nbrN2, nbrN3 - 1)
	   + FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, nbrN1, nbrN2, nbrN3));
}

// evaluate Hilbert space dimension using previously generated Hilbert space dimension files (or compute them if they don't exist)
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

// long EvaluateHilbertSpaceDimensionWithDiskStorage(int nbrParticles, int nbrFluxQuanta, bool statistics,
// 						  long* lzDimensions, long* lDimensions, int lzMin, int lzMax)
// {
//   OutputFileName = new char[256];
//   if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
//     sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
//   else
//     sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
  
//   return FermionShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
// }

// save dimensions in a given file
//
// outputFileName = output file name
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true for bosons, false for fermions
// lzDimensions = array taht contains dimension of each Lz sector
// lDimensions = array that contains dimension of each L sector (without taking into account the Lz degeneracy)
// lzMin = twice the minimum Lz value
// lzMax = twice the maximum Lz value
// totalDimension = total Hilbert space dimension

bool WriteDimensionToDisk(char* outputFileName, int nbrParticles, int nbrFluxQuanta, bool statistics,
			  long* lzDimensions, long* lDimensions, int lzMin, int lzMax, long totalDimension)
{
  ofstream File;
  File.open(outputFileName, ios::binary | ios::out);
  File << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " ";
  if (statistics == true)
    File << "bosons";
  else
    File << "femions";
  File << " on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  File << "# total Hilbert space dimension = " << totalDimension << endl << endl 
       << "N = " << nbrParticles << endl
       << "2S = " << nbrFluxQuanta << endl << endl
       << "#  dimensions for the Lz subspaces (starting from 2Lz = " << lzMin << " " << (lzMin + 2) << " ..." << endl
       << "Lz =";
  for (int x = lzMin; x <= lzMax; x += 2)
    File << " " << lzDimensions[(x - lzMin) >> 1];
  File << endl
       << "#  dimensions for the L subspaces (starting from 2L = " << lzMin << " " << (lzMin + 2) << " ..." << endl
       << "L =";
  for (int x = lzMin; x <= lzMax; x += 2)
    File << " " << lDimensions[(x - lzMin) >> 1];
  File << endl;
  File.close();
  return true;
}

// save dimensions in a given output stream for fermions with SU(2) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& FermionSU2WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " fermions" << endl;
  output << "# with SU(2) spin on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Sz value) (twice the total Lz/L value) (dimension of the subspace with fixed Lz, Sz) (dimension of the subspace with fixed L, Lz=L, Sz)" << endl << endl;
  for (int Sz = nbrParticles & 1; Sz <= nbrParticles; Sz += 2)
    {
      int NUp = nbrParticles + Sz;
      int NDown = nbrParticles - Sz;
      if ((NUp >= 0) && (NDown >=0) && ((NUp & 0x1) == 0) && ((NDown & 0x1) == 0))
	{
	  NUp >>= 1;
	  NDown >>= 1;
	  int Min = (nbrParticles * nbrFluxQuanta) & 1;
	  int Max  = ((((nbrFluxQuanta - NUp + 1) * NUp) + ((nbrFluxQuanta - NDown + 1) * NDown)));
	  if ((Max >=  Min) && (NUp <= (nbrFluxQuanta + 1)) && (NDown <= (nbrFluxQuanta + 1)))
	    {
	      long* LzDimension = new long [((Max - Min) >> 1) + 1];
	      for (int Lz = Min; Lz <= Max; Lz += 2)
		LzDimension[(Lz - Min) >> 1] = FermionSU2ShiftedEvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, 
											      (Lz + (nbrParticles * nbrFluxQuanta)) >> 1, NUp);
	      for (int Lz = Min; Lz < Max; Lz += 2)
		output << Sz << " " << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " 
		       << (LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1]) << endl;
	      output << Sz << " " << Max << " " << LzDimension[(Max - Min) >> 1] << " " 
		     << LzDimension[(Max - Min) >> 1] << endl;
	      delete[] LzDimension;	      
	    }
	}
    }
  return output;
}

// save dimensions in a given output stream for fermions with SU(3) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& FermionSU3WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " fermions" << endl;
  output << "# with SU(3) spin on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Tz value) (three time the total Y value) (twice the total Lz/L value) (dimension of the subspace with fixed Lz, Tz, Y) (dimension of the subspace with fixed L, Lz=L, Tz, Y)" << endl << endl;
  for (int Tz = 0; Tz <= nbrParticles; ++Tz)
    for (int Y = - 2 * nbrParticles; Y <= nbrParticles; Y += 3)
      {
	int N1 = (2 * nbrParticles) + Y + (3 * Tz);
	int N2 = (2 * nbrParticles) + Y - (3 * Tz);
	int N3 = nbrParticles - Y;
	if ((N1 >= 0) && (N2 >= 0) && (N3 >= 0) && ((N1 % 6) == 0) && ((N2 % 6) == 0) && ((N3 % 3) == 0))
	  {
	    N1 /= 6;
	    N2 /= 6;
	    N3 /= 3;
	    int Min = (nbrParticles * nbrFluxQuanta) & 1;
	    int Max  = ((nbrFluxQuanta - N1 + 1) * N1) + ((nbrFluxQuanta - N2 + 1) * N2) + ((nbrFluxQuanta - N3 + 1) * N3);
	    if ((Max >=  Min) && (N1 <= (nbrFluxQuanta + 1)) && (N2 <= (nbrFluxQuanta + 1)) && (N3 <= (nbrFluxQuanta + 1)))
	      {
		long* LzDimension = new long [((Max - Min) >> 1) + 1];
		for (int Lz = Min; Lz <= Max; Lz += 2)
		  LzDimension[(Lz - Min) >> 1] = FermionSU3ShiftedEvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, 
												(Lz + (nbrParticles * nbrFluxQuanta)) >> 1, N1, N2, N3);
		for (int Lz = Min; Lz < Max; Lz += 2)
		  output << Tz << " " << Y << " " << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " 
			 << (LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1]) << endl;
		output << Tz << " " << Y << " " << Max << " " << LzDimension[(Max - Min) >> 1] << " " 
		       << LzDimension[(Max - Min) >> 1] << endl;
		delete[] LzDimension;	      
	      }
	  }
      }
  return output;
}

// save dimensions in a given output stream for fermions with SU(4) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& FermionSU4WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " fermions" << endl;
  output << "# with SU(4) spin on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Sz value) (twice the total Iz value) (twice the total Pz value) (twice the total Lz/L value) (dimension of the subspace with fixed Lz, Sz, Iz and Pz) (dimension of the subspace with fixed L, Lz=L, Sz, Iz and Pz)" << endl << endl;
  for (int Sz = nbrParticles & 1; Sz <= nbrParticles; Sz += 2)
    for (int Iz = nbrParticles & 1; Iz <= Sz; Iz += 2)
      for (int Pz = nbrParticles & 1; Pz <= Iz; Pz += 2)
	{
	  int NUpPlus = nbrParticles + Sz + Iz + Pz;
	  int NUpMinus = nbrParticles + Sz - Iz - Pz;
	  int NDownPlus = nbrParticles - Sz + Iz - Pz;
	  int NDownMinus = nbrParticles - Sz - Iz + Pz;
	  if ((NUpPlus >= 0) && ((NUpPlus & 0x3) == 0) && (NUpMinus >= 0) && ((NUpMinus & 0x3) == 0) &&
	      (NDownPlus >= 0) && ((NDownPlus & 0x3) == 0) && (NDownMinus >= 0) && ((NDownMinus & 0x3) == 0))
	    {
	      NUpPlus >>= 2;
	      NUpMinus >>= 2;
	      NDownPlus >>= 2;
	      NDownMinus >>= 2;
	      int Min = (nbrParticles * nbrFluxQuanta) & 1;
	      int Max  = (((nbrFluxQuanta - NUpPlus + 1) * NUpPlus) + ((nbrFluxQuanta - NUpMinus + 1) * NUpMinus) + 
			  ((nbrFluxQuanta - NDownPlus+ 1) * NDownPlus) + ((nbrFluxQuanta - NDownMinus + 1) * NDownMinus));
	      if ((Max >=  Min) && (NUpPlus <= (nbrFluxQuanta + 1)) && (NUpMinus <= (nbrFluxQuanta + 1)) && 
		  (NDownPlus <= (nbrFluxQuanta + 1))   && (NDownMinus <= (nbrFluxQuanta + 1))) 
		{
		  long* LzDimension = new long [((Max - Min) >> 1) + 1];
		  for (int Lz = Min; Lz <= Max; Lz += 2)
		    LzDimension[(Lz - Min) >> 1] = FermionSU4ShiftedEvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, 
												  (Lz + (nbrParticles * nbrFluxQuanta)) >> 1,
												  NUpPlus + NUpMinus,
												  NUpPlus + NDownPlus, NUpPlus + NDownMinus);
		  for (int Lz = Min; Lz < Max; Lz += 2)
		    output << Sz << " " << Iz << " " << Pz << " " << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " 
			   << (LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1]) << endl;
		  output << Sz << " " << Iz << " " << Pz << " " << Max << " " << LzDimension[(Max - Min) >> 1] << " " 
			 << LzDimension[(Max - Min) >> 1] << endl;
		  delete[] LzDimension;	      
		}
	    }
	}
  return output;
}

