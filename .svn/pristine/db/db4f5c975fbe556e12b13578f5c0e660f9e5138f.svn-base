#include "HilbertSpace/BosonOnSphere.h"

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

// evaluate Hilbert space dimension for bosons with SU(2) spin
//
// nbrBosons = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value (with shift nbrBosons * lzMax)
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension
long BosonSU2ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin);

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

// evaluate Hilbert space dimension for bosons with SU(3) spin
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension
long BosonSU3ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3);

// evaluate Hilbert space dimension for fermions in two Landau levels
//
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz);

// evaluate Hilbert space dimension for fermions in two Landau levels
//
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz);

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

// save dimensions in a given output stream for bosons with SU(2) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& BosonSU2WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for fermions with SU(3) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& FermionSU3WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for bosons with SU(3) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& BosonSU3WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for fermions with SU(4) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& FermionSU4WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for fermions in 2 Landau levels
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& Fermion2LLWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// save dimensions in a given output stream for bosons in 2 Landau levels
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& Boson2LLWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// LzMaxUp = maximum lz value on SLL 
// LzMaxDown = maximum lz value on LLL
// return value = Hilbert space dimension
long Boson2LLShiftedEvaluateFullHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int LzMaxUp, int LzMaxDown);

// save dimensions in a given output stream for fermions in 2 Landau levels
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& Fermion3LLWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

//evaluate Hilbert space dimension for bosons on the 4D sphere
//
//nbrBosons = number of bosons
//nbrFluxQuanta = number of flux quanta
//shiftedTotalJz = momentum jz total value (shifted by nbrBosons*nbrFluxQuanta)
//totalShiftedKz = momentum kz total value (shifted by nbrBosons*nbrFluxQuanta)
//return value = Hilbert space dimension
long Boson4DSphereEvaluateHilbertSpaceDimension(int nbrBosons, int nbrFluxQuanta, int shiftedTotalJz, int shiftedTotalKz, int currentJ, int currentJz, int currentKz, int currentTotalJz, int currentTotalKz);

// save dimensions in a given output stream for bosons on the 4D sphere
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& Boson4DSphereWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);

// evaluate Hilbert space dimension for fermions on the CP2 geometry
//
// nbrFermions = number of fermions
// currentTz = current value of Tz for a single particle
// currentTzMax = current maiximum value of Tz that can be reached with the currentY Y value
// currentY = current value of Y for a single particle
// currentTotalTz = current total value of Tz
// currentTotalY = current total value of Y
// totalTz = required total value of Tz
// totalTz = required total value of Y
// minY = three times the minimum value for Y
// return value = Hilbert space dimension
long FermionCP2EvaluateHilbertSpaceDimension(int nbrFermions, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, int totalTz, int totalY, int minY);

// save dimensions in a given output stream for fermions on the CP geometry
// 
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream
ostream& FermionCP2SphereWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta);


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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-lz", "Lz sector", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-sz", "spin projection for SU(2) case", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "4-D", "consider particles on the 4D sphere (only available for bosons)");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(3) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su2su2-spin", "consider particles with SU(2)xSU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "2-ll", "consider particles within two Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "3-ll", "consider particles within three Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "cp2", "consider particles on the CP2 ");
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
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
  int LzMin = Manager.GetInteger("total-lz"); 
  int JzMin = 0;
  int KzMin = 0;
  
  if (((NbrParticles * NbrFluxQuanta) & 1) != 0) //if odd then lzmin is 1/2
    {
    LzMin = 1;
    JzMin = 1;
    }
  
  if ((Manager.GetBoolean("su4-spin") == false) && (Manager.GetBoolean("su2-spin") == false) && 
      (Manager.GetBoolean("su3-spin") == false) && (Manager.GetBoolean("su2su2-spin") == false) && 
      (Manager.GetBoolean("2-ll") == false) && (Manager.GetBoolean("3-ll") == false) && (Manager.GetBoolean("4-D") == false) &&
      (Manager.GetBoolean("cp2") == false))
    {
      if (Manager.GetBoolean("ground-only") == true)
	{
	  if (Manager.GetBoolean("boson") == true)
	    cout << BosonEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
	  else
	    cout << FermionEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, LzMin) << endl;
	  return 0;
	}
      else
	{
	  int LzMax = 0;
	  if (Manager.GetBoolean("boson") == true)
	    LzMax = (NbrParticles * NbrFluxQuanta);
	  else
	    LzMax = ((NbrFluxQuanta - NbrParticles + 1) * NbrParticles);
	  long* LzDimensions = new long [1 + ((LzMax - LzMin) >> 1)];
	  long* LDimensions = new long [1 + ((LzMax - LzMin) >> 1)];
	  if (Manager.GetBoolean("boson") == true)
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
	  
	  if (Manager.GetBoolean("save-disk") == true)
	    {
	      char* OutputFileName = 0;
	      if (Manager.GetString("output-file") == 0)
		{
		  OutputFileName = new char[256];
		  if (Manager.GetBoolean("boson") == true)
		    sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		  else
		    sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		}
	      else
		{
		  OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		  strcpy (OutputFileName, Manager.GetString("output-file"));
		}
	      WriteDimensionToDisk (OutputFileName, NbrParticles, NbrParticles, Manager.GetBoolean("boson"),
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
	  return 0;
	}
      return 0;
    }
  if (Manager.GetBoolean("su4-spin") == true)
    {
      int Sz = 0;
      if (NbrParticles & 1)
	Sz = 1;
      if (Manager.GetBoolean("boson") == true)
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
	  if (Manager.GetBoolean("ground-only") == true)
	    cout << FermionSU4ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								   (Sz + NbrParticles) >> 1, (Sz + NbrParticles) >> 1, (Sz + NbrParticles) >> 1) << endl;
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "fermions_sphere_su4_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
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
      return 0;
    }
  if (Manager.GetBoolean("su2-spin") == true)
    {
      int Sz = 0;
      if (NbrParticles & 1)
	Sz = 1;
      if (Manager.GetInteger("total-sz") > 0)
        Sz = Manager.GetInteger("total-sz");

      if (Manager.GetBoolean("boson") == true)
	{

	  if (Manager.GetBoolean("ground-only") == true)
	    cout << BosonSU2ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								 (Sz + NbrParticles) >> 1) << endl;
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "bosons_sphere_su2_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
		    }		  
		  ofstream File;
		  File.open(OutputFileName, ios::binary | ios::out);
		  BosonSU2WriteDimension(File, NbrParticles, NbrFluxQuanta);
		  File.close();
		  delete[] OutputFileName;
		}
	      else
		BosonSU2WriteDimension (cout, NbrParticles, NbrFluxQuanta);
	    }	      
	}
      else
	{
	  if (NbrParticles > (((NbrFluxQuanta + 1) << 1)))
	    {
	      cout << "error : number of flux quanta is too low" << endl;
	      return -1;
	    }
	  if (Manager.GetBoolean("ground-only") == true)
	    cout << FermionSU2ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								   (Sz + NbrParticles) >> 1) << endl;
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "fermions_sphere_su2_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
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
      return 0;
    }
  if (Manager.GetBoolean("su3-spin") == true)
    {
      if (Manager.GetBoolean("boson") == true)
	{
	  if (Manager.GetBoolean("ground-only") == true)
	    {
	      int MeanNbrParticles = NbrParticles / 3;
	      cout << BosonSU3ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								   MeanNbrParticles, MeanNbrParticles, (NbrParticles - (2 * MeanNbrParticles))) << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "bosons_sphere_su3_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
		    }		  
		  ofstream File;
		  File.open(OutputFileName, ios::binary | ios::out);
		  BosonSU3WriteDimension(File, NbrParticles, NbrFluxQuanta);
		  File.close();
		  delete[] OutputFileName;
		}
	      else
		BosonSU3WriteDimension (cout, NbrParticles, NbrFluxQuanta);
	    }	      
	}
      else
	{
	  if (NbrParticles > (((NbrFluxQuanta + 1) * 3)))
	    {
	      cout << "error : number of flux quanta is too low" << endl;
	      return -1;
	    }
	  if (Manager.GetBoolean("ground-only") == true)
	    {
	      int MeanNbrParticles = NbrParticles / 3;
	      cout << FermionSU3ShiftedEvaluateHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, (LzMin + (NbrFluxQuanta * NbrParticles)) >> 1, 
								     MeanNbrParticles, MeanNbrParticles, (NbrParticles - (2 * MeanNbrParticles))) << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "fermions_sphere_su3_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
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
      return 0;
     }
  if (Manager.GetBoolean("2-ll") == true)
    {
      if (Manager.GetBoolean("boson") == true)
	{
	  if (Manager.GetBoolean("save-disk") == true)
	    {
	      char* OutputFileName = 0;
	      if (Manager.GetString("output-file") == 0)
		{
		  OutputFileName = new char[256];
		  sprintf (OutputFileName, "bosons_sphere_2ll_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		}
	      else
		{
		  OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		  strcpy (OutputFileName, Manager.GetString("output-file"));
		}		  
	      ofstream File;
	      File.open(OutputFileName, ios::binary | ios::out);
	      Boson2LLWriteDimension(File, NbrParticles, NbrFluxQuanta);
	      File.close();
	      delete[] OutputFileName;
	    }
	  else
	      Boson2LLWriteDimension(cout, NbrParticles, NbrFluxQuanta);
	}
      else
	{
	  if (NbrParticles > ((2 * NbrFluxQuanta) + 4))
	    {
	      cout << "error : number of flux quanta is too low" << endl;
	      return -1;
	    }
	  if (Manager.GetBoolean("ground-only") == true)
	    {
	      cout << Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, NbrFluxQuanta + 2, (LzMin + ((NbrFluxQuanta + 2) * NbrParticles)) >> 1) << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "fermions_sphere_2ll_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
		    }		  
		  ofstream File;
		  File.open(OutputFileName, ios::binary | ios::out);
		  Fermion2LLWriteDimension(File, NbrParticles, NbrFluxQuanta);
		  File.close();
		  delete[] OutputFileName;
		}
	      else
		Fermion2LLWriteDimension (cout, NbrParticles, NbrFluxQuanta);
	    }
	}
      return 0;
    }
  if (Manager.GetBoolean("3-ll") == true)
    {
      if (Manager.GetBoolean("boson") == true)
	{
	  cout << "3 Landau level mode not yet available" << endl;	
	  return -1;
	}
      else
	{
	  if (NbrParticles > ((3 * NbrFluxQuanta) + 9))
	    {
	      cout << "error : number of flux quanta is too low" << endl;
	      return -1;
	    }
	  if (Manager.GetBoolean("ground-only") == true)
	    {
	      cout << Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(NbrParticles, NbrFluxQuanta, NbrFluxQuanta + 4, (LzMin + ((NbrFluxQuanta + 4) * NbrParticles)) >> 1) << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  char* OutputFileName = 0;
		  if (Manager.GetString("output-file") == 0)
		    {
		      OutputFileName = new char[256];
		      sprintf (OutputFileName, "fermions_sphere_3ll_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		    }
		  else
		    {
		      OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
		      strcpy (OutputFileName, Manager.GetString("output-file"));
		    }		  
		  ofstream File;
		  File.open(OutputFileName, ios::binary | ios::out);
		  Fermion3LLWriteDimension(File, NbrParticles, NbrFluxQuanta);
		  File.close();
		  delete[] OutputFileName;
		}
	      else
		Fermion3LLWriteDimension (cout, NbrParticles, NbrFluxQuanta);
	    }
	}
      return 0;
    }
  if (Manager.GetBoolean("4-D") == true)
    {
      if (Manager.GetBoolean("boson") == false)
	{
	  cout << "Fermionic mode not implemented for FQHE on the 4D sphere" << endl;
	}
      else
	{
	  if (Manager.GetBoolean("save-disk") == true)
	    {
	      char* OutputFileName = 0;
	      if (Manager.GetString("output-file") == 0)
		{
		  OutputFileName = new char[256];
		  sprintf (OutputFileName, "bosons_sphere4d_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		}
	      ofstream File;
	      File.open(OutputFileName, ios::binary | ios::out);
	      Boson4DSphereWriteDimension(File, NbrParticles, NbrFluxQuanta);
	      File.close();
	      delete[] OutputFileName;
	    }
	  else
	    Boson4DSphereWriteDimension(cout, NbrParticles, NbrFluxQuanta);
	  
	}
      return 0;
    }
  if (Manager.GetBoolean("cp2") == true)
    {
      if (Manager.GetBoolean("boson") == false)
	{
	  if (Manager.GetBoolean("save-disk") == true)
	    {
	      char* OutputFileName = 0;
	      if (Manager.GetString("output-file") == 0)
		{
		  OutputFileName = new char[256];
		  sprintf (OutputFileName, "fermions_cp2_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
		}
	      ofstream File;
	      File.open(OutputFileName, ios::binary | ios::out);
	      FermionCP2SphereWriteDimension(File, NbrParticles, NbrFluxQuanta);
	      File.close();
	      delete[] OutputFileName;
	    }
	  else
	    {
	      FermionCP2SphereWriteDimension(cout, NbrParticles, NbrFluxQuanta);
	    }
	}
      else
	{
// 	  if (Manager.GetBoolean("save-disk") == true)
// 	    {
// 	      char* OutputFileName = 0;
// 	      if (Manager.GetString("output-file") == 0)
// 		{
// 		  OutputFileName = new char[256];
// 		  sprintf (OutputFileName, "bosons_sphere4d_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
// 		}
// 	      ofstream File;
// 	      File.open(OutputFileName, ios::binary | ios::out);
// 	      Boson4DSphereWriteDimension(File, NbrParticles, NbrFluxQuanta);
// 	      File.close();
// 	      delete[] OutputFileName;
// 	    }
// 	  else
// 	    Boson4DSphereWriteDimension(cout, NbrParticles, NbrFluxQuanta);	  
	}
      return 0;
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
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }

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
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
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
    {
      if ((lzMax >= totalLz) && (totalEntanglement != (totalSpin ^ totalIsospin)))
	return 1l;
      else
	return 0l;
    }

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
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
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

// evaluate Hilbert space dimension for bosons with SU(3) spin
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long BosonSU3ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrN1; i >= 0; --i)
    for (int j = nbrN2; j >= 0; --j)
      for (int k = nbrN3; k >= 0; --k)
	Tmp += BosonSU3ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j + k), lzMax - 1, totalLz - (lzMax * (i + j + k)), 
							    nbrN1 - i, nbrN2 - j, nbrN3 - k);
  return  Tmp;
}

// evaluate Hilbert space dimension for fermions in two Landau levels
//
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((nbrFluxQuanta + 2) >= totalLz) && (totalLz >= 0))
	    ++Tmp;
	  if (((nbrFluxQuanta + 1) >= totalLz) && (totalLz >= 1))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 1)) && (lzMax >= 1))
	Tmp += Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 1)) && (lzMax >= 1))    
    Tmp += Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
  Tmp += Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions, nbrFluxQuanta, lzMax - 1, totalLz);
  return Tmp;
}

// evaluate Hilbert space dimension for fermions in three Landau levels
//
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 0))
	    ++Tmp;
	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 1))
	    ++Tmp;
	  if (((nbrFluxQuanta + 2) >= totalLz) && (totalLz >= 2))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
	{
	  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	    Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	  Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
    {
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))    
    Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
  Tmp += Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrFermions, nbrFluxQuanta, lzMax - 1, totalLz);
  return Tmp;
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
//   if (Manager.GetBoolean("boson") == true)
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

// save dimensions in a given output stream for bosons with SU(2) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& BosonSU2WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " bosons" << endl;
  output << "# with SU(2) spin on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Sz value) (twice the total Lz/L value) (dimension of the subspace with fixed Lz, Sz) (dimension of the subspace with fixed L, Lz=L, Sz)" << endl << endl;
  for (int Sz = nbrParticles & 1; Sz <= nbrParticles; Sz += 2)
    {
      int NUp = nbrParticles + Sz;
      int NDown = nbrParticles - Sz;
      if ((NUp >= 0) && (NDown >=0))// && ((NUp & 0x1) == 0) && ((NDown & 0x1) == 0))
	{
	  NUp >>= 1;
	  NDown >>= 1;
	  int Min = (nbrParticles * nbrFluxQuanta) & 1;
	  int Max  = ((((nbrFluxQuanta + 1) * NUp) + ((nbrFluxQuanta + 1) * NDown)));
	  if (Max >=  Min)
	    {
	      long* LzDimension = new long [((Max - Min) >> 1) + 1];
	      for (int Lz = Min; Lz <= Max; Lz += 2)
		LzDimension[(Lz - Min) >> 1] = BosonSU2ShiftedEvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, 
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

// save dimensions in a given output stream for bosons with SU(3) spin
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& BosonSU3WriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " bosons" << endl;
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
	    int Max  = nbrFluxQuanta * nbrParticles;
	    long* LzDimension = new long [((Max - Min) >> 1) + 1];
	    for (int Lz = Min; Lz <= Max; Lz += 2)
	      LzDimension[(Lz - Min) >> 1] = BosonSU3ShiftedEvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, 
											  (Lz + (nbrParticles * nbrFluxQuanta)) >> 1, N1, N2, N3);
	    for (int Lz = Min; Lz < Max; Lz += 2)
	      output << Tz << " " << Y << " " << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " 
		     << (LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1]) << endl;
	    output << Tz << " " << Y << " " << Max << " " << LzDimension[(Max - Min) >> 1] << " " 
		   << LzDimension[(Max - Min) >> 1] << endl;
	    delete[] LzDimension;	      
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

// save dimensions in a given output stream for fermions in 2 Landau levels
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& Fermion2LLWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " fermions" << endl;
  output << "# in two Landau levels on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Lz/L value) (dimension of the subspace with fixed Lz) (dimension of the subspace with fixed L, Lz=L)" << endl << endl;
  int Min = (nbrParticles * nbrFluxQuanta) & 1;
  int TmpNbrParticleDown = (nbrParticles - 1) / 2;
  int Max = (((nbrFluxQuanta - TmpNbrParticleDown + 1) * TmpNbrParticleDown) +
	     ((nbrFluxQuanta + 2 - (nbrParticles - TmpNbrParticleDown) + 1) * (nbrParticles - TmpNbrParticleDown)));
  long* LzDimension = new long [((Max - Min) >> 1) + 1];
  for (int Lz = Min; Lz <= Max; Lz += 2)
    LzDimension[(Lz - Min) >> 1] = Fermion2LLShiftedEvaluateFullHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, nbrFluxQuanta + 2, (Lz + (nbrParticles * (nbrFluxQuanta + 2))) >> 1);
  for (int Lz = Min; Lz < Max; Lz += 2)
    output << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " 
	   << (LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1]) << endl;
    output << Max << " " << LzDimension[(Max - Min) >> 1] << " " 
	   << LzDimension[(Max - Min) >> 1] << endl;
  delete[] LzDimension;
  return output;
}

// save dimensions in a given output stream for fermions in 3 Landau levels
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& Fermion3LLWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " fermions" << endl;
  output << "# in three Landau levels on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Lz/L value) (dimension of the subspace with fixed Lz) (dimension of the subspace with fixed L, Lz=L)" << endl << endl;
  int Min = (nbrParticles * nbrFluxQuanta) & 1;
  int TmpNbrParticle1 = (nbrParticles + 3) / 3;
  int TmpNbrParticle2 = (nbrParticles - TmpNbrParticle1 + 1) / 2;
  int TmpNbrParticle3 = nbrParticles - TmpNbrParticle1 - TmpNbrParticle2;
  while (TmpNbrParticle3 > (nbrFluxQuanta + 1))
    {
      --TmpNbrParticle3;
      ++TmpNbrParticle2;
    }
  while (TmpNbrParticle2 > (nbrFluxQuanta + 3))
    {
      --TmpNbrParticle2;
      ++TmpNbrParticle1;
    }
  int Max = (((nbrFluxQuanta - TmpNbrParticle3 + 1) * TmpNbrParticle3)
	     + ((nbrFluxQuanta + 2 - TmpNbrParticle2 + 1) * TmpNbrParticle2)
	     + ((nbrFluxQuanta + 4 - TmpNbrParticle1 + 1) * TmpNbrParticle1));
  long* LzDimension = new long [((Max - Min) >> 1) + 1];
  for (int Lz = Min; Lz <= Max; Lz += 2)
    LzDimension[(Lz - Min) >> 1] = Fermion3LLShiftedEvaluateFullHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, nbrFluxQuanta + 4, (Lz + (nbrParticles * (nbrFluxQuanta + 4))) >> 1);
  for (int Lz = Min; Lz < Max; Lz += 2)
    output << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " 
	   << (LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1]) << endl;
    output << Max << " " << LzDimension[(Max - Min) >> 1] << " " 
	   << LzDimension[(Max - Min) >> 1] << endl;
  delete[] LzDimension;
  return output;
}

// save dimensions in a given output stream for bosons in 2 Landau levels
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& Boson2LLWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each L and Lz sector for " << nbrParticles << " bosons" << endl;
  output << "# in two Landau levels on the sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Lz/L value) (dimension of the subspace with fixed Lz) (dimension of the subspace with fixed L, Lz=L)" << endl << endl;
  int Min = (nbrParticles * nbrFluxQuanta) & 1;
  //int Min = -nbrParticles * (nbrFluxQuanta+2);
  int TmpNbrParticleDown = (nbrParticles - 1) / 2;
  int Max = (nbrFluxQuanta + 2) * nbrParticles; 
  long* LzDimension = new long [((Max - Min) >> 1) + 1];
  for (int Lz = Min; Lz <= Max; Lz += 2)
    {
      /*BosonOnSphereTwoLandauLevels* space = new BosonOnSphereTwoLandauLevels(nbrParticles, Lz, nbrFluxQuanta+2, nbrFluxQuanta);
      LzDimension[(Lz - Min) >> 1] = space->GetTargetHilbertSpaceDimension();
      delete space; */
      LzDimension[(Lz - Min) >> 1] = Boson2LLShiftedEvaluateFullHilbertSpaceDimension(nbrParticles, 0, (Lz + (nbrParticles * (nbrFluxQuanta+2))) >> 1 , nbrFluxQuanta+2, nbrFluxQuanta);
    }
  for (int Lz = Min; Lz < Max; Lz += 2)
    {
      output << Lz << " " << LzDimension[(Lz - Min) >> 1] << " " << LzDimension[(Lz - Min) >> 1] - LzDimension[((Lz - Min) >> 1) + 1] << endl;      
    }
    output << Max << " " << LzDimension[(Max - Min) >> 1] << " " << LzDimension[(Max - Min) >> 1] << endl;
  delete[] LzDimension;
  return output;
}

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// LzMaxUp = maximum lz value on SLL 
// LzMaxDown = maximum lz value on LLL
// return value = Hilbert space dimension

long Boson2LLShiftedEvaluateFullHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int LzMaxUp, int LzMaxDown)
{
  if ((nbrBosons == 0) && (totalLz == 0)) //all bosons gone and correct total totalLz been realised. Fill in correct entries on way back up.
    {
      return 1l;
    }
    
  if (nbrBosons < 0 || totalLz < 0  )//|| (this->MaxLzLeft(nbrBosons,pos) < totalLz) ) //if state not working out. 
    return 0l;
    
  if (lzMax < 0) //if the position is negative. This should never happen.
    return 0l;
  
  if ((lzMax == (LzMaxUp + LzMaxDown + 2)) && (totalLz != 0)) //if at the final position and still hav not satisfied the totallz.
    return 0l;
  
  int currentLz = 0; // this is the lz value of the current position.
  if ( lzMax <=  LzMaxUp )  //if its on the Up (second) Landau level.
    {
      currentLz = LzMaxUp - lzMax; 
    }
  else //if its on the down (lowest) Landau level
    {
      currentLz = LzMaxDown + 1 - (lzMax - (LzMaxUp+1));
    }
    
  long Tmp = 0;
  if ( lzMax < (LzMaxUp + LzMaxDown + 2) ) //if the position has not reached the end.
    {
      //Can place between 0 and nbrBosons on this spot and then move on.
      for ( int i = nbrBosons ; i >= 1 ; i-- ) 
        {
	  //place i bosons
	  Tmp += Boson2LLShiftedEvaluateFullHilbertSpaceDimension(nbrBosons - i, lzMax+1, totalLz - (currentLz*i), LzMaxUp, LzMaxDown);
	}
	Tmp += Boson2LLShiftedEvaluateFullHilbertSpaceDimension(nbrBosons, lzMax+1, totalLz, LzMaxUp, LzMaxDown);
    }
    return Tmp;
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// totalSpin = twice the total spin value
// return value = Hilbert space dimension

long BosonSU2EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin)
{
  return BosonSU2ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, 
						      (totalLz + lzMax * nbrBosons) >> 1, (totalSpin + nbrBosons) >> 1);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension      

long BosonSU2ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrBosons < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrBosons) || ((lzMax * nbrBosons) < totalLz))
    return 0l;
    
  if (nbrBosons == 1) 
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  if (totalLz == 0)
    return 1l;

  unsigned long Tmp = 0l;  
  for (int i = totalSpin; i >= 0; --i)
    for (int j = (nbrBosons - totalSpin); j >= 0; --j)
      Tmp += BosonSU2ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j), lzMax - 1, totalLz - (lzMax * (i + j)), totalSpin - i);
  return Tmp;  
}


//evaluate Hilbert space dimension for bosons on the 4D sphere
//
//nbrBosons = number of bosons
//return value = Hilbert space dimension

long Boson4DSphereEvaluateHilbertSpaceDimension(int nbrBosons, int nbrFluxQuanta, int shiftedTotalJz, int shiftedTotalKz, int currentJ, int currentJz, int currentKz, int currentTotalJz, int currentTotalKz)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentTotalKz > shiftedTotalKz)
    return 0l;
  if (currentTotalJz > shiftedTotalJz)
    return 0l;
  if (currentTotalJz + nbrBosons*(currentJ + nbrFluxQuanta) < shiftedTotalJz)
    return 0l;
  if (currentTotalKz + nbrBosons*(2*nbrFluxQuanta + currentJ) < shiftedTotalKz)
    return 0l;
  
  
  if (currentKz < 0)
   {
     currentKz = nbrFluxQuanta - currentJ;
     currentJz--;
   }
    
  if (currentJz < 0)
  {
   currentJ--;
   currentJz = currentJ;
   currentKz = nbrFluxQuanta - currentJ;
  }
  
  if (nbrBosons == 0)
    {
      if ((currentTotalJz == shiftedTotalJz) && (currentTotalKz == shiftedTotalKz))
      {
	return 1l;
      }
      else	
	return 0l;
    }
    
  if (currentJ < 0)
    return 0l;
  
  long Count = 0;
  for (int k = nbrBosons; k >= 0; --k)
    Count += Boson4DSphereEvaluateHilbertSpaceDimension(nbrBosons - k, nbrFluxQuanta, shiftedTotalJz, shiftedTotalKz, currentJ, currentJz, currentKz - 1, currentTotalJz + (k * (2*currentJz - currentJ + nbrFluxQuanta)), currentTotalKz + (k * (2*currentKz + currentJ)));
  return Count;
}

// save dimensions in a given output stream for bosons on the 4D sphere
//
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& Boson4DSphereWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each Jz and Kz sector for " << nbrParticles << " bosons" << endl;
  output << "# on the 4D sphere geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (twice the total Jz value) (twice the total Kz value) (dimension of the subspace with fixed Jz, Kz)" << endl << endl;
	 
  for (int jz = 0; jz <= nbrParticles*nbrFluxQuanta; jz ++)
    {
      for (int kz = 0; kz <= nbrParticles*nbrFluxQuanta - jz; kz ++)
	 {
	  if ((kz <= jz) and (((nbrParticles*nbrFluxQuanta) & 1) == ((kz + jz) & 1)) )
	    {
	      int shiftedTotalJz = jz + nbrParticles*nbrFluxQuanta;
	      int shiftedTotalKz = kz + nbrParticles*nbrFluxQuanta;
	      output << jz << " " << kz << " " << Boson4DSphereEvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, shiftedTotalJz, shiftedTotalKz, nbrFluxQuanta, nbrFluxQuanta, 0, 0, 0) << endl;
	    }
	 }
    }
  return output;
}	 

// evaluate Hilbert space dimension for fermions on the CP2 geometry
//
// nbrFermions = number of fermions
// currentTz = current value of Tz for a single particle
// currentTzMax = current maiximum value of Tz that can be reached with the currentY Y value
// currentY = current value of Y for a single particle
// currentTotalTz = current total value of Tz
// currentTotalY = current total value of Y
// totalTz = required total value of Tz
// totalTz = required total value of Y
// minY = three times the minimum value for Y
// return value = Hilbert space dimension

long FermionCP2EvaluateHilbertSpaceDimension(int nbrFermions, int currentTz, int currentTzMax, int currentY, int currentTotalTz, int currentTotalY, int totalTz, int totalY, int minY)
{
  if (nbrFermions < 0)
    return 0l;
  if ((currentTotalTz + (currentTzMax * nbrFermions)) < totalTz)
    return 0l;
  if ((currentTotalY + (currentY * nbrFermions)) < totalY)
    return 0l;
  
  if (currentTz < -currentTzMax)
   {
     --currentTzMax;
     currentTz = currentTzMax;
     currentY = currentY - 3;
   }
    
  if (nbrFermions == 0)
    {
      if ((currentTotalTz == totalTz) && (currentTotalY == totalY))
      {
	return 1l;
      }
      else	
	return 0l;
    }

  if (currentY < minY)
    return 0l;
      
  long Count = 0;  
  Count += FermionCP2EvaluateHilbertSpaceDimension(nbrFermions - 1, currentTz - 2, currentTzMax, currentY, currentTotalTz + currentTz, currentTotalY + currentY, totalTz, totalY, minY);
  Count += FermionCP2EvaluateHilbertSpaceDimension(nbrFermions, currentTz - 2, currentTzMax, currentY, currentTotalTz, currentTotalY, totalTz, totalY, minY);
  return Count;
}

// save dimensions in a given output stream for fermions on the CP geometry
// 
// output = reference on the output stream
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = reference on the output stream

ostream& FermionCP2SphereWriteDimension(ostream& output, int nbrParticles, int nbrFluxQuanta)
{
  output << "# Hilbert space dimension in each Jz and Kz sector for " << nbrParticles << " fermions" << endl;
  output << "# on the C2 geometry with " << nbrFluxQuanta << " flux quanta" << endl;
  output << "#" << endl << "#  dimensions for each subspaces with the following convention " << endl 
	 << "# (three time the total Y value) (twice the total Tz value) (dimension of the subspace with fixed Y and Tz)" << endl;
	 
  int MinR = 0;
  int MaxR = nbrFluxQuanta * nbrParticles;
  for (int r = MinR; r <= MaxR; ++r)
    {
      int MinS = 0;
      int MaxS = nbrFluxQuanta * nbrParticles - r;
      if (MaxS > r)
	MaxS = r;
      for (int s = MinS; s <= MaxS ; ++s)
	{
	  int Tz = r - s;
	  int Y = 3 * (r + s) - (2 * nbrParticles * nbrFluxQuanta);
	  long TmpDimension = FermionCP2EvaluateHilbertSpaceDimension(nbrParticles, nbrFluxQuanta, nbrFluxQuanta, nbrFluxQuanta, 
								      0, 0, Tz, Y, -2 * nbrFluxQuanta);
	  if (TmpDimension > 0l)
	    {
	      output << Y << " " << Tz << " " << TmpDimension << endl;
	    }
	}
    }
  return output;
}	 

