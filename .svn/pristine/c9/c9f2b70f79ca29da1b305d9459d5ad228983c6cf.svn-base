#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"

#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunctionOneOverR.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEFermionsOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "momentum", "maximum single particle momentum to study", 10, true, 1);
  (*SystemGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction", "laughlin");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for energy calculation", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskFermionMonteCarloEnergy -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((BooleanOption*) Manager["list-wavefunctions"])->GetBoolean() == true)
    {
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MMax = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();


  Abstract1DComplexFunction* WaveFunction = new LaughlinOnDiskWaveFunctionOneOverR(NbrFermions, 3);
  RealVector Location(2 * NbrFermions, true);

  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);

  double Energy = 0.0;
  double EnergyError = 0.0;
  double Normalization = 0.0;
  double NormalizationError = 0.0;
  double Tmp = 0.0;

  Complex Tmp3;
  double Tmp2;
  double Tmp2bis;
  int NextCoordinates = 0;
  double Radius;
  double Theta;
  for (int j = 0; j < NbrFermions; ++j)
    {
      Radius = RandomNumber->GetRealRandomNumber();      
      Radius = sqrt (- 2.0 * log(1.0 - Radius));
      Theta = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      Location[j << 1] = Radius * cos (Theta);
      Location[1 + (j << 1)] = Radius * sin (Theta);
    }
  double PreviousProbabilities = SqrNorm((*WaveFunction)(Location));
  double CurrentProbabilities = PreviousProbabilities;
  double PreviousCoordinates1;
  double PreviousCoordinates2;
  for (int i = 1; i < NbrIter; ++i)
    {
      PreviousCoordinates1 = Location[NextCoordinates << 1];
      PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Radius = RandomNumber->GetRealRandomNumber();      
      Radius = sqrt (- 2.0 * log(1.0 - Radius));
      Theta = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      Location[NextCoordinates << 1] = Radius * cos (Theta);
      Location[1 + (NextCoordinates << 1)] = Radius * sin (Theta);
      CurrentProbabilities = SqrNorm((*WaveFunction)(Location));
      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	PreviousProbabilities = CurrentProbabilities;
      else
 	{
 	  Location[NextCoordinates << 1] = PreviousCoordinates1;
 	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
 	  CurrentProbabilities = PreviousProbabilities;
 	}
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrFermions)
	--NextCoordinates;
      
      Tmp = 0.0;
      for (int j = 0; j < NbrFermions; ++j)
	for (int k = j + 1; k < NbrFermions; ++k)
	  Tmp += 1.0 / sqrt (((Location[k << 1] - Location[j << 1]) * (Location[k << 1] - Location[j << 1])) + 
			      ((Location[1 + (k << 1)] - Location[1 + (j << 1)]) * (Location[1 + (k << 1)] - Location[1 + (j << 1)])));
      Energy += Tmp;
      EnergyError += Tmp * Tmp;
      if ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0)
 	{
	  cout << " i = " << i << endl;
	  double TmpEnergy = Energy / ((double) i);
	  double TmpEnergyError = sqrt (((EnergyError  / ((double) i)) - (TmpEnergy * TmpEnergy)) / ((double) i));
 	  cout << TmpEnergy << " +/- " << TmpEnergyError << endl;
 	  cout << "-----------------------------------------------" << endl;
 	}
    } 
  cout << " final results :" << endl;
  double TmpEnergy = Energy / ((double) NbrIter);
  double TmpEnergyError = sqrt (((EnergyError  / ((double) NbrIter)) - (TmpEnergy * TmpEnergy)) / ((double) NbrIter));
  cout << TmpEnergy << " +/- " << TmpEnergyError << endl;
  cout << ((TmpEnergy - ((double) (NbrFermions * NbrFermions))) / ((double) NbrFermions)) << " +/- " << (TmpEnergyError / ((double) NbrFermions)) << endl;
  cout << "-----------------------------------------------" << endl;
 return 0;
}


