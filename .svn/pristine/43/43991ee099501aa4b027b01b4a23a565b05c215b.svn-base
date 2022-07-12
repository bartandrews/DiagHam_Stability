#include "Vector/RealVector.h"

#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/SU4HalperinOnSphereWaveFunction.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "GeneralTools/ConfigurationParser.h"

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
  OptionManager Manager ("FQHESphereFermionMonteCarloEnergy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "config-file", "name of the file that describes the wave function");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for energy calculation", false);
  (*MiscGroup) += new BooleanOption  ('h', "test", "test Monte Carlo observable");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionMonteCarloEnergy -h" << endl;
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
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();

  Abstract1DComplexFunction* WaveFunction = 0;

//  Abstract1DComplexFunction* WaveFunction = new LaughlinOnSphereWaveFunction(NbrFermions, 3);
//  Abstract1DComplexFunction* WaveFunction = new HalperinOnSphereWaveFunction(NbrFermions - 8, 8, 3, 3, 2);
//  Abstract1DComplexFunction* WaveFunction = new SU4HalperinOnSphereWaveFunction(NbrFermions / 4, NbrFermions / 4, NbrFermions / 4, NbrFermions / 4, 
//										3, 3, 3, 3, 3, 3, 3);
  if (((SingleStringOption*) Manager["config-file"])->GetString() == 0)
    {
      WaveFunction = new SU4HalperinOnSphereWaveFunction(8, 5, 4, 3, 3, 3, 3, 3, 3, 2, 2);
    }
  else
    {
      ConfigurationParser WaveFunctionDefinition ;
      if (WaveFunctionDefinition.Parse(((SingleStringOption*) Manager["config-file"])->GetString()) == false)
	{
	  WaveFunctionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      bool ErrorFlag;
      WaveFunction = new SU4HalperinOnSphereWaveFunction(WaveFunctionDefinition, ErrorFlag, NbrFermions, LzMax);
      if (ErrorFlag == false)
	{
	  cout << "error while parsing " << ((SingleStringOption*) Manager["config-file"])->GetString() << endl;
	  return -1;
	}
    }
  int TwiceNbrFermions = NbrFermions << 1;

  RealVector Location(TwiceNbrFermions, true);

  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);

  double Energy = 0.0;
  double EnergyError = 0.0;
  double PositiveBackgroundEnergy = ((0.5 * ((double) (NbrFermions * NbrFermions))) / sqrt(0.5 * ((double) LzMax)));
  double Tmp = 0.0;
  double Dist;
  
  int RecordStep = ((SingleIntegerOption*) Manager["record-step"])->GetInteger();
  double* RecordedEnergy = 0;
  if (RecordStep != 0)
    RecordedEnergy = new double [(NbrIter / RecordStep) + 1];
  int RecordIndex = 0;

  Complex Tmp3;
  int NextCoordinates = 0;
  for (int j = 0; j < TwiceNbrFermions; j += 2)
    {
      Location[j] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));
      Location[1 + j] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
    }
  double PreviousProbabilities = SqrNorm((*WaveFunction)(Location));
  double CurrentProbabilities = PreviousProbabilities;
  double MonteCarloDistance;
  double MinGaussianVariance = 0.0;
  double MaxGaussianVariance = 2.0;
  double GaussianVariance = 1.0;
  double GaussianVarianceStep = 0.01;
  double GaussianFactor =  pow(2.0 * M_PI * GaussianVariance * (1.0 - exp (-2.0 / GaussianVariance)), -0.5 * ((double) NbrFermions)); 
  
  double PreviousCoordinates1;
  double PreviousCoordinates2;
  double CosTheta;
  double SinTheta;
  double Phi;
  int Acceptance = 0;
  double AcceptanceRate = 1.0;
  double DistFactor = 1.0 / sqrt(((double) LzMax));
  int DisplayStep = ((SingleIntegerOption*) Manager["display-step"])->GetInteger();

  if (NbrWarmUpIter > 0)
    cout << "starting warm-up sequence" << endl;
  for (int i = 1; i < NbrWarmUpIter; ++i)
    {      
      PreviousCoordinates1 = Location[NextCoordinates << 1];
      PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      CurrentProbabilities = SqrNorm((*WaveFunction)(Location));
      MonteCarloDistance = 1.0;
//       MonteCarloDistance = exp (- 2.0 * ((sin(0.5 * (Location[NextCoordinates << 1])) * sin(0.5 * (Location[NextCoordinates << 1]))) - 
// 					 (sin(0.5 * PreviousCoordinates1) * sin(0.5 * PreviousCoordinates1))) / GaussianVariance);
//      cout << CurrentProbabilities << " " << PreviousProbabilities << endl;
      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < (CurrentProbabilities * MonteCarloDistance)))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  ++Acceptance;
	}
      else
 	{
 	  Location[NextCoordinates << 1] = PreviousCoordinates1;
 	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
 	  CurrentProbabilities = PreviousProbabilities;
 	}
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if ((i % 1000) == 0)
	{
	  AcceptanceRate = ((double) Acceptance) / ((double) i);
	  cout << Acceptance << " / " << i << " = " << ((100.0 * ((double) Acceptance)) / ((double) i)) << "%" << " (" << GaussianVariance << ")" << endl;
// 	  if (NewAcceptanceRate < AcceptanceRate)
// 	    {
// 	      GaussianVariance -= GaussianVarianceStep;
// 	      GaussianVarianceStep *= 0.5;
// 	    }
	  if (AcceptanceRate < 0.5)
	    MaxGaussianVariance = GaussianVariance;
	  else
	    MinGaussianVariance = GaussianVariance;
	  GaussianVariance = 0.5 * (MaxGaussianVariance + MinGaussianVariance);
//	  if (NewAcceptanceRate > 0.6)
	    
// 	  GaussianVariance += GaussianVarianceStep;	      
// 	  AcceptanceRate = NewAcceptanceRate;
//	  GaussianVariance = 0.5 * (MaxGaussianVariance + MinGaussianVariance);
//	  GaussianFactor =  pow(2.0 * M_PI * GaussianVariance * (1.0 - exp (-2.0 / GaussianVariance)), -0.5 * ((double) NbrFermions)); 	  
	}
    }
  if (NbrWarmUpIter > 0)
    cout << "warm-up sequence is over" << endl;
  Acceptance = 0;
  GaussianVariance = 0.0018;

  for (int i = 1; i < NbrIter; ++i)
    {      
      PreviousCoordinates1 = Location[NextCoordinates << 1];
      PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      CurrentProbabilities = SqrNorm((*WaveFunction)(Location));
      MonteCarloDistance = 1.0;
//       MonteCarloDistance = exp (- 2.0 * ((sin(0.5 * (Location[NextCoordinates << 1])) * sin(0.5 * (Location[NextCoordinates << 1]))) - 
// 					 (sin(0.5 * PreviousCoordinates1) * sin(0.5 * PreviousCoordinates1))) / GaussianVariance);
      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < (CurrentProbabilities * MonteCarloDistance)))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  ++Acceptance;
	}
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
      for (int j = 0; j < TwiceNbrFermions; j += 2)
	{
	  Phi = Location[j + 1];
	  CosTheta = Location[j];
	  SinTheta = sin(CosTheta);
	  CosTheta = cos(CosTheta);
	  for (int k = j + 2; k < TwiceNbrFermions; k += 2)
	    {
	      Dist = sqrt ((1.0 - (CosTheta * cos(Location[k])) - (SinTheta * sin(Location[k]) * cos(Phi - Location[1 + k]))));
	      if (Dist > 1e-6)
		Tmp += 1.0 / Dist;
	    }
	}
      Tmp *= DistFactor;
      Energy += Tmp;
      EnergyError += Tmp * Tmp;
      if ((RecordStep != 0) && ((i % RecordStep) == 0))
	{
	  RecordedEnergy[RecordIndex] = ((Energy / ((double) i))  - PositiveBackgroundEnergy) / ((double) NbrFermions);
	  ++RecordIndex;
	}
      if ((i % DisplayStep) == 0)
 	{
	  cout << " i = " << i << endl;
	  double TmpEnergy = Energy / ((double) i);
	  double TmpEnergyError = sqrt (((EnergyError  / ((double) i)) - (TmpEnergy * TmpEnergy)) / ((double) i));
	  cout << ((TmpEnergy - ((0.5 * ((double) (NbrFermions * NbrFermions))) / sqrt(0.5 * ((double) LzMax)))) / ((double) NbrFermions)) << " +/- " << (TmpEnergyError / ((double) NbrFermions)) << endl;
 	  cout << "raw : " << Energy << " " << (TmpEnergy * TmpEnergy) << " " << (EnergyError  / ((double) i)) << " " << TmpEnergy << " +/- " << TmpEnergyError << endl;
	  cout << "acceptance rate = " << (((double) Acceptance) / (((double) i))) << endl;
 	  cout << "-----------------------------------------------" << endl;
 	}
    }
  cout << " final results :" << endl;
  double TmpEnergy = Energy / ((double) NbrIter);
  double TmpEnergyError = sqrt (((EnergyError  / ((double) NbrIter)) - (TmpEnergy * TmpEnergy)) / ((double) NbrIter));
  cout << TmpEnergy << " +/- " << TmpEnergyError << endl;
  cout << ((TmpEnergy - PositiveBackgroundEnergy) / ((double) NbrFermions)) << " +/- " << (TmpEnergyError / ((double) NbrFermions)) << endl;
  cout << "-----------------------------------------------" << endl;

  if (RecordStep != 0)
    {
      RecordedEnergy[RecordIndex] = ((Energy / ((double) NbrIter))  - PositiveBackgroundEnergy) / ((double) NbrFermions);
      ofstream EnergyRecordFile;
      EnergyRecordFile.precision(14);
      EnergyRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary);
      int NbrRecords = NbrIter / RecordStep;
      for (int i = 0; i < NbrRecords; ++i)
	EnergyRecordFile << i << " " << RecordedEnergy[i] << endl;
      EnergyRecordFile.close();
    }

 return 0;
}


