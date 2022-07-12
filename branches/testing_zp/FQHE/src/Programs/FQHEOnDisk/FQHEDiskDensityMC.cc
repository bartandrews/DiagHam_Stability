#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHEDiskLaughlinOneQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHEDiskLaughlinOneJainQuasielectronWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasielectronWaveFunction.h"

#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;


void RandomZ (RealVector& positions, double scale, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomZOneCoordinate(RealVector& positions, double scale, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

bool RandomZOneCoordinateWithJump(RealVector& positions, double scale, double jump, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

bool RandomZOneCoordinateWithJumpRadial(RealVector& positions, double scale, double jump, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

// flip two one-body coordinates
//
// coordinates = reference on the n-body coordinate vector
// i = index of the first one body coordinate
// j = index of the second one body coordinate
void FlipCoordinates (RealVector& coordinates, int i, int j);

// find the index of a radius in a given array of radius
//
// radius = radius to find
// radiusArray = array of radius
// nbrRadius = number of radius in the array
// return value = radius index
int GetRadiusCoordinate (double radius, double* radiusArray, int nbrRadius);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEDiskDensityMC" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "laughlin", "do the calculation for the Laughlin state instead of the Moore-Read state");  
  (*SystemGroup) += new BooleanOption ('\n', "test-symmetry", "check the test wave function is symmetric/antisymmetric");  
  (*SystemGroup) += new SingleStringOption  ('\n', "load-permutations", "read all the permutations needed to compute the reference wave function from a file");  
  (*SystemGroup) += new SingleStringOption  ('\n', "save-permutations", "file name where all the permutations needed to compute the reference wave function have to be stored");  
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
  (*SystemGroup) += new BooleanOption ('\n', "quasielectron", "plot quasiparticles instead of the ground state");  
  (*SystemGroup) += new BooleanOption ('\n', "quasihole", "plot quasiholes instead of the ground state");  
  (*SystemGroup) += new SingleDoubleOption ('\n', "excitation-position", "position of the excitation from the center (in grid length unit)", 0.15);  
  (*SystemGroup) += new BooleanOption ('\n', "force-xsymmetry", "assume the wave function is invariant under the x <->-x symmetry");  
  (*SystemGroup) += new BooleanOption ('\n', "force-ysymmetry", "assume the wave function is invariant under the y <->-y symmetry");  
  (*SystemGroup) += new BooleanOption ('\n', "force-rsymmetry", "assume the wave function is invariant under rotation");  
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name", "density.dat"); 

  (*MonteCarloGroup) += new SingleIntegerOption  ('s', "nbr-step", "number of steps for the density profil in each direction", 20);
  (*MonteCarloGroup) += new SingleDoubleOption  ('\n', "grid-length", "size of the sampling grid (in unit of the magnetic length", 20);
  (*MonteCarloGroup) += new SingleDoubleOption  ('\n', "jump", "length of the jump used for the metropolis algorithm", 0.3);
  
  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new BooleanOption  ('r', "resume", "resume from a previous run");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);
  (*MonteCarloGroup) += new  SingleStringOption ('\n', "record-wavefunctions", "optional file where each wavefunctions will be tabulated and recorded");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskDensityMC -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  bool OverlapFlag = true;
  if (((BooleanOption*) Manager["test-symmetry"])->GetBoolean() == true)
    {
      OverlapFlag = false;
    }
  long NbrWarmUpIter = (long) ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  long NbrIter = (long) ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrSteps = ((SingleIntegerOption*) Manager["nbr-step"])->GetInteger();
  bool ResumeFlag = Manager.GetBoolean("resume");
  bool StatisticFlag = !(((BooleanOption*) Manager["boson"])->GetBoolean());
  bool QuasielectronFlag = (((BooleanOption*) Manager["quasielectron"])->GetBoolean());
  bool QuasiholeFlag = (((BooleanOption*) Manager["quasihole"])->GetBoolean());
  bool LaughlinFlag = (((BooleanOption*) Manager["laughlin"])->GetBoolean());
  char* RecordFileName = ((SingleStringOption*) Manager["record-file"])->GetString();
  long RecordStep = (long) Manager.GetInteger("record-step");
  if ((RecordFileName != 0) && (RecordStep > 0))
    {
      ofstream RecordFile;
      RecordFile.open(RecordFileName, ios::out | ios::binary);
      RecordFile.close();
    }
  else
    if (RecordFileName != 0)
      RecordFileName = 0;

  bool XSymmetryFlag = ((BooleanOption*) Manager["force-xsymmetry"])->GetBoolean();
  bool YSymmetryFlag = ((BooleanOption*) Manager["force-ysymmetry"])->GetBoolean();
  if ((XSymmetryFlag == true) || (YSymmetryFlag == true))
    {
      if ((NbrSteps  & 1) != 0)
	{
	  --NbrSteps; 
	  cout << "number of grid steps has to be even when using the --force-xsymmetry or --force-ysymmetry options, it will reduced to " << NbrSteps << endl;
	}
    }

  double GridScale =  ((SingleDoubleOption*) Manager["grid-length"])->GetDouble();
  double InvGridStep = ((double) NbrSteps) / GridScale;
  double GridStep = 1.0 / InvGridStep;
  double GridShift = -0.5 * GridScale;
  double MCJump = ((SingleDoubleOption*) Manager["jump"])->GetDouble();
  bool WaveFunctionMemory = false;

  bool RSymmetryFlag = ((BooleanOption*) Manager["force-rsymmetry"])->GetBoolean();
  if (RSymmetryFlag == true)
    {
      XSymmetryFlag = false;
      YSymmetryFlag = false;
    }
  double* RadiusArray = new double [NbrSteps + 1];
  double AreaStep = 2.0 * GridScale * GridScale / ((double) NbrSteps); 
  RadiusArray[0] = 0.0;
  for (int i = 1; i <= NbrSteps; ++i)
    {
//      RadiusArray[i] = sqrt ((RadiusArray[i - 1] * RadiusArray[i - 1]) + AreaStep); 
      RadiusArray[i] = RadiusArray[i - 1] + (0.5 * M_SQRT2 * GridScale / ((double) NbrSteps)); 
    }
  double SquareScale = RadiusArray[NbrSteps];
  double ExcitationPosition = ((SingleDoubleOption*) Manager["excitation-position"])->GetDouble();

  Abstract1DComplexFunction* SymmetrizedFunction = 0;
  if (QuasielectronFlag == true)
    {
      if (LaughlinFlag == true)
	{
	  if (StatisticFlag == true)
	    SymmetrizedFunction = new FQHEDiskLaughlinOneJainQuasielectronWaveFunction(NbrParticles, 0.0, GridScale, 3);
	  else
	    SymmetrizedFunction = new FQHEDiskLaughlinOneJainQuasielectronWaveFunction(NbrParticles, 0.0, GridScale, 2);
	}
      else
 	{
 	  if (((SingleStringOption*) Manager["load-permutations"])->GetString() == 0)
 	    SymmetrizedFunction = new PfaffianOnDiskTwoQuasielectronWaveFunction(NbrParticles,  - ExcitationPosition * GridScale, ExcitationPosition * GridScale, StatisticFlag);
 	  else
 	    SymmetrizedFunction = new PfaffianOnDiskTwoQuasielectronWaveFunction(((SingleStringOption*) Manager["load-permutations"])->GetString(),
										 -ExcitationPosition * GridScale, ExcitationPosition * GridScale, StatisticFlag);
 	  if (((SingleStringOption*) Manager["save-permutations"])->GetString() != 0)
 	    {
 	      ((PfaffianOnDiskTwoQuasielectronWaveFunction*) SymmetrizedFunction)->WritePermutations(((SingleStringOption*) Manager["save-permutations"])->GetString());
 	      return 0;
 	    }
	  WaveFunctionMemory = false;
 	}
    }
  else
    if (QuasiholeFlag == true)
      {
	if (LaughlinFlag == true)
	  if (StatisticFlag == true)
	    SymmetrizedFunction = new FQHEDiskLaughlinOneQuasiholeWaveFunction(NbrParticles, 0.0, 3);
	  else
	    SymmetrizedFunction = new FQHEDiskLaughlinOneQuasiholeWaveFunction(NbrParticles, 0.0, 2);
	else
	  SymmetrizedFunction = new PfaffianOnDiskTwoQuasiholeWaveFunction(NbrParticles, -ExcitationPosition * GridScale, ExcitationPosition * GridScale, StatisticFlag);
      }
    else
      {
	if (LaughlinFlag == true)
	  if (StatisticFlag == true)
	    SymmetrizedFunction = new LaughlinOnDiskWaveFunction(NbrParticles, 3, 0.5 * GridScale); 
	  else
	    SymmetrizedFunction = new LaughlinOnDiskWaveFunction(NbrParticles, 2, 0.5 * GridScale);
	else
	  SymmetrizedFunction = new PfaffianOnDiskWaveFunction(NbrParticles);
      }

  AbstractRandomNumberGenerator* RandomNumber = 0;
  if (((SingleStringOption*) Manager["random-file"])->GetString() != 0)
    {
      if (ResumeFlag == true)
	{
	  ifstream MCState;
	  MCState.open("mcstate.dat", ios::in | ios::binary);
	  int TmpNbrIter;
	  ReadLittleEndian(MCState, TmpNbrIter);
	  unsigned long TmpNumber;
	  ReadLittleEndian(MCState, TmpNumber);
	  MCState.close();
	  RandomNumber = new FileRandomNumberGenerator(((SingleStringOption*) Manager["random-file"])->GetString(), TmpNumber, 
						       ((SingleIntegerOption*) Manager["random-seek"])->GetInteger());	  
	}
      else
	RandomNumber = new FileRandomNumberGenerator(((SingleStringOption*) Manager["random-file"])->GetString(), (NbrWarmUpIter * 4) + (NbrIter * 4) + 2000, 
						     ((SingleIntegerOption*) Manager["random-seek"])->GetInteger());
    }
  else
    {
      //      RandomNumber = new StdlibRandomNumberGenerator (29457);
      RandomNumber = new NumRecRandomGenerator(29457);
    }

  if (OverlapFlag == true)
    {
      int TwiceNbrParticles = NbrParticles * 2;
      RealVector TmpZ (TwiceNbrParticles, true);
      double PreviousTmpZRe;
      double PreviousTmpZIm;
      double** FunctionBasisDecomposition = new double* [NbrSteps];
      for (int k = 0; k < NbrSteps; ++k)
	{
	  FunctionBasisDecomposition[k] = new double [NbrSteps];
	  for (int j = 0; j < NbrSteps; ++j)
	    FunctionBasisDecomposition[k][j] = 0.0;
	}
      double TmpExp = 0.0;

      int* GridLocations = new int [TwiceNbrParticles];

      double TotalProbability = 0.0;
      double TotalProbabilityError = 0.0;
      long CurrentPercent = 0;

      RandomZ (TmpZ, GridScale, NbrParticles, RandomNumber);
      int NextCoordinate = 0;
      double PreviousProbabilities = 0.0;
      double CurrentProbabilities = 0.0;
      TotalProbability = 0.0;
      int Acceptance = 0;	  
      
      Complex Tmp = (*SymmetrizedFunction)(TmpZ);
      CurrentProbabilities = SqrNorm(Tmp);
      PreviousProbabilities = CurrentProbabilities;

      cout << "starting warm-up sequence" << endl;
      for (long i = 1; i < NbrWarmUpIter; ++i)
	{
	  PreviousTmpZRe = TmpZ[NextCoordinate << 1];
	  PreviousTmpZIm = TmpZ[(NextCoordinate << 1) + 1];
	  bool TmpFlag = false;
	  if (RSymmetryFlag == false)
	    TmpFlag = RandomZOneCoordinateWithJump(TmpZ, GridScale, MCJump, NextCoordinate, RandomNumber);
	  else
	    TmpFlag = RandomZOneCoordinateWithJumpRadial(TmpZ, SquareScale, MCJump, NextCoordinate, RandomNumber);
	  if (TmpFlag == true)
	    {
	      if (WaveFunctionMemory == true)
		((PfaffianOnDiskTwoQuasielectronWaveFunction*) SymmetrizedFunction)->SetNextCoordinate(NextCoordinate);
	      Complex TmpMetropolis = (*SymmetrizedFunction)(TmpZ);
	      CurrentProbabilities = SqrNorm(TmpMetropolis);
	      TmpExp = exp(-0.5 * (((TmpZ[(NextCoordinate << 1)] * TmpZ[(NextCoordinate << 1)]) - (PreviousTmpZRe * PreviousTmpZRe))
				   + ((TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])) - (PreviousTmpZIm * PreviousTmpZIm)));
	      CurrentProbabilities *= TmpExp;
	      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
		{
		  PreviousProbabilities = CurrentProbabilities;
		  ++Acceptance;
		}
	      else
		{
		  TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
		  TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
		  CurrentProbabilities = PreviousProbabilities;
		  if (WaveFunctionMemory == true)
		    ((PfaffianOnDiskTwoQuasielectronWaveFunction*) SymmetrizedFunction)->RestorePreviousData();
		}
	    }
	  else
	    {
	      TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
	      TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
	    }
	  NextCoordinate = (int) (RandomNumber->GetRealRandomNumber() * (double) NbrParticles);
	  if ((( i * 20l) / NbrWarmUpIter) != CurrentPercent)
	    {
	      CurrentPercent = (i * 20l) / NbrWarmUpIter;
	      cout << (CurrentPercent * 5l) << "% " << flush;
	    }
	} 
      cout << endl << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrWarmUpIter)) * 100.0) << "%" <<endl;
      
      CurrentPercent = 0l;
      for (int i = 0; i < TwiceNbrParticles; ++i)
	GridLocations[i] = (int) ((TmpZ[i] - GridShift) * InvGridStep);

      cout << "starting MC sequence" << endl;
      Acceptance = 0;

      for (long i = 0; i < NbrIter; ++i)
	{
	  PreviousTmpZRe = TmpZ[NextCoordinate << 1];
	  PreviousTmpZIm = TmpZ[(NextCoordinate << 1) + 1];
	  bool TmpFlag = false;
	  if (RSymmetryFlag == false)
	    TmpFlag = RandomZOneCoordinateWithJump(TmpZ, GridScale, MCJump, NextCoordinate, RandomNumber);
	  else
	    TmpFlag = RandomZOneCoordinateWithJumpRadial(TmpZ, SquareScale, MCJump, NextCoordinate, RandomNumber);
	  if (TmpFlag == true)
	    {
	      if (WaveFunctionMemory == true)
		((PfaffianOnDiskTwoQuasielectronWaveFunction*) SymmetrizedFunction)->SetNextCoordinate(NextCoordinate);
	      Complex TmpMetropolis = (*SymmetrizedFunction)(TmpZ);
	      CurrentProbabilities = SqrNorm(TmpMetropolis);
	      TmpExp = exp(-0.5 * (((TmpZ[(NextCoordinate << 1)] * TmpZ[(NextCoordinate << 1)]) - (PreviousTmpZRe * PreviousTmpZRe))
				   + ((TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])) - (PreviousTmpZIm * PreviousTmpZIm)));
	      CurrentProbabilities *= TmpExp;
	      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
		{
		  PreviousProbabilities = CurrentProbabilities;
		  if (RSymmetryFlag == false)
		    {
		      GridLocations[NextCoordinate << 1] = (int) ((TmpZ[NextCoordinate << 1] - GridShift) * InvGridStep);
		      GridLocations[(NextCoordinate << 1) + 1] = (int) ((TmpZ[(NextCoordinate << 1) + 1] - GridShift) * InvGridStep);
		    }
		  else
		    {
		      GridLocations[NextCoordinate] = GetRadiusCoordinate(sqrt((TmpZ[NextCoordinate << 1] * TmpZ[NextCoordinate << 1]) + (TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])), RadiusArray, NbrSteps);
		    }
		  ++Acceptance;	      
		}
	      else
		{
		  TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
		  TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
		  CurrentProbabilities = PreviousProbabilities;
		  if (WaveFunctionMemory == true)
		    ((PfaffianOnDiskTwoQuasielectronWaveFunction*) SymmetrizedFunction)->RestorePreviousData();
		}
	    }
	  else
	    {
	      TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
	      TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
	    }

	  if (RSymmetryFlag == false)
	    {
	      for (int j = 0; j < TwiceNbrParticles; j += 2)
		FunctionBasisDecomposition[GridLocations[j]][GridLocations[j + 1]] += 1.0;
	    }
	  else
	    {
	      for (int j = 0; j < NbrParticles; ++j)
		FunctionBasisDecomposition[0][GridLocations[j]] += sqrt((TmpZ[j << 1] * TmpZ[j << 1]) + (TmpZ[(j << 1) + 1] * TmpZ[(j << 1) + 1]));
	    }
	  TotalProbability++;	  
	  TotalProbabilityError++;

	  NextCoordinate = (int) (RandomNumber->GetRealRandomNumber() * (double) NbrParticles);

	  if (((i * 20l) / NbrIter) != CurrentPercent)
	    {
	      CurrentPercent = (i * 20l) / NbrIter;
	      cout << (CurrentPercent * 5l) << "% " << flush;
	    }

	  if ((RecordFileName != 0) && ((i % RecordStep) == 0l))
	    {
	      ofstream RecordFile;
	      RecordFile.open(RecordFileName, ios::out | ios::binary | ios::app);
	      RecordFile.precision(14);
	      RecordFile << i << " " << TotalProbability;
	      if (RSymmetryFlag == false)
		{
		  for (int j = 0; j < NbrSteps; ++j)
		    for (int k = 0; k < NbrSteps; ++k)
		      RecordFile << " " << FunctionBasisDecomposition[j][k];
		}
	      else
		{
		  for (int k = 0; k < NbrSteps; ++k)
		    RecordFile << " " << FunctionBasisDecomposition[0][k];
		}
	      RecordFile << endl;
	      RecordFile.close();
	    }
	}
      cout << endl << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrIter)) * 100.0) << "%" <<endl;

      TotalProbabilityError /= (double) NbrIter;
      TotalProbabilityError -= (TotalProbability * TotalProbability) / (((double) NbrIter) * ((double) NbrIter));
      TotalProbabilityError = sqrt(TotalProbabilityError) / (TotalProbability / ((double) NbrIter));
      TotalProbabilityError /= sqrt ((double) NbrIter);

      TotalProbability =  1.0 / TotalProbability;
     if (RSymmetryFlag == false)
       {
	 if (XSymmetryFlag == true) 
	   {
	     int HalfNbrSteps = NbrSteps >> 1;
	     if (YSymmetryFlag == true) 
	       {
		 for (int i = 0; i < HalfNbrSteps; ++i)
		   for (int j = 0; j < HalfNbrSteps; ++j)	  
		     {
		       double TmpTotal = FunctionBasisDecomposition[i][j];
		       TmpTotal += FunctionBasisDecomposition[NbrSteps - i - 1][NbrSteps - j - 1];
		       TmpTotal += FunctionBasisDecomposition[i][NbrSteps - j - 1];
		       TmpTotal += FunctionBasisDecomposition[NbrSteps - i - 1][j];
		       TmpTotal *= 0.25;
		       TmpTotal *= TotalProbability;
		       FunctionBasisDecomposition[i][j] = TmpTotal;
		       FunctionBasisDecomposition[NbrSteps - i - 1][NbrSteps - j - 1] = TmpTotal;
		       FunctionBasisDecomposition[i][NbrSteps - j - 1] = TmpTotal;
		       FunctionBasisDecomposition[NbrSteps - i - 1][j] = TmpTotal;
		     }
	       }
	     else
	       {
		 for (int i = 0; i < HalfNbrSteps; ++i)
		   for (int j = 0; j < NbrSteps; ++j)	  
		     {
		       double TmpTotal = FunctionBasisDecomposition[i][j];
		       TmpTotal += FunctionBasisDecomposition[NbrSteps - i - 1][j];
		       TmpTotal *= 0.5;
		       TmpTotal *= TotalProbability;
		       FunctionBasisDecomposition[i][j] = TmpTotal;
		       FunctionBasisDecomposition[NbrSteps - i - 1][j] = TmpTotal;
		     }
	       }
	   }
	 else
	   if (YSymmetryFlag == true) 
	     {
	       int HalfNbrSteps = NbrSteps >> 1;
	       for (int i = 0; i < NbrSteps; ++i)
		 for (int j = 0; j < HalfNbrSteps; ++j)	  
		   {
		     double TmpTotal = FunctionBasisDecomposition[i][j];
		     TmpTotal += FunctionBasisDecomposition[i][NbrSteps - j - 1];
		     TmpTotal *= 0.5;
		     TmpTotal *= TotalProbability;
		     FunctionBasisDecomposition[i][j] = TmpTotal;
		     FunctionBasisDecomposition[i][NbrSteps - j - 1] = TmpTotal;
		   }
	     }
	   else
	     {
	       for (int i = 0; i < NbrSteps; ++i)
		 for (int j = 0; j < NbrSteps; ++j)	  
		   FunctionBasisDecomposition[i][j] *= TotalProbability;
	     }
       }
     else
       {
	 for (int i = 0; i < NbrSteps; ++i)
	   {
	     FunctionBasisDecomposition[0][i] *= TotalProbability;
	     FunctionBasisDecomposition[0][i] /= (M_PI * ((RadiusArray[i + 1] * RadiusArray[i + 1]) - (RadiusArray[i] * RadiusArray[i])));
	   }
       }
     
     ofstream DensityRecordFile;
      DensityRecordFile.precision(14);
      DensityRecordFile.open(((SingleStringOption*) Manager["output"])->GetString(), ios::out);

      Manager.DisplayOption(DensityRecordFile, true, '#');
      DensityRecordFile << "#" << endl;
      double Sum = 0.0;

      if (RSymmetryFlag == false)
	{
	  DensityRecordFile << "#" << endl << "# density wave function " << endl << "# x y  density density_error" << endl;
	  double TmpX = GridShift + (0.5 * GridStep);
	  for (int i = 0; i < NbrSteps; ++i)
	    {
	      double TmpY = GridShift + (0.5 * GridStep);
	      for (int j = 0; j < NbrSteps; ++j)
		{
		  DensityRecordFile << TmpX << " " << TmpY << " " << FunctionBasisDecomposition[i][j] << endl;
		  Sum += FunctionBasisDecomposition[i][j];
		  TmpY += GridStep;
		}
	      DensityRecordFile << endl;
	      TmpX += GridStep;
	    }
	}
      else 
	{
	  DensityRecordFile << "#" << endl << "# density wave function " << endl << "# radius  density density_error" << endl;
	  for (int j = 0; j < NbrSteps; ++j)
	    {
	      DensityRecordFile << "# " << RadiusArray[j] << " " << FunctionBasisDecomposition[0][j] << endl;
	    }
	  DensityRecordFile << "#" << endl << "# density wave function " << endl << "# x y  density density_error" << endl;
	  double TmpX = GridShift + (0.5 * GridStep);
	  for (int i = 0; i < NbrSteps; ++i)
	    {
	      double TmpY = GridShift + (0.5 * GridStep);
	      for (int j = 0; j < NbrSteps; ++j)
		{
		  int TmpCoordinate = GetRadiusCoordinate(sqrt((TmpX * TmpX) + (TmpY * TmpY)), RadiusArray, NbrSteps);
		  DensityRecordFile << TmpX << " " << TmpY << " " << FunctionBasisDecomposition[0][TmpCoordinate] << endl;
		  //		  FunctionBasisDecomposition[0][i] /= (M_PI * ((RadiusArray[i + 1] * RadiusArray[i + 1]) - (RadiusArray[i] * RadiusArray[i])));
//		  Sum += FunctionBasisDecomposition[0][TmpCoordinate] * GridStep * GridStep;
		  TmpY += GridStep;
		}
	      DensityRecordFile << endl;
	      TmpX += GridStep;
	    }	  
	}
      DensityRecordFile.close();
      cout << "Sum = " << Sum << endl;
    }
   else
     {
       RealVector TmpZ (2 * NbrParticles, true);
       RandomZ (TmpZ, GridScale, NbrParticles, RandomNumber);
       cout << (*SymmetrizedFunction)(TmpZ) << endl;;
       for (int i = 0; i < NbrParticles; ++i)
	 {
	   for (int j = i + 1; j < NbrParticles; ++j)
	     {
	       cout << i << "<->" << j << " : ";
	       cout << (*SymmetrizedFunction)(TmpZ) << " | " ;
	       FlipCoordinates(TmpZ, i, j);
	       cout << (*SymmetrizedFunction)(TmpZ) << endl;
	       FlipCoordinates(TmpZ, i, j);	       
	     }
	 }
       for (int i = 1; i < NbrParticles; ++i)
	 {
	   TmpZ[(i << 1)] = TmpZ[0];
	   TmpZ[(i << 1) + 1] = TmpZ[1];
	   cout << (i  + 1) << " body cancellation : " << (*SymmetrizedFunction)(TmpZ) << endl;
	 }
     }
  return 0;
}

void RandomZ (RealVector& positions, double scale, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      positions[2 * j] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
      positions[(2 * j) + 1] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
    }
}

void RandomZOneCoordinate(RealVector& positions, double scale, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  positions[coordinate] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
  ++coordinate;
  positions[coordinate] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
}

bool RandomZOneCoordinateWithJump(RealVector& positions, double scale, double jump, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double OldCoordinate =  positions[coordinate];
  double TmpJump = scale * jump * (0.5 - randomNumberGenerator->GetRealRandomNumber());
  if (((OldCoordinate + TmpJump) > (0.5 * scale)) || ((OldCoordinate + TmpJump) < -(0.5 * scale)))
    return false;
  positions[coordinate] += TmpJump;
  ++coordinate;
  OldCoordinate =  positions[coordinate];
  TmpJump = scale * jump * (0.5 - randomNumberGenerator->GetRealRandomNumber()); 
  if (((OldCoordinate + TmpJump) > (0.5 * scale)) || ((OldCoordinate + TmpJump) < -(0.5 * scale)))
    return false;
  positions[coordinate] += TmpJump;
  return true;
}

bool RandomZOneCoordinateWithJumpRadial(RealVector& positions, double scale, double jump, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double& OldCoordinate1 =  positions[coordinate];
  double TmpJump1 = scale * jump * (0.5 - randomNumberGenerator->GetRealRandomNumber());
  double& OldCoordinate2 =  positions[coordinate + 1];
  double TmpJump2 = scale * jump * (0.5 - randomNumberGenerator->GetRealRandomNumber());
  if ((((OldCoordinate1 + TmpJump1) * (OldCoordinate1 + TmpJump1)) + ((OldCoordinate2 + TmpJump2) * (OldCoordinate2 + TmpJump2))) > (scale * scale))
    return false;
  OldCoordinate1 += TmpJump1;
  OldCoordinate2 += TmpJump2;
  return true;
}

// flip two one-body coordinates
//
// coordinates = reference on the n-body coordinate vector
// i = index of the first one body coordinate
// j = index of the second one body coordinate

void FlipCoordinates (RealVector& coordinates, int i, int j)
{
  double Tmp = coordinates[i << 1];  
  coordinates[i << 1] = coordinates[j << 1];
  coordinates[j << 1] = Tmp;
  Tmp = coordinates[(i << 1) + 1];  
  coordinates[(i << 1) + 1] = coordinates[(j << 1) + 1];
  coordinates[(j << 1) + 1] = Tmp;
}


// get locations on the grid from a position vector
//
// coordinates = reference on the n-body coordinate vector
// locations = array of integer location on the grid
// invGridStep = sinvert of the grid cell size
// gridShift = X (or Y) coordinate of the grid upper left corner

void GetGridCoordinates (RealVector& coordinates, int* locations, double invGridStep, double gridShift)
{
  for (int i = 0; i < coordinates.GetVectorDimension(); ++i)
    locations[i] = (int) ((coordinates[i] - gridShift) * invGridStep);
}

// find the index of a radius in a given array of radius
//
// radius = radius to find
// radiusArray = array of radius
// nbrRadius = number of radius in the array
// return value = radius index

int GetRadiusCoordinate (double radius, double* radiusArray, int nbrRadius)
{
  --nbrRadius;
  if (radius >= radiusArray[nbrRadius])
    return nbrRadius;
  int Pos = 0;
  while (Pos < (nbrRadius - 1))
    {
      int Mid = (Pos + nbrRadius) >> 1;
      if (radius < radiusArray[Mid])
	nbrRadius = Mid;
      else 
	Pos = Mid;
    }
  if (radius < radiusArray[nbrRadius])
    return Pos;
  else
    return nbrRadius;
}
