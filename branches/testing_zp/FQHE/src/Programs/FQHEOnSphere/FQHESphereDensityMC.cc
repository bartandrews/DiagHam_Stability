#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
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

#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU2HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU3HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU3GeneralizedGaffnianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU4HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/HundRuleCFStates.h"
#include "Tools/FQHEWaveFunction/SU3HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasielectronWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESphereLaughlinOneQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESphereLaughlinOneQuasielectronWaveFunction.h"

#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

#define M_2PI 6.283185307179586477

using std::ios;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;


void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, double jumpDistance, 
			   AbstractRandomNumberGenerator* randomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);

void FlipLzMinusLz (ComplexVector& uv);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereDensityMC" , "0.01");
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
  (*SystemGroup) += new BooleanOption ('\n', "force-symmetry", "assume the wave function is invariant under the Lz <->-Lz symmetry");  
  (*SystemGroup) += new SingleStringOption  ('\n', "load-permutations", "read all the permutations needed to compute the reference wave function from a file");  
  (*SystemGroup) += new SingleStringOption  ('\n', "save-permutations", "file name where all the permutations needed to compute the reference wave function have to be stored");  
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
  (*SystemGroup) += new BooleanOption ('\n', "quasielectron", "plot quasiparticles instead of quasiholes");  
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name", "density.dat"); 

  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals used to sample the sphere geometry", 100);
  (*MonteCarloGroup) += new SingleIntegerOption  ('s', "nbr-step", "number of steps for the density profil", 100);
  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleDoubleOption  ('\n', "jump", "length of the jump used for the metropolis algorithm", 0.3);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new SingleDoubleOption  ('\n', "jump", "length of the jump used for the metropolis algorithm", 0.3);
  (*MonteCarloGroup) += new BooleanOption  ('r', "resume", "resume from a previous run");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "weight-symmetrized" , "use the norm of the symmetrized wave fonction as probalbility density instead of the exact state");
  (*MonteCarloGroup) += new  SingleStringOption ('\n', "record-wavefunctions", "optional file where each wavefunctions will be tabulated and recorded");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereDensityMC -h" << endl;
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
  bool SymmetryFlag = false;
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrSteps = ((SingleIntegerOption*) Manager["nbr-step"])->GetInteger();
  int NbrOrbitals = ((SingleIntegerOption*) Manager["nbr-orbitals"])->GetInteger();
  bool ResumeFlag = Manager.GetBoolean("resume");
  bool StatisticFlag = !(((BooleanOption*) Manager["boson"])->GetBoolean());
  bool QuasielectronFlag = (((BooleanOption*) Manager["quasielectron"])->GetBoolean());
  bool LaughlinFlag = (((BooleanOption*) Manager["laughlin"])->GetBoolean());
  double JumpDistance = ((SingleDoubleOption*) Manager["jump"])->GetDouble();
  double StepSize = 0.25;
  if (((BooleanOption*) Manager["force-symmetry"])->GetBoolean() == true)
    {
      SymmetryFlag = true;
      if ((NbrOrbitals  & 1) != 0)
	{
	  --NbrOrbitals; 
	  cout << "number of orbitals has to be even when using the --force-symmetry option, it will reduced to " << NbrOrbitals << endl;
	}
    }
  char* RecordWaveFunctions = ((SingleStringOption*) Manager["record-wavefunctions"])->GetString();
  if (RecordWaveFunctions != 0)
    {
      ofstream RecordFile;
      RecordFile.open(RecordWaveFunctions, ios::out | ios::binary);
      RecordFile.close();
    }
  int RecordStep = Manager.GetInteger("record-step");

  double GridScale =  M_PI;
  double InvGridStep = ((double) NbrOrbitals) / GridScale;
  double GridStep = 1.0 / InvGridStep;
  double GridShift = -0.5 * GridScale;



//   int LzMax = 2 * NbrParticles - 2;
//   if (QuasielectronFlag == true)
//     LzMax -= 2;
//   AbstractQHEParticle* ExactSpace = new FermionOnSphere (NbrParticles, 0, LzMax);
//   RealVector ExactState;
//   ExactState.ReadVector ("fermions_hardcore_nbody_3_n_6_2s_9_lz_0.0.vec");

  Abstract1DComplexFunctionOnSphere* SymmetrizedFunction = 0;
  if (QuasielectronFlag == true)
    {
      if (LaughlinFlag == true)
	SymmetrizedFunction = new FQHESphereLaughlinOneQuasielectronWaveFunction(NbrParticles, 0.0, 0.0, StatisticFlag);
      else
	{
	  if (((SingleStringOption*) Manager["load-permutations"])->GetString() == 0)
	    SymmetrizedFunction = new PfaffianOnSphereTwoQuasielectronWaveFunction(NbrParticles, 0.0, 0.0, M_PI, 0.0, StatisticFlag);
	  else
	    SymmetrizedFunction = new PfaffianOnSphereTwoQuasielectronWaveFunction(((SingleStringOption*) Manager["load-permutations"])->GetString(),
										   0.0, 0.0, M_PI, 0.0, StatisticFlag);
	  if (((SingleStringOption*) Manager["save-permutations"])->GetString() != 0)
	    {
	      ((PfaffianOnSphereTwoQuasielectronWaveFunction*) SymmetrizedFunction)->WritePermutations(((SingleStringOption*) Manager["save-permutations"])->GetString());
	      return 0;
	    }
	}
    }
  else
    {
      if (LaughlinFlag == true)
	SymmetrizedFunction = new FQHESphereLaughlinOneQuasiholeWaveFunction(NbrParticles, 0.0, 0.0, StatisticFlag);
      else
	SymmetrizedFunction = new PfaffianOnSphereTwoQuasiholeWaveFunction(NbrParticles, 0.0, 0.0, M_PI, 0.0, StatisticFlag);
    }
  Abstract1DComplexFunctionOnSphere* TestFunction;
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
      RealVector TmpPositions (NbrParticles * 2, true);
      RealVector PreviousTmpPositions (NbrParticles * 2, true);
      ComplexVector TmpUV (NbrParticles * 2, true);
      ComplexVector PreviousTmpUV (NbrParticles * 2, true);
      Complex PreviousTmpU;
      Complex PreviousTmpV;
      double PreviousTmpPositionTheta;
      double  PreviousTmpPositionPhi;
      ParticleOnSphereFunctionBasis FunctionBasis(NbrOrbitals - 1);
      Complex** FunctionBasisEvaluation = new Complex* [NbrParticles];
      for (int k = 0; k < NbrParticles; ++k)
	FunctionBasisEvaluation[k] = new Complex [NbrOrbitals];
      double* FunctionBasisDecomposition = new double [NbrOrbitals];
      double* TmpFunctionBasisDecomposition = new double [NbrOrbitals];
      double* FunctionBasisDecompositionError  = new double [NbrOrbitals];
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  FunctionBasisDecomposition[k] = 0.0;
	  FunctionBasisDecompositionError[k] = 0.0;
	}
      double* SinTable = new double [NbrParticles];
      double* FunctionBasisDecompositionGrid = new double [NbrOrbitals];
      for (int k = 0; k < NbrOrbitals; ++k)
	FunctionBasisDecompositionGrid[k] = 0.0;
      int* GridLocations = new int [NbrParticles];

      double TotalProbability = 0.0;
      double TotalProbabilityError = 0.0;
      double ThetaStep = M_PI / NbrSteps;
      int CurrentPercent = 0;
      double Theta = 0.0;
      double* Density = new double [NbrSteps];
      for (int i = 0; i < NbrSteps; ++i)
	Density[i] = 0.0;

      RandomUV (TmpUV, TmpPositions, NbrParticles, RandomNumber);
      int NextCoordinate = 0;
      double PreviousProbabilities = 0.0;
      double CurrentProbabilities = 0.0;
      TotalProbability = 0.0;
      int Acceptance = 0;	  
      
      Complex Tmp = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
      CurrentProbabilities = SqrNorm(Tmp);
      for (int k = 0; k < NbrParticles; ++k)
	{
	  SinTable[k] = sin(TmpPositions[k << 1]);
	  CurrentProbabilities *= SinTable[k];
	}
      PreviousProbabilities = CurrentProbabilities;
      
      cout << "starting warm-up sequence" << endl;
      for (int i = 1; i < NbrWarmUpIter; ++i)
	{
	  PreviousTmpU = TmpUV[(NextCoordinate << 1)];
	  PreviousTmpV = TmpUV[(NextCoordinate << 1) + 1];
	  PreviousTmpPositionTheta = TmpPositions[(NextCoordinate << 1)];
	  PreviousTmpPositionPhi = TmpPositions[(NextCoordinate << 1) + 1];
	  RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinate, JumpDistance, RandomNumber);
	  Complex TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	  CurrentProbabilities = SqrNorm(TmpMetropolis);
	  SinTable[NextCoordinate] = sin(TmpPositions[NextCoordinate << 1]);
	  for (int k = 0; k < NbrParticles; ++k)
	    CurrentProbabilities *= SinTable[k];
	  if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	    {
	      PreviousProbabilities = CurrentProbabilities;
	      ++Acceptance;
	    }
	  else
	    {
	      TmpUV[(NextCoordinate << 1)] = PreviousTmpU;
	      TmpUV[(NextCoordinate << 1) + 1] = PreviousTmpV;
	      TmpPositions[(NextCoordinate << 1)] = PreviousTmpPositionTheta;
	      TmpPositions[(NextCoordinate << 1) + 1] = PreviousTmpPositionPhi;
	      SinTable[NextCoordinate] = sin(TmpPositions[NextCoordinate << 1]);
	      CurrentProbabilities = PreviousProbabilities;
	    }
	  NextCoordinate = (int) (RandomNumber->GetRealRandomNumber() * (double) NbrParticles);
	  if (((i * 20) / NbrWarmUpIter) != CurrentPercent)
	    {
	      CurrentPercent = (i * 20) / NbrWarmUpIter;
	      cout << (CurrentPercent * 5) << "% " << flush;
	    }
	} 
      cout << endl << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrWarmUpIter)) * 100.0) << "%" <<endl;
      
      CurrentPercent = 0;
      for (int k = 0; k < NbrParticles; ++k)
	FunctionBasis.GetAllFunctionValues(TmpUV[k << 1], TmpUV[(k << 1) + 1], FunctionBasisEvaluation[k]);
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  double Tmp1 = 0.0;
	  for (int j = 0; j < NbrParticles; ++j)
	    Tmp1 +=  SqrNorm(FunctionBasisEvaluation[j][k]);
	  TmpFunctionBasisDecomposition[k] = Tmp1;
	}
      for (int i = 0; i < NbrParticles; ++i)
	GridLocations[i] = (int) (TmpPositions[i << 1] * InvGridStep);

      cout << "starting MC sequence" << endl;
      Acceptance = 0;

      for (int i = 0; i < NbrIter; ++i)
	{
	  PreviousTmpU = TmpUV[(NextCoordinate << 1)];
	  PreviousTmpV = TmpUV[(NextCoordinate << 1) + 1];
	  PreviousTmpPositionTheta = TmpPositions[(NextCoordinate << 1)];
	  PreviousTmpPositionPhi = TmpPositions[(NextCoordinate << 1) + 1];
	  RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinate, JumpDistance, RandomNumber);
	  Complex TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	  CurrentProbabilities = SqrNorm(TmpMetropolis);
	  SinTable[NextCoordinate] = sin(TmpPositions[NextCoordinate << 1]);
	  for (int k = 0; k < NbrParticles; ++k)
	    CurrentProbabilities *= SinTable[k];
	  if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	    {
	      PreviousProbabilities = CurrentProbabilities;
	      ++Acceptance;	      
	      FunctionBasis.GetAllFunctionValues(TmpUV[NextCoordinate << 1], TmpUV[(NextCoordinate << 1) + 1], FunctionBasisEvaluation[NextCoordinate]);
	      for (int k = 0; k < NbrOrbitals; ++k)
		{
		  double Tmp1 = 0.0;
		  for (int j = 0; j < NbrParticles; ++j)
		    Tmp1 +=  SqrNorm(FunctionBasisEvaluation[j][k]);
		  TmpFunctionBasisDecomposition[k] = Tmp1;
		}
//	      GridLocations[NextCoordinate] = (int) (TmpPositions[NextCoordinate << 1] * InvGridStep);
	    }
	  else
	    {
	      TmpUV[(NextCoordinate << 1)] = PreviousTmpU;
	      TmpUV[(NextCoordinate << 1) + 1] = PreviousTmpV;
	      TmpPositions[(NextCoordinate << 1)] = PreviousTmpPositionTheta;
	      TmpPositions[(NextCoordinate << 1) + 1] = PreviousTmpPositionPhi;
	      SinTable[NextCoordinate] = sin(TmpPositions[NextCoordinate << 1]);
	      CurrentProbabilities = PreviousProbabilities;
	    }
	  if (SymmetryFlag == false)
	    for (int k = 0; k < NbrOrbitals; ++k)
	      {
		FunctionBasisDecomposition[k] += TmpFunctionBasisDecomposition[k];
		FunctionBasisDecompositionError[k] += TmpFunctionBasisDecomposition[k] * TmpFunctionBasisDecomposition[k];
	      }
	  else
	    for (int k = 0; k < NbrOrbitals; ++k)
	      {
		FunctionBasisDecomposition[k] += 0.5 * (TmpFunctionBasisDecomposition[k] + TmpFunctionBasisDecomposition[NbrOrbitals - 1 - k]);
		FunctionBasisDecompositionError[k] += 0.25 * ((TmpFunctionBasisDecomposition[k] + TmpFunctionBasisDecomposition[NbrOrbitals - 1 - k])
							      * (TmpFunctionBasisDecomposition[k] + TmpFunctionBasisDecomposition[NbrOrbitals - 1 - k]));
	      }
// 	  for (int j = 0; j < NbrParticles; ++j)
// 	    FunctionBasisDecompositionGrid[GridLocations[j]] += 1.0 / SinTable[j];
	  TotalProbability++;	  
	  TotalProbabilityError++;

	  NextCoordinate = (int) (RandomNumber->GetRealRandomNumber() * (double) NbrParticles);

	  if ((RecordWaveFunctions != 0) && ((i % RecordStep) == 0l))
	    {
	      ofstream RecordFile;
	      RecordFile.open(RecordWaveFunctions, ios::out | ios::binary | ios::app);
	      RecordFile.precision(14);
	      RecordFile << i << " " << TotalProbability;
	      for (int j = 0; j < NbrOrbitals; ++j)
		RecordFile << " " << FunctionBasisDecomposition[j] << " " << FunctionBasisDecompositionError[j];
	      RecordFile << endl;
	      RecordFile.close();
	    }

	  if (((i * 20) / NbrIter) != CurrentPercent)
	    {
	      CurrentPercent = (i * 20) / NbrIter;
	      cout << (CurrentPercent * 5) << "% " << flush;
	    }
	}
      cout << endl << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrIter)) * 100.0) << "%" <<endl;

      TotalProbabilityError /= (double) NbrIter;
      TotalProbabilityError -= (TotalProbability * TotalProbability) / (((double) NbrIter) * ((double) NbrIter));
      TotalProbabilityError = sqrt(TotalProbabilityError) / (TotalProbability / ((double) NbrIter));
      TotalProbabilityError /= sqrt ((double) NbrIter);
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  FunctionBasisDecompositionError[k] *= ((double) NbrIter);
	  FunctionBasisDecompositionError[k] -= FunctionBasisDecomposition[k] * FunctionBasisDecomposition[k];
	  FunctionBasisDecompositionError[k] = sqrt(FunctionBasisDecompositionError[k] / ((double) NbrIter)) / FunctionBasisDecomposition[k];
	  FunctionBasisDecompositionError[k] += TotalProbabilityError;
	}

      TotalProbability =  1.0 / TotalProbability;
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  FunctionBasisDecomposition[k] *= TotalProbability;
	  FunctionBasisDecompositionError[k] *= FunctionBasisDecomposition[k];
	}
      for (int k = 0; k < NbrOrbitals; ++k)
	FunctionBasisDecompositionGrid[k] *= TotalProbability;

      ofstream DensityRecordFile;
      DensityRecordFile.precision(14);
      DensityRecordFile.open(((SingleStringOption*) Manager["output"])->GetString(), ios::out);

      RealVector TmpPos(2, true);
      Theta = 0.0;
      double Sum = 0.0; 
      double TmpFactor = 4.0 * M_PI / ((double) NbrOrbitals);
      Manager.DisplayOption(DensityRecordFile, true, '#');
      DensityRecordFile << "# decomposition of the wave on " << NbrOrbitals << " orbitals" << endl; 
      DensityRecordFile << "# orbital_index     component    error" << endl;
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  DensityRecordFile << "# " << k << " " << FunctionBasisDecomposition[k] << " " <<  FunctionBasisDecompositionError [k] << endl;
	  cout << k << " " << FunctionBasisDecomposition[k] << " " <<  FunctionBasisDecompositionError [k] << endl;
//	  cout << k << " " << FunctionBasisDecompositionGrid[k] << endl;
	}
      DensityRecordFile << "#" << endl << "# density wave function " << endl << "# theta density" << endl;
      for (int i = 0; i <= NbrSteps; ++i)
	{
	  TmpPos[0] = Theta;
	  double Tmp = 0.0;
	  Complex Tmp2;
	  for (int k = 0; k < NbrOrbitals; ++k)
	    {
	      FunctionBasis.GetFunctionValue(TmpPos, Tmp2, k);
	      Tmp += FunctionBasisDecomposition[k] * SqrNorm(Tmp2);
	    }
	  Tmp *= TmpFactor;
	  Sum += sin(Theta) * Tmp;
	  DensityRecordFile << Theta << " " << Tmp << endl;
	  Theta += ThetaStep;
	}
//       for (int i = 0; i < NbrOrbitals; ++i)
// 	{
// 	  Sum += sin(Theta) *  FunctionBasisDecompositionGrid[i];
// 	  DensityRecordFile << Theta << " " <<  FunctionBasisDecompositionGrid[i] << endl;
// 	  Theta += GridStep;
// 	}
//       Sum *= 2.0 * M_PI;
      DensityRecordFile.close();
      cout << "Sum = " << (Sum * ThetaStep * 2.0 * M_PI) << endl;
    }
   else
     {
       
       ComplexVector UV (NbrParticles * 2, true);
       RealVector TmpPositions (NbrParticles * 2, true);
       RandomUV (UV, TmpPositions, NbrParticles, RandomNumber);
       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;;
       for (int i = 0; i < NbrParticles; ++i)
	 {
	   for (int j = i + 1; j < NbrParticles; ++j)
	     {
	       cout << i << "<->" << j << " : ";
	       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << " | " ;
	       FlipCoordinates(UV, i, j);
	       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;
	       FlipCoordinates(UV, i, j);	       
	     }
	 }
       cout << "Lz <-> -Lz : " << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << " | " ;
       FlipLzMinusLz(UV);
       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;
       FlipLzMinusLz(UV);
       for (int i = 1; i < NbrParticles; ++i)
	 {
	   UV[(i << 1)] = UV[0];
	   UV[(i << 1) + 1] = UV[1];
	   cout << (i  + 1) << " body cancellation : " << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;
	 }
     }
  return 0;
}

void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
      double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
      positions[2 * j] = x;
      positions[(2 * j) + 1] = y;
      uv.Re(2 * j) = cos(0.5 * x);
      uv.Im(2 * j) = uv.Re(2 * j) * sin(0.5 * y);
      uv.Re(2 * j) *= cos(0.5 * y);
      uv.Re(2 * j + 1) = sin(0.5 * x);
      uv.Im(2 * j + 1) = - uv.Re(2 * j + 1) * sin(0.5 * y);
      uv.Re(2 * j + 1) *= cos(0.5 * y);      
    }
}

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, double jumpDistance,
			   AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double x = positions[coordinate];
  double y = positions[coordinate + 1];
  x += jumpDistance *  (1.0 - 2.0 * randomNumberGenerator->GetRealRandomNumber());
  if (x < 0.0) 
    {
      x *= -1.0;
      y += M_PI;
    }
  else
    if (x > M_PI) 
      {
	x = M_2PI - x;
	y += M_PI;
      }
  y += 2.0 * jumpDistance *  (1.0 - 2.0 * randomNumberGenerator->GetRealRandomNumber());
  if (y < 0.0)
    y += M_2PI;
  else
    if (y > M_2PI)
    y -= M_2PI;

  positions[coordinate] = x;
  uv.Re(coordinate) = cos(0.5 * x);
  uv.Im(coordinate) = uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);
  ++coordinate;
  positions[coordinate] = y;
  uv.Re(coordinate) = sin(0.5 * x);
  uv.Im(coordinate) = - uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);      
}


void FlipCoordinates (ComplexVector& uv, int i, int j)
{
  Complex Tmp = uv[2 * i];
  uv.Re(2 * i) = uv.Re(2 * j);
  uv.Re(2 * j) = Tmp.Re;
  uv.Im(2 * i) = uv.Im(2 * j);
  uv.Im(2 * j) = Tmp.Im;
  Tmp = uv[2 * i + 1];
  uv.Re(2 * i + 1) = uv.Re(2 * j + 1);
  uv.Re(2 * j + 1) = Tmp.Re;
  uv.Im(2 * i + 1) = uv.Im(2 * j + 1);
  uv.Im(2 * j + 1) = Tmp.Im;
}

void FlipLzMinusLz (ComplexVector& uv)
{
  for (int i = 0; i < uv.GetVectorDimension(); i += 2)
    {
      Complex Tmp = uv[i + 1];
      uv[i + 1] = uv[i];
      uv[i] = Tmp;
    }
}


