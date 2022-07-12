#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"

#include "Options/Options.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"

#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"

#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"

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
using std::ofstream;


void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);

// evaluate the wave function of a given fermionic Fock state
//
// indices = array of indices describing which orbitals are filled
// nbrFermions = number of fermions
// functions = reference on the matrix where orbitals are evaluted at particle position
// slater = reference on a temporary matrix ued to compute Slater determinant
// return value = evaluated wave function 
Complex EvaluateFermionicComponent (int* indices, int nbrFermions, ComplexMatrix& functions, ComplexMatrix& slater);

// evaluate the wave function of a given bosonic Fock state
//
// indices = array of indices describing which orbitals are filled
// nbrBosons = number of bosons
// functions = reference on the matrix where orbitals are evaluted at particle position
// slater = reference on a temporary matrix ued to compute Slater determinant
// return value = evaluated wave function 
Complex EvaluateBosonicComponent (int* indices, int nbrBosons, double symmetryFactor, ComplexMatrix& functions, 
				  ComplexMatrix& slater, int* changeBit, int* changeBitSign);

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereSUKToU1MCOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-ll", "number of pseudo Landau-level", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
 
  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new BooleanOption  ('r', "resume", "resume from a previous run");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereJainCentralCharge -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int NbrPseudoLandauLevels = Manager.GetInteger("nbr-ll");
  if ((NbrParticles % NbrPseudoLandauLevels) != 0)
    {
      cout << "the number of particles has to be a multiple of the number of pseudo Landau levels" << endl;
      return -1;
    }
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  bool ResumeFlag = Manager.GetBoolean("resume");
  bool StatisticFlag = !(((BooleanOption*) Manager["boson"])->GetBoolean());

  int* HighestWeightDescription = new int [NbrParticles];
  int* GaffnianHighestWeightDescription = new int [NbrParticles];
  int LzMax = 0;
  double SymmetryFactor1 = 1.0;
  double SymmetryFactor2 = 1.0;
  if (StatisticFlag == true)
    {
      LzMax = (5 * (NbrParticles / 2)) - 4;
      for (int i = 0; i < (NbrParticles >> 1); ++i)
	{
	  GaffnianHighestWeightDescription[(i << 1)] = (5 * i);
	  GaffnianHighestWeightDescription[(i << 1) + 1] = (5 * i) + 1;
	}
      HighestWeightDescription[0] = 0;
      HighestWeightDescription[1] = 1;
      for (int i = 1; i < ((NbrParticles >> 1) - 1); ++i)
	{
	  HighestWeightDescription[(i << 1)] = (5 * i) - 1;
	  HighestWeightDescription[(i << 1) + 1] = (5 * i) + 2;
	}
      HighestWeightDescription[NbrParticles - 2] = LzMax - 1;
      HighestWeightDescription[NbrParticles - 1] = LzMax;
     }
  else
    {
      LzMax = (3 * (NbrParticles / 2)) - 3;
      for (int i = 0; i < (NbrParticles >> 1); ++i)
	{
	  GaffnianHighestWeightDescription[(i << 1)] = (3 * i);
	  GaffnianHighestWeightDescription[(i << 1) + 1] = (3 * i);
	}
      SymmetryFactor2 = pow(2.0, 0.25 * ((double) NbrParticles));
      HighestWeightDescription[0] = 0;
      HighestWeightDescription[1] = 0;
      for (int i = 1; i < ((NbrParticles >> 1) - 1); ++i)
	{
	  HighestWeightDescription[(i << 1)] = (3 * i) - 1;
	  HighestWeightDescription[(i << 1) + 1] = (3 * i) + 1;
	}
      HighestWeightDescription[NbrParticles - 2] = LzMax;
      HighestWeightDescription[NbrParticles - 1] = LzMax;
      SymmetryFactor1 = 2.0;
    }
   for (int i = 0; i < NbrParticles; ++i)
     cout << HighestWeightDescription[i] << " ";
   cout << endl;
  for (int i = 0; i < NbrParticles; ++i)
    cout << GaffnianHighestWeightDescription[i] << " ";
  cout << endl;

  ParticleOnSphereFunctionBasis ExactBasis(LzMax);
//   AbstractQHEParticle* ExactSpace = new FermionOnSphere (NbrParticles, 0, LzMax);
//   RealVector ExactState;
  //  ExactState.ReadVector("fermions_gaffnianhighestweight_n_6_2s_11_lz_0.0.vec");
  //   ExactState.ReadVector("bosons_jainhighestweight_n_6_2s_6_lz_0.0.vec");

  Abstract1DComplexFunctionOnSphere* SymmetrizedFunction = 0;
  if (StatisticFlag == true)
    {
      SymmetrizedFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, NbrPseudoLandauLevels, 2);
    }
  else
    {
      SymmetrizedFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, NbrPseudoLandauLevels, 1);
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
      RandomNumber = new StdlibRandomNumberGenerator (29457);
    }

  
  int RecordStep = Manager.GetInteger("record-step");
  
  Complex Tmp;
  Complex Tmp2;
  Complex Component1 = 0.0;
  Complex Component2 = 0.0;
  Complex ErrorComponent1 = 0.0;
  Complex ErrorComponent2 = 0.0;
  double NormalizationCF = 0.0;
  double NormalizationExact1 = 0.0;
  double NormalizationExact2 = 0.0;
  double TmpNormalizationCF = 0.0;
  double TmpNormalizationExact1 = 0.0;
  double TmpNormalizationExact2 = 0.0;
  Complex TmpComponent1;
  Complex TmpComponent2;
  Complex TmpExact;
  int NextCoordinates = 0;
  ComplexVector TmpUV (NbrParticles * 2, true);
  RealVector TmpPositions (NbrParticles * 2, true);
  double PreviousProbabilities = 0.0;
  double CurrentProbabilities = 0.0;
  double TotalProbability = 0.0;
  int Acceptance = 0;
  double AcceptanceRate = 1.0;
  int InitialNbrIter = 0;
// #ifdef __LAPACK__
//   ComplexLapackDeterminant Slater;
// #else
  ComplexMatrix Slater(NbrParticles, NbrParticles);
  //#endif
  ComplexMatrix Functions(LzMax + 1, NbrParticles);
  RealVector TmpCoordinates(2);
  int* ChangeBitSign;
  int* ChangeBit;
  Slater.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign);

  if (ResumeFlag == false)
    {
      RandomUV (TmpUV, TmpPositions, NbrParticles, RandomNumber);
      TmpExact = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
      Tmp = Conj(TmpExact);
      for (int i = 0; i < NbrParticles; ++i)
	{
	  TmpCoordinates[0] = TmpPositions[i << 1];
	  TmpCoordinates[1] = TmpPositions[1 + (i << 1)];
	  for (int j = 0; j <= LzMax; ++j)
	    {
	      ExactBasis.GetFunctionValue(TmpCoordinates, Tmp2, j);
	      Functions[i].Re(j) = Tmp2.Re;
	      Functions[i].Im(j) = Tmp2.Im;
	    }
	}
      if (StatisticFlag == true)
	Tmp *= EvaluateFermionicComponent(HighestWeightDescription, NbrParticles, Functions, Slater);
      else
	Tmp *= EvaluateBosonicComponent(HighestWeightDescription, NbrParticles, SymmetryFactor1, Functions, Slater, 
					ChangeBit, ChangeBitSign);
      //      PreviousProbabilities = Norm(Tmp);
      PreviousProbabilities = fabs(Tmp.Re);
      CurrentProbabilities = PreviousProbabilities;
      TotalProbability = PreviousProbabilities;
      Acceptance = 0;
      AcceptanceRate = 1.0;
      if (NbrWarmUpIter > 0)
	cout << "starting warm-up sequence" << endl;
      for (int i = 1; i < NbrWarmUpIter; ++i)
	{      
	  Complex PreviousCoordinatesU = TmpUV[NextCoordinates << 1];
	  Complex PreviousCoordinatesV = TmpUV[1 + (NextCoordinates << 1)];
	  double PreviousCoordinates1 = TmpPositions[NextCoordinates << 1];
	  double PreviousCoordinates2 = TmpPositions[1 + (NextCoordinates << 1)];
	  Complex TmpMetropolis;
	  RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinates, RandomNumber);
	  TmpExact = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	  TmpMetropolis = Conj(TmpExact);
	  for (int l = 0; l < NbrParticles; ++l)
	    {
	      TmpCoordinates[0] = TmpPositions[l << 1];
	      TmpCoordinates[1] = TmpPositions[1 + (l << 1)];
	      for (int j = 0; j <= LzMax; ++j)
		{
		  ExactBasis.GetFunctionValue(TmpCoordinates, Tmp2, j);
		  Functions[l].Re(j) = Tmp2.Re;
		  Functions[l].Im(j) = Tmp2.Im;
		}
	    }
	  if (StatisticFlag == true)
	    TmpMetropolis *= EvaluateFermionicComponent(HighestWeightDescription, NbrParticles, Functions, Slater);
	  else
	    TmpMetropolis *= EvaluateBosonicComponent(HighestWeightDescription, NbrParticles, SymmetryFactor1, Functions, Slater, 
						      ChangeBit, ChangeBitSign);
	  CurrentProbabilities = fabs(TmpMetropolis.Re);
	  //	  CurrentProbabilities = Norm(TmpMetropolis);
	  if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	    {
	      PreviousProbabilities = CurrentProbabilities;
	      ++Acceptance;
	    }
	  else
	    {
	      TmpUV.Re(NextCoordinates << 1) = PreviousCoordinatesU.Re;
	      TmpUV.Im(NextCoordinates << 1) = PreviousCoordinatesU.Im;
	      TmpUV.Re(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Re;
	      TmpUV.Im(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Im;
	      TmpPositions[NextCoordinates << 1] = PreviousCoordinates1;
	      TmpPositions[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	      CurrentProbabilities = PreviousProbabilities;
	    }
	  NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
	  if ((i % 1000) == 0)
	    {
	      AcceptanceRate = ((double) Acceptance) / ((double) i);
	      cout << Acceptance << " / " << i << " = " << ((100.0 * ((double) Acceptance)) / ((double) i)) << "%" << endl;
	    }
	}
      
      if (NbrWarmUpIter > 0)
	cout << "warm-up sequence is over" << endl;
      Acceptance = 0;
      if ((RecordStep != 0) && (ResumeFlag == false))
	{
	  ofstream OverlapRecordFile;
	  OverlapRecordFile.precision(14);
	  OverlapRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary);
	  OverlapRecordFile << "# Monte Carlo overlap calculation" << endl
			    << "# step overlap.Re overlap.Im error.Re error.Im [(scalar_product.Re scalar_product.Im error_scalar_product.Re error_scalar_product.Im normalization_exact error_normalization_exact) per state] normalization error_normalization" << endl;
	}
    }
  else
    {
      ifstream MCState;
      MCState.open("mcstate.dat", ios::in | ios::binary);
      ReadLittleEndian(MCState, InitialNbrIter);
      unsigned long TmpNumber;
      ReadLittleEndian(MCState, TmpNumber);
      ReadLittleEndian(MCState, Acceptance);
      ReadLittleEndian(MCState, PreviousProbabilities);
      ReadLittleEndian(MCState, CurrentProbabilities);
      ReadLittleEndian(MCState, TotalProbability);
      ReadLittleEndian(MCState, NextCoordinates);
      for (int j = 0; j < (NbrParticles >> 1); ++j)
	{		   	       
	  ReadLittleEndian(MCState, TmpUV.Re(j));
	  ReadLittleEndian(MCState, TmpUV.Im(j));
	}
      ReadLittleEndian(MCState, Component1);
      ReadLittleEndian(MCState, Component2);
      ReadLittleEndian(MCState, ErrorComponent1);
      ReadLittleEndian(MCState, ErrorComponent2);
      MCState.close();	     
    }
  for (int i = 0; i < NbrParticles; ++i)
    {
      TmpCoordinates[0] = TmpPositions[i << 1];
      TmpCoordinates[1] = TmpPositions[1 + (i << 1)];
      for (int j = 0; j <= LzMax; ++j)
	{
	  ExactBasis.GetFunctionValue(TmpCoordinates, Tmp2, j);
	  Functions[i].Re(j) = Tmp2.Re;
	  Functions[i].Im(j) = Tmp2.Im;
	}
    }
  if (StatisticFlag == true)
    {
      TmpComponent1 = EvaluateFermionicComponent(HighestWeightDescription, NbrParticles, Functions, Slater);
      TmpComponent2 = EvaluateFermionicComponent(GaffnianHighestWeightDescription, NbrParticles, Functions, Slater);
    }
  else
    {
      TmpComponent1 = EvaluateBosonicComponent(HighestWeightDescription, NbrParticles, SymmetryFactor1, Functions, Slater, 
					       ChangeBit, ChangeBitSign);
      TmpComponent2 = EvaluateBosonicComponent(GaffnianHighestWeightDescription, NbrParticles, SymmetryFactor2, Functions, Slater, 
					       ChangeBit, ChangeBitSign);
    }
  TmpExact = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
  TmpNormalizationExact1 = SqrNorm(TmpComponent1);
  TmpNormalizationExact2 = SqrNorm(TmpComponent2);
  TmpNormalizationCF = SqrNorm(TmpExact);

  for (int i = InitialNbrIter; i < NbrIter; ++i)
    {
      Complex PreviousCoordinatesU = TmpUV[NextCoordinates << 1];
      Complex PreviousCoordinatesV = TmpUV[1 + (NextCoordinates << 1)];
      double PreviousCoordinates1 = TmpPositions[NextCoordinates << 1];
      double PreviousCoordinates2 = TmpPositions[1 + (NextCoordinates << 1)];
      RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinates, RandomNumber);
      Complex TmpMetropolis;
      Complex TmpExactOld = TmpExact;
      TmpExact = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
      TmpMetropolis = Conj(TmpExact);
      for (int l = 0; l < NbrParticles; ++l)
	{
	  TmpCoordinates[0] = TmpPositions[l << 1];
	  TmpCoordinates[1] = TmpPositions[1 + (l << 1)];
	  for (int j = 0; j <= LzMax; ++j)
	    {
	      ExactBasis.GetFunctionValue(TmpCoordinates, Tmp2, j);
	      Functions[l].Re(j) = Tmp2.Re;
	      Functions[l].Im(j) = Tmp2.Im;
	    }
	}
      if (StatisticFlag == true)
	TmpMetropolis *= EvaluateFermionicComponent(HighestWeightDescription, NbrParticles, Functions, Slater);
      else
	TmpMetropolis *= EvaluateBosonicComponent(HighestWeightDescription, NbrParticles, SymmetryFactor1, Functions, Slater, 
						  ChangeBit, ChangeBitSign);
      CurrentProbabilities = 1.0;
      //     CurrentProbabilities = fabs(TmpMetropolis.Re);
      //      CurrentProbabilities = Norm(TmpMetropolis);
     //      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  Tmp = TmpMetropolis;

	  TmpCoordinates[0] = TmpPositions[NextCoordinates << 1];
	  TmpCoordinates[1] = TmpPositions[1 + (NextCoordinates << 1)];
	  for (int j = 0; j <= LzMax; ++j)
	    {
	      ExactBasis.GetFunctionValue(TmpCoordinates, Tmp2, j);
	      Functions[NextCoordinates].Re(j) = Tmp2.Re;
	      Functions[NextCoordinates].Im(j) = Tmp2.Im;
	    }
	  if (StatisticFlag == true)
	    {
	      TmpComponent1 = EvaluateFermionicComponent(HighestWeightDescription, NbrParticles, Functions, Slater);
	      TmpComponent2 = EvaluateFermionicComponent(GaffnianHighestWeightDescription, NbrParticles, Functions, Slater);
	    }
	  else
	    {
	      TmpComponent1 = EvaluateBosonicComponent(HighestWeightDescription, NbrParticles, SymmetryFactor1, Functions, Slater, 
						       ChangeBit, ChangeBitSign);
	      TmpComponent2 = EvaluateBosonicComponent(GaffnianHighestWeightDescription, NbrParticles, SymmetryFactor2, Functions, Slater, 
						       ChangeBit, ChangeBitSign);
	    }
	  TmpNormalizationExact1 = SqrNorm(TmpComponent1);
	  TmpNormalizationExact2 = SqrNorm(TmpComponent2);
	  TmpNormalizationCF = SqrNorm(TmpExact);

	  ++Acceptance;
	}
//       else
// 	{
// 	  TmpUV.Re(NextCoordinates << 1) = PreviousCoordinatesU.Re;
// 	  TmpUV.Im(NextCoordinates << 1) = PreviousCoordinatesU.Im;
// 	  TmpUV.Re(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Re;
// 	  TmpUV.Im(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Im;
// 	  TmpPositions[NextCoordinates << 1] = PreviousCoordinates1;
// 	  TmpPositions[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
// 	  CurrentProbabilities = PreviousProbabilities;
// 	  TmpExact = TmpExactOld;
// 	}
      TotalProbability += CurrentProbabilities;
      NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrParticles)
	--NextCoordinates;
      
      //      cout << Tmp << " " << TmpComponent1 << " " << TmpComponent2 << " " << Component1 << " " << Component2 << endl;
      Complex Tmp3 = Conj(TmpExact) * TmpComponent1;
      Tmp3 /= CurrentProbabilities;      
      Component1 += Tmp3;
      ErrorComponent1.Re += Tmp3.Re * Tmp3.Re;
      ErrorComponent1.Im += Tmp3.Im * Tmp3.Im;
      Tmp3 = Conj(TmpExact) * TmpComponent2;
      Tmp3 /= CurrentProbabilities;      
      Component2 += Tmp3;
      ErrorComponent2.Re += Tmp3.Re * Tmp3.Re;
      ErrorComponent2.Im += Tmp3.Im * Tmp3.Im;
      NormalizationExact1 += TmpNormalizationExact1 / CurrentProbabilities;
      NormalizationExact2 += TmpNormalizationExact2 / CurrentProbabilities;
      NormalizationCF += TmpNormalizationCF / CurrentProbabilities;

      if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	{
	  ofstream OverlapRecordFile;
	  OverlapRecordFile.precision(14);
	  OverlapRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary | ios::app);
	  OverlapRecordFile << i ;
	  OverlapRecordFile << endl;
	  OverlapRecordFile.close();

	  ofstream MCState;
	  MCState.open("mcstate.dat", ios::out | ios::binary);
	  WriteLittleEndian(MCState, i);
	  unsigned long TmpNumber = RandomNumber->GetNbrGeneratedNumbers();
	  WriteLittleEndian(MCState, TmpNumber);
	  WriteLittleEndian(MCState, Acceptance);
	  WriteLittleEndian(MCState, PreviousProbabilities);
	  WriteLittleEndian(MCState, CurrentProbabilities);
	  WriteLittleEndian(MCState, TotalProbability);
	  WriteLittleEndian(MCState, NextCoordinates);
	  for (int j = 0; j < (NbrParticles >> 1); ++j)
	    {		   	       
	      WriteLittleEndian(MCState, TmpUV.Re(j));
	      WriteLittleEndian(MCState, TmpUV.Im(j));
	    }
	  WriteLittleEndian(MCState, Component1);
	  WriteLittleEndian(MCState, Component2);
	  WriteLittleEndian(MCState, ErrorComponent1);
	  WriteLittleEndian(MCState, ErrorComponent2);
	  MCState.close();	     
	}
      if ((i > 0) && ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0))
	{
	  cout << " i = " << i << endl;
	  Complex Tmp1 = Component1  / ((double) i);
	  Complex Tmp2 = Component2  / ((double) i);
	  Complex TmpError1 (sqrt( ((ErrorComponent1.Re / ((double) i)) - (Component1.Re * Component1.Re)) / ((double) i) ),
			     sqrt( ((ErrorComponent1.Im / ((double) i)) - (Component1.Im * Component1.Im)) / ((double) i) ));
	  Complex TmpError2 (sqrt( ((ErrorComponent2.Re / ((double) i)) - (Component2.Re * Component2.Re)) / ((double) i) ),
			     sqrt( ((ErrorComponent2.Im / ((double) i)) - (Component2.Im * Component2.Im)) / ((double) i) ));
	  //	  cout << Tmp1 << " +/- " << ErrorComponent1  << "   " << Tmp2 <<  " +/- " << ErrorComponent2  << endl;
	  cout << (Component1 / sqrt(NormalizationExact1 * NormalizationCF)) << " " 
	       << (Component2 / sqrt(NormalizationExact2 * NormalizationCF)) << endl;
	  cout << ((Component2 / sqrt(NormalizationExact2 * NormalizationCF)) / (Component1 / sqrt(NormalizationExact1 * NormalizationCF))) << " " << ((Norm(Component2) / sqrt(NormalizationExact2 * NormalizationCF)) / (Norm(Component1) / sqrt(NormalizationExact1 * NormalizationCF))) <<endl;
	  cout << "x = " << (Tmp2 / Tmp1) << " +/- " << TmpError1 << endl;
	  cout << "x.Re = " << (Tmp2.Re / Tmp1.Re) << " +/- " << TmpError1 << endl;
	  cout << "acceptance rate = " << (((double) Acceptance) / (((double) i))) << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    }
  if (((RecordStep != 0) && ((NbrIter % RecordStep) == 0)))
    {
      ofstream OverlapRecordFile;
      OverlapRecordFile.precision(14);
      OverlapRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary | ios::app);
      OverlapRecordFile << NbrIter ;
      OverlapRecordFile << endl;
      OverlapRecordFile.close();
      
      cout << " final results :" << endl;
      Complex Tmp1 = Component1  / ((double) NbrIter);
      Complex Tmp2 = Component2  / ((double) NbrIter);
      Complex TmpError1 (sqrt( ((ErrorComponent1.Re / ((double) NbrIter)) - (Component1.Re * Component1.Re)) / ((double) NbrIter) ),
			 sqrt( ((ErrorComponent1.Im / ((double) NbrIter)) - (Component1.Im * Component1.Im)) / ((double) NbrIter) ));
      Complex TmpError2 (sqrt( ((ErrorComponent2.Re / ((double) NbrIter)) - (Component2.Re * Component2.Re)) / ((double) NbrIter) ),
			 sqrt( ((ErrorComponent2.Im / ((double) NbrIter)) - (Component2.Im * Component2.Im)) / ((double) NbrIter) ));
      cout << (Tmp2 / Tmp1) << " +/- " << TmpError1 << endl;
      cout << "-----------------------------------------------" << endl;
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

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
  double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
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

// evaluate the wave function of a given fermionic Fock state
//
// indices = array of indices describing which orbitals are filled
// nbrFermions = number of fermions
// functions = reference on the matrix where orbitals are evaluted at particle position
// slater = reference on a temporary matrix ued to compute Slater determinant
// return value = evaluated wave function 
  
Complex EvaluateFermionicComponent (int* indices, int nbrFermions, ComplexMatrix& functions, ComplexMatrix& slater)
{
  double Factor = 1.0;
  for (int i = 2; i <= nbrFermions; ++i)
    Factor *= sqrt((double) i);
  Factor = 1.0 / Factor;

  for (int i = 0; i < nbrFermions; ++i)
    {
      ComplexVector& TmpColum2 = functions[i];	  
      for (int j = 0; j < nbrFermions; ++j)
	{
// #ifdef __LAPACK__
// 	  slater.SetMatrixElement(i,j,TmpColum2.Re(indices[j]), TmpColum2.Im(indices[j]));
// #else
	  slater[i].Re(j) = TmpColum2.Re(indices[j]);
	  slater[i].Im(j) = TmpColum2.Im(indices[j]);
	  //#endif
	}
    }
  
  return (slater.Determinant() * Factor);
}

// evaluate the wave function of a given bosonic Fock state
//
// indices = array of indices describing which orbitals are filled
// nbrBosons = number of bosons
// functions = reference on the matrix where orbitals are evaluted at particle position
// slater = reference on a temporary matrix ued to compute Slater determinant
// return value = evaluated wave function 
  
Complex EvaluateBosonicComponent (int* indices, int nbrBosons, double symmetryFactor, ComplexMatrix& functions, 
				  ComplexMatrix& slater, int* changeBit, int* changeBitSign)
{
  double Factor = 1.0;
  for (int i = 2; i <= nbrBosons; ++i)
    Factor *= sqrt((double) i);
  Factor = 1.0 / (Factor * symmetryFactor);

  for (int i = 0; i < nbrBosons; ++i)
    {
      ComplexVector& TmpColum2 = functions[i];	  
      for (int j = 0; j < nbrBosons; ++j)
	{
	  slater[i].Re(j) = TmpColum2.Re(indices[j]);
	  slater[i].Im(j) = TmpColum2.Im(indices[j]);
	}
    }
  
  return (slater.FastPermanent(changeBit, changeBitSign) * Factor);
}

