#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

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
using std::ofstream;


void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);


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
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "k index of the SU(k) symmetry group", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "intra-corr", "power of the intra-color correlations", 3);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "inter-corr", "power of the inter-color correlations", 2);  
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "test-symmetry", "check the test wave function is symmetric/antisymmetric");  
  (*SystemGroup) += new BooleanOption ('\n', "reverse-flux", "use reverse flux attachment for the composite fermions");  
  (*SystemGroup) += new BooleanOption ('\n', "jain-cf", "use composite fermion state instead of the symetrized state");  
  (*SystemGroup) += new BooleanOption ('\n', "cf-symmetrized", "compute ovelaps between exact states and both the symetrized state and Jain CF state (norm of the CF state is used as probalbility density)");  
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function");  
  (*SystemGroup) += new SingleStringOption  ('\n', "use-set", "file name of a file describing a set of exact states that have to be used as test wave functions");  
  (*SystemGroup) += new SingleStringOption  ('\n', "load-permutations", "read all the permutations needed to compute the reference wave function from a file");  
  (*SystemGroup) += new SingleStringOption  ('\n', "save-permutations", "file name where all the permutations needed to compute the reference wave function have to be stored");  
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
 
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
  (*MonteCarloGroup) += new BooleanOption  ('\n', "weight-symmetrized" , "use the norm of the symmetrized wave fonction as probalbility density instead of the exact state");
  (*MonteCarloGroup) += new  SingleStringOption ('\n', "record-wavefunctions", "optional file where each wavefunctions will be tabulated and recorded");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSUKToU1MCOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int KValue = Manager.GetInteger("k-value");
  if ((NbrParticles % KValue) != 0)
    {
      cout << "the number of particles has to be a multiple of the number of colors (i.e. k)" << endl;
      return -1;
    }
  int NbrParticlePerColor = NbrParticles / KValue;
  int IntraCorrelation = Manager.GetInteger("intra-corr");
  int InterCorrelation = Manager.GetInteger("inter-corr");
  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool OverlapFlag = true;
  if (((BooleanOption*) Manager["test-symmetry"])->GetBoolean() == true)
    {
      OverlapFlag = false;
    }
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  bool InvertFlag = Manager.GetBoolean("reverse-flux");
  bool ResumeFlag = Manager.GetBoolean("resume");
  int LzMax = NbrParticlePerColor * (((KValue - 1) * InterCorrelation) + IntraCorrelation) - IntraCorrelation - KValue + 1;
  bool UseExactFlag = false;
  bool StatisticFlag = !(((BooleanOption*) Manager["boson"])->GetBoolean());
  bool UseBaseAsWeightFlag = ((BooleanOption*) Manager["weight-symmetrized"])->GetBoolean();
  bool JainAndSymmetrizedFlag = ((BooleanOption*) Manager["cf-symmetrized"])->GetBoolean();
  bool JainFlag =  ((BooleanOption*) Manager["jain-cf"])->GetBoolean();
  if (JainAndSymmetrizedFlag == true)
    {
      UseBaseAsWeightFlag = true;
      JainFlag = false;
    }
  char* RecordWaveFunctions = ((SingleStringOption*) Manager["record-wavefunctions"])->GetString();
  if (RecordWaveFunctions != 0)
    {
      ofstream RecordFile;
      RecordFile.open(RecordWaveFunctions, ios::out | ios::binary);
      RecordFile.close();
    }

  AbstractQHEParticle* ExactSpace = 0;
  RealVector* ExactState = 0;
  AbstractFunctionBasis* ExactBasis = 0;
  int NbrExactStates = 1;
  if (((((SingleStringOption*) Manager["use-exact"])->GetString() != 0) || (((SingleStringOption*) Manager["use-set"])->GetString() != 0))
      && (OverlapFlag == true))
    {
      UseExactFlag = true;
      if (StatisticFlag == true)
	{
	  if (HaldaneBasisFlag == false)
	    ExactSpace = new FermionOnSphere (NbrParticles, 0, LzMax);
	  else
	    {
	      int* ReferenceState = 0;
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
		  return -1;     
		}
	      if (MaxNbrLz != (LzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;     
		}
	      int TotalLz = 0;
	      ExactSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState);
	    }
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    {
#ifdef  __64_BITS__
	  if ((LzMax + NbrParticles - 1) < 63)
#else
	    if ((LzMax + NbrParticles - 1) < 31)	
#endif
	      ExactSpace = new BosonOnSphereShort (NbrParticles, 0, LzMax);
	    else
	      ExactSpace = new BosonOnSphere (NbrParticles, 0, LzMax);
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
		  return -1;     
		}
	      if (MaxNbrLz != (LzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;     
		}
	      int TotalLz = 0;
#ifdef  __64_BITS__
	      if ((LzMax + NbrParticles - 1) < 63)
#else
		if ((LzMax + NbrParticles - 1) < 31)	
#endif
		  ExactSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);
	    }
	}
      if (((SingleStringOption*) Manager["use-set"])->GetString() != 0)
	{
	  ConfigurationParser ExactSetDefinition;
	  if (ExactSetDefinition.Parse(((SingleStringOption*) Manager["use-set"])->GetString()) == false)
	    {
	      ExactSetDefinition.DumpErrors(cout) << endl;
	      return -1;
	    }
	  char** ExactStateFilenames;
	  if (ExactSetDefinition.GetAsStringArray("Vectors", ' ', ExactStateFilenames, NbrExactStates) == false)
	    {
	      cout << "error while parsing Vectors in " << ((SingleStringOption*) Manager["use-set"])->GetString() << endl;
	      return -1;     
	    }
	  ExactState = new RealVector [NbrExactStates];
	  for (int i = 0; i < NbrExactStates; ++i)
	    {
	      if (ExactState[i].ReadVector (ExactStateFilenames[i]) == false)
		{
		  cout << "can't open vector file " << ExactStateFilenames[i] << endl;
		  return -1;      
		}	      
	    }
	}
      else
	{
	  NbrExactStates = 1;
	  ExactState = new RealVector [1];
	  if (ExactState[0].ReadVector (((SingleStringOption*) Manager["use-exact"])->GetString()) == false)
	    {
	      cout << "can't open vector file " << ((SingleStringOption*) Manager["use-exact"])->GetString() << endl;
	      return -1;      
	    }
	}
      if (ExactSpace->GetHilbertSpaceDimension() != ExactState[0].GetVectorDimension())
	{
	  cout << "dimension mismatch : hilbert space = " << ExactSpace->GetHilbertSpaceDimension() << ", exact state = " << ExactState[0].GetVectorDimension() << endl;
	  return -1;
	}
      ExactBasis = new ParticleOnSphereFunctionBasis (LzMax);	
    }

  Abstract1DComplexFunctionOnSphere* BaseFunction = 0;
  switch (KValue)
    {
    case 2:
      {
	BaseFunction = new FQHESU2HalperinPermanentOnSphereWaveFunction (NbrParticlePerColor, NbrParticlePerColor, 
									 IntraCorrelation - 1, IntraCorrelation - 1, InterCorrelation - 1, InvertFlag);
      }
      break;
    case 3:
      {
	BaseFunction = new FQHESU3HalperinPermanentOnSphereWaveFunction(NbrParticlePerColor, NbrParticlePerColor, NbrParticlePerColor, 
									IntraCorrelation - 1, IntraCorrelation - 1, IntraCorrelation - 1, 
									InterCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1, InvertFlag);
      }
      break;
    case 4:
      {
	BaseFunction = new FQHESU4HalperinPermanentOnSphereWaveFunction(NbrParticlePerColor, NbrParticlePerColor, NbrParticlePerColor, 
									NbrParticlePerColor,
									IntraCorrelation - 1, IntraCorrelation - 1, IntraCorrelation - 1,
                                                                        IntraCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1,
									InterCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1,
									InterCorrelation - 1, InvertFlag);
      }
      break;
    default:
      {
	cout << "invalid or unsupported number of colors (i.e. k)" << endl;
	return -1;
      }
    }
  Abstract1DComplexFunctionOnSphere* SymmetrizedFunction = 0;
  Abstract1DComplexFunctionOnSphere* SymmetrizedFunction2 = 0;
  if (JainFlag == true)
    {
      SymmetrizedFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, KValue, InterCorrelation);
    }
  else
    {
      if (((SingleStringOption*) Manager["load-permutations"])->GetString() == 0)
	SymmetrizedFunction = new FQHESphereSymmetrizedSUKToU1WaveFunction (NbrParticles, KValue, BaseFunction, true);      
      else
	SymmetrizedFunction = new FQHESphereSymmetrizedSUKToU1WaveFunction (((SingleStringOption*) Manager["load-permutations"])->GetString(), BaseFunction, true);  
      if (((SingleStringOption*) Manager["save-permutations"])->GetString() != 0)
	{
	  ((FQHESphereSymmetrizedSUKToU1WaveFunction*) SymmetrizedFunction)->WritePermutations(((SingleStringOption*) Manager["save-permutations"])->GetString());
	  return 0;
	}
      if (JainAndSymmetrizedFlag == true)
	{
	  SymmetrizedFunction2 = SymmetrizedFunction;
	  SymmetrizedFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, KValue, InterCorrelation);
	}
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
      RandomNumber = new StdlibRandomNumberGenerator (29457);
    }

  if ((UseExactFlag == true) ||  (OverlapFlag == false))
    {
      TestFunction = 0;
    }
  else
    {
      if (InvertFlag == false)
	{
	  TestFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, KValue, 2);
//	  TestFunction = new FQHESU3GeneralizedGaffnianOnSphereWaveFunction(NbrParticles, 2, 1);
	}
      else
	{
	  TestFunction = new HundRuleCFStates (NbrParticles,  - (NbrParticles / 2) + 2, 1);
	}
    }

  StdlibRandomNumberGenerator RandomNumberGenerator(29457);
  
  if (OverlapFlag == true)
    {
       int RecordStep = Manager.GetInteger("record-step");
       
       double Normalization = 0.0;
       double ErrorNormalization = 0.0;
       double Normalization2 = 0.0;
       double ErrorNormalization2 = 0.0;
       Complex Tmp;
       Complex Value2;
       double Tmp2;
       Complex* Overlap = new Complex[NbrExactStates];
       Complex* Overlap2 = new Complex[NbrExactStates];
       Complex** ErrorOverlap = new Complex*[NbrExactStates];
       Complex** ErrorOverlap2 = new Complex*[NbrExactStates];
       for (int i = 0; i < NbrExactStates; ++i)
	 {
	   ErrorOverlap[i] = new  Complex[NbrExactStates];
	   ErrorOverlap2[i] = new  Complex[NbrExactStates];
	 }
       double* NormalizationExact = new double[NbrExactStates];
       double* ErrorNormalizationExact = new double[NbrExactStates];
       Complex* Tmp3 = new Complex[NbrExactStates];
       double* Tmp2bis = new double[NbrExactStates];
       int NextCoordinates = 0;
       ComplexVector TmpUV (NbrParticles * 2, true);
       RealVector TmpPositions (NbrParticles * 2, true);
       Complex* ValueExact = new Complex[NbrExactStates];
       double PreviousProbabilities = 0.0;
       double CurrentProbabilities = 0.0;
       double TotalProbability = 0.0;
       int Acceptance = 0;
       double AcceptanceRate = 1.0;
       int InitialNbrIter = 0;
       if (ResumeFlag == false)
	 {
	   for (int j = 0; j < NbrExactStates; ++j)
	     {
	       Overlap[j] = 0.0;
	       for (int k = 0; k <= j; ++k)
		 {
		   ErrorOverlap[j][k] = 0.0;
		   ErrorOverlap2[j][k] = 0.0;
		 }
	       Overlap2[j] = 0.0;
	       NormalizationExact[j] = 0.0;	   
	       ErrorNormalizationExact[j] = 0.0;	 
	       Tmp3[j] = 0.0;	   
	       Tmp2bis[j] = 0.0;	 
	     }
	   RandomUV (TmpUV, TmpPositions, NbrParticles, RandomNumber);
	   if (UseExactFlag == false)
	     Tmp = TestFunction->CalculateFromSpinorVariables(TmpUV);
	   else
	     {
	       if (UseBaseAsWeightFlag == false)
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   Tmp = Operation.GetScalar();
		 }
	       else
		 Tmp = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	     }
	   if (UseBaseAsWeightFlag == false)
	     ValueExact[0] = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	   else
	     {
	       if (NbrExactStates == 1)
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   ValueExact[0] = Operation.GetScalar();
		 }
	       else
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, NbrExactStates, &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   for (int j = 0; j < NbrExactStates; ++j)
		     {
		       ValueExact[j] = Operation.GetScalar(j);
		     }
		 }
	       if (JainAndSymmetrizedFlag == true)
		 Value2 = SymmetrizedFunction2->CalculateFromSpinorVariables(TmpUV);
	     }
	   PreviousProbabilities = Norm(Tmp);
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
	       if (UseExactFlag == false)
		 TmpMetropolis = TestFunction->CalculateFromSpinorVariables(TmpUV);
	       else
		 {
		   if (UseBaseAsWeightFlag == false)
		     {
		       QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
		       Operation.ApplyOperation(Architecture.GetArchitecture());      
		       TmpMetropolis = Operation.GetScalar();
		     }
		   else
		     TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
		 }
	       CurrentProbabilities = Norm(TmpMetropolis);
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
	   ReadLittleEndian(MCState, Normalization);
	   ReadLittleEndian(MCState, ErrorNormalization);
	   for (int j = 0; j < NbrExactStates; ++j)
	     {		   	       
	       ReadLittleEndian(MCState, Overlap[j].Re);
	       ReadLittleEndian(MCState, Overlap[j].Im);
	       for (int k = 0; k <=j; ++k)
		 {
		   ReadLittleEndian(MCState, ErrorOverlap[j][k].Re);
		   ReadLittleEndian(MCState, ErrorOverlap[j][k].Im);
		 }
	       ReadLittleEndian(MCState, NormalizationExact[j]);
	       ReadLittleEndian(MCState, ErrorNormalizationExact[j]);
	     }
	   if (JainAndSymmetrizedFlag == true)
	     {
	       ReadLittleEndian(MCState, Normalization2);
	       ReadLittleEndian(MCState, ErrorNormalization2);
	       for (int j = 0; j < NbrExactStates; ++j)
		 {		   	       
		   ReadLittleEndian(MCState, Overlap2[j].Re);
		   ReadLittleEndian(MCState, Overlap2[j].Im);
		   for (int k = 0; k <=j; ++k)
		     {
		       ReadLittleEndian(MCState, ErrorOverlap2[j][k].Re);
		       ReadLittleEndian(MCState, ErrorOverlap2[j][k].Im);
		     }
		 }
	     }
	   MCState.close();	     
	 }
       for (int i = InitialNbrIter; i < NbrIter; ++i)
	 {
	   Complex PreviousCoordinatesU = TmpUV[NextCoordinates << 1];
	   Complex PreviousCoordinatesV = TmpUV[1 + (NextCoordinates << 1)];
	   double PreviousCoordinates1 = TmpPositions[NextCoordinates << 1];
	   double PreviousCoordinates2 = TmpPositions[1 + (NextCoordinates << 1)];
	   RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinates, RandomNumber);
	   Complex TmpMetropolis;
	   if (UseExactFlag == false)
	     TmpMetropolis = TestFunction->CalculateFromSpinorVariables(TmpUV);
	   else
	     {
	       if (UseBaseAsWeightFlag == false)
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   TmpMetropolis = Operation.GetScalar();
		 }
	       else
		 TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	     }
	   CurrentProbabilities = Norm(TmpMetropolis);
	   if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	     {
	       PreviousProbabilities = CurrentProbabilities;
	       Tmp = TmpMetropolis;
	       if (UseBaseAsWeightFlag == false)
		 ValueExact[0] = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	       else
		 {
		   if (NbrExactStates == 1)
		     {
		       QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
		       Operation.ApplyOperation(Architecture.GetArchitecture());      
		       ValueExact[0] = Operation.GetScalar();
		     }
		   else
		     {
		       QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, NbrExactStates, &TmpPositions, ExactBasis);
		       Operation.ApplyOperation(Architecture.GetArchitecture());      
		       for (int j = 0; j < NbrExactStates; ++j)
			 {
			   ValueExact[j] = Operation.GetScalar(j);
			 }
		     }
		   if (JainAndSymmetrizedFlag == true)
		     Value2 = SymmetrizedFunction2->CalculateFromSpinorVariables(TmpUV);
		 }
	       if (RecordWaveFunctions != 0)
		 {
		   ofstream RecordFile;
		   RecordFile.precision(14);
		   RecordFile.open(RecordWaveFunctions, ios::out | ios::binary | ios::app);
		   for (int j = 0; j < (NbrParticles << 1); ++j)
		     RecordFile << TmpUV[j] << "|a";
		   RecordFile << TmpMetropolis;
		   if (UseBaseAsWeightFlag == false)
		     RecordFile << "|" << ValueExact[0];
		   else
		     for (int j = 0; j < NbrExactStates; ++j)
		       RecordFile << "|c" << ValueExact[j];
		   RecordFile << endl;
		   RecordFile.close();
		 }
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
	   TotalProbability += CurrentProbabilities;
	   NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
	   if (NextCoordinates == NbrParticles)
	     --NextCoordinates;

	   Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
	   Tmp2 /= CurrentProbabilities;
	   for (int j = 0; j < NbrExactStates; ++j)
	     {
	       Tmp2bis[j] = (ValueExact[j].Re * ValueExact[j].Re) + (ValueExact[j].Im * ValueExact[j].Im);
	       Tmp3[j] = (Conj(Tmp) * ValueExact[j]);
	       Tmp3[j] /= CurrentProbabilities;      
	       Tmp2bis[j] /= CurrentProbabilities;  
	       Overlap[j] += Tmp3[j];
	       for (int k = 0; k <= j; ++k)
		 {
		   ErrorOverlap[j][k].Re += Tmp3[j].Re * Tmp3[k].Re;
		   ErrorOverlap[j][k].Im += Tmp3[j].Im * Tmp3[k].Im;
		 }
	       NormalizationExact[j] += Tmp2bis[j];
	       ErrorNormalizationExact[j] += Tmp2bis[j] * Tmp2bis[j];
	     }
	   Normalization += Tmp2;
	   ErrorNormalization += Tmp2 * Tmp2;

	   if (JainAndSymmetrizedFlag == true)
	     {
	       Tmp2 = (Value2.Re * Value2.Re) + (Value2.Im * Value2.Im);
	       Tmp2 /= CurrentProbabilities;
	       for (int j = 0; j < NbrExactStates; ++j)
		 {
		   Tmp2bis[j] = (ValueExact[j].Re * ValueExact[j].Re) + (ValueExact[j].Im * ValueExact[j].Im);
		   Tmp3[j] = (Conj(Value2) * ValueExact[j]);
		   Tmp3[j] /= CurrentProbabilities;      
		   Tmp2bis[j] /= CurrentProbabilities;  
		   Overlap2[j] += Tmp3[j];
		   for (int k = 0; k <= j; ++k)
		     {
		       ErrorOverlap2[j][k].Re += Tmp3[j].Re * Tmp3[k].Re;
		       ErrorOverlap2[j][k].Im += Tmp3[j].Im * Tmp3[k].Im;
		     }
		 }
	       Normalization2 += Tmp2;
	       ErrorNormalization2 += Tmp2 * Tmp2;	   
	     }

	   if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	     {
	       ofstream OverlapRecordFile;
	       OverlapRecordFile.precision(14);
	       OverlapRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary | ios::app);
	       OverlapRecordFile << i ;
	       double Tmp6 = Normalization  / ((double) i);
	       double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	       for (int j = 0; j < NbrExactStates; ++j)
		 {
		   Complex Tmp4 = Overlap[j] / ((double) i);
		   Complex Tmp5 (sqrt( ((ErrorOverlap[j][j].Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			     sqrt( ((ErrorOverlap[j][j].Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
		   double Tmp8 = NormalizationExact[j]  / ((double) i);
		   double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
		   Tmp5.Re /= Tmp4.Re;
		   Tmp5.Im /= Tmp4.Im;
		   Tmp5.Re = fabs(Tmp5.Re);
		   Tmp5.Im = fabs(Tmp5.Im);
		   Tmp5.Re += (Tmp7 / Tmp6);
		   Tmp5.Im += (Tmp7 / Tmp6);
		   Tmp5.Re += (Tmp9 / Tmp8);
		   Tmp5.Im += (Tmp9 / Tmp8);
		   Tmp4 /= sqrt(Tmp6 * Tmp8);	  
		   Tmp5.Re *= Tmp4.Re;
		   Tmp5.Im *= Tmp4.Im;
		   OverlapRecordFile << " " << Tmp4.Re << " " << Tmp4.Im << " " << Tmp5.Re << " " << Tmp5.Im << " " << Overlap[j].Re << " " << Overlap[j].Im << " ";
		   for (int k = 0; k <= j; ++k)
		     OverlapRecordFile << ErrorOverlap[j][k].Re << " " << ErrorOverlap[j][k].Im << " ";
		   OverlapRecordFile << NormalizationExact[j] << " " << ErrorNormalizationExact[j];
		 }
	       OverlapRecordFile << " " << Normalization << " " << ErrorNormalization;
 	       if (JainAndSymmetrizedFlag == true)
		 {
		   Tmp6 = Normalization2 / ((double) i);
		   Tmp7 = sqrt( ((ErrorNormalization2 / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
		   for (int j = 0; j < NbrExactStates; ++j)
		     {
		       Complex Tmp4 = Overlap2[j] / ((double) i);
		       Complex Tmp5 (sqrt( ((ErrorOverlap2[j][j].Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
				     sqrt( ((ErrorOverlap2[j][j].Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
		       double Tmp8 = NormalizationExact[j]  / ((double) i);
		       double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
		       Tmp5.Re /= Tmp4.Re;
		       Tmp5.Im /= Tmp4.Im;
		       Tmp5.Re = fabs(Tmp5.Re);
		       Tmp5.Im = fabs(Tmp5.Im);
		       Tmp5.Re += (Tmp7 / Tmp6);
		       Tmp5.Im += (Tmp7 / Tmp6);
		       Tmp5.Re += (Tmp9 / Tmp8);
		       Tmp5.Im += (Tmp9 / Tmp8);
		       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
		       Tmp5.Re *= Tmp4.Re;
		       Tmp5.Im *= Tmp4.Im;
		       OverlapRecordFile << " " << Tmp4.Re << " " << Tmp4.Im << " " << Tmp5.Re << " " << Tmp5.Im << " " << Overlap2[j].Re << " " << Overlap2[j].Im << " ";
		       for (int k = 0; k <= j; ++k)
			 OverlapRecordFile << ErrorOverlap2[j][k].Re << " " << ErrorOverlap2[j][k].Im << " ";
		     }
		   OverlapRecordFile << Normalization2 << " " << ErrorNormalization2;
 		 }
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
	       WriteLittleEndian(MCState, Normalization);
	       WriteLittleEndian(MCState, ErrorNormalization);
	       for (int j = 0; j < NbrExactStates; ++j)
		 {		   	       
		   WriteLittleEndian(MCState, Overlap[j].Re);
		   WriteLittleEndian(MCState, Overlap[j].Im);
		   for (int k = 0; k <=j; ++k)
		     {
		       WriteLittleEndian(MCState, ErrorOverlap[j][k].Re);
		       WriteLittleEndian(MCState, ErrorOverlap[j][k].Im);
		     }
		   WriteLittleEndian(MCState, NormalizationExact[j]);
		   WriteLittleEndian(MCState, ErrorNormalizationExact[j]);
		 }
	       if (JainAndSymmetrizedFlag == true)
		 {
		   WriteLittleEndian(MCState, Normalization2);
		   WriteLittleEndian(MCState, ErrorNormalization2);
		   for (int j = 0; j < NbrExactStates; ++j)
		     {		   	       
		       WriteLittleEndian(MCState, Overlap2[j].Re);
		       WriteLittleEndian(MCState, Overlap2[j].Im);
		       for (int k = 0; k <=j; ++k)
			 {
			   WriteLittleEndian(MCState, ErrorOverlap2[j][k].Re);
			   WriteLittleEndian(MCState, ErrorOverlap2[j][k].Im);
			 }
		     }
		 }
	       MCState.close();	     
	     }
	   if ((i > 0) && ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0))
	     {
	       cout << " i = " << i << endl;
	       if (JainAndSymmetrizedFlag == true)
		 cout << "  overlap with the Jain CF state :" << endl;
	       double Tmp6 = Normalization  / ((double) i);
	       double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	       for (int j = 0; j < NbrExactStates; ++j)
		 {		   
		   Complex Tmp4 = Overlap[j] / ((double) i);
		   Complex Tmp5 (sqrt( ((ErrorOverlap[j][j].Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
				 sqrt( ((ErrorOverlap[j][j].Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
		   double Tmp8 = NormalizationExact[j]  / ((double) i);
		   double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
		   if (NbrExactStates > 1)
		     cout << "overlap " << j << " : ";
		   if (((BooleanOption*) Manager["show-details"])->GetBoolean() == true)
		     {
		       cout << Tmp4;
		       cout << " +/- " << Tmp5 << endl;
		       cout << Tmp6;
		       cout << " +/- " << Tmp7 << endl;	  
		       cout << Tmp8;
		       cout << " +/- " << Tmp9 << endl;	  
		     }
		   Tmp5.Re /= Tmp4.Re;
		   Tmp5.Im /= Tmp4.Im;
		   Tmp5.Re = fabs(Tmp5.Re);
		   Tmp5.Im = fabs(Tmp5.Im);
		   Tmp5.Re += (Tmp7 / Tmp6);
		   Tmp5.Im += (Tmp7 / Tmp6);
		   Tmp5.Re += (Tmp9 / Tmp8);
		   Tmp5.Im += (Tmp9 / Tmp8);
		   Tmp4 /= sqrt(Tmp6 * Tmp8);	  
		   Tmp5.Re *= Tmp4.Re;
		   Tmp5.Im *= Tmp4.Im;
		   cout << Tmp4 << " +/- " << Tmp5 << endl;
		 }
	       if (JainAndSymmetrizedFlag == true)
		 {
		   cout << "  overlap with the symmetrized state :" << endl;
		   Tmp6 = Normalization2  / ((double) i);
		   Tmp7 = sqrt( ((ErrorNormalization2 / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
		   for (int j = 0; j < NbrExactStates; ++j)
		     {		   
		       Complex Tmp4 = Overlap2[j] / ((double) i);
		       Complex Tmp5 (sqrt( ((ErrorOverlap2[j][j].Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
				     sqrt( ((ErrorOverlap2[j][j].Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
		       double Tmp8 = NormalizationExact[j]  / ((double) i);
		       double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
		       if (NbrExactStates > 1)
			 cout << "overlap " << j << " : ";
		       if (((BooleanOption*) Manager["show-details"])->GetBoolean() == true)
			 {
			   cout << Tmp4;
			   cout << " +/- " << Tmp5 << endl;
			   cout << Tmp6;
			   cout << " +/- " << Tmp7 << endl;	  
			   cout << Tmp8;
			   cout << " +/- " << Tmp9 << endl;	  
			 }
		       Tmp5.Re /= Tmp4.Re;
		       Tmp5.Im /= Tmp4.Im;
		       Tmp5.Re = fabs(Tmp5.Re);
		       Tmp5.Im = fabs(Tmp5.Im);
		       Tmp5.Re += (Tmp7 / Tmp6);
		       Tmp5.Im += (Tmp7 / Tmp6);
		       Tmp5.Re += (Tmp9 / Tmp8);
		       Tmp5.Im += (Tmp9 / Tmp8);
		       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
		       Tmp5.Re *= Tmp4.Re;
		       Tmp5.Im *= Tmp4.Im;
		       cout << Tmp4 << " +/- " << Tmp5 << endl;
		     }
		 }
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
	   double Tmp6 = Normalization  / ((double) NbrIter);
	   double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
	   for (int j = 0; j < NbrExactStates; ++j)
	     {
	       Complex Tmp4 = Overlap[j] / ((double) NbrIter);
	       Complex Tmp5 (sqrt( ((ErrorOverlap[j][j].Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
			     sqrt( ((ErrorOverlap[j][j].Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
	       double Tmp8 = NormalizationExact[j]  / ((double) NbrIter);
	       double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
	       Tmp5.Re /= Tmp4.Re;
	       Tmp5.Im /= Tmp4.Im;
	       Tmp5.Re = fabs(Tmp5.Re);
	       Tmp5.Im = fabs(Tmp5.Im);
	       Tmp5.Re += (Tmp7 / Tmp6);
	       Tmp5.Im += (Tmp7 / Tmp6);
	       Tmp5.Re += (Tmp9 / Tmp8);
	       Tmp5.Im += (Tmp9 / Tmp8);
	       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
	       Tmp5.Re *= Tmp4.Re;
	       Tmp5.Im *= Tmp4.Im;
	       OverlapRecordFile << " " << Tmp4.Re << " " << Tmp4.Im << " " << Tmp5.Re << " " << Tmp5.Im << " " << Overlap[j].Re << " " << Overlap[j].Im << " ";
	       for (int k = 0; k <= j; ++k)
		 OverlapRecordFile << ErrorOverlap[j][k].Re << " " << ErrorOverlap[j][k].Im << " ";
	       OverlapRecordFile << NormalizationExact[j] << " " << ErrorNormalizationExact[j];
	     }
           OverlapRecordFile << " " << Normalization << " " << ErrorNormalization;
	   if (JainAndSymmetrizedFlag == true)
	     {
	       Tmp6 = Normalization2 / ((double) NbrIter);
	       Tmp7 = sqrt( ((ErrorNormalization2 / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
	       for (int j = 0; j < NbrExactStates; ++j)
		 {
		   Complex Tmp4 = Overlap2[j] / ((double) NbrIter);
		   Complex Tmp5 (sqrt( ((ErrorOverlap2[j][j].Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
				 sqrt( ((ErrorOverlap2[j][j].Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
		   double Tmp8 = NormalizationExact[j]  / ((double) NbrIter);
		   double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
		   Tmp5.Re /= Tmp4.Re;
		   Tmp5.Im /= Tmp4.Im;
		   Tmp5.Re = fabs(Tmp5.Re);
		   Tmp5.Im = fabs(Tmp5.Im);
		   Tmp5.Re += (Tmp7 / Tmp6);
		   Tmp5.Im += (Tmp7 / Tmp6);
		   Tmp5.Re += (Tmp9 / Tmp8);
		   Tmp5.Im += (Tmp9 / Tmp8);
		   Tmp4 /= sqrt(Tmp6 * Tmp8);	  
		   Tmp5.Re *= Tmp4.Re;
		   Tmp5.Im *= Tmp4.Im;
		   OverlapRecordFile << " " << Tmp4.Re << " " << Tmp4.Im << " " << Tmp5.Re << " " << Tmp5.Im << " " << Overlap2[j].Re << " " << Overlap2[j].Im << " ";
		   for (int k = 0; k <= j; ++k)
		     OverlapRecordFile << ErrorOverlap2[j][k].Re << " " << ErrorOverlap2[j][k].Im << " ";
		 }
	       OverlapRecordFile << Normalization2 << " " << ErrorNormalization2;
	     }
	   OverlapRecordFile << endl;
	   OverlapRecordFile.close();
	 }
       cout << " final results :" << endl;
       if (JainAndSymmetrizedFlag == true)
	 cout << "  overlap with the Jain CF state :" << endl;
       double Tmp6 = Normalization  / ((double) NbrIter);
       double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
       double TotalNorm = 0.0;
       double TotalNormError = 0.0;
       for (int j = 0; j < NbrExactStates; ++j)
	 {
	   Complex Tmp4 = Overlap[j] / ((double) NbrIter);
	   Complex Tmp5 (sqrt( ((ErrorOverlap[j][j].Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
			 sqrt( ((ErrorOverlap[j][j].Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
	   double Tmp8 = NormalizationExact[j]  / ((double) NbrIter);
	   double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
	   
	   Tmp5.Re /= Tmp4.Re;
	   Tmp5.Im /= Tmp4.Im;
	   Tmp5.Re = fabs(Tmp5.Re);
	   Tmp5.Im = fabs(Tmp5.Im);
	   Tmp5.Re += (Tmp7 / Tmp6);
	   Tmp5.Im += (Tmp7 / Tmp6);
	   Tmp5.Re += (Tmp9 / Tmp8);
	   Tmp5.Im += (Tmp9 / Tmp8);
	   Tmp4 /= sqrt(Tmp6 * Tmp8);	  
	   Tmp5.Re *= Tmp4.Re;
	   Tmp5.Im *= Tmp4.Im;
	   if (NbrExactStates > 1)
	     cout << "overlap " << j << " : ";
	   cout << Tmp4 << " +/- " << Tmp5 << endl;
	   cout << Norm(Tmp4) << " +/- " << ((fabs(Tmp4.Re * Tmp5.Re) + fabs(Tmp4.Im * Tmp5.Im))  / Norm(Tmp4)) << endl;
	   TotalNorm += SqrNorm(Tmp4);
	   TotalNormError += 2.0 * fabs(Tmp4.Re * Tmp5.Re) + fabs(Tmp4.Im * Tmp5.Im);
	 }
       cout << "total overlap = " << TotalNorm << " +/- " << TotalNormError << endl;
       if (JainAndSymmetrizedFlag == true)
	 {
	   cout << "  overlap with the symmetrized state :" << endl;
	   Tmp6 = Normalization2  / ((double) NbrIter);
	   Tmp7 = sqrt( ((ErrorNormalization2 / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
	   TotalNorm = 0.0;
	   for (int j = 0; j < NbrExactStates; ++j)
	     {		   
	       Complex Tmp4 = Overlap2[j] / ((double) NbrIter);
	       Complex Tmp5 (sqrt( ((ErrorOverlap2[j][j].Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
			     sqrt( ((ErrorOverlap2[j][j].Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
	       double Tmp8 = NormalizationExact[j]  / ((double) NbrIter);
	       double Tmp9 = sqrt( ((ErrorNormalizationExact[j] / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
	       if (NbrExactStates > 1)
		 cout << "overlap " << j << " : ";
	       Tmp5.Re /= Tmp4.Re;
	       Tmp5.Im /= Tmp4.Im;
	       Tmp5.Re = fabs(Tmp5.Re);
	       Tmp5.Im = fabs(Tmp5.Im);
	       Tmp5.Re += (Tmp7 / Tmp6);
	       Tmp5.Im += (Tmp7 / Tmp6);
	       Tmp5.Re += (Tmp9 / Tmp8);
	       Tmp5.Im += (Tmp9 / Tmp8);
	       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
	       Tmp5.Re *= Tmp4.Re;
	       Tmp5.Im *= Tmp4.Im;
	       cout << Tmp4 << " +/- " << Tmp5 << endl;
	       cout << Norm(Tmp4) << " +/- " << ((fabs(Tmp4.Re * Tmp5.Re) + fabs(Tmp4.Im * Tmp5.Im))  / Norm(Tmp4)) << endl;
	       TotalNorm += SqrNorm(Tmp4);
	       TotalNormError += 2.0 * fabs(Tmp4.Re * Tmp5.Re) + fabs(Tmp4.Im * Tmp5.Im);
	     }
	   cout << "total overlap = " << TotalNorm << " +/- " << TotalNormError << endl;	   
	 }
      cout << "-----------------------------------------------" << endl;
       
       
     }
   else
     {
//       SymmetrizedFunction = TestFunction;
       
       ComplexVector UV (NbrParticles * 2, true);
       RealVector TmpPositions (NbrParticles * 2, true);
       RandomUV (UV, TmpPositions, NbrParticles, RandomNumber);
       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;;
       
//        UV.Re(2) = UV.Re(0);
//        UV.Im(2) = UV.Im(0);
//        UV.Re(3) = UV.Re(1);
//        UV.Im(3) = UV.Im(1);
//        UV.Re(4) = UV.Re(0);
//        UV.Im(4) = UV.Im(0);
//        UV.Re(5) = UV.Re(1);
//        UV.Im(5) = UV.Im(1);
//        UV.Re(6) = UV.Re(0);
//        UV.Im(6) = UV.Im(0);
//        UV.Re(7) = UV.Re(1);
//        UV.Im(7) = UV.Im(1);
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
