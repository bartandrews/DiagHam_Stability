#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/ParticleOnSphereDeltaHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereDeltaModifiedHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereCoulombDeltaHamiltonian.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"

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
  OptionManager Manager ("QHEBosonsOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;
  QHEWaveFunctionManager WaveFunctionManager;

  Manager += SystemGroup;
  WaveFunctionManager.AddOptionGroup(&Manager);
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lz", "twice the momentum projection", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function");

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDeltaOverlap -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((BooleanOption*) Manager["list-wavefunctions"])->GetBoolean() == true)
    {
      WaveFunctionManager.ShowAvalaibleWaveFunctions(cout);
      return 0;
    }

  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz"])->GetInteger();

  if (((SingleStringOption*) Manager["exact-state"])->GetString() == 0)
    {
      cout << "QHEBosonsDeltaOverlap requires an exact state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["exact-state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["exact-state"])->GetString() << endl;
      return -1;      
    }
  if (((SingleStringOption*) Manager["use-exact"])->GetString() != 0)
    {
      RealVector TestState;
      if (TestState.ReadVector (((SingleStringOption*) Manager["use-exact"])->GetString()) == false)
	{
	  cout << "can't open vector file " << ((SingleStringOption*) Manager["use-exact"])->GetString() << endl;
	  return -1;      
	}
      if (State.GetVectorDimension() != TestState.GetVectorDimension())
	{
	  cout << "dimension mismatch" << endl;
	  return -1;      
	}
      cout << "overlap = " << (TestState * State) << endl;
      return 0;
    }

  Abstract1DComplexFunction* TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  if (TestWaveFunction == 0)
    {
      cout << "no or unknown analytical wave function" << endl;
      return -1;
    }
  BosonOnSphere Space (NbrBosons, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);
  
  RealVector Location(2 * NbrBosons, true);

  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);

  
  int RecordStep = Manager.GetInteger("record-step");

  Complex* RecordedOverlap = 0;
  Complex* RecordedOverlapError = 0;
  if (RecordStep != 0)
    {
      RecordedOverlap = new Complex [(NbrIter / RecordStep) + 1];
      RecordedOverlapError = new Complex [(NbrIter / RecordStep) + 1];
    }
  int RecordIndex = 0;
  double Factor = 1.0;
  for (int j = 0; j < NbrBosons; ++j)
    {
      Factor *= 4.0 * M_PI;
    }
  Complex Overlap;
  Complex ErrorOverlap;
  double Normalization = 0.0;
  double ErrorNormalization = 0.0;
  double NormalizationExact = 0.0;
  double ErrorNormalizationExact = 0.0;
  Complex Tmp;
  Complex Tmp3;
  double Tmp2;
  double Tmp2bis;
  int NextCoordinates = 0;
  for (int j = 0; j < NbrBosons; ++j)
    {
      Location[j << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));
      Location[1 + (j << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
    }
  Tmp = (*TestWaveFunction)(Location);
  double PreviousProbabilities = Norm(Tmp);
  double CurrentProbabilities = PreviousProbabilities;
  double TotalProbability = PreviousProbabilities;
  for (int i = 0; i < NbrIter; ++i)
    {
      double PreviousCoordinates1 = Location[NextCoordinates << 1];
      double PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      Complex TmpMetropolis = (*TestWaveFunction)(Location);
//      Complex TmpMetropolis2 = (*TestWaveFunction2)(Location);
//      cout << TmpMetropolis << " " << (8 * TmpMetropolis2) << endl;
      CurrentProbabilities = Norm(TmpMetropolis);
      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  Tmp = TmpMetropolis;
	}
      else
	{
	  Location[NextCoordinates << 1] = PreviousCoordinates1;
	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	  CurrentProbabilities = PreviousProbabilities;
	}
      TotalProbability += CurrentProbabilities;
      NextCoordinates = (int) (((double) NbrBosons) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrBosons)
	--NextCoordinates;

      int TimeCoherence = NextCoordinates;
      if (((BooleanOption*) Manager["with-timecoherence"])->GetBoolean() == false)
	TimeCoherence = -1;
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis, TimeCoherence);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      Complex ValueExact (Operation.GetScalar());
      Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
      Tmp2bis = (ValueExact.Re * ValueExact.Re) + (ValueExact.Im * ValueExact.Im);
      Tmp3 = (Conj(Tmp) * ValueExact);
      Tmp2 /= CurrentProbabilities;
      Tmp3 /= CurrentProbabilities;      
      Tmp2bis /= CurrentProbabilities;  
      Overlap += Tmp3;
      ErrorOverlap.Re += Tmp3.Re * Tmp3.Re;
      ErrorOverlap.Im += Tmp3.Im * Tmp3.Im;
      Normalization += Tmp2;
      ErrorNormalization += Tmp2 * Tmp2;
      NormalizationExact += Tmp2bis;
      ErrorNormalizationExact += Tmp2bis * Tmp2bis;
      if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	{
	  Complex Tmp4 = Overlap / ((double) i);
	  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	  double Tmp6 = Normalization  / ((double) i);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  double Tmp8 = NormalizationExact  / ((double) i);
	  double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
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
	  RecordedOverlap[RecordIndex] = Tmp4;
	  RecordedOverlapError[RecordIndex] = Tmp5;
	  ++RecordIndex;
	}
      if ((i > 0) && ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0))
	{
	  cout << " i = " << i << endl;
	  Complex Tmp4 = Overlap / ((double) i);
	  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	  double Tmp6 = Normalization  / ((double) i);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  double Tmp8 = NormalizationExact  / ((double) i);
	  double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  

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
	  cout << "-----------------------------------------------" << endl;
	}
    } 
  cout << " final results :" << endl;
  Complex Tmp4 = Overlap / ((double) NbrIter);
  cout << Tmp4;
  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
		sqrt( ((ErrorOverlap.Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
  cout << " +/- " << Tmp5 << endl;
  double Tmp6 = Normalization  / ((double) NbrIter);
  cout << Tmp6;
  double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
  cout << " +/- " << Tmp7 << endl;	  
  double Tmp8 = NormalizationExact  / ((double) NbrIter);
  cout << Tmp8;
  double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
  cout << " +/- " << Tmp9 << endl;	  
  
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
  cout << "-----------------------------------------------" << endl;


  if (RecordStep != 0)
    {
      RecordedOverlap[RecordIndex] = Tmp4;
      RecordedOverlapError[RecordIndex] = Tmp5;
      ofstream OverlapRecordFile;
      OverlapRecordFile.precision(14);
      OverlapRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary);
      int NbrRecords = NbrIter / RecordStep;
      OverlapRecordFile << "# Monte Carlo overlap calculation" << endl
		       << "# step overlap.Re overlap.Im overlap2 error.Re error.Im error2" << endl;
      for (int i = 0; i < NbrRecords; ++i)
	OverlapRecordFile << i << " " << RecordedOverlap[i].Re << " " << RecordedOverlap[i].Im << " " << SqrNorm(RecordedOverlap[i])
			  << " " << RecordedOverlapError[i].Re << " " <<  RecordedOverlapError[i].Im << " " 
			  << (2 * ((RecordedOverlap[i].Re * RecordedOverlapError[i].Re) + (RecordedOverlap[i].Im * RecordedOverlapError[i].Im))) << endl;
      OverlapRecordFile.close();
    }

  return 0;
}


