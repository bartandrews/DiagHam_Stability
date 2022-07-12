#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEWaveFunction/PolarizedProductWavefunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/WaveFunctionOverlapOptimizer.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Tools/FQHEWaveFunction/AdvancedReadRezayiOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/ExplicitMooreReadOnSphereWaveFunction.h"

#include "MCObservables/RealObservable.h"
#include "MCObservables/ComplexObservable.h"
#include "MCObservables/WeightedRealObservable.h"
#include "MCObservables/WeightedComplexObservable.h"
#include "MCObservables/MCHistoryRecord.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <cstring>

double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
double dsqrarg1;
#define DSQR1(a) ((dsqrarg1=(a)) == 0.0 ? 0.0 : dsqrarg1*dsqrarg1)
double dsqrarg2;
#define DSQR2(a) ((dsqrarg2=(a)) == 0.0 ? 0.0 : dsqrarg2*dsqrarg2)
double dsqrarg3;
#define DSQR3(a) ((dsqrarg3=(a)) == 0.0 ? 0.0 : dsqrarg3*dsqrarg3)

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

Complex OverlapValue(ComplexObservable &ScalarProduct, RealObservable &Norm);
double OverlapError(ComplexObservable &ScalarProduct, RealObservable &Norm);

Complex OverlapValue(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &Norm1, WeightedRealObservable &Norm2);
double OverlapError(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &Norm1, WeightedRealObservable &Norm2);

void FillRecursiveGrid(int NbrCoordinate, int MaxGrid, RealVector &GridPositions, Abstract1DComplexFunction* TestWaveFunction, MCHistoryRecord *History);


double MinDist (RealVector &Positions)
{
  int NbrFermions = Positions.GetVectorDimension()/2;
  double sp, spmax=-1.0;
  for ( int i=0; i<NbrFermions; ++i)
    for ( int j=i+1; j<NbrFermions; ++j)
      {
	sp=sin(Positions[2*i])*sin(Positions[2*j])*cos(Positions[2*i+1]-Positions[2*j+1]) +
	  cos(Positions[2*i])*cos(Positions[2*j]);
	if (sp>spmax) spmax=sp;
      }
  return acos(spmax);
}

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEFermionsOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;
  QHEWaveFunctionManager WaveFunctionManager;

  Manager += SystemGroup;
  WaveFunctionManager.AddOptionGroup(&Manager);
  PolarizedProductWavefunction::AddPolarizedProductStateOptionGroup(Manager);
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lz", "twice the momentum projection", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption  ('\n', "pfaffian-as-exact", "use Moore-Read Pfaffian instead of exact state");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function");

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('t', "nbr-warmup-iter", "number of steps for thermalization", 500);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "randomSeed", "value of random seed to be used", 29457);
  (*MonteCarloGroup) += new SingleIntegerOption ('H', "history-mode", "use on-file history: (0=off, 1=generate new, 2=read history, 3=optimize with history, 4=continue to generate given history)", 1);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "history-file", "name of the file where overlap recording has to be done", NULL);
  (*MonteCarloGroup) += new BooleanOption ('\n', "varyMR", "vary coefficient of 1/z in pair wavefunction");  
  (*MonteCarloGroup) += new SingleIntegerOption ('d', "sample-density", "spacing of samples to be saved in History-mode", 1);
  (*MonteCarloGroup) += new SingleIntegerOption ('\n', "linearPoints", "number of function evaluations along the gradient in optimising History mode ", 25);
  (*MonteCarloGroup) += new SingleIntegerOption ('\n', "randomPoints", "number of random function evaluations in optimising History mode ", 0);
  (*MonteCarloGroup) += new SingleIntegerOption ('\n', "limitSamples", "maximal number of samples to be used from history-file", 10000000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iterations between two consecutive result recording the overlap value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where overlap recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MiscGroup) += new BooleanOption ('\n', "write-grid", "output wave function on grid points");
  (*MiscGroup) += new BooleanOption ('\n', "write-random", "output wave function on random points");
  (*MiscGroup) += new SingleIntegerOption ('\n', "max-grid", "number of grid points in each direction",10);
  (*MiscGroup) += new SingleIntegerOption  ('\n', "test", "test Monte Carlo observable", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionOverlap -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  if (Manager.GetInteger("test")>0)
    {
      RealObservable testObservable(16,1);
      for (int i=0; i<Manager.GetInteger("test"); ++i)
	testObservable << i;
      cout << "Mean="<<testObservable.Average()<<endl;
      cout << "Variance="<<testObservable.Variance()<<endl;
      cout << "BinVariance="<<testObservable.VarianceOfBins()<<endl;
      cout << "ErrorEstimate="<<testObservable.ErrorEstimate()<<endl;
      exit(0);
    }


  if (Manager.GetBoolean("list-wavefunctions") == true)
    {
      WaveFunctionManager.ShowAvalaibleWaveFunctions(cout);
      return 0;
    }

  bool UsePfaffian = Manager.GetBoolean("pfaffian-as-exact");
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int NbrIter = Manager.GetInteger("nbr-iter");
  int NbrWarmUpIter = Manager.GetInteger("nbr-warmup-iter");
  int Lz = Manager.GetInteger("lz");
  RealVector State;

  if ((Manager.GetInteger("history-mode")<2) && (! UsePfaffian))
    {
      if (Manager.GetString("exact-state") == 0)
	{
	  cout << "QHEFermionOverlap requires an exact state" << endl;
	  return -1;
	}  
      if (State.ReadVector (Manager.GetString("exact-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("exact-state") << endl;
	  return -1;      
	}
      if (Manager.GetString("use-exact") != 0)
	{
	  RealVector TestState;
	  if (TestState.ReadVector (Manager.GetString("use-exact")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("use-exact") << endl;
	      return -1;      
	    }
	  if (State.GetVectorDimension() != TestState.GetVectorDimension())
	    {
	      cout << "dimension mismatch" << endl;
	      return -1;      
	    }
	  Complex ov = (TestState * State);
	  cout << "overlap = " << ov << endl;
	  cout << "Sqrovrl = " << SqrNorm(ov) << endl;
	  return 0;
	}
    }

  if (Manager.GetInteger("history-mode")==2)
    {
      if (Manager.GetString("exact-state") != 0)
	{
	  if (State.ReadVector (Manager.GetString("exact-state")))
	    {
	      cout << "Using exact state to check eventual outliers" << endl;
	    }
	  else
	    {
	      cout << "can't open vector file " << Manager.GetString("exact-state") << endl;
	      return -1;      
	    }
	}      
    }

  Abstract1DComplexFunction* TestWaveFunction = NULL;
  if (Manager.GetBoolean("product-state"))
    TestWaveFunction = new PolarizedProductWavefunction(Architecture.GetArchitecture(), Manager, NbrFermions, LzMax, Lz, &WaveFunctionManager);
  else
    TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  
  
  if (TestWaveFunction == 0)
    {
      cout << "no or unknown analytical wave function: try to add option --test-wavefunction" << endl;
      return -1;
    }

  Abstract1DComplexFunction* PfaffianWaveFunction = NULL;
  if (UsePfaffian)
    {
      double *Coefficients = new double[1];
      Coefficients[0]=0.0;
      PfaffianWaveFunction = new PairedCFOnSphereWaveFunction(NbrFermions, 1, -1, 1.0, Coefficients, false, 2);
      ((PairedCFOnSphereWaveFunction*) PfaffianWaveFunction)->AdaptAverageMCNorm();
      delete [] Coefficients;
      LzMax = 2*NbrFermions-3;
      Lz = 0;
    }

  AbstractRandomNumberGenerator* RandomNumber = 0;
  if (Manager.GetString("random-file") != 0)
    {
      RandomNumber = new FileRandomNumberGenerator(Manager.GetString("random-file"), (unsigned long)((NbrWarmUpIter + NbrIter) * 4.33) + 2000, 
						     Manager.GetInteger("random-seek"));
    }
  else
    {
      RandomNumber = new StdlibRandomNumberGenerator (Manager.GetInteger("randomSeed"));
    }

  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(NbrFermions, RandomNumber);
  Particles->MultiplyStepLength(sqrt((double)LzMax/NbrFermions/3.0));
  FermionOnSphere *Space=NULL;
  Complex ValueExact;
  Complex TrialValue;
  double CurrentSamplingAmplitude;
  
  int HistoryMode = Manager.GetInteger("history-mode");
  int SampleDensity = Manager.GetInteger("sample-density");
  MCHistoryRecord *History=NULL;
  char *HistoryFileName=NULL;
  if (HistoryMode>0)
    {
      HistoryFileName=Manager.GetString("history-file");      
      if (HistoryMode==1)
	{
	  if (HistoryFileName==NULL)
	    {
	      if (UsePfaffian)
		{
		  HistoryFileName = new char[30];
		  sprintf(HistoryFileName,"pfaffian_n%d.samp",NbrFermions);
		}
	      else
		{
		  // default filename: add extension to exact vector
		  HistoryFileName = new char[strlen(Manager.GetString("exact-state"))+6];
		  sprintf(HistoryFileName,"%s.samp",Manager.GetString("exact-state"));
		  char * tmpC = new char[strlen(HistoryFileName)+5];
		  sprintf(tmpC,"%s",HistoryFileName);
		  std::ifstream testExistant(tmpC,std::ios::in);
		  int count=1;
		  while (testExistant.is_open())
		    {
		      testExistant.close();
		      sprintf(tmpC,"%s%d",HistoryFileName,count++);
		      testExistant.open(tmpC,std::ios::in);
		    }
		  delete [] HistoryFileName;
		  HistoryFileName = tmpC;
		}
	    }
	  char *tmpC = WaveFunctionManager.GetDescription();
	  History=new MCHistoryRecord(NbrIter, 2*NbrFermions, Manager.GetString("exact-state"), tmpC, HistoryFileName
				      /* could add additional observables here */);

	  /* testing code: write coordinates on cluster of grid-points */
	  if (Manager.GetBoolean("write-grid"))
	    {
	      RealVector GridPositions(2*NbrFermions,true);

	      cout << "Writing grid"<<endl;
	      FillRecursiveGrid(2*NbrFermions-1, Manager.GetInteger("max-grid"), GridPositions, TestWaveFunction, History);

	      cout << "Grid-points written"<<endl;
	      
	      return 0;
	    }
	  /* testing code: write coordinates on cluster of grid-points */
	  if (Manager.GetBoolean("write-random"))
	    {
	      cout << "Writing random sampling"<<endl;
	      for (int i=0; i<NbrIter; ++i)
		{
		  Particles->Randomize();
		  Complex Value = (*TestWaveFunction)(Particles->GetPositions());
		  double SA=1.0;
		  History->RecordAcceptedStep( SA, Particles->GetPositions(), Value);
		}
	      cout << "Random points written"<<endl;
	      return 0;
	    }
	  

	  delete [] tmpC;
	}
      else if ((HistoryFileName==NULL)&&(HistoryMode>1))
	{
	  cout << "History mode "<<HistoryMode<<" requires a history file!" << endl;
	  return -1;
	}
    }
  
  if (HistoryMode ==2)
    {
      // code to process available samples
      int sampleCount;
      int totalSampleCount=0;
      RealVector Positions(2*NbrFermions);
      WeightedRealObservable NormTrialObs(100);
      WeightedRealObservable NormExactObs(100);
      WeightedComplexObservable OverlapObs(100);
      int NbrCoordinates = 2*NbrFermions;
      History = new MCHistoryRecord(HistoryFileName, NbrCoordinates /* could add additional observables here */);
      double typicalSA=0.0, typicalWF=0.0, typicalTV=0.0;
      int averageTypical=50;
      int initialSkip=10;
      for (int i=0; i<initialSkip; ++i)
	History->GetMonteCarloStep(sampleCount, CurrentSamplingAmplitude, &(Positions[0]), ValueExact);
      bool MoreData=true;
      for (int i=0; (i<averageTypical)&MoreData; ++i)
	{
	  MoreData=History->GetMonteCarloStep(sampleCount, CurrentSamplingAmplitude, &(Positions[0]), ValueExact);
	  typicalSA+=CurrentSamplingAmplitude;
	  typicalWF+=Norm(ValueExact);
	  typicalTV+=Norm((*TestWaveFunction)(Positions));
	}
      typicalSA/=averageTypical;
      typicalWF/=averageTypical;
      typicalTV/=averageTypical;
      History->RewindHistory();
      cout << "typicalSA= " << typicalSA << ", typicalWF="<<typicalWF<<", typicalTV="<<typicalTV<<endl;
      int i=0;
      Complex rawExact;
      while ( (i++<NbrIter) && (History->GetMonteCarloStep(sampleCount, CurrentSamplingAmplitude, &(Positions[0]), ValueExact)))
	{
	  rawExact=ValueExact;
	  totalSampleCount+=sampleCount;
	  TrialValue = (*TestWaveFunction)(Positions)/typicalTV;	  
	  CurrentSamplingAmplitude /= typicalSA;
	  ValueExact /= typicalWF;
	  if (Norm(ValueExact)>100.0)
	    {
	      cout << i << ": excluding large Psi: " << ValueExact << endl;
	      if (State.GetVectorDimension()>0) // have exact state!
		{
		  cout << "Checking if abnormal result is reproduced!" << endl;
		  if (Space==NULL) Space = new FermionOnSphere(NbrFermions, Lz, LzMax);
		  ParticleOnSphereFunctionBasis Basis(LzMax,ParticleOnSphereFunctionBasis::LeftHanded);  
		  QHEParticleWaveFunctionOperation Operation(Space, &State, &Positions, &Basis, /* TimeCoherence */ -1);
		  Operation.ApplyOperation(Architecture.GetArchitecture());      
		  ValueExact = Operation.GetScalar();
		  cout << "Comparing: " << ValueExact << " (new) to "<< rawExact << " (old) (ratio " << ValueExact/rawExact << ")" <<endl;
		  cout << "Minimal distance: " << MinDist(Positions);
		}
	    }
	  else
	    {
	      // could change here to add multiple observations in order to get errors right in sampling density of one
	      // (see FQHESphereFermionsWithSpinOverlap.cc)
	      NormTrialObs.Observe(SqrNorm(TrialValue)/CurrentSamplingAmplitude,(double)sampleCount);
	      NormExactObs.Observe(SqrNorm(ValueExact)/CurrentSamplingAmplitude,(double)sampleCount);
	      OverlapObs.Observe(Conj(TrialValue)*ValueExact/CurrentSamplingAmplitude,(double)sampleCount);
	    }
//              if (sampleCount>10)
// 	    {
// 	      cout << "Total "<<sampleCount<<" samples in these coordinates: " << endl;
// 	      cout << "Psi^2=" << SqrNorm(ValueExact) <<" trial^2="<<CurrentSamplingAmplitude<<endl;
// 	      cout << "Minimal distance: " << MinDist(Positions);
// 	      //cout << Positions<< endl;
// 	    }		
	}
      if (--i>NbrIter) cout << "Attention, step number limited by NbrIter!" << endl;
      History->RewindHistory();      
      
      cout << " final results :" << endl;
      cout << "overlap:   " << OverlapValue(OverlapObs, NormTrialObs, NormExactObs) << " +/- "
	   << OverlapError(OverlapObs, NormTrialObs, NormExactObs) << endl;
      cout << "overlap^2: " << SqrNorm(OverlapValue(OverlapObs, NormTrialObs, NormExactObs)) << " +/- " <<
	2.0*Norm(OverlapValue(OverlapObs, NormTrialObs, NormExactObs))*OverlapError(OverlapObs, NormTrialObs, NormExactObs) << endl;
      cout << "Processed a total of "<<totalSampleCount<<" MC samples" << endl;
      delete History;
      return 0;
    }


  if (HistoryMode ==3)
    {
      if (!(WaveFunctionManager.GetWaveFunctionType() & QHEWaveFunctionManager::TrialWaveFunction))
	{
	  cout << "This type of wavefunction cannot be optimized." << endl;
	  return -1;
	}
      // code to optimize paired state with available samples:
      bool varyMR = Manager.GetBoolean("varyMR");
      WaveFunctionOverlapOptimizer *Optimizer =
	new WaveFunctionOverlapOptimizer( TestWaveFunction, HistoryFileName, NbrFermions,
					  /* excludeLastParameter */ !varyMR,
					  Manager.GetInteger("linearPoints"), Manager.GetInteger("randomPoints"),
					  Manager.GetInteger("limitSamples"));
      RealVector optimalParameters( ((Abstract1DComplexTrialFunction*) TestWaveFunction)->GetNbrParameters());
      Complex optimalOverlap;
      Optimizer->GetMaximumSqrOverlap(optimalParameters, optimalOverlap);
      return 0;
    }
      
    
  if (!UsePfaffian)
    {
      Space = new FermionOnSphere(NbrFermions, Lz, LzMax);
      cout << "Hilbertspace dimension="<<Space->GetHilbertSpaceDimension()<<endl;
    }
  ParticleOnSphereFunctionBasis Basis(LzMax,ParticleOnSphereFunctionBasis::LeftHanded);  
  
  double PreviousSamplingAmplitude;
  
  int RecordStep = Manager.GetInteger("record-step");

  Complex* RecordedOverlap = 0;
  double* RecordedOverlapError = 0;
  if (RecordStep != 0)
    {
      RecordedOverlap = new Complex [(NbrIter / RecordStep) + 1];
      RecordedOverlapError = new double [(NbrIter / RecordStep) + 1];
    }
  int RecordIndex = 0;
  double Factor = 1.0;
  for (int j = 0; j < NbrFermions; ++j)
    {
      Factor *= 4.0 * M_PI;
    }

  // MC observables:
  Complex Overlap(0.0);
  ComplexObservable ScalarProduct(100);    
  RealObservable NormExact(100);
  // stuff
  Complex TmpMetropolis;
  int NextCoordinates = 0;
  int Accepted = 0;
  bool NoTimeCoherence=!(Manager.GetBoolean("with-timecoherence"));
  int TimeCoherence;

  if (Manager.GetBoolean("product-state"))
    ((PolarizedProductWavefunction*)TestWaveFunction)->TestSymmetries(Particles);


  if (HistoryMode == 4) // continuing to work on old History
    {
      if (Manager.GetString("exact-state") == 0)
	{
	  cout << "QHEFermionOverlap requires an exact state" << endl;
	  return -1;
	}  
      if (State.ReadVector (Manager.GetString("exact-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("exact-state") << endl;
	  return -1;      
	}
      History = new MCHistoryRecord (HistoryFileName, 2*NbrFermions, PreviousSamplingAmplitude, &(Particles->GetPositions()[0]), ValueExact, NULL);
      // initialize function values at initial positions: - trial function
      TrialValue = (*TestWaveFunction)(Particles->GetPositions());  

      // - exact function
      TimeCoherence = -1;
      QHEParticleWaveFunctionOperation Operation(Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      if (SqrNorm(ValueExact - Operation.GetScalar()) > 1e-10)
	{
	  cout << "Exact wavefunction from record not exactly reproduced: "
	       << ValueExact << " vs " <<Operation.GetScalar()<<endl;
	  
	  cout << "Exact state may have changed!" << endl;
	}
      
      CurrentSamplingAmplitude = SqrNorm(TrialValue);;
      if ( fabs(PreviousSamplingAmplitude -  CurrentSamplingAmplitude) > 1e-10)
	{
	  cout << "SamplingAmplitude from record not exactly reproduced: "
	       << PreviousSamplingAmplitude << " vs " <<CurrentSamplingAmplitude<<endl;
	  cout << "Sampling function may have changed" << endl;
	}
    }
  else
    {
      // - exact function
      if (UsePfaffian)
	ValueExact = (*PfaffianWaveFunction)(Particles->GetPositions());
      else
	{
	  TimeCoherence = -1;
	  QHEParticleWaveFunctionOperation Operation(Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  ValueExact = Operation.GetScalar();
	}
      // initialize function values at initial positions: - trial function
      TrialValue = (*TestWaveFunction)(Particles->GetPositions());
      PreviousSamplingAmplitude = SqrNorm(TrialValue);
      CurrentSamplingAmplitude = PreviousSamplingAmplitude;

      History->RecordInitialPositions( CurrentSamplingAmplitude, Particles->GetPositions(), ValueExact);
    } 


// testing code: comparing wavefunctions

//   AdvancedReadRezayiOnSphereWaveFunction *MR1 = new AdvancedReadRezayiOnSphereWaveFunction(NbrFermions/2, 2, true);
//   ExplicitMooreReadOnSphereWaveFunction *MR2 = new ExplicitMooreReadOnSphereWaveFunction(NbrFermions/2, 2, true);

//    Complex T1, T2;

  // thermalization loop
  for (int i = 0; i < NbrWarmUpIter; ++i)
    {
      // make a random move of particle NextCoordinates
      Particles->Move(NextCoordinates);
      // evaluate new trial wavefunction value
      TmpMetropolis = (*TestWaveFunction)(Particles->GetPositions());

      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((RandomNumber->GetRealRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      // determine next particle to move
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrFermions) --NextCoordinates;
    } // end thermalization
  
  for (int i = 0; i < NbrIter; ++i)
    {
      // make a random move of particle NextCoordinates
      Particles->Move(NextCoordinates);
      // evaluate new trial wavefunction value
      TmpMetropolis = (*TestWaveFunction)(Particles->GetPositions());

      //testing:
//        T1=(*MR1)(Particles->GetPositions());
//        T2=(*MR2)(Particles->GetPositions());
//        cout << "Values: " << T1 << " vs "<< T2 << endl;

      
      // accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((RandomNumber->GetRealRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	  ++Accepted;	  
	  if (SampleDensity==1) // do not recalculate Exact or count steps here if using spread out samples
	    {
	      // recalculate exact function value:
	      TimeCoherence = NextCoordinates;
	      if (UsePfaffian)
		ValueExact = (*PfaffianWaveFunction)(Particles->GetPositions());
	      else
		{
		  if (NoTimeCoherence) TimeCoherence = -1;      
		  QHEParticleWaveFunctionOperation Operation(Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
		  Operation.ApplyOperation(Architecture.GetArchitecture());      
		  ValueExact = Operation.GetScalar();
		}
	      if (History)
		History->RecordAcceptedStep( CurrentSamplingAmplitude, Particles->GetPositions(), ValueExact);
	    }
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	  if ((History)&&(SampleDensity==1)) // do not count steps here if using spread out samples
	    History->RecordRejectedStep();
	}
      if ((SampleDensity>1)&&( (i % SampleDensity)==0 )) // recalculate and record spread out samples here if SampleDensity>1
	{
	  // recalculate exact function value:
	  TimeCoherence = NextCoordinates;
	  if (UsePfaffian)
	    ValueExact = (*PfaffianWaveFunction)(Particles->GetPositions());
	  else
	    {
	      if (NoTimeCoherence) TimeCoherence = -1;      
	      QHEParticleWaveFunctionOperation Operation(Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
	      Operation.ApplyOperation(Architecture.GetArchitecture());      
	      ValueExact = Operation.GetScalar();
	    }
	  if (History) History->RecordAcceptedStep( CurrentSamplingAmplitude, Particles->GetPositions(), ValueExact);
	}
      // determine next particle to move
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrFermions) --NextCoordinates;
      
      // note observations:
      if ( (i % SampleDensity)==0 )
	{
	  double norm = SqrNorm(ValueExact)/CurrentSamplingAmplitude;
	  NormExact << norm;
	  Overlap = Conj(TrialValue)*ValueExact/CurrentSamplingAmplitude;
	  ScalarProduct << Overlap;
	}
      if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	{
	  RecordedOverlap[RecordIndex] = OverlapValue(ScalarProduct, NormExact);
	  RecordedOverlapError[RecordIndex] = OverlapError(ScalarProduct, NormExact);
	  ++RecordIndex;
	}
      if ((i > 0) && ((i % (Manager.GetInteger("display-step"))) == 0))
	{
	  cout << " i = " << i << endl;
	  cout << OverlapValue(ScalarProduct, NormExact) << " +/- " << OverlapError(ScalarProduct, NormExact) << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    }
  delete History;
  cout << " final results :" << endl;
  cout << "overlap:   " << OverlapValue(ScalarProduct, NormExact) << " +/- "
       << OverlapError(ScalarProduct, NormExact) << endl;
  cout << "overlap^2: " << SqrNorm(OverlapValue(ScalarProduct, NormExact)) << " +/- " <<
    2.0*Norm(OverlapValue(ScalarProduct, NormExact))*OverlapError(ScalarProduct, NormExact) << endl;
  cout << "Fraction of moves accepted: " << (double) Accepted / NbrIter *100.0<< " %"<< endl;
  cout << "-----------------------------------------------" << endl;  
  if (RecordStep != 0)
    {
      RecordedOverlap[RecordIndex] = OverlapValue(ScalarProduct, NormExact);
      RecordedOverlapError[RecordIndex] = OverlapError(ScalarProduct, NormExact);
      ofstream OverlapRecordFile;
      OverlapRecordFile.precision(14);
      OverlapRecordFile.open(Manager.GetString("record-file"), ios::out | ios::binary);
      int NbrRecords = NbrIter / RecordStep;
      OverlapRecordFile << "# Monte Carlo overlap calculation" << endl
		       << "# step overlap.Re overlap.Im overlap2 error error2" << endl;
      for (int i = 0; i < NbrRecords; ++i)
	OverlapRecordFile << i << " " << RecordedOverlap[i].Re << " " << RecordedOverlap[i].Im << " " << SqrNorm(RecordedOverlap[i])
			  << " " << RecordedOverlapError[i] << " " <<  DSQR(RecordedOverlapError[i]) << endl;
      OverlapRecordFile.close();
    }
  delete Particles;
  if (Space!=NULL)
    delete Space;
  return 0;
}


Complex OverlapValue(ComplexObservable &ScalarProduct, RealObservable &NormObs)
{
  return Complex(ScalarProduct.Average()/sqrt(NormObs.Average()));
}
  
double OverlapError(ComplexObservable &ScalarProduct, RealObservable &NormObs)
{
  double norm = NormObs.Average();
  double prod = Norm(ScalarProduct.Average());
  return sqrt( DSQR1(ScalarProduct.ErrorEstimate())/norm + DSQR2(prod*NormObs.ErrorEstimate()/2.0/norm)/norm);
}

Complex OverlapValue(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &NormObs1, WeightedRealObservable &NormObs2)
//Complex OverlapValue(ComplexObservable &ScalarProduct, RealObservable &NormObs1, RealObservable &NormObs2)
{
  return Complex(ScalarProduct.Average()/sqrt(NormObs1.Average()*NormObs2.Average()));
}

double OverlapError(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &NormObs1, WeightedRealObservable &NormObs2)
//double OverlapError(ComplexObservable &ScalarProduct, RealObservable &NormObs1, RealObservable &NormObs2)
{
  double norm1 = NormObs1.Average();
  double norm2 = NormObs2.Average();
  double norm = sqrt(norm1*norm2);
  double prod = Norm(ScalarProduct.Average());
  return sqrt( DSQR1(ScalarProduct.ErrorEstimate()/norm) + DSQR2(prod*NormObs1.ErrorEstimate()*NormObs2.Average()/2.0/norm/norm/norm) + DSQR3(prod*NormObs2.ErrorEstimate()*NormObs1.Average()/2.0/norm/norm/norm));
}
  

void FillRecursiveGrid(int NbrCoordinate, int MaxGrid, RealVector &GridPositions, Abstract1DComplexFunction* TestWaveFunction, MCHistoryRecord *History)
{
  if (NbrCoordinate<0)
    {
      Complex Value = (*TestWaveFunction)(GridPositions);
      double SA=1.0;
      History->RecordAcceptedStep( SA, GridPositions, Value);
    }
  else
    {
      for (int NbrGrid=0; NbrGrid<MaxGrid; ++NbrGrid)
	{
	  if ((NbrCoordinate & 1) != 0)
	    GridPositions[NbrCoordinate] = NbrGrid * 2.0*M_PI/(double)MaxGrid;
	  else
	    GridPositions[NbrCoordinate] = NbrGrid * M_PI/(double)(MaxGrid-1);
	  FillRecursiveGrid(NbrCoordinate-1,MaxGrid, GridPositions, TestWaveFunction, History);
	}
    }
}

