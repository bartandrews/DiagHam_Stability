#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/ClebschGordanCoefficients.h"

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

#include <iostream>
#include <cstring>
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
  OptionManager Manager ("FQHESphereHistoryCheck" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lz", "twice the momentum projection", 0);
  (*SystemGroup) += new SingleStringOption  ('r', "reference-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new SingleStringOption ('c', "history-to-check", "name of the file where MC samples were recorded", NULL);
  (*SystemGroup) += new SingleStringOption ('o', "history-output", "name of the fixed output file", NULL);
  (*SystemGroup) += new SingleDoubleOption  ('t', "threshold", "threshold for testing exact values (0.0=test all)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('f', "force-discard", "dump all samples above this threshold", 1e4);
  (*SystemGroup) += new SingleDoubleOption  ('d', "discard-limit", "deviation before discarding a sample", 1e-10);
  (*SystemGroup) += new SingleIntegerOption  ('s', "step-limit", "maximum number of samples to be checked", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereHistoryCheck -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz");
  RealVector State;
  
  if (Manager.GetString("reference-state") == 0)
    {
      cout << "FQHESphereHistoryCheck requires the exact state used to generate the History" << endl;
      return -1;
    }
  if (State.ReadVector (Manager.GetString("reference-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("reference-state") << endl;
      return -1;      
    }

  MCHistoryRecord *History;
  if (Manager.GetString("history-to-check") == 0)
    {
      cout << "FQHESphereHistoryCheck requires a History-file to be checking" << endl;
      return -1;
    }
  int NbrCoordinates = 2*NbrFermions;
  History = new MCHistoryRecord(Manager.GetString("history-to-check"), NbrCoordinates
				/* could add additional observables here */);

  
  // get an idea of the average amplitude of the different numbers in HistoryFile:
  int sampleCount;
  double SamplingAmplitude;
  Complex ExactValue;
  int toCheck = History->GetProjectedSamples();
  if (toCheck > 1000) toCheck=1000;
  int averageTypical=toCheck/2;
  int initialSkip=toCheck/2;
  bool HaveMoreHistory=true;
  RealVector Positions(2*NbrFermions);
  double typicalSA=0.0, typicalWF=0.0;
  for (int i=0; (HaveMoreHistory&&(i<initialSkip)); ++i)
    HaveMoreHistory=History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), ExactValue);
  if (HaveMoreHistory)
    {
      HaveMoreHistory=true;
      for (int i=0;(HaveMoreHistory&&(i<averageTypical)); ++i)
	{
	  HaveMoreHistory=History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), ExactValue);
	  typicalSA+=SamplingAmplitude;
	  typicalWF+=Norm(ExactValue);
	}
      typicalSA/=averageTypical;
      typicalWF/=averageTypical;
    }
  else
    {
      typicalSA=1.0;
      typicalWF=1.0;
    }

  cout << "typicalSA= " << typicalSA << ", typicalWF="<<typicalWF<<endl;

  History->RewindHistory();

  FermionOnSphere Space (NbrFermions, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax,ParticleOnSphereFunctionBasis::LeftHanded);


  char *HistoryOutputName;

  double threshold = Manager.GetDouble("threshold");
  double forceDiscard = Manager.GetDouble("force-discard");
  double limit = Manager.GetDouble("discard-limit");
  int StepLimit  = Manager.GetInteger("step-limit");
  if (StepLimit==0) StepLimit=History->GetProjectedSamples();
  
  if (Manager.GetString("history-output") == 0)
    {
      HistoryOutputName = new char[strlen(Manager.GetString("history-to-check"))+40];
      sprintf(HistoryOutputName,"%s.fix_l%g_t%g",Manager.GetString("history-to-check"),limit,threshold);
    }
  else HistoryOutputName =  Manager.GetString("history-output");
  
  MCHistoryRecord *HistoryOutput=new MCHistoryRecord(History->GetProjectedSamples(), 2*NbrFermions, Manager.GetString("reference-state"), NULL, HistoryOutputName);

  int totalSampleCount=0;
  int conservedSampleCount=0;
  int checkedConfigurations=0;
  int keptFromChecked=0;
  double maxKeep=0.0;
  double minDiscard=1e300;
  double norm;
  double factor=1.0;
  Complex ratio, readValue, recalculatedValue;

  // accomodate possibility that some arbitrary norm has been introduced
 nextSample:
  History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), readValue);
  QHEParticleWaveFunctionOperation Operation(&Space, &State, &Positions, &Basis, /* TimeCoherence */ -1);
  Operation.ApplyOperation(Architecture.GetArchitecture());      
  recalculatedValue = Operation.GetScalar();
  ratio = recalculatedValue/readValue;
  if (fabs(Real(ratio)-1.0) > limit)
    {
      if (fabs(Imag(ratio)/Real(ratio)) < limit)
	{
	nextSample2:
	  factor = Real(ratio);
	  History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), readValue);
	  QHEParticleWaveFunctionOperation Operation(&Space, &State, &Positions, &Basis, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(Architecture.GetArchitecture());      
	  recalculatedValue = Operation.GetScalar();
	  ratio = recalculatedValue/readValue;
	  if (fabs(Imag(ratio)/Real(ratio)) < limit)
	    {
	      if( fabs(Real(ratio)/factor-1.0) < limit)
		{
		  factor = 0.5*(factor+Real(ratio));
		  cout << "Attention: Found a change of normalization by a factor of " << factor << endl;
		}
	      else
		{
		  cout << "No consistent change in normalization apparent" << endl;
		  exit(1);
		}
	    }
	  else goto nextSample2;
	}
      else goto nextSample;
    }
  
  while ((totalSampleCount<StepLimit)&&(History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), readValue)))
    {
      SamplingAmplitude /= typicalSA;
      ExactValue = readValue/typicalWF;
      totalSampleCount+=sampleCount;
      if ((norm=Norm(ExactValue))>threshold)
	{
	  cout << "testing " << totalSampleCount << ": " << norm<<" ";
	  if (norm < forceDiscard)
	    {
	      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Positions, &Basis, /* TimeCoherence */ -1);
	      Operation.ApplyOperation(Architecture.GetArchitecture());      
	      recalculatedValue = Operation.GetScalar();
	      ratio = recalculatedValue/readValue/factor;
	      cout  << recalculatedValue << "\t"<< readValue << "\t" << ratio;
	      checkedConfigurations++;
	      if ((fabs(Real(ratio)-1.0) < limit) && (fabs(Imag(ratio)) < limit) && (norm < forceDiscard))
		{	      
		  cout << " -> keep" << endl;
		  keptFromChecked++;
		  conservedSampleCount+=sampleCount;	      
		  HistoryOutput->RecordAcceptedStep( SamplingAmplitude, Positions, readValue);
		  for (int i=1; i<sampleCount; ++i) HistoryOutput->RecordRejectedStep();
		  if (norm>maxKeep) maxKeep = norm;
		}
	      else
		{	      
		  cout << " -> discard";
		  if (norm > forceDiscard) cout << " forced" << endl;
		  else cout  << endl;
		  if (norm<minDiscard) minDiscard = norm;
		}
	    }
	  else
	    cout << " -> discard forced" << endl;	      
	}
      else
	{
	  conservedSampleCount+=sampleCount;	  
	  HistoryOutput->RecordAcceptedStep( SamplingAmplitude, Positions, readValue);
	  for (int i=1; i<sampleCount; ++i) HistoryOutput->RecordRejectedStep();
	}	
    }
  cout << "Kept " << conservedSampleCount<< " / "<<totalSampleCount << " samples."<< endl;
  cout << "Among examined configurations, kept " << keptFromChecked << " / " << checkedConfigurations<< endl;
  if ((checkedConfigurations==0) || (keptFromChecked==0))
    cout << "You might want to redo this procedure with a lower value of the threshold." << endl;
  cout << "Maximal Norm of kept sample: " << maxKeep << " Minimal Norm of discarded sample: "<< minDiscard << endl;
  cout << "Output saved to " << HistoryOutputName << endl;
  delete HistoryOutput;

}
