#include "MultiFileMCHistoryRecord.h"

using std::ios;
using std::cout;
using std::endl;

MultiFileMCHistoryRecord::MultiFileMCHistoryRecord()
{
  InputHistories=NULL;
}

// to merge a number of records into a single file, (remain open in reading mode)
// outputFileName: new history file to be generated (NULL => only open for reading)
// nbrInputFiles: number of input files
// inputFileNames: history files to be merged
// additionalData: description of additional fields
MultiFileMCHistoryRecord::MultiFileMCHistoryRecord(char* outputFileName, int nbrInputFiles, char **inputFileNames, int &nbrPositions, List<AbstractMCHistoryData> *additionalData)
{
  this->NbrHistoryFiles = nbrInputFiles;
  this->InputHistories = new MCHistoryRecord*[NbrHistoryFiles];
  this->RecordMode = AbstractMCHistoryRecord::Reading;
  this->ProjectedStepNum = 0;
  for (int i=0; i<NbrHistoryFiles; ++i)
    {
      this->InputHistories[i] = new MCHistoryRecord(inputFileNames[i], nbrPositions, additionalData);
      this->ProjectedStepNum += this->InputHistories[i]->GetProjectedSamples();
    }
  this->NbrPositions = nbrPositions; // value will be read from first input file! 
  
  if (outputFileName!=NULL)
    {
      int SampleCount;
      double SamplingAmplitude;
      RealVector Positions(NbrPositions);      
      Complex ValueExact;
      
      MCHistoryRecord *OutputHistory = new MCHistoryRecord(this->ProjectedStepNum, this->NbrPositions, "as in source files",
							   "merged", outputFileName);
      
      for (NbrPresentHistory=0; NbrPresentHistory<NbrHistoryFiles; ++ NbrPresentHistory)
	{
	  while (InputHistories[NbrPresentHistory]->GetMonteCarloStep(SampleCount, SamplingAmplitude, &(Positions[0]), ValueExact))
	    {
	      // record accepted step - to be called for each accepted step, or at every step to be written to file
	      OutputHistory->RecordAcceptedStep( SamplingAmplitude, Positions, ValueExact);
	    }
	  InputHistories[NbrPresentHistory]->RewindHistory();
	}
      delete OutputHistory;
    }
  this->NbrPresentHistory=0; //effectively rewind to start!

}



// destructor -> automatically closes LogFile
MultiFileMCHistoryRecord::~MultiFileMCHistoryRecord()
{
  if (this->InputHistories!=0)
    {
      for (NbrPresentHistory=0; NbrPresentHistory<NbrHistoryFiles; ++ NbrPresentHistory)
	delete this->InputHistories[NbrPresentHistory];
      delete [] this->InputHistories;
    }
}
    
  // read one MC sample back from file, gives back the parameters in call of RecordAcceptedStep
  // sampleCount additionally gives the multiplicity of each record
bool MultiFileMCHistoryRecord::GetMonteCarloStep( int &sampleCount, double & samplingAmplitude, double *positions, Complex &valueExact)
{
  bool success = InputHistories[NbrPresentHistory]->GetMonteCarloStep(sampleCount, samplingAmplitude, positions, valueExact);
  if (success) return success;
  else
    {
      if (NbrPresentHistory<NbrHistoryFiles-1)
	{
	  ++NbrPresentHistory;
	  success = InputHistories[NbrPresentHistory]->GetMonteCarloStep(sampleCount, samplingAmplitude, positions, valueExact);
	}
      return success;
    }
}


// rewind in reading mode:
void MultiFileMCHistoryRecord::RewindHistory()
{
  for (NbrPresentHistory=0; NbrPresentHistory<NbrHistoryFiles; ++ NbrPresentHistory)
    InputHistories[NbrPresentHistory]->RewindHistory();
  this->NbrPresentHistory=0; //effectively rewind to start!
}


