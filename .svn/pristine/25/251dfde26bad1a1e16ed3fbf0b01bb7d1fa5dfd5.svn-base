#include "MCHistoryRecord.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ListIterator.h"

#include <iostream>
#include <cstdlib>

using std::ios;
using std::cout;
using std::endl;

// dis/activate testing output
//#define MCHistoryRecord_TESTING

// constructors
// for recording mode
// projectedSamples : number of samples to be expected
// nbrPositions : number of (individual 1d) coordinates contained within each record
// exactFile : name of exact file that was sampled
// samplingDescriptor : description of sampling function
// fileName : name of outputFile to be created
// additionalData : description of (optional) additional data to be recorded
MCHistoryRecord::MCHistoryRecord(int projectedSamples, int nbrPositions, const char* exactFile, const char* samplingDescriptor, char* fileName, List<AbstractMCHistoryData> *additionalData)
{
  this->RecordMode = AbstractMCHistoryRecord::Recording;
  this->LastSampleCount=0;
  this->TotalSampleCount=0;
  this->TotalRecordCount=0;
  this->ProjectedStepNum=projectedSamples;
  this->NbrPositions = nbrPositions;
  LogFile.open(fileName, ios::binary | ios::out);  
  WriteLittleEndian(LogFile, this->ProjectedStepNum);
  WriteLittleEndian(LogFile, this->NbrPositions);  
  if (additionalData != NULL)
    {
      this->NumAdditionalData=additionalData->GetNbrElement();
      WriteLittleEndian(LogFile, this->NumAdditionalData);
      this->AdditionalData=new AbstractMCHistoryData*[NumAdditionalData];
      AbstractMCHistoryData *Data;
      int i=0;
      unsigned tmpU;
      for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;++i)
	{
	  tmpU=Data->GetHistoryDataType();
	  WriteLittleEndian(LogFile,tmpU);
	  tmpU=Data->GetSizeOnDisk();
	  WriteLittleEndian(LogFile,tmpU);
	  this->AdditionalData[i]=Data;
	}
    }
  else
    {
      this->NumAdditionalData=0;
      this->AdditionalData=NULL;
      WriteLittleEndian(LogFile, this->NumAdditionalData);
    }
}


// to continue aborted calculation
// fileName: file to continue
// double & samplingAmplitude, double *positions, Complex &valueExact: describing last recorded positions
MCHistoryRecord::MCHistoryRecord(char* fileName, int nbrPositions, double & samplingAmplitude, double *positions, Complex &valueExact, List<AbstractMCHistoryData> *additionalData)
{
  this->RecordMode = AbstractMCHistoryRecord::Recording;
  HistoryFile.open(fileName, ios::binary | ios::in);
  this->TotalSampleCount=0;
  ReadLittleEndian(HistoryFile, this->ProjectedStepNum);
  ReadLittleEndian(HistoryFile, this->NbrPositions);
  if (this->NbrPositions != nbrPositions)
    {
      cout << "Mismatch of nbrPositions in reading MCHistoryRecord" << endl;
      exit (-1);
    }
  int sampleCount=0;
  int numAdditional;
  ReadLittleEndian(HistoryFile,numAdditional);
  this->AdditionalData=NULL;
  if (numAdditional > 0)
    {
      cout << "Additional data present in MCHistoryRecord "<< fileName << endl;
      if (additionalData != NULL)
	{
	  if (numAdditional != additionalData->GetNbrElement())
	    {
	      cout << "Mismatch of the number of additional MCHistoryRecord's you wish to read" << endl;
	      exit (1);
	    }	  
	  this->NumAdditionalData=numAdditional;
	  this->AdditionalData=new AbstractMCHistoryData*[NumAdditionalData];
	  this->SkipAdditional=0;
	  AbstractMCHistoryData *Data;
	  int i=0;
	  for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;++i)
	    {
	      unsigned check, size;
	      ReadLittleEndian(HistoryFile,check);
	      ReadLittleEndian(HistoryFile,size);
	      if (check != Data->GetHistoryDataType())
		{
		  cout << "Mismatch of field number "<<i+1<<" for the additional MCHistoryRecord's you wish to read" << endl;
		  exit (1);
		}
	      this->AdditionalData[i]=Data;
	    }
	}
      else  // data present, but not requested.
	{
	  cout << "Additional data present in file is being discarded." << endl;
	  this->NumAdditionalData=0;
	  this->SkipAdditional=0;
	  unsigned check,skip;
	  AbstractMCHistoryData *Data;
	  for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;)
	    {
	      ReadLittleEndian(HistoryFile,check);
	      ReadLittleEndian(HistoryFile,skip);
	      this->SkipAdditional+=skip;
	    }
	}
    }
  else
    { // no additional data:
      this->NumAdditionalData=0;
      this->SkipAdditional=0;
    }
  // first entry from writing is lastSampleCount that was zero...
  ReadLittleEndian(HistoryFile, this->LastSampleCount);
//   if (this->LastSampleCount!=0)
//     {
//       cout << "Problem with header of History record " <<fileName << endl;
//       exit(2);
//     }

  // find last block:
  int skip = ((unsigned)NumAdditionalData>SkipAdditional ? NumAdditionalData : SkipAdditional);
  std::streampos secondLastBlock, lastBlock=0;
  char signature; 
  while ( ! (HistoryFile.eof()))
    {
      ReadLittleEndian(HistoryFile,signature);
      if (signature == 'b') // recognized a new block
	{
	  secondLastBlock=lastBlock;
	  lastBlock=HistoryFile.tellg();
	  ReadLittleEndian(HistoryFile,samplingAmplitude);
	  for (int i = 0; i < this->NbrPositions; ++i)
	    ReadLittleEndian(HistoryFile, positions[i]);
	  ReadLittleEndian(HistoryFile,valueExact);
	  if (skip>0)
	    HistoryFile.seekg(skip,ios::cur);
	  ReadLittleEndian(HistoryFile,this->LastSampleCount);
	  if ( ! (HistoryFile.eof()))
	    {
	      this->TotalSampleCount+=this->LastSampleCount;
	      sampleCount=this->LastSampleCount;
	    }
	}      
    }
  if (HistoryFile.eof()) HistoryFile.clear();
  if (signature == 'e')
    {
      cout << "Regular end of file detected." << endl;	        
      HistoryFile.seekg(lastBlock);
    }
  else if (signature == 'b')
    {
      cout << "File ends with half finished block..." << endl;
      HistoryFile.seekg(secondLastBlock);
    }     
  ReadLittleEndian(HistoryFile,samplingAmplitude);
  for (int i = 0; i < this->NbrPositions; ++i)
    ReadLittleEndian(HistoryFile, positions[i]);
  ReadLittleEndian(HistoryFile,valueExact);
  if (skip>0)
    HistoryFile.seekg(skip,ios::cur);
  std::streampos WritePos = HistoryFile.tellg();
  if (signature == 'e')
    {
      ReadLittleEndian(HistoryFile,this->LastSampleCount);
      if (sampleCount!=this->LastSampleCount) cout << "problem with reading last block at second time" << endl;
    }
  else this->LastSampleCount=1;
  sampleCount=this->LastSampleCount;
  HistoryFile.close();
  LogFile.open(fileName, ios::binary | ios::in | ios::out );
  LogFile.seekp(WritePos, ios::beg);

}


// for reading mode
// Input = name of MCHistoryRecord file to be read
// nbrPositions = requested number of positions in each record (0=return number of positions in record)
// additionalData = descriptor of eventual additional fields
MCHistoryRecord::MCHistoryRecord(char *Input, int &nbrPositions, List<AbstractMCHistoryData> *additionalData)
{
  this->RecordMode = AbstractMCHistoryRecord::Reading;
  HistoryFile.open(Input, ios::binary | ios::in);
  this->TotalSampleCount=0;
  ReadLittleEndian(HistoryFile, this->ProjectedStepNum);
  ReadLittleEndian(HistoryFile, this->NbrPositions);
  if ((nbrPositions>0)&&(this->NbrPositions != nbrPositions))
    {
      cout << "Mismatch of nbrPositions in reading MCHistoryRecord" << endl;
      exit (-1);
    }
  else
    {
      nbrPositions = this->NbrPositions;
    }
  cout << "MCHistoryRecord "<< Input << " for system with "<<NbrPositions/2<< " particles"<<endl;
  int numAdditional;
  ReadLittleEndian(HistoryFile,numAdditional);
  this->AdditionalData=NULL;
  if (numAdditional > 0)
    {
      cout << "Additional data present in MCHistoryRecord "<< Input << endl;
      if (additionalData != NULL)
	{
	  if (numAdditional != additionalData->GetNbrElement())
	    {
	      cout << "Mismatch of the number of additional MCHistoryRecord's you wish to read" << endl;
	      exit (1);
	    }	  
	  this->NumAdditionalData=numAdditional;
	  this->AdditionalData=new AbstractMCHistoryData*[NumAdditionalData];
	  this->SkipAdditional=0;
	  AbstractMCHistoryData *Data;
	  int i=0;
	  for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;++i)
	    {
	      unsigned check, size;
	      ReadLittleEndian(HistoryFile,check);
	      ReadLittleEndian(HistoryFile,size);
	      if (check != Data->GetHistoryDataType())
		{
		  cout << "Mismatch of field number "<<i+1<<" for the additional MCHistoryRecord's you wish to read" << endl;
		  exit (1);
		}
	      this->AdditionalData[i]=Data;
	    }
	}
      else  // data present, but not requested.
	{
	  cout << "Additional data present in file is being discarded." << endl;
	  this->NumAdditionalData=0;
	  this->SkipAdditional=0;
	  unsigned check,skip;
	  AbstractMCHistoryData *Data;
	  for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;)
	    {
	      ReadLittleEndian(HistoryFile,check);
	      ReadLittleEndian(HistoryFile,skip);
	      this->SkipAdditional+=skip;
	    }
	}
    }
  else
    { // no additional data:
      cout << "No additional data in MCHistoryRecord "<< Input << endl;
      this->NumAdditionalData=0;
      this->SkipAdditional=0;
    }
  // first entry from writing is lastSampleCount that was zero...
  ReadLittleEndian(HistoryFile, this->LastSampleCount);
   if (this->LastSampleCount!=0)
     {
       cout << "Problem with header of History record " <<Input << endl;
       cout << "LastSampleCount = "<<LastSampleCount<<endl;
       exit(2);
     }
  this->StartPos=HistoryFile.tellg();
  // cout << "Constructor: StartPos is: " << StartPos << " peeking: " <<HistoryFile.peek() << endl;
}


MCHistoryRecord::~MCHistoryRecord()
{
  if (this->RecordMode & MCHistoryRecord::Recording)
    {
      // write last multiplicity
      WriteLittleEndian(LogFile,this->LastSampleCount);
#ifdef MCHistoryRecord_TESTING
      cout << "Recorded final count: "<<this->LastSampleCount<<endl;
#endif
      cout << "Total " << TotalRecordCount << " History records written." << endl;
      char signature = 'e'; // signature for the END
      WriteLittleEndian(LogFile,signature);
      LogFile.flush();
      LogFile.close();
      if (AdditionalData != 0) delete [] AdditionalData;
    }
  else if (this->RecordMode & AbstractMCHistoryRecord::Reading)
    {
      HistoryFile.close();
      if (AdditionalData != 0) delete [] AdditionalData;
    }
}

// record rejected step - to be called for rejected Microsteps
void MCHistoryRecord::RecordRejectedStep()
{
  this->LastSampleCount++;
  this->TotalSampleCount++;
#ifdef MCHistoryRecord_TESTING
  cout << "Step rejected"<<endl;
#endif
}

// record accepted step - to be called for each accepted step, or at every step to be written to file
bool MCHistoryRecord::RecordAcceptedStep( double samplingAmplitude, RealVector &positions, Complex &valueExact)
{
  // fix to assure that first recorded step is one that was accepted:
  //if (LastSampleCount==TotalSampleCount) LastSampleCount=0;
  //
  this->TotalSampleCount++;
  this->TotalRecordCount++;
  // first: write multiplicity of last configuration (is part of prior record)  
  WriteLittleEndian(LogFile,this->LastSampleCount);
  // then: write info about new positions (starting new record)
  char signature = 'b'; // signature for a new Block
  WriteLittleEndian(LogFile,signature);
  WriteLittleEndian(LogFile,samplingAmplitude);
  for (int i = 0; i < this->NbrPositions; ++i)
    WriteLittleEndian(LogFile, positions[i]);
  WriteLittleEndian(LogFile,valueExact);
  if (NumAdditionalData>0)
    {
      for (int i=0; i<NumAdditionalData; ++i) 
	this->AdditionalData[i]->WriteMCHistoryData();
    }
#ifdef MCHistoryRecord_TESTING
  cout << "Recorded previous count: "<<this->LastSampleCount<<", S: "<<samplingAmplitude<<", E: "<<valueExact<<endl;
#endif
  this->LastSampleCount=1; // reset counter for present position
  return true;
}

// record accepted step - to be called for each accepted step, or at every step to be written to file
bool MCHistoryRecord::RecordInitialPositions( double samplingAmplitude, RealVector &positions, Complex &valueExact)
{
  // first: write zero multiplicity by convention
  int tmp=0;
  WriteLittleEndian(LogFile,tmp);

  // start writing info about new positions (starting new record)
  char signature = 'b'; // signature for a new Block
  WriteLittleEndian(LogFile,signature);
  WriteLittleEndian(LogFile,samplingAmplitude);
  for (int i = 0; i < this->NbrPositions; ++i)
    WriteLittleEndian(LogFile, positions[i]);
  WriteLittleEndian(LogFile,valueExact);
  if (NumAdditionalData>0)
    {
      for (int i=0; i<NumAdditionalData; ++i) 
	this->AdditionalData[i]->WriteMCHistoryData();
    }
#ifdef MCHistoryRecord_TESTING
  cout << "Recorded start at S: "<<samplingAmplitude<<", E: "<<valueExact<<endl;
#endif
  this->LastSampleCount=0; // set counter to zero for a start...
  return true;
}

// read one MC sample back from file, gives back the parameters in call of RecordAcceptedStep
// sampleCount additionally gives the multiplicity of each record
// sampleCount : number of microsteps that the simulation remained at this position (for information only)
// samplingAmplitude : effective samplingAmplitude, including the factor for the multiplicity
// positions : particle positions for this sample
// valueExact : value of the exact wavefunction for this sample
bool MCHistoryRecord::GetMonteCarloStep( int &sampleCount, double &samplingAmplitude, double *positions, Complex &valueExact)
{
  if ( ! (HistoryFile.eof()))
    {
      char signature; 
      ReadLittleEndian(HistoryFile,signature);
      if (signature == 'b') // recognized a new block
	{
	  ReadLittleEndian(HistoryFile,samplingAmplitude);
	  for (int i = 0; i < this->NbrPositions; ++i)
	    ReadLittleEndian(HistoryFile, positions[i]);
	  ReadLittleEndian(HistoryFile,valueExact);
	  if (NumAdditionalData>0)
	    {
	      for (int i=0; i<NumAdditionalData; ++i) 
		this->AdditionalData[i]->ReadMCHistoryData();
	    }
	  else if (SkipAdditional>0)
	    {	      
	      HistoryFile.seekg(SkipAdditional,ios::cur);
	    }
	  ReadLittleEndian(HistoryFile,this->LastSampleCount);
	  this->TotalSampleCount+=this->LastSampleCount;
	  sampleCount=this->LastSampleCount;
#ifdef MCHistoryRecord_TESTING
	  cout << "Read S: "<<samplingAmplitude<<", E: "<<valueExact<<", count: "<<sampleCount<<endl;
#endif
	  return true;
	}
      else if (signature == 'e')
	{
	  // cout << "End of file detected!" << endl;
	  return false;
	}      
    }
  // cout << "Reached end of History file!"<<endl;
  return false;
}

void MCHistoryRecord::RewindHistory()
{
  if (this->RecordMode & AbstractMCHistoryRecord::Reading)
    {
      HistoryFile.seekg(StartPos);
      char c = HistoryFile.peek();
      if (c == 'b') return;
      else
	{
	  if (HistoryFile.eof())
	    {
	      cout << "need to clear eof-bit, here!!!" << endl;
	      HistoryFile.clear();
	    }
	  HistoryFile.seekg(StartPos);
	  c = HistoryFile.peek();
	  if (c == 'b') return;
	  else
	    {
	      cout << "Rewind failed, character was " <<c << endl;
	      exit(-2);
	    }
	}
    }
}
