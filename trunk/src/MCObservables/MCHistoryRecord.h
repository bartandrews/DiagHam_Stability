////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                      Copyright (C) 2007 Gunnar Moeller                     //
//                                                                            //
//                                                                            //
//     class for storing the history of a Monte-Carlo overlap calculation     //
//                                                                            //
//                       last modification : 28/05/2007                       //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef MCHISTORYRECORD_H
#define MCHISTORYRECORD_H

#include "config.h"

#include <fstream>
#include "AbstractMCHistoryData.h"
#include "AbstractMCHistoryRecord.h"
#include "GeneralTools/List.h"
#include "Vector/RealVector.h"
#include "MathTools/Complex.h"

using std::fstream;

class MCHistoryRecord : public AbstractMCHistoryRecord
{
 protected:
  int LastSampleCount;
  int TotalSampleCount;
  int TotalRecordCount;
  int ProjectedStepNum;

  int NbrPositions;
  
  int NumAdditionalData;
  AbstractMCHistoryData **AdditionalData;
  unsigned SkipAdditional;

  ofstream LogFile;
  ifstream HistoryFile;

  int RecordMode;
  std::streampos StartPos;
 public:

  // constructors
  // for recording mode
  // projectedSamplex: expected number of MC steps
  // nbrPositions: Number of coordinates for each call of wavefunction (2*N for a two-D system)
  MCHistoryRecord(int projectedSamples, int nbrPositions, const char* exactFile, const char* samplingDescriptor, char* fileName, List<AbstractMCHistoryData> *additionalData=NULL);

  // to continue aborted calculation
  // fileName: file to continue
  // double & samplingAmplitude, double *positions, Complex &valueExact: describing last recorded positions
  MCHistoryRecord(char* fileName, int nbrPositions, double & samplingAmplitude, double *positions, Complex &valueExact, List<AbstractMCHistoryData> *additionalData=NULL);
  
  // for reading mode
  // Input = name of MCHistoryRecord file to be read
  // nbrPositions = requested number of positions in each record (call with value 0 -> return number of positions in record)
  // additionalData = descriptor of eventual additional fields
  MCHistoryRecord(char *Input, int &nbrPositions, List<AbstractMCHistoryData> *additionalData=NULL);  

  // destructor -> automatically closes LogFile
  virtual ~MCHistoryRecord();
  
  // record rejected step - to be called for rejected Microsteps
  void RecordRejectedStep();

  // record accepted step - to be called for each accepted step, or at every step to be written to file
  bool RecordAcceptedStep(double samplingAmplitude, RealVector &positions, Complex &valueExact);

  // record accepted step - to be called for each accepted step, or at every step to be written to file
  bool RecordInitialPositions( double samplingAmplitude, RealVector &positions, Complex &valueExact);


  // read one MC sample back from file, gives back the parameters in call of RecordAcceptedStep
  // sampleCount additionally gives the multiplicity of each record
  // sampleCount : number of microsteps that the simulation remained at this position (for information only)
  // samplingAmplitude : effective samplingAmplitude, including the factor for the multiplicity
  // positions : particle positions for this sample
  // valueExact : value of the exact wavefunction for this sample
  virtual bool GetMonteCarloStep( int &sampleCount, double & samplingAmplitude, double *positions, Complex &valueExact);

  // rewind in reading mode:
  virtual void RewindHistory();

  // get projected Samples
  virtual int GetProjectedSamples() {return ProjectedStepNum;}
};

#endif // MCHISTORYRECORD_H
