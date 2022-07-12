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


#ifndef MULTIFILEMCHISTORYRECORD_H
#define MULTIFILEMCHISTORYRECORD_H

#include "config.h"

#include <fstream>
#include "AbstractMCHistoryData.h"
#include "AbstractMCHistoryRecord.h"
#include "MCHistoryRecord.h"
#include "GeneralTools/List.h"
#include "Vector/RealVector.h"
#include "MathTools/Complex.h"

using std::fstream;

class MultiFileMCHistoryRecord : public AbstractMCHistoryRecord
{
 protected:
  int ProjectedStepNum;

  int NbrPositions;
  
  int NbrHistoryFiles;
  int NbrPresentHistory;
  
  MCHistoryRecord **InputHistories;

  int RecordMode;

 public:

  enum
    {
      Recording = 0x01,
      Reading = 0x02
    };

  // phantom standard constructor
  MultiFileMCHistoryRecord();

  // to merge a number of records into a single file, (remain open in reading mode)
  // outputFileName: new history file to be generated (NULL => only open for reading)
  // nbrInputFiles: number of input files
  // inputFileNames: history files to be merged
  // nbrPositions: requested number of particle positions (call with 0 to return value given in input files)
  // additionalData: description of additional fields
  MultiFileMCHistoryRecord(char* outputFileName, int nbrInputFiles, char **inputFileNames, int &nbrPositions, List<AbstractMCHistoryData> *additionalData=NULL, bool verbose=false);
    
  // destructor -> automatically closes LogFile
  virtual ~MultiFileMCHistoryRecord();
    
  // read one MC sample back from file, gives back the parameters in call of RecordAcceptedStep
  // sampleCount additionally gives the multiplicity of each record
  virtual bool GetMonteCarloStep( int &sampleCount, double & samplingAmplitude, double *positions, Complex &valueExact);

  // rewind in reading mode:
  virtual void RewindHistory();

  // get projected Samples
  virtual int GetProjectedSamples() {return ProjectedStepNum;}
};

#endif // MULTIFILEMCHISTORYRECORD_H
