////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                      Copyright (C) 2007 Gunnar Moeller                     //
//                                                                            //
//                                                                            //
//   interface to storing the history of a Monte-Carlo overlap calculation    //
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


#ifndef ABSTRACTMCHISTORYRECORD_H
#define ABSTRACTMCHISTORYRECORD_H

#include "MathTools/Complex.h"

class AbstractMCHistoryRecord
{
 public:
  
  enum
    {
      Recording = 0x01,
      Reading = 0x02
    };

  // destructor -> automatically closes LogFile
  virtual ~AbstractMCHistoryRecord();  

  // read one MC sample back from file, gives back the parameters in call of RecordAcceptedStep
  // sampleCount additionally gives the multiplicity of each record
  virtual bool GetMonteCarloStep( int &sampleCount, double & samplingAmplitude, double *positions, Complex &valueExact) = 0;

  // rewind in reading mode:
  virtual void RewindHistory() = 0;

  // get projected Samples
  virtual int GetProjectedSamples() = 0;
};

#endif // MCHISTORYRECORD_H
