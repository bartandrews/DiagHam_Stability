////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract scalar sum  operation                  //
//                                                                            //
//                        last modification : 29/07/2003                      //
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


#ifndef ABSTRACTSCALARSUMOPERATION_H
#define ABSTRACTSCALARSUMOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "MathTools/Complex.h"


class AbstractScalarSumOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  // complex scalar used to store the results
  Complex Scalar;
  // array complex scalars used to store the results in case of multiple scalar sum evaluation
  Complex* Scalars;
  // number of complex scalars in the Scalars array
  int NbrScalars;

 public:
  
  // destructor
  //
  virtual ~AbstractScalarSumOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // get dimension (i.e. Hilbert space dimension, nbr of subdivisions,...)
  // 
  // return value = dimension  
  virtual int GetDimension () = 0;

  // return scalar corresponding to the result of the operation
  //
  // index = optional index in case of multiple scalar sum evaluation
  // return value = reference ont the scalar used to store the operation
  virtual Complex& GetScalar(int index = 0);

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};


// return scalar corresponding to the result of the operation
//
// index = optional index in case of multiple scalar sum evaluation
// return value = reference ont the scalar used to store the operation

inline Complex& AbstractScalarSumOperation::GetScalar(int index)
{
  if (this->Scalars == 0)
    return this->Scalar;
  else
    return this->Scalars[index];
}                                                                                                                                                                                  

#endif
