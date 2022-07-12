////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of FQHE on disk quasihole propagator operation          //
//                                                                            //
//                        last modification : 08/03/2009                      //
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


#ifndef FQHEDISKQUASIHOLEPROPAGATOROPERATION_H
#define FQHEDISKQUASIHOLEPROPAGATOROPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"


class AbstractFunctionBasis;
class RealVector;


class FQHEDiskQuasiholePropagatorOperation: public AbstractScalarSumOperation
{

 protected:

  // pointer to the Hilbert space
  ParticleOnSphere* HilbertSpace;

  // vector corresponding to the Jack polynomial  
  RealVector* JackPolynomial;

 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space to use
  // jackPolynomial = vector corresponding to the Jack polynomial
  FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, RealVector* jackPolynomial);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHEDiskQuasiholePropagatorOperation(const FQHEDiskQuasiholePropagatorOperation& operation);
  
  // destructor
  //
  ~FQHEDiskQuasiholePropagatorOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetLargeIndicesRange (const long& firstComponent, const long& nbrComponent);

  // get dimension (i.e. Hilbert space dimension)
  //
  // return value = dimension
  long GetLargeDimension ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
 protected:

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get dimension (i.e. Hilbert space dimension)
//
// return value = dimension

inline long FQHEDiskQuasiholePropagatorOperation::GetLargeDimension ()
{
  return this->HilbertSpace->GetLargeHilbertSpaceDimension();
}

#endif
