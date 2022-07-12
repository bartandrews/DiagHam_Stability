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
class LongRationalVector;


class FQHEDiskQuasiholePropagatorOperation: public AbstractScalarSumOperation
{

 protected:

  // pointer to the Hilbert space
  ParticleOnSphere* HilbertSpace;

  // vector corresponding to the Jack polynomial  
  RealVector* JackPolynomial;
  // vector corresponding to the Jack polynomial, rational version  
  LongRationalVector* RationalJackPolynomial;

  // vector corresponding to the second Jack polynomial  (to compute scalar product instead of norm)
  RealVector* JackPolynomial2;
  // vector corresponding to the Jack polynomial, rational version   (to compute scalar product instead of norm)
  LongRationalVector* RationalJackPolynomial2;

 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space to use
  // jackPolynomial = vector corresponding to the Jack polynomial
  FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, RealVector* jackPolynomial);

  // constructor 
  //
  // space = pointer to the Hilbert space to use
  // jackPolynomial = vector corresponding to the Jack polynomial
  FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, LongRationalVector* jackPolynomial);

  // constructor for a scalar product
  //
  // space = pointer to the Hilbert space to use
  // jackPolynomial1 = vector corresponding to the first Jack polynomial
  // jackPolynomial2 = vector corresponding to the second Jack polynomial
  FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, RealVector* jackPolynomial1, RealVector* jackPolynomial2);

  // constructor  for a scalar product
  //
  // space = pointer to the Hilbert space to use
  // jackPolynomial1 = vector corresponding to the first Jack polynomial
  // jackPolynomial2 = vector corresponding to the second Jack polynomial
  FQHEDiskQuasiholePropagatorOperation(ParticleOnSphere* space, LongRationalVector* jackPolynomial1, LongRationalVector* jackPolynomial2);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHEDiskQuasiholePropagatorOperation(const FQHEDiskQuasiholePropagatorOperation& operation);
  
  // destructor
  //
  ~FQHEDiskQuasiholePropagatorOperation();
  
  // get dimension (i.e. Hilbert space dimension, nbr of subdivisions,...), return 0 if large number are required
  // 
  // return value = dimension  
  virtual int GetDimension ();

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
