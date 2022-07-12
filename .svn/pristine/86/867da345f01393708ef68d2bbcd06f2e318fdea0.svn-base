////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of operator martrix element evaluation operation         //
//                                                                            //
//                        last modification : 12/02/2010                      //
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


#ifndef OPERATORMATRIXELEMENTOPERATION_H
#define OPERATORMATRIXELEMENTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


class AbstractOperator;


class OperatorMatrixElementOperation: public AbstractScalarSumOperation
{

 protected:

  // array of vectors to use for the right hand side of the matrix element
  RealVector RealRightVector;
  // pointer to the vector to use for the left hand side of the matrix element
  RealVector RealLeftVector;  

  // array of vectors to use for the right hand side of the matrix element
  ComplexVector ComplexRightVector;
  // pointer to the vector to use for the left hand side of the matrix element
  ComplexVector ComplexLeftVector;  

  // operator which has to be evaluated
  AbstractOperator* Operator;


 public:

  // constructor 
  //
  // oper = operator which has to be evaluated  
  // leftVector = pointer to the vector to use for the left hand side of the matrix element
  // rightVector = pointer to the vector to use for the right hand side of the matrix element
  // nbrComponents = if positive, override the number of components that have to be computed (usually deduced from the vector dimensions)
  OperatorMatrixElementOperation(AbstractOperator* oper, RealVector& leftVector, RealVector& rightVectors, long nbrComponents = -1l);

  // constructor 
  //
  // oper = operator which has to be evaluated  
  // leftVector = pointer to the vector to use for the left hand side of the matrix element
  // rightVector = pointer to the vector to use for the right hand side of the matrix element
  // nbrComponents = if positive, override the number of components that have to be computed (usually deduced from the vector dimensions)
  OperatorMatrixElementOperation(AbstractOperator* oper, ComplexVector& leftVector, ComplexVector& rightVectors, long nbrComponents = -1l);

  // copy constructor 
  //
  // operation = reference on operation to copy
  OperatorMatrixElementOperation(const OperatorMatrixElementOperation& operation);
  
  // destructor
  //
  ~OperatorMatrixElementOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get dimension (i.e. Hilbert space dimension, nbr of subdivisions,...), return 0 if large number are required
  // 
  // return value = dimension  
  virtual int GetDimension ();

  // get dimension (i.e. Hilbert space dimension, nbr of subdivisions,...) when large number are required
  // 
  // return value = dimension  
  virtual long GetLargeDimension ();

 protected:

  // apply operation(architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

#endif
