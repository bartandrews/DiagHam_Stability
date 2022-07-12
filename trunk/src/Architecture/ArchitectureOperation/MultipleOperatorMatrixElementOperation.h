////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of mulitple operators matrix element evaluation operation     //
//                                                                            //
//                        last modification : 22/02/2013                      //
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


#ifndef MULTIPLEOPERATORMATRIXELEMENTOPERATION_H
#define MULTIPLEOPERATORMATRIXELEMENTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


class AbstractOperator;


class MultipleOperatorMatrixElementOperation: public AbstractArchitectureOperation
{

 protected:

  // vector to use for the right hand side of the matrix element
  RealVector RealRightVector;
  // vector to use for the left hand side of the matrix element
  RealVector RealLeftVector;  

  // vector to use for the right hand side of the matrix element
  ComplexVector ComplexRightVector;
  // vector to use for the left hand side of the matrix element
  ComplexVector ComplexLeftVector;  

  // array of operators which have to be evaluated
  AbstractOperator** Operators;
  // number of operators
  int NbrOperators;

  // index of the first operator to evaluate 
  int FirstOperator;
  // number of operators that have to be evaluated bu the local architecture
  int LocalNbrOperators;

  // array where all matrix elements will be stored
  Complex* MatrixElements;

 public:

  // constructor 
  //
  // oper = array of operators which have to be evaluated
  // nbrOpartors = number of operators
  // leftVector = reference on the vector to use for the left hand side of the matrix elements
  // rightVector = reference on the vector to use for the right hand side of the matrix elements
  MultipleOperatorMatrixElementOperation(AbstractOperator** oper, int nbrOpartors, RealVector& leftVector, RealVector& rightVectors);

  // constructor 
  //
  // oper = operator which has to be evaluated  
  // nbrOpartors = number of operators
  // leftVector = reference on the vector to use for the left hand side of the matrix element
  // rightVector = reference on the vector to use for the right hand side of the matrix element
  MultipleOperatorMatrixElementOperation(AbstractOperator** oper, int nbrOpartors, ComplexVector& leftVector, ComplexVector& rightVectors);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MultipleOperatorMatrixElementOperation(const MultipleOperatorMatrixElementOperation& operation);
  
  // destructor
  //
  ~MultipleOperatorMatrixElementOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get the matrix element produced by the i-th operator
  //
  // index = operator index
  // return value = matrix element
  Complex GetMatrixElement(int index);

  // set the number of operators that have to be locally evaluated
  // 
  // firstOperator = index of the first operator to evaluate 
  // nbrOperators =number of operators that have to be evaluated bu the local architecture
  void SetOperatorRange(int firstOperator, int nbrOperators);

 protected:

  // apply operation(architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation();
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};

// get the matrix element produced by the i-th operator
//
// index = operator index
// return value = matrix element

inline Complex MultipleOperatorMatrixElementOperation::GetMatrixElement(int index)
{
  return this->MatrixElements[index];
}

#endif
