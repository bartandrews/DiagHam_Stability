////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
// class of single operator and mulitple matrix element evaluation operation  //
//                                                                            //
//                        last modification : 10/06/2016                      //
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


#ifndef OPERATORMULTIPLEMATRIXELEMENTOPERATION_H
#define OPERATORMULTIPLEMATRIXELEMENTOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


class AbstractOperator;


class OperatorMultipleMatrixElementOperation: public AbstractArchitectureOperation
{

 protected:

  // vectors to use for the right hand side of the matrix element
  RealVector* RealRightVectors;
  // vectors to use for the left hand side of the matrix element
  RealVector* RealLeftVectors;  

  // vectors to use for the right hand side of the matrix element
  ComplexVector* ComplexRightVectors;
  // vectors to use for the left hand side of the matrix element
  ComplexVector* ComplexLeftVectors;  

  // number of of states in RealLeftVectors/ComplexLeftVectors
  int NbrStates;

  // operator which has to be evaluated
  AbstractOperator* Operator;

  // index of the first state to evaluate 
  int FirstState;
  // number of states that have to be evaluated bu the local architecture
  int LocalNbrStates;

  // array where all matrix elements will be stored
  Complex* MatrixElements;

 public:

  // constructor 
  //
  // oper = operators which has to be evaluated
  // leftVectors = array of vectors to use for the left hand side of the matrix elements 
  // rightVectors = array of vectors to use for the right hand side of the matrix elements (warning, crossed terms are not evaluated)
  // nbrStates = number of states in leftVector (should be equal to the number of states in rightVector)
  OperatorMultipleMatrixElementOperation(AbstractOperator* oper, RealVector* leftVectors, RealVector* rightVectors, int nbrStates);

  // constructor 
  //
  // oper = operators which has to be evaluated
  // leftVectors = array of vectors to use for the left hand side of the matrix elements 
  // rightVectors = array of vectors to use for the right hand side of the matrix elements (warning, crossed terms are not evaluated)
  // nbrStates = number of states in leftVector (should be equal to the number of states in rightVector)
  OperatorMultipleMatrixElementOperation(AbstractOperator* oper, ComplexVector* leftVectors, ComplexVector* rightVectors, int nbrStates);

  // copy constructor 
  //
  // operation = reference on operation to copy
  OperatorMultipleMatrixElementOperation(const OperatorMultipleMatrixElementOperation& operation);
  
  // destructor
  //
  ~OperatorMultipleMatrixElementOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get the matrix element produced by the i-th operator
  //
  // index = operator index
  // return value = matrix element
  Complex GetMatrixElement(int index);

  // set the number of states that have to be locally evaluated
  // 
  // firstState = index of the first state to evaluate 
  // nbrStates =number of states that have to be evaluated bu the local architecture
  void SetStateRange(int firstState, int nbrStates);

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

// get the matrix element produced by the i-th state
//
// index = state index
// return value = matrix element

inline Complex OperatorMultipleMatrixElementOperation::GetMatrixElement(int index)
{
  return this->MatrixElements[index];
}

#endif
