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


#include "config.h"
#include "Architecture/ArchitectureOperation/MultipleOperatorMatrixElementOperation.h"
#include "Operator/AbstractOperator.h"


// constructor 
//
// oper = array of operators which have to be evaluated
// nbrOpartors = number of operators
// leftVector = reference on the vector to use for the left hand side of the matrix elements
// rightVector = reference on the vector to use for the right hand side of the matrix elements

MultipleOperatorMatrixElementOperation::MultipleOperatorMatrixElementOperation(AbstractOperator** oper, int nbrOpartors, RealVector& leftVector, RealVector& rightVector)
{
  this->RealLeftVector = leftVector;
  this->RealRightVector = rightVector;
  this->NbrOperators = nbrOpartors;
  this->LocalNbrOperators = this->NbrOperators;
  this->Operators = oper;
  this->OperationType = AbstractArchitectureOperation::MultipleOperatorMatrixElement;  
  this->MatrixElements = new Complex [this->NbrOperators];
  for (int i = 0; i < this->NbrOperators; ++i)
    this->MatrixElements[i] = 0.0;
  this->FirstOperator = 0;
}

// constructor 
//
// oper = operator which has to be evaluated  
// nbrOpartors = number of operators
// leftVector = reference on the vector to use for the left hand side of the matrix element
// rightVector = reference on the vector to use for the right hand side of the matrix element

MultipleOperatorMatrixElementOperation::MultipleOperatorMatrixElementOperation(AbstractOperator** oper, int nbrOpartors, ComplexVector& leftVector, ComplexVector& rightVector)
{
  this->ComplexLeftVector = leftVector;
  this->ComplexRightVector = rightVector;
  this->NbrOperators = nbrOpartors;
  this->LocalNbrOperators = this->NbrOperators;
  this->Operators = oper;
  this->OperationType = AbstractArchitectureOperation::MultipleOperatorMatrixElement;
  this->MatrixElements = new Complex [this->NbrOperators];
  for (int i = 0; i < this->NbrOperators; ++i)
    this->MatrixElements[i] = 0.0;
  this->FirstOperator = 0;
}
 
// copy constructor 
//
// operation = reference on operation to copy

MultipleOperatorMatrixElementOperation::MultipleOperatorMatrixElementOperation(const MultipleOperatorMatrixElementOperation& operation)
{
  this->RealRightVector = operation.RealRightVector;
  this->RealLeftVector = operation.RealLeftVector;
  this->ComplexRightVector = operation.ComplexRightVector;
  this->ComplexLeftVector = operation.ComplexLeftVector;
  this->OperationType = operation.OperationType;
  this->NbrOperators = operation.NbrOperators;
  this->LocalNbrOperators = operation.LocalNbrOperators;
  this->Operators = operation.Operators;
  this->MatrixElements = new Complex [this->NbrOperators];
  for (int i = 0; i < this->NbrOperators; ++i)
    this->MatrixElements[i] = operation.MatrixElements[i];
  this->FirstOperator = operation.FirstOperator;
}
  
// destructor
//

MultipleOperatorMatrixElementOperation::~MultipleOperatorMatrixElementOperation()
{
  if (this->MatrixElements != 0)
    delete[] this->MatrixElements;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MultipleOperatorMatrixElementOperation::Clone()
{
  return new MultipleOperatorMatrixElementOperation (*this);
}
  

// set the number of operators that have to be locally evaluated
// 
// firstOperator = index of the first operator to evaluate 
// nbrOperators = number of operators that have to be evaluated bu the local architecture

void MultipleOperatorMatrixElementOperation::SetOperatorRange(int firstOperator, int nbrOperators)
{
  this->FirstOperator = firstOperator;
  this->LocalNbrOperators = nbrOperators;
}

// apply operation(architecture independent)
//
// return value = true if no error occurs

bool MultipleOperatorMatrixElementOperation::RawApplyOperation()
{
  int LastOperator = this->FirstOperator + this->LocalNbrOperators;
  if (this->RealRightVector.GetVectorDimension() > 0)
    {
      for (int i = this->FirstOperator; i < LastOperator; ++i)
	this->MatrixElements[i] = this->Operators[i]->MatrixElement(this->RealLeftVector, this->RealRightVector);
    }
  else
    {
      for (int i = this->FirstOperator; i < LastOperator; ++i)
	this->MatrixElements[i] = this->Operators[i]->MatrixElement(this->ComplexLeftVector, this->ComplexRightVector);
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool MultipleOperatorMatrixElementOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrOperators / architecture->GetNbrThreads();
  if (Step == 0)
    Step = 1;
  int TmpFirstOperator = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  MultipleOperatorMatrixElementOperation** TmpOperations = new MultipleOperatorMatrixElementOperation* [architecture->GetNbrThreads()];
  for (int OperatorIndex = 0; OperatorIndex < ReducedNbrThreads; ++OperatorIndex)
    {
      TmpOperations[OperatorIndex] = (MultipleOperatorMatrixElementOperation*) this->Clone();
      if ((TmpFirstOperator + Step) >= this->NbrOperators)
	{
	  TmpOperations[OperatorIndex]->SetOperatorRange(TmpFirstOperator, this->NbrOperators - TmpFirstOperator);
	  Step = 0;
	  TmpFirstOperator = 0;
	}
      else
	{
	  TmpOperations[OperatorIndex]->SetOperatorRange(TmpFirstOperator,Step);
	  TmpFirstOperator += Step;
	}
      architecture->SetThreadOperation(TmpOperations[OperatorIndex], OperatorIndex);
    }
  TmpOperations[ReducedNbrThreads] = (MultipleOperatorMatrixElementOperation*) this->Clone();
  if (Step > 0)
    TmpOperations[ReducedNbrThreads]->SetOperatorRange(TmpFirstOperator, this->NbrOperators - TmpFirstOperator);
  else
    TmpOperations[ReducedNbrThreads]->SetOperatorRange(TmpFirstOperator, Step);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      for (int j = 0; j < this->NbrOperators; ++j)
	this->MatrixElements[j] += TmpOperations[i]->GetMatrixElement(j);
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}
  
