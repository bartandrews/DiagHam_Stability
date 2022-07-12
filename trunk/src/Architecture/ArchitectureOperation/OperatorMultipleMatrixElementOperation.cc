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


#include "config.h"
#include "Architecture/ArchitectureOperation/OperatorMultipleMatrixElementOperation.h"
#include "Operator/AbstractOperator.h"


// constructor 
//
// oper = operators which has to be evaluated
// leftVector = array of vectors to use for the left hand side of the matrix elements 
// rightVector = array of vectors to use for the right hand side of the matrix elements (warning, crossed terms are not evaluated)
// nbrStates = number of states in leftVector (should be equal to the number of states in rightVector)

OperatorMultipleMatrixElementOperation::OperatorMultipleMatrixElementOperation(AbstractOperator* oper, RealVector* leftVectors, RealVector* rightVectors, int nbrStates)
{
  this->RealLeftVectors = leftVectors;
  this->RealRightVectors = rightVectors;
  this->ComplexLeftVectors = 0;
  this->ComplexRightVectors = 0;
  this->NbrStates = nbrStates;
  this->LocalNbrStates = this->NbrStates;
  this->Operator = oper;
  this->OperationType = AbstractArchitectureOperation::OperatorMultipleMatrixElement;  
  this->MatrixElements = new Complex [this->NbrStates];
  for (int i = 0; i < this->NbrStates; ++i)
    this->MatrixElements[i] = 0.0;
  this->FirstState = 0;
}

// constructor 
//
// oper = operators which has to be evaluated
// leftVectors = array of vectors to use for the left hand side of the matrix elements 
// rightVectors = array of vectors to use for the right hand side of the matrix elements (warning, crossed terms are not evaluated)
// nbrStates = number of states in leftVector (should be equal to the number of states in rightVector)

OperatorMultipleMatrixElementOperation::OperatorMultipleMatrixElementOperation(AbstractOperator* oper, ComplexVector* leftVectors, ComplexVector* rightVectors, int nbrStates)
{
  this->RealLeftVectors = 0;
  this->RealRightVectors = 0;
  this->ComplexLeftVectors = leftVectors;
  this->ComplexRightVectors = rightVectors;
  this->NbrStates = nbrStates;
  this->LocalNbrStates = this->NbrStates;
  this->Operator = oper;
  this->OperationType = AbstractArchitectureOperation::OperatorMultipleMatrixElement;  
  this->MatrixElements = new Complex [this->NbrStates];
  for (int i = 0; i < this->NbrStates; ++i)
    this->MatrixElements[i] = 0.0;
  this->FirstState = 0;
}
 
// copy constructor 
//
// operation = reference on operation to copy

OperatorMultipleMatrixElementOperation::OperatorMultipleMatrixElementOperation(const OperatorMultipleMatrixElementOperation& operation)
{
  this->RealRightVectors = operation.RealRightVectors;
  this->RealLeftVectors = operation.RealLeftVectors;
  this->ComplexRightVectors = operation.ComplexRightVectors;
  this->ComplexLeftVectors = operation.ComplexLeftVectors;
  this->OperationType = operation.OperationType;
  this->NbrStates = operation.NbrStates;
  this->LocalNbrStates = operation.LocalNbrStates;
  this->Operator = operation.Operator;
  this->MatrixElements = new Complex [this->NbrStates];
  for (int i = 0; i < this->NbrStates; ++i)
    this->MatrixElements[i] = operation.MatrixElements[i];
  this->FirstState = operation.FirstState;
}
  
// destructor
//

OperatorMultipleMatrixElementOperation::~OperatorMultipleMatrixElementOperation()
{
  if (this->MatrixElements != 0)
    delete[] this->MatrixElements;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* OperatorMultipleMatrixElementOperation::Clone()
{
  return new OperatorMultipleMatrixElementOperation (*this);
}
  

// set the number of operators that have to be locally evaluated
// 
// firstState = index of the first state to evaluate 
// nbrStates =number of states that have to be evaluated bu the local architecture

void OperatorMultipleMatrixElementOperation::SetStateRange(int firstState, int nbrStates)
{
  this->FirstState = firstState;
  this->LocalNbrStates = nbrStates;
}

// apply operation(architecture independent)
//
// return value = true if no error occurs

bool OperatorMultipleMatrixElementOperation::RawApplyOperation()
{
  int LastOperator = this->FirstState + this->LocalNbrStates;
  if (this->RealRightVectors != 0)
    {
      for (int i = this->FirstState; i < LastOperator; ++i)
	this->MatrixElements[i] = this->Operator->MatrixElement(this->RealLeftVectors[i], this->RealRightVectors[i]);
    }
  else
    {
      for (int i = this->FirstState; i < LastOperator; ++i)
	this->MatrixElements[i] = this->Operator->MatrixElement(this->ComplexLeftVectors[i], this->ComplexRightVectors[i]);
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool OperatorMultipleMatrixElementOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrStates / architecture->GetNbrThreads();
  if (Step == 0)
    Step = 1;
  int TmpFirstState = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  OperatorMultipleMatrixElementOperation** TmpOperations = new OperatorMultipleMatrixElementOperation* [architecture->GetNbrThreads()];
  for (int OperatorIndex = 0; OperatorIndex < ReducedNbrThreads; ++OperatorIndex)
    {
      TmpOperations[OperatorIndex] = (OperatorMultipleMatrixElementOperation*) this->Clone();
      if ((TmpFirstState + Step) >= this->NbrStates)
	{
	  TmpOperations[OperatorIndex]->SetStateRange(TmpFirstState, this->NbrStates - TmpFirstState);
	  Step = 0;
	  TmpFirstState = 0;
	}
      else
	{
	  TmpOperations[OperatorIndex]->SetStateRange(TmpFirstState,Step);
	  TmpFirstState += Step;
	}
      architecture->SetThreadOperation(TmpOperations[OperatorIndex], OperatorIndex);
    }
  TmpOperations[ReducedNbrThreads] = (OperatorMultipleMatrixElementOperation*) this->Clone();
  if (Step > 0)
    TmpOperations[ReducedNbrThreads]->SetStateRange(TmpFirstState, this->NbrStates - TmpFirstState);
  else
    TmpOperations[ReducedNbrThreads]->SetStateRange(TmpFirstState, Step);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      for (int j = 0; j < this->NbrStates; ++j)
	this->MatrixElements[j] += TmpOperations[i]->GetMatrixElement(j);
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}
  
