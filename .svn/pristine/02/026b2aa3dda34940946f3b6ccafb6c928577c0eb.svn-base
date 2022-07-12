////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hermitian matrix built from Q^t Q  operation         //
//                                                                            //
//                        last modification : 22/12/2016                      //
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
#include "Architecture/ArchitectureOperation/HermitianMatrixFromMatrixOperation.h"
#include "Operator/AbstractOperator.h"


// constructor 
//
// initialMatrix = reference on the matrix Q to use
// destinationMatrix = reference on the matrix where the result Q^t Q will be stored

HermitianMatrixFromMatrixOperation::HermitianMatrixFromMatrixOperation(RealMatrix& initialMatrix, RealSymmetricMatrix& destinationMatrix)
{
  this->RealInitialMatrix = initialMatrix;
  this->RealDestinationMatrix = destinationMatrix;
  this->NbrElements = (((long) initialMatrix.GetNbrRow()) * (((long) initialMatrix.GetNbrRow()) + 1l)) / 2l;
  this->LocalNbrElements = this->NbrElements;
  this->OperationType = AbstractArchitectureOperation::HermitianMatrixFromMatrix;  
  this->FirstElement = 0;
}

// constructor 
//
// initialMatrix = reference on the matrix Q to use
// destinationMatrix = reference on the matrix where the result Q^t Q will be stored

HermitianMatrixFromMatrixOperation::HermitianMatrixFromMatrixOperation(ComplexMatrix& initialMatrix, HermitianMatrix& destinationMatrix)
{
  this->ComplexInitialMatrix = initialMatrix;
  this->ComplexDestinationMatrix = destinationMatrix;
  this->NbrElements = (((long) initialMatrix.GetNbrRow()) * (((long) initialMatrix.GetNbrRow()) + 1l)) / 2l;
  this->LocalNbrElements = this->NbrElements;
  this->OperationType = AbstractArchitectureOperation::HermitianMatrixFromMatrix;  
  this->FirstElement = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

HermitianMatrixFromMatrixOperation::HermitianMatrixFromMatrixOperation(const HermitianMatrixFromMatrixOperation& operation)
{
  this->RealInitialMatrix = operation.RealInitialMatrix;
  this->ComplexInitialMatrix = operation.ComplexInitialMatrix;
  this->RealDestinationMatrix = operation.RealDestinationMatrix;
  this->ComplexDestinationMatrix = operation.ComplexDestinationMatrix;
  this->OperationType = operation.OperationType;
  this->NbrElements = operation.NbrElements;
  this->LocalNbrElements = operation.LocalNbrElements;
  this->FirstElement = operation.FirstElement;
}
  
// destructor
//

HermitianMatrixFromMatrixOperation::~HermitianMatrixFromMatrixOperation()
{
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* HermitianMatrixFromMatrixOperation::Clone()
{
  return new HermitianMatrixFromMatrixOperation (*this);
}
  

// set the number of elements that have to be locally evaluated
// 
// firstElement = index of the first element to evaluate 
// nbrElements = number of elements that have to be evaluated bu the local architecture

void HermitianMatrixFromMatrixOperation::SetElementRange(long firstElement, long nbrElements)
{
  this->FirstElement = firstElement;
  this->LocalNbrElements = nbrElements;
}

// apply operation(architecture independent)
//
// return value = true if no error occurs

bool HermitianMatrixFromMatrixOperation::RawApplyOperation()
{
  int TmpNbrRow = 0;
  if (this->RealDestinationMatrix.GetNbrRow() > 0)
    {
      TmpNbrRow = this->RealDestinationMatrix.GetNbrRow();
    }
  else
    {
      TmpNbrRow = this->ComplexDestinationMatrix.GetNbrRow();
    }
  int CurrentRowIndex = 0;
  int CurrentColumnIndex = 0;
  long TmpLocalNbrElements = this->FirstElement;
  while (TmpLocalNbrElements > 0l)
    {
      ++CurrentRowIndex;
      if (CurrentRowIndex == TmpNbrRow)
	{
	  ++CurrentColumnIndex;
	  CurrentRowIndex = CurrentColumnIndex;
	}
      --TmpLocalNbrElements;
    }
  TmpLocalNbrElements = this->LocalNbrElements;
  if (this->RealDestinationMatrix.GetNbrRow() > 0)
    {
      while (TmpLocalNbrElements > 0l)
	{
	  this->RealDestinationMatrix.SetMatrixElement(CurrentRowIndex, CurrentColumnIndex, 
						       this->RealInitialMatrix[CurrentRowIndex] * this->RealInitialMatrix[CurrentColumnIndex]);
	  ++CurrentRowIndex;
	  if (CurrentRowIndex == TmpNbrRow)
	    {
	      ++CurrentColumnIndex;
	      CurrentRowIndex = CurrentColumnIndex;
	    }
	  --TmpLocalNbrElements;
	}
    }
  else
    {
      while (TmpLocalNbrElements > 0l)
	{
	  this->ComplexDestinationMatrix.SetMatrixElement(CurrentRowIndex, CurrentColumnIndex, 
							  this->ComplexInitialMatrix[CurrentRowIndex] * this->ComplexInitialMatrix[CurrentColumnIndex]);
	  ++CurrentRowIndex;
	  if (CurrentRowIndex == TmpNbrRow)
	    {
	      ++CurrentColumnIndex;
	      CurrentRowIndex = CurrentColumnIndex;
	    }
	  --TmpLocalNbrElements;
	}
    }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool HermitianMatrixFromMatrixOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->NbrElements / ((long) architecture->GetNbrThreads());
  if (Step == 0l)
    Step = 1l;
  long TmpFirstElement = 0l;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  HermitianMatrixFromMatrixOperation** TmpOperations = new HermitianMatrixFromMatrixOperation* [architecture->GetNbrThreads()];
  for (int ElementIndex = 0; ElementIndex < ReducedNbrThreads; ++ElementIndex)
    {
      TmpOperations[ElementIndex] = (HermitianMatrixFromMatrixOperation*) this->Clone();
      if ((TmpFirstElement + Step) >= this->NbrElements)
	{
	  TmpOperations[ElementIndex]->SetElementRange(TmpFirstElement, this->NbrElements - TmpFirstElement);
	  Step = 0l;
	  TmpFirstElement = 0l;
	}
      else
	{
	  TmpOperations[ElementIndex]->SetElementRange(TmpFirstElement,Step);
	  TmpFirstElement += Step;
	}
      architecture->SetThreadOperation(TmpOperations[ElementIndex], ElementIndex);
    }
  TmpOperations[ReducedNbrThreads] = (HermitianMatrixFromMatrixOperation*) this->Clone();
  if (Step > 0l)
    TmpOperations[ReducedNbrThreads]->SetElementRange(TmpFirstElement, this->NbrElements - TmpFirstElement);
  else
    TmpOperations[ReducedNbrThreads]->SetElementRange(TmpFirstElement, Step);
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  delete[] TmpOperations;
  return true;
}
  
