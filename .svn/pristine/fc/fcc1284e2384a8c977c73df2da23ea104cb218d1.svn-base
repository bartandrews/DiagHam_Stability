////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of add linear combination operation                //
//                                                                            //
//                        last modification : 26/05/2003                      //
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
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "Architecture/SMPArchitecture.h"

#include <sys/time.h>
#include <stdlib.h>


// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = array containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, ComplexVector* sourceVector, 
									    int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = sourceVector;
  this->SourceVectorByPointers = 0;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = array containing pointers to the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, ComplexVector** sourceVector, 
									    int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = sourceVector;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = matrix containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexMatrix& sourceVector, 
									   int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = 0;
  this->SourceVectorMatrix = sourceVector;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = array containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, ComplexVector* sourceVector, 
									    int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = sourceVector;
  this->SourceVectorByPointers = 0;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = 0;
  this->RealCoefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = array containing poiinters to the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, ComplexVector** sourceVector, 
									    int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = sourceVector;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = 0;
  this->RealCoefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = matrix containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(ComplexVector* destinationVector, ComplexMatrix& sourceVector, 
									   int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = 0;
  this->SourceVectorMatrix = sourceVector;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = 0;
  this->RealCoefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = array containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, RealVector* sourceVector, 
									    int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = 0;
  this->RealSourceVector = sourceVector;
  this->RealSourceVectorByPointers = 0;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = array containing pointers to the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation (ComplexVector* destinationVector, RealVector** sourceVector, 
									    int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = 0;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = matrix containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(ComplexVector* destinationVector, RealMatrix& sourceVector, 
									   int nbrVector, Complex* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorByPointers = 0;
  this->RealSourceVector = 0;
  this->RealSourceVectorByPointers = 0;
  this->RealSourceVectorMatrix = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->RealCoefficients = 0;
  this->DestinationVector = destinationVector;  
  this->ExecutionTime = 0.0;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
}

// copy constructor 
//
// operation = reference on operation to copy

AddComplexLinearCombinationOperation::AddComplexLinearCombinationOperation(const AddComplexLinearCombinationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->NbrVector = operation.NbrVector;
  this->Coefficients = operation.Coefficients;
  this->RealCoefficients = operation.RealCoefficients;
  this->OperationType = AbstractArchitectureOperation::AddComplexLinearCombination;
  this->SourceVector = operation.SourceVector;
  this->SourceVectorByPointers = operation.SourceVectorByPointers;
  this->RealSourceVector = operation.RealSourceVector;
  this->RealSourceVectorByPointers = operation.RealSourceVectorByPointers;
  this->DestinationVector = operation.DestinationVector;  
  this->SourceVectorMatrix = operation.SourceVectorMatrix;
  this->RealSourceVectorMatrix = operation.RealSourceVectorMatrix;
  this->ExecutionTime = operation.ExecutionTime;
}
  
// destructor
//

AddComplexLinearCombinationOperation::~AddComplexLinearCombinationOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void AddComplexLinearCombinationOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* AddComplexLinearCombinationOperation::Clone()
{
  return new AddComplexLinearCombinationOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool AddComplexLinearCombinationOperation::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);
  if ((this->RealSourceVector != 0) || (this->RealSourceVectorByPointers != 0) || (this->RealSourceVectorMatrix.GetNbrRow() != 0))
    {
      if (this->RealSourceVector != 0)
	for (int i = 0; i < this->NbrVector; ++i)
	  {
	    this->DestinationVector->AddLinearCombination(Coefficients[i], (this->RealSourceVector[i]), this->FirstComponent, 
							  this->NbrComponent);
	  }
      else
	if (this->RealSourceVectorByPointers != 0)
	  for (int i = 0; i < this->NbrVector; ++i)
	    {
	      this->DestinationVector->AddLinearCombination(Coefficients[i], *(this->RealSourceVectorByPointers[i]), this->FirstComponent, 
							    this->NbrComponent);
	    }
	else
	  for (int i = 0; i < this->NbrVector; ++i)
	    {
	      this->DestinationVector->AddLinearCombination(Coefficients[i], (this->RealSourceVectorMatrix[i]), this->FirstComponent, 
							    this->NbrComponent);
	    }     
    }
  else
    {
      if (this->RealCoefficients == 0)
	{
	  if (this->SourceVector != 0)
	    for (int i = 0; i < this->NbrVector; ++i)
	      {
		this->DestinationVector->AddLinearCombination(Coefficients[i], (this->SourceVector[i]), this->FirstComponent, 
							      this->NbrComponent);
	      }
	  else
	    if (this->SourceVectorByPointers != 0)
	      for (int i = 0; i < this->NbrVector; ++i)
		{
		  this->DestinationVector->AddLinearCombination(Coefficients[i], *(this->SourceVectorByPointers[i]), this->FirstComponent, 
								this->NbrComponent);
		}
	    else
	      for (int i = 0; i < this->NbrVector; ++i)
		{
		  this->DestinationVector->AddLinearCombination(Coefficients[i], (this->SourceVectorMatrix[i]), this->FirstComponent, 
								this->NbrComponent);
		}
	}
      else
	{
	  if (this->SourceVector != 0)
	    for (int i = 0; i < this->NbrVector; ++i)
	      {
		this->DestinationVector->AddLinearCombination(RealCoefficients[i], (this->SourceVector[i]), this->FirstComponent, 
							      this->NbrComponent);
	      }
	  else
	    if (this->SourceVectorByPointers != 0)
	      for (int i = 0; i < this->NbrVector; ++i)
		{
		  this->DestinationVector->AddLinearCombination(RealCoefficients[i], *(this->SourceVectorByPointers[i]), this->FirstComponent, 
								this->NbrComponent);
		}
	    else
	      for (int i = 0; i < this->NbrVector; ++i)
		{
		  this->DestinationVector->AddLinearCombination(RealCoefficients[i], (this->SourceVectorMatrix[i]), this->FirstComponent, 
								this->NbrComponent);
		}
	}
    }
  gettimeofday (&(TotalEndingTime2), 0);
  this->ExecutionTime = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool AddComplexLinearCombinationOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->DestinationVector->GetVectorDimension() / architecture->GetNbrThreads();
  int FirstComponent = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  AddComplexLinearCombinationOperation** TmpOperations = new AddComplexLinearCombinationOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (AddComplexLinearCombinationOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(FirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      FirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (AddComplexLinearCombinationOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(FirstComponent, this->DestinationVector->GetVectorDimension() - FirstComponent);  
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  if (architecture->VerboseMode() == true)
    {
      char TmpString[512];
      double MinTime = 1e14;
      double MaxTime = 0.0;
      for (int i = 0; i < architecture->GetNbrThreads(); ++i)
	{
	  sprintf (TmpString, "AddComplexLinearCombinationOperation core operation on SMP id %d done in %.3f seconds", i, TmpOperations[i]->ExecutionTime);
	  if (MinTime > TmpOperations[i]->ExecutionTime)
	    {
	      MinTime = TmpOperations[i]->ExecutionTime;
	    }
	  if (MaxTime < TmpOperations[i]->ExecutionTime)
	    {
	      MaxTime = TmpOperations[i]->ExecutionTime;
	    }
	  architecture->AddToLog(TmpString);
	}
      sprintf (TmpString, "AddComplexLinearCombinationOperation core operation min time=%.3f sec, max time=%.3f sec", MinTime, MaxTime);
      architecture->AddToLog(TmpString);
   }
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}

