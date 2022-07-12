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
//                        last modification : 24/10/2002                      //
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
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Vector/RealVector.h"


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVector = vector to be multiplied by the hamiltonian
// destinationVector = vector where the result has to be stored

AddRealLinearCombinationOperation::AddRealLinearCombinationOperation (RealVector* destinationVector, RealVector* sourceVector, 
								      int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::AddRealLinearCombination;
}

// constructor 
//
// destinationVector = vector to which the linear combination has to be added
// sourceVector = matrix containing the vectors that have to be added
// nbrVector = number of vector that have to be added
// coefficients = coefficient of the linear combination

AddRealLinearCombinationOperation::AddRealLinearCombinationOperation(RealVector* destinationVector, RealMatrix& sourceVector, int nbrVector, double* coefficients)
{
  this->FirstComponent = 0;
  this->NbrComponent = destinationVector->GetVectorDimension();
  this->SourceVector = 0;
  this->SourceVectorMatrix = sourceVector;
  this->NbrVector = nbrVector;
  this->Coefficients = coefficients;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::AddRealLinearCombination;  
}

// copy constructor 
//
// operation = reference on operation to copy

AddRealLinearCombinationOperation::AddRealLinearCombinationOperation(const AddRealLinearCombinationOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->NbrVector = operation.NbrVector;
  this->Coefficients = operation.Coefficients;
  this->OperationType = AbstractArchitectureOperation::AddRealLinearCombination;
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;  
  this->SourceVectorMatrix = operation.SourceVectorMatrix;
}
  
// destructor
//

AddRealLinearCombinationOperation::~AddRealLinearCombinationOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void AddRealLinearCombinationOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* AddRealLinearCombinationOperation::Clone()
{
  return new AddRealLinearCombinationOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool AddRealLinearCombinationOperation::ApplyOperation()
{
  if (this->SourceVector != 0)
    for (int i = 0; i < this->NbrVector; ++i)
      {
	this->DestinationVector->AddLinearCombination(Coefficients[i], (this->SourceVector[i]), this->FirstComponent, 
						      this->NbrComponent);
      }
  else
    for (int i = 0; i < this->NbrVector; ++i)
      {
	this->DestinationVector->AddLinearCombination(Coefficients[i], (this->SourceVectorMatrix[i]), this->FirstComponent, 
						      this->NbrComponent);
      }
  return true;
}

