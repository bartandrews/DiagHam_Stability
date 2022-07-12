////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of hamiltonian vector multiplication operation           //
//                                                                            //
//                        last modification : 23/10/2002                      //
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
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// sourceVector = vector to be multiplied by the hamiltonian
// destinationVector = vector where the result has to be stored

VectorHamiltonianMultiplyOperation::VectorHamiltonianMultiplyOperation (AbstractHamiltonian* hamiltonian, Vector* sourceVector, Vector* destinationVector)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceVector->GetVectorDimension();
  this->Hamiltonian = hamiltonian;
  this->SourceVector = sourceVector;
  this->DestinationVector = destinationVector;  
  this->OperationType = AbstractArchitectureOperation::VectorHamiltonianMultiply;
}

// copy constructor 
//
// operation = reference on operation to copy

VectorHamiltonianMultiplyOperation::VectorHamiltonianMultiplyOperation(const VectorHamiltonianMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->Hamiltonian = operation.Hamiltonian;
  this->SourceVector = operation.SourceVector;
  this->DestinationVector = operation.DestinationVector;  
  this->OperationType = AbstractArchitectureOperation::VectorHamiltonianMultiply;
}
  
// destructor
//

VectorHamiltonianMultiplyOperation::~VectorHamiltonianMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void VectorHamiltonianMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void VectorHamiltonianMultiplyOperation::SetDestinationVector (Vector* DestinationVector)
{
  this->DestinationVector = DestinationVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* VectorHamiltonianMultiplyOperation::Clone()
{
  return new VectorHamiltonianMultiplyOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool VectorHamiltonianMultiplyOperation::ApplyOperation()
{
  this->Hamiltonian->Multiply((*(this->SourceVector)), (*(this->DestinationVector)), this->FirstComponent, 
			      this->NbrComponent);
  return true;
}

