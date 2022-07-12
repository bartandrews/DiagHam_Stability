////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of operation for the core part of                  //
//                         vector tensor multplication                        //
//                                                                            //
//                        last modification : 20/07/2013                      //
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
#include "Architecture/ArchitectureOperation/VectorTensorMultiplicationCoreOperation.h"
#include "Architecture/SMPArchitecture.h"
#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"


// constructor
//
// hamiltonian = pointer to the tensor hamiltonian
// tensorIndex = index of tensor to consider
// vSource = vector to be multiplied

VectorTensorMultiplicationCoreOperation::VectorTensorMultiplicationCoreOperation(TensorProductSparseMatrixHamiltonian* tensorHamiltonian, 
										 int tensorIndex, RealVector& vSource)
{
  this->TensorHamiltonian = tensorHamiltonian;
  this->TensorIndex = tensorIndex;
  this->FirstComponent = 0;
  this->NbrComponent = this->TensorHamiltonian->RightMatrices[this->TensorIndex].GetNbrRow();
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = (long) this->NbrComponent;
  this->OperationType = AbstractArchitectureOperation::GenericHamiltonianParticlePrecalculation;
  this->VectorSource = vSource;
  this->ComplexVectorSource = ComplexVector();
}
  
// constructor
//
// hamiltonian = pointer to the tensor hamiltonian
// tensorIndex = index of tensor to consider
// vSource = vector to be multiplied

VectorTensorMultiplicationCoreOperation::VectorTensorMultiplicationCoreOperation(TensorProductSparseMatrixHamiltonian* tensorHamiltonian, 
										 int tensorIndex, ComplexVector& vSource)
{
  this->TensorHamiltonian = tensorHamiltonian;
  this->TensorIndex = tensorIndex;
  this->FirstComponent = 0;
  this->NbrComponent = this->TensorHamiltonian->RightMatrices[this->TensorIndex].GetNbrRow();
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = (long) this->NbrComponent;
  this->OperationType = AbstractArchitectureOperation::GenericHamiltonianParticlePrecalculation;
  this->ComplexVectorSource = vSource;
  this->VectorSource = RealVector();
}
  
// copy constructor
//
// operation = reference on operation to copy

VectorTensorMultiplicationCoreOperation::VectorTensorMultiplicationCoreOperation(const VectorTensorMultiplicationCoreOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::GenericHamiltonianParticlePrecalculation;
  this->TensorHamiltonian = operation.TensorHamiltonian;
  this->TensorIndex = operation.TensorIndex;
  this->VectorSource = operation.VectorSource;
  this->ComplexVectorSource = operation.ComplexVectorSource;
}
  
// destructor
//

VectorTensorMultiplicationCoreOperation::~VectorTensorMultiplicationCoreOperation()
{
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* VectorTensorMultiplicationCoreOperation::Clone()
{
  return new VectorTensorMultiplicationCoreOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool VectorTensorMultiplicationCoreOperation::RawApplyOperation()
{
  if (this->VectorSource.GetVectorDimension() > 0)
    {
      this->TensorHamiltonian->LowLevelAddMultiplyTensorCore(this->TensorIndex, this->TensorHamiltonian->TemporaryArray, this->VectorSource, 
							     this->FirstComponent, this->NbrComponent);
    }
  else
    {
      this->TensorHamiltonian->LowLevelAddMultiplyTensorCore(this->TensorIndex, this->TensorHamiltonian->ComplexTemporaryArray, 
							     this->ComplexVectorSource, 
							     this->FirstComponent, this->NbrComponent);
    }
  return true;
}
