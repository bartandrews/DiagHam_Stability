////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of hamiltonian-sparse tensor multiplication operation      //
//                                                                            //
//                        last modification : 21/07/2013                      //
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


#ifndef VECTORSPARSETENSORMULTIPLYOPERATION_H
#define VECTORSPARSETENSORMULTIPLYOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"


class AbstractSparseTensor;
class Vector;


class VectorSparseTensorMultiplyOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the hamiltonian
  TensorProductSparseMatrixHamiltonian* TensorHamiltonian;

  // index of tensor to consider
  int TensorIndex;
 
  // vector to be multiplied by the hamiltonian
  Vector* SourceVector;
  // vector where the result has to be stored
  Vector* DestinationVector;  

  // execution time measured in RawApply
  double ExecutionTime;
  
 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // tensorIndex = index of tensor to consider
  // destinationVector = vector where the result has to be stored
  // fullHilbertSpace = split the workload with respect to the full Hilbert space dimension instead of the auxillary space dimension
  VectorSparseTensorMultiplyOperation(TensorProductSparseMatrixHamiltonian* hamiltonian, int tensorIndex,  
				      Vector* sourceVector, bool fullHilbertSpace = false);

  // copy constructor 
  //
  // operation = reference on operation to copy
  VectorSparseTensorMultiplyOperation(const VectorSparseTensorMultiplyOperation& operation);

  // destructor
  //
  ~VectorSparseTensorMultiplyOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (Vector* DestinationVector);

  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* VectorSparseTensorMultiplyOperation::GetDestinationVector ()
{
  return this->DestinationVector;
}

#endif
