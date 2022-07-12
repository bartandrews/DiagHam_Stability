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


#ifndef VECTORTENSORMULTIPLICATIONCOREOPERATION_H
#define VECTORTENSORMULTIPLICATIONCOREOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"



class VectorTensorMultiplicationCoreOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the tensor hamiltonian
  TensorProductSparseMatrixHamiltonian* TensorHamiltonian;

  // vector to be multiplied
  RealVector VectorSource;
  // vector to be multiplied (complex version)
  ComplexVector ComplexVectorSource;

  // index of tensor to consider
  int TensorIndex;

 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the tensor hamiltonian
  // tensorIndex = index of tensor to consider
  // vSource = vector to be multiplied
  VectorTensorMultiplicationCoreOperation(TensorProductSparseMatrixHamiltonian* tensorHamiltonian, 
					  int tensorIndex, RealVector& vSource);

  // constructor 
  //
  // hamiltonian = pointer to the tensor hamiltonian
  // tensorIndex = index of tensor to consider
  // vSource = vector to be multiplied
  VectorTensorMultiplicationCoreOperation(TensorProductSparseMatrixHamiltonian* tensorHamiltonian, 
					  int tensorIndex, ComplexVector& vSource);

  // copy constructor 
  //
  // operation = reference on operation to copy
  VectorTensorMultiplicationCoreOperation(const VectorTensorMultiplicationCoreOperation& operation);
  
  // destructor
  //
  ~VectorTensorMultiplicationCoreOperation();
  
  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  virtual int GetHilbertSpaceDimension ();

  // clone operation
  //
  // return value = pointer to cloned operation
  virtual AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation();

  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int VectorTensorMultiplicationCoreOperation::GetHilbertSpaceDimension ()
{
  return this->TensorHamiltonian->RightMatrices[this->TensorIndex].GetNbrRow();
}

#endif
