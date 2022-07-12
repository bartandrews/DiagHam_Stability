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


#ifndef HERMITIANMATRIXFROMMATRIXOPERATION_H
#define HERMITIANMATRIXFROMMATRIXOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"


class AbstractOperator;


class HermitianMatrixFromMatrixOperation: public AbstractArchitectureOperation
{

 protected:


  // matrix Q to use (for real numbers)
  RealMatrix RealInitialMatrix;
  // matrix Q to use (for complex numbers)  
  ComplexMatrix ComplexInitialMatrix;

  // matrix where the result Q^t Q will be stored(for real numbers)
  RealSymmetricMatrix RealDestinationMatrix;
  // matrix where the result Q^t Q will be stored(for complex numbers)
  HermitianMatrix ComplexDestinationMatrix;
  

  // index of the first elements to evaluate 
  long FirstElement;
  // number of elements that have to be evaluated by the local architecture
  long LocalNbrElements;
  // total number of elements that have to be evaluated
  long NbrElements;

 public:

  // constructor 
  //
  // initialMatrix = reference on the matrix Q to use
  // destinationMatrix = reference on the matrix where the result Q^t Q will be stored
  HermitianMatrixFromMatrixOperation(RealMatrix& initialMatrix, RealSymmetricMatrix& destinationMatrix);

  // constructor 
  //
  // initialMatrix = reference on the matrix Q to use
  // destinationMatrix = reference on the matrix where the result Q^t Q will be stored
  HermitianMatrixFromMatrixOperation(ComplexMatrix& initialMatrix, HermitianMatrix& destinationMatrix);

  // copy constructor 
  //
  // operation = reference on operation to copy
  HermitianMatrixFromMatrixOperation(const HermitianMatrixFromMatrixOperation& operation);
  
  // destructor
  //
  ~HermitianMatrixFromMatrixOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // set the number of elements that have to be locally evaluated
  // 
  // firstElement = index of the first elements to evaluate 
  // nbrElements = number of elements that have to be evaluated bu the local architecture
  void SetElementRange(long firstElement, long nbrElements);

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

#endif
