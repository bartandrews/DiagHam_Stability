////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of matrix matrix multiplication operation             //
//                                                                            //
//                        last modification : 15/01/2003                      //
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


#ifndef MATRIXMATRIXMULTIPLYOPERATION_H
#define MATRIXMATRIXMULTIPLYOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class Matrix;


class MatrixMatrixMultiplyOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the matrix used as left matrix for the multiplication and where the result will be stored
  Matrix* LeftMatrix;
  // pointer to the matrix used as right matrix for the multiplication
  Matrix* RightMatrix;

 public:
  
  // constructor 
  //
  // leftMatrix = pointer to the matrix used as left matrix for the multiplication and where the result will be stored
  // rightMatrix = pointer to the matrix used as right matrix for the multiplication
  MatrixMatrixMultiplyOperation(Matrix* leftMatrix, Matrix* rightMatrix);

  // copy constructor 
  //
  // operation = reference on operation to copy
  MatrixMatrixMultiplyOperation(const MatrixMatrixMultiplyOperation& operation);
  
  // destructor
  //
  ~MatrixMatrixMultiplyOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get destination matrix
  // 
  // return value = pointer to destination matrix
  Matrix* GetDestinationMatrix ();

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

// get destination matrix
// 
// return value = pointer to destination matrix

inline Matrix* MatrixMatrixMultiplyOperation::GetDestinationMatrix ()
{
  return this->LeftMatrix;
}

#endif
