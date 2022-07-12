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


#include "config.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealMatrix.h"


// constructor 
//
// leftMatrix = pointer to the matrix used as left matrix for the multiplication and where the result will be stored
// rightMatrix = pointer to the matrix used as right matrix for the multiplication

MatrixMatrixMultiplyOperation::MatrixMatrixMultiplyOperation (Matrix* leftMatrix, Matrix* rightMatrix)
{
  this->FirstComponent = 0;
  this->NbrComponent = leftMatrix->GetNbrRow();
  this->LeftMatrix = leftMatrix;
  this->RightMatrix = rightMatrix;
  this->OperationType = AbstractArchitectureOperation::MatrixMatrixMultiply;
}

// copy constructor 
//
// operation = reference on operation to copy

MatrixMatrixMultiplyOperation::MatrixMatrixMultiplyOperation(const MatrixMatrixMultiplyOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LeftMatrix = operation.LeftMatrix;
  this->RightMatrix = operation.RightMatrix;
  this->OperationType = AbstractArchitectureOperation::MatrixMatrixMultiply;
}
  
// destructor
//

MatrixMatrixMultiplyOperation::~MatrixMatrixMultiplyOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void MatrixMatrixMultiplyOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* MatrixMatrixMultiplyOperation::Clone()
{
  return new MatrixMatrixMultiplyOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool MatrixMatrixMultiplyOperation::ApplyOperation()
{
  if ((this->LeftMatrix->GetMatrixType() == Matrix::RealElements) && (this->RightMatrix->GetMatrixType() == Matrix::RealElements))
    {
      ((RealMatrix*) this->LeftMatrix)->Multiply(*((RealMatrix*) this->RightMatrix), this->FirstComponent, this->NbrComponent);
      return true;
    }
  return true;
}

