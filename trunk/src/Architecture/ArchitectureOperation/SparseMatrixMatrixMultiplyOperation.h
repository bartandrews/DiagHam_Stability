////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of sparse matrix matrix multiplication operation         //
//                                                                            //
//                        last modification : 04/11/2012                      //
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


#ifndef SPARSEMATRIXMATRIXMULTIPLYOPERATION_H
#define SPARSEMATRIXMATRIXMULTIPLYOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Matrix/SparseRealMatrix.h"


class SparseMatrixMatrixMultiplyOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the matrix used as left matrix for the multiplication 
  SparseRealMatrix* LeftMatrix;
  // pointer to the matrix used as right matrix for the multiplication
  SparseRealMatrix* RightMatrix;
  // pointer to the matrix used as the middle matrix for the 3 matrix multiplication
  SparseRealMatrix* MiddleMatrix;
  // matrix where the product is stored
  SparseRealMatrix DestinationMatrix;

  // number of matrix elements that is computed by the operation
  long LocalNbrMatrixElements;

  // temporary arrays that are used druing the matrix product
  double* LocalTmpMatrixElements;
  int* LocalTmpColumnIndices;
  // maximum number of matrix elements that can ba stored in LocalTmpMatrixElements;
  long MaxTmpMatrixElements;


 public:
  
  // constructor 
  //
  // leftMatrix = pointer to the matrix used as left matrix for the multiplication
  // rightMatrix = pointer to the matrix used as right matrix for the multiplication
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
  SparseMatrixMatrixMultiplyOperation(SparseRealMatrix* leftMatrix, SparseRealMatrix* rightMatrix,
				      double* tmpMatrixElements, int* tmpColumnIndices, 
				      long nbrTmpMatrixElements);

  // constructor for a product of three matrices 
  //
  // leftMatrix = pointer to the matrix used as left matrix for the multiplication
  // middleMatrix = pointer to the matrix used as the middle matrix for the multiplication
  // rightMatrix = pointer to the matrix used as right matrix for the multiplication
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
  SparseMatrixMatrixMultiplyOperation (SparseRealMatrix* leftMatrix, SparseRealMatrix* middleMatrix, SparseRealMatrix* rightMatrix,
				       double* tmpMatrixElements, int* tmpColumnIndices, 
				       long nbrTmpMatrixElements);

  // copy constructor 
  //
  // operation = reference on operation to copy
  SparseMatrixMatrixMultiplyOperation(const SparseMatrixMatrixMultiplyOperation& operation);
  
  // destructor
  //
  ~SparseMatrixMatrixMultiplyOperation();
  
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
  SparseRealMatrix& GetDestinationMatrix ();

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
 protected:

  // set the temporary arrays required during the matrix product
  //
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  void SetTemporaryArrays (double* tmpMatrixElements, int* tmpColumnIndices);

  // apply operation for mono processor architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture);
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
};

// get destination matrix
// 
// return value = pointer to destination matrix

inline SparseRealMatrix& SparseMatrixMatrixMultiplyOperation::GetDestinationMatrix ()
{
  return this->DestinationMatrix;
}

#endif
