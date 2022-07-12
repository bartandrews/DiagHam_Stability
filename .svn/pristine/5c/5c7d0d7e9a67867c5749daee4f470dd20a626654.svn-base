////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of complex matrix with sparse storage                //
//                                                                            //
//                        last modification : 04/10/2012                      //
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


#ifndef SPARSECOMPLEXMATRIX_H
#define SPARSECOMPLEXMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/SparseRealMatrix.h"

#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class SparseComplexMatrix : public Matrix
{

  friend class ComplexMatrix;
  friend class RealVector;
  friend class ComplexVector;

 protected:

  // number of non-zero mtarix elements
  long NbrMatrixElements; 

  // maximum number of matrix elements that can be stored
  long MaximumNbrMatrixElements;
  //number of matrix elements that have to be added if the size of MatrixElements has to be increaed
  long NbrMatrixElementPacketSize;

  // array that contains the matrix elements
  Complex* MatrixElements; 

  // array that give the column index of each matrix element
  int* ColumnIndices;

  // array that gives the entry point of each row in MatrixElements
  long* RowPointers;
  // array that gives the entry point last element of each row in MatrixElements
  long* RowLastPointers;
  
  // garbage collector flag
  GarbageFlag Flag;

 public:

  // default constructor
  //
  SparseComplexMatrix();

  // constructor for a sparse matrix without any specific struture
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  SparseComplexMatrix(int nbrRow, int nbrColumn);

  // constructor for a sparse matrix without any specific struture but a given number of non-zero matrix elements
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // nbrMatrixElements = number of non-zero matrix elements
  // zero = true if matrix elements have to be set to zero
  SparseComplexMatrix(int nbrRow, int nbrColumn, long nbrMatrixElements, bool zero = false);

  // constructor for a sparse matrix knowing how many non-zero elements per row will be required
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // nbrElementPerRow = number of non-zero matrix elements per row
  SparseComplexMatrix(int nbrRow, int nbrColumn, int* nbrElementPerRow);

  // copy constructor from a complex matrix
  //
  // M = matrix to copy
  SparseComplexMatrix(const SparseComplexMatrix& M);

  // copy constructor from a sparse real matrix (duplicating all data)
  //
  // M = matrix to copy
  SparseComplexMatrix(const SparseRealMatrix& M);

  // copy constructor (duplicating all data)
  //
  // M = matrix to copy
  // accuracy = value below which a matrix element is considered to be zero
  SparseComplexMatrix(Matrix& M, double accuracy = 0.0);

  // destructor
  //
  ~SparseComplexMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  SparseComplexMatrix& operator = (const SparseComplexMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();  

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  SparseComplexMatrix& Copy (SparseComplexMatrix& matrix);

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  SparseComplexMatrix& Copy (SparseRealMatrix& matrix);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, double x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const Complex& x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // real = new real value for matrix element
  // imag = new imaginary value for matrix element
  void SetMatrixElement(int i, int j, double real, double imag);

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, double& x) const;

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, Complex& x) const;

  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, double x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, const Complex& x);

  // Resize matrix
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void Resize (int nbrRow, int nbrColumn);

  // Resize matrix and set to zero all elements that have been added
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  void ResizeAndClean (int nbrRow, int nbrColumn);

  // Set all entries in matrix to zero
  //
  void ClearMatrix ();

  // set matrix to identity 
  //
  void SetToIdentity();

  // add two matrices
  //
  // matrix1 = first matrix
  // matrix2 = second matrix
  // return value = sum of the two matrices
  friend SparseComplexMatrix operator + (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // difference of two matrices
  //
  // matrix1 = first matrix
  // matrix2 = second matrix
  // return value = difference of the two matrices
  friend SparseComplexMatrix operator - (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // create the linear combination of two matrices
  //
  // x1 = prefactor of the first matrix
  // matrix1 = first matrix
  // x2 = prefactor of the second matrix
  // matrix2 = second matrix
   // return value = linear combination
  friend SparseComplexMatrix SparseComplexMatrixLinearCombination(const Complex& x1, const SparseComplexMatrix& matrix1, const Complex& x2, const SparseComplexMatrix& matrix2);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  SparseComplexMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  SparseComplexMatrix& operator /= (double x);
  
  // multiply a matrix to the right by another matrix
  //
  // matrix = matrix used as multiplicator
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseComplexMatrix& matrix);

  // multiply two matrices, providing all the required temporary arrays
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseComplexMatrix Multiply (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2, 
				       Complex* tmpMatrixElements, int* tmpColumnIndices, Complex* tmpElements);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = temporary array of complex numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of complex numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseComplexMatrix& matrix, Complex* tmpMatrixElements, 
				 int* tmpColumnIndices, Complex* tmpElements);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays, extend their capacity if needed
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = reference on the temporary array of complex numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = reference on the temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrElements = reference ont the number of elements in tmpMatrixElements and tmpColumnIndices
  // tmpElements = temporary array of complex numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseComplexMatrix& matrix, Complex*& tmpMatrixElements, int*& tmpColumnIndices, 
				 long& nbrElements, Complex* tmpElements);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = temporary array of complex numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of complex numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseComplexMatrix& Multiply (const SparseRealMatrix& matrix, Complex* tmpMatrixElements, 
				 int* tmpColumnIndices, Complex* tmpElements);

  // multiply two matrices, providing all the required temporary arrays
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseComplexMatrix Multiply (const SparseComplexMatrix& matrix1, const SparseRealMatrix& matrix2, 
				       Complex* tmpMatrixElements, int* tmpColumnIndices, Complex* tmpElements);

  // multiply two matrices, providing all the required temporary arrays
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseComplexMatrix Multiply (const SparseRealMatrix& matrix1, const SparseComplexMatrix& matrix2, 
				       Complex* tmpMatrixElements, int* tmpColumnIndices, Complex* tmpElements);

  // multiply two matrices, minimizing the amount of temporary storage
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // return value = reference on current matrix
  friend SparseComplexMatrix MemoryEfficientMultiply (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // multiply two matrices, minimizing the amount of temporary storage
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // return value = reference on current matrix
  friend SparseComplexMatrix MemoryEfficientMultiply (const SparseRealMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // multiply two matrices, minimizing the amount of temporary storage
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // return value = reference on current matrix
  friend SparseComplexMatrix MemoryEfficientMultiply (const SparseComplexMatrix& matrix1, const SparseRealMatrix& matrix2);

  // conjugate the current sparse matrix (M1^+ A M2), assuming A is hermitian
  //
  // matrix1 = left matrix used for the conjugation
  // matrix2 = left matrix used for the conjugation
  // return value = conjugated hermitian matrix
  HermitianMatrix Conjugate (ComplexMatrix& matrix1, ComplexMatrix& matrix2);

  // conjugate a matrix
  //
  // matrix1 = left matrix
  // matrix2 = matrix to conjugate
  // matrix3 = right matrix
  // return value = reference on conjugated matrix
  friend SparseComplexMatrix Conjugate (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2, 
					const SparseComplexMatrix& matrix3);

  // multiply three matrices, providing all the required temporary arrays
  //
  // matrix1 = left matrix
  // matrix2 = matrix to conjugate
  // matrix3 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseComplexMatrix Conjugate (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2, const SparseComplexMatrix& matrix3, 
					Complex* tmpMatrixElements, int* tmpColumnIndices, Complex* tmpElements);

  // conjugate a matrix
  //
  // matrix1 = left matrix
  // matrix2 = matrix to conjugate
  // matrix3 = right matrix
  // return value = reference on conjugated matrix
  friend SparseComplexMatrix Conjugate (const SparseComplexMatrix& matrix1, const SparseRealMatrix& matrix2, 
					const SparseComplexMatrix& matrix3);

  // multiply three matrices, providing all the required temporary arrays
  //
  // matrix1 = left matrix
  // matrix2 = matrix to conjugate
  // matrix3 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseComplexMatrix Conjugate (const SparseComplexMatrix& matrix1, const SparseRealMatrix& matrix2, const SparseComplexMatrix& matrix3, 
					Complex* tmpMatrixElements, int* tmpColumnIndices, Complex* tmpElements);


  // matrix-vector multiplication action to the right (i.e. v^h M)
  //
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be stored
  void RightMultiply (ComplexVector& inputVector, ComplexVector& outputVector);

  // matrix-vector multiplication action to the right including a global scaling factor (i.e. alpha v^h M)
  //
  // coefficient = global multiplicative coefficient 
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be stored
  void RightMultiply (double coefficient, ComplexVector& inputVector, ComplexVector& outputVector);

  // matrix-vector multiplication action to the right (i.e. v^h M), adding the result to another vector
  //
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be added
  void RightAddMultiply (ComplexVector& inputVector, ComplexVector& outputVector);

  // matrix-vector multiplication action to the right including a global scaling factor (i.e. alpha v^h M), adding the result to another vector
  //
  // coefficient = global multiplicative coefficient 
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be added
  void RightAddMultiply (double coefficient, ComplexVector& inputVector, ComplexVector& outputVector);

  // compute the number of non-zero matrix elements (zero having strictly zero square norm)
  //
  // return value = number of non-zero matrix elements
  virtual long ComputeNbrNonZeroMatrixElements();

  // compute the total amount of memory needed to store the sparse matrix
  //
  // return value = amount of memory (in bytes)
  unsigned long GetAllocatedMemory();

  // evaluate the real part of the matrix trace
  //
  // return value = real part of the matrix trace 
  double Tr ();

  // evaluate the matrix trace
  //
  // return value = matrix trace 
  Complex ComplexTr ();

  // compute the tensor product of two sparse matrices (matrix1 x matrix2), and store the result in a sparse matrix
  //
  // matrix1 = reference on the left matrix
  // matrix2 = reference on the right matrix
  // return value = tensor product
  friend SparseComplexMatrix TensorProduct (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // compute the tensor product of two sparse matrices, applying conjugation on the left one (conj(matrix1) x matrix2), and store the result in a sparse matrix
  //
  // matrix1 = reference on the left matrix
  // matrix2 = reference on the right matrix
  // return value = tensor product
  friend SparseComplexMatrix TensorProductWithConjugation (const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);

  // compute the hermitian transpose of the current matrix
  //
  // return value = hermitian transposed matrix
  SparseComplexMatrix HermitianTranspose ();

  // create a block diagonal matrix from two matrices 
  //
  // matrix1 = first matrix (i.e. the one at starting from the first row, first column)
  // matrix2 = second matrix
  // return value = sparse block diagonal matrix
  friend SparseComplexMatrix CreateBlockDiagonalMatrix(const SparseComplexMatrix& matrix1, const SparseComplexMatrix& matrix2);
  
  // output the matrix in a sparse display (column formatted output)
  //
  // str = reference on output stream
  // error = numerical accuracy below which a matrix element is considered to be equal to zero
  // return value = reference on output stream
  ostream& PrintNonZero (ostream& str, double error = MACHINE_PRECISION);

  // output the matrix in a sparse display (column formatted output), using labels for the row and column indices
  //
  // str = reference on output stream
  // rowLabels = array of labels for the row indices
  // columnLabels = array of labels for the column indices
  // error = numerical accuracy below which a matrix element is considered to be equal to zero
  // return value = reference on output stream  
  virtual ostream& PrintNonZero (ostream& str, char** rowLabels, char** columnLabels, double error = MACHINE_PRECISION);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const SparseComplexMatrix& P);

#ifdef USE_OUTPUT

 // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const SparseComplexMatrix& P);

#endif

 protected:

  // find the position of a given column index in a row
  //
  // index = column index
  // minPosition = minimum position of the row in the compressed row storage
  // maxPosition = maximum position of the row in the compressed row storage
  // return value = position of the column index (-1 if it does not exist)
  long FindColumnIndexPosition(int index, long minPosition, long maxPosition) const;

  // increase the number of matrix elements
  //
  // nbrElements = number of elements to add
  void IncreaseNbrMatrixElements(long nbrElements = 1l);

};

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void SparseComplexMatrix::GetMatrixElement(int i, int j, double& x) const
{
  if (this->RowPointers[i] == -1l)
    {
      x = 0.0;
      return;
    }
  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
  if (TmpIndex == -1l)
    {
      x = 0.0;
      return;      
    }
  x = this->MatrixElements[TmpIndex].Re;
  return;
}

// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void SparseComplexMatrix::GetMatrixElement(int i, int j, Complex& x) const
{
  if (this->RowPointers[i] == -1l)
    {
      x = 0.0;
      return;
    }
  long TmpIndex = this->FindColumnIndexPosition(j, this->RowPointers[i], this->RowLastPointers[i]);
  if (TmpIndex == -1l)
    {
      x = 0.0;
      return;      
    }
  x = this->MatrixElements[TmpIndex];
  return;
}

// find the position of a given column index in a row
//
// index = column index
// minPosition = minimum position of the row in the compressed row storage
// maxPosition = maximum position of the row in the compressed row storage
// return value = position of the column index (-1 if it does not exist)

inline long SparseComplexMatrix::FindColumnIndexPosition(int index, long minPosition, long maxPosition) const
{
  while (minPosition <= maxPosition)
    {
      if (this->ColumnIndices[minPosition] == index)
	return minPosition;
      ++minPosition;
    }
  return -1l;
}

#endif
