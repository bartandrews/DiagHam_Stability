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


#ifndef SPARSEREALMATRIX_H
#define SPARSEREALMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;
class SparseComplexMatrix;


class SparseRealMatrix;

// create a block diagonal matrix from two matrices 
//
// matrix1 = first matrix (i.e. the one at starting from the first row, first column)
// matrix2 = second matrix
// coefficient1 = optional multiplicative coefficient in front of matrix1
// coefficient2 = optional multiplicative coefficient in front of matrix2
// return value = sparse block diagonal matrix
SparseRealMatrix CreateBlockDiagonalMatrix(const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, double coefficient1 = 1.0, double coefficient2 = 1.0);

// create an  block off-diagonal matrix from two matrices 
//
// matrix1 = first matrix (i.e. the one at starting from the first row)
// matrix2 = second matrix
// coefficient1 = optional multiplicative coefficient in front of matrix1
// coefficient2 = optional multiplicative coefficient in front of matrix2
// return value = sparse block diagonal matrix
SparseRealMatrix CreateBlockOffDiagonalMatrix(const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, double coefficient1 = 1.0, double coefficient2 = 1.0);
  
class SparseRealMatrix : public Matrix
{

  friend class RealMatrix;
  friend class RealVector;
  friend class SparseMatrixMatrixMultiplyOperation;
  friend class TensorProductSparseMatrixHamiltonian;
  friend class TensorProductSparseMatrixSelectedBlockHamiltonian;
  friend class TripleTensorProductSparseMatrixHamiltonian;
  friend class MultipleTensorProductSparseMatrixHamiltonian;
  friend class SparseComplexMatrix;

 protected:

  // number of non-zero matrix elements
  long NbrMatrixElements; 

  // maximum number of matrix elements that can be stored
  long MaximumNbrMatrixElements;
  //number of matrix elements that have to be added if the size of MatrixElements has to be increaed
  long NbrMatrixElementPacketSize;

  // array that contains the matrix elements
  double* MatrixElements; 

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
  SparseRealMatrix();

  // constructor for a sparse matrix without any specific struture
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  SparseRealMatrix(int nbrRow, int nbrColumn);

  // constructor for a sparse matrix without any specific struture but a given number of non-zero matrix elements
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // nbrMatrixElements = number of non-zero matrix elements
  // zero = true if matrix elements have to be set to zero
  SparseRealMatrix(int nbrRow, int nbrColumn, long nbrMatrixElements, bool zero = false);

  // constructor for a sparse matrix knowing how many non-zero elements per row will be required
  //
  // nbrRow = number of rows
  // nbrColumn = number of columns
  // nbrElementPerRow = number of non-zero matrix elements per row
  SparseRealMatrix(int nbrRow, int nbrColumn, int* nbrElementPerRow);

  // copy constructor from a real matrix
  //
  // M = matrix to copy
  SparseRealMatrix(const SparseRealMatrix& M);

  // copy constructor (duplicating all datas)
  //
  // M = matrix to copy
  // accuracy = value below which a matrix element is considered to be zero
  SparseRealMatrix(Matrix& M, double accuracy = 0.0);

  // destructor
  //
  ~SparseRealMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  SparseRealMatrix& operator = (const SparseRealMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();  

  // copy a matrix into another (duplicating data)
  //
  // matrix = matrix to copy
  // return value = reference on current matrix
  SparseRealMatrix& Copy (SparseRealMatrix& matrix);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, double x);

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, double& x) const;

  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  void AddToMatrixElement(int i, int j, double x);

  // lock the sparse matrix such that no additional element can be added (still SetMatrixElement/AddToMatrixElement can still be used to alter the existing matrix elements)
  //
  void LockMatrix();

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

  // check if a sparse matrix is correctly formed
  //
  // return value = true if the sparse matrix is correctly formed
  bool CheckSanity();

  // add two matrices
  //
  // matrix1 = first matrix
  // matrix2 = second matrix
  // return value = sum of the two matrices
  friend SparseRealMatrix operator + (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2);

  // difference of two matrices
  //
  // matrix1 = first matrix
  // matrix2 = second matrix
  // return value = difference of the two matrices
  friend SparseRealMatrix operator - (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend SparseRealMatrix operator * (SparseRealMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend SparseRealMatrix operator * (double x,SparseRealMatrix& M);

  // create the linear combination of two matrices
  //
  // x1 = prefactor of the first matrix
  // matrix1 = first matrix
  // x2 = prefactor of the second matrix
  // matrix2 = second matrix
   // return value = linear combination
  friend SparseRealMatrix SparseRealMatrixLinearCombination(const double& x1, const SparseRealMatrix& matrix1, const double& x2, const SparseRealMatrix& matrix2);

  // create the linear combination of several matrices
  //
  // nbrMatrices = number of matrices that should be added
  // prefactors = array of prefactors
  // matrices = array of matrices
  // return value = linear combination
  friend SparseRealMatrix SparseRealMatrixLinearCombination(int nbrMatrices, double* prefactors, SparseRealMatrix* matrices);
  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  SparseRealMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  SparseRealMatrix& operator /= (double x);
  
  // multiply a matrix to the right by another matrix
  //
  // matrix = matrix used as multiplicator
  // return value = reference on current matrix
  SparseRealMatrix& Multiply (const SparseRealMatrix& matrix);
  
  // multiply two matrices
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // return value = reference on current matrix
  friend SparseRealMatrix Multiply (const SparseRealMatrix& matrix, const SparseRealMatrix& matrix2);

  // multiply two matrices, minimizing the amount of temporary storage
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // return value = reference on current matrix
  friend SparseRealMatrix MemoryEfficientMultiply (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2);

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

  // multiply two matrices, providing all the required temporary arrays
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseRealMatrix Multiply (const SparseRealMatrix& matrix, const SparseRealMatrix& matrix2, 
				    double* tmpMatrixElements, int* tmpColumnIndices, double* tmpElements);

  // multiply two matrices, providing all the required temporary arrays and using architecture optimisation
  //
  // matrix1 = pointer to the left matrix
  // matrix2 = pointer to the right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
  // architecture = pointer to the architecture
  // return value = reference on current matrix
  friend SparseRealMatrix Multiply (SparseRealMatrix* matrix1, SparseRealMatrix* matrix2, 
				    double* tmpMatrixElements, int* tmpColumnIndices, 
				    long nbrTmpMatrixElements, AbstractArchitecture* architecture);

  // multiply two matrices, providing all the required temporary arrays and using architecture optimisation
  //
  // matrix1 = left matrix
  // matrix2 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
  // architecture = pointer to the architecture
  // return value = reference on current matrix
  friend  SparseRealMatrix Multiply (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, 
				     double* tmpMatrixElements, int* tmpColumnIndices, 
				     long nbrTmpMatrixElements, AbstractArchitecture* architecture);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseRealMatrix& Multiply (const SparseRealMatrix& matrix, double* tmpMatrixElements, 
			      int* tmpColumnIndices, double* tmpElements);

  // multiply a matrix to the right by another matrix, providing all the required temporary arrays, extend their capacity if needed
  //
  // matrix = matrix used as multiplicator
  // tmpMatrixElements = reference on the temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = reference on the temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrElements = reference ont the number of elements in tmpMatrixElements and tmpColumnIndices
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  SparseRealMatrix& Multiply (const SparseRealMatrix& matrix, double*& tmpMatrixElements, int*& tmpColumnIndices, 
			      long& nbrElements, double* tmpElements);

  // conjugate the current sparse matrix (M1^+ A M2), assuming A is symmetric
  //
  // matrix1 = left matrix used for the conjugation
  // matrix2 = left matrix used for the conjugation
  // return value = conjugated symmetric matrix
  RealSymmetricMatrix Conjugate (RealMatrix& matrix1, RealMatrix& matrix2);

  // conjugate a matrix
  //
  // matrix1 = left matrix
  // matrix2 = matrix to conjugate
  // matrix3 = right matrix
  // return value = reference on conjugated matrix
  friend SparseRealMatrix Conjugate (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, 
				     const SparseRealMatrix& matrix3);

  // conjugate a matrix, providing all the required temporary arrays and using architecture optimisation
  //
  // matrix1 = left matrix
  // matrix2 = matrix to conjugate
  // matrix3 = right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpElements = temporary array of real numbers, the dimension should be equal to the "matrix" number of rows 
  // return value = reference on current matrix
  friend SparseRealMatrix Conjugate (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, const SparseRealMatrix& matrix3, 
				     double* tmpMatrixElements, int* tmpColumnIndices, double* tmpElements);

  // multiply three matrices, providing all the required temporary arrays
  //
  // matrix1 = pointer to the left matrix
  // matrix2 = pointer to the matrix to conjugate
  // matrix3 = pointer to the right matrix
  // tmpMatrixElements = temporary array of real numbers, the dimension should be equal or higher to the resulting number of non zero elements
  // tmpColumnIndices = temporary array of integers, the dimension should be equal or higher to the resulting number of non zero elements
  // nbrTmpMatrixElements = maximum number of elements available in tmpMatrixElements
  // architecture = pointer to the architecture
  // return value = reference on current matrix
  friend SparseRealMatrix Conjugate (SparseRealMatrix* matrix1, SparseRealMatrix* matrix2, SparseRealMatrix* matrix3, 
				     double* tmpMatrixElements, int* tmpColumnIndices, 
				     long nbrTmpMatrixElements, AbstractArchitecture* architecture);

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

  // matrix-vector multiplication action to the right (i.e. v^t M)
  //
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be stored
  void RightMultiply (RealVector& inputVector, RealVector& outputVector);

  // matrix-vector multiplication action to the right including a global scaling factor (i.e. alpha v^t M)
  //
  // coefficient = global multiplicative coefficient 
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be stored
  void RightMultiply (double coefficient, RealVector& inputVector, RealVector& outputVector);

  // matrix-vector multiplication action to the right (i.e. v^t M), adding the result to another vector
  //
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be added
  void RightAddMultiply (RealVector& inputVector, RealVector& outputVector);

  // matrix-vector multiplication action to the right including a global scaling factor (i.e. alpha v^t M), adding the result to another vector
  //
  // coefficient = global multiplicative coefficient 
  // inputVector = vector that will be multiplied
  // outputVector = vector where the result will be added
  void RightAddMultiply (double coefficient, RealVector& inputVector, RealVector& outputVector);

  // compute the number of non-zero matrix elements (zero having strictly zero square norm)
  //
  // return value = number of non-zero matrix elements
  virtual long ComputeNbrNonZeroMatrixElements();

  // compute the total amount of memory needed to store the sparse matrix
  //
  // return value = amount of memory (in bytes)
  unsigned long GetAllocatedMemory();

  // get the matrix element
  // i = position

  double GetMatrixElement(int i);

  // get the column index
  // i = position

  int GetColumnIndex(int i);

  // get the nbr of matrix elements
  // i = position

  long GetNbrMatrixElements();

  //returns the array with indices of rows

  void GetRowIndices(int* RowIndices);

  // evaluate the real part of the matrix trace
  //
  // return value = real part of the matrix trace 
  double Tr ();

  // evaluate the real part of the matrix partial trace 
  //
  // indices = array of indices that describes the partial trace 
  // nbrIndices = number of indices
  // return value = real part of the matrix partial trace 
  double PartialTr(int* indices, int nbrIndices);

  // compute the tensor product of two sparse matrices (matrix1 x matrix2), and store the result in a sparse matrix
  //
  // matrix1 = reference on the left matrix
  // matrix2 = reference on the right matrix
  // return value = tensor product
  friend SparseRealMatrix TensorProduct (const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2);

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

  // compute the transpose of the current matrix
  //
  // return value = hermitian transposed matrix
  SparseRealMatrix Transpose ();

  // create a block diagonal matrix from two matrices 
  //
  // matrix1 = first matrix (i.e. the one at starting from the first row, first column)
  // matrix2 = second matrix
  // coefficient1 = optional multiplicative coefficient in front of matrix1
  // coefficient2 = optional multiplicative coefficient in front of matrix2
  // return value = sparse block diagonal matrix
  friend SparseRealMatrix CreateBlockDiagonalMatrix(const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, double coefficient1, double coefficient2);
  
  // create an  block off-diagonal matrix from two matrices 
  //
  // matrix1 = first matrix (i.e. the one at starting from the first row)
  // matrix2 = second matrix
  // coefficient1 = optional multiplicative coefficient in front of matrix1
  // coefficient2 = optional multiplicative coefficient in front of matrix2
  // return value = sparse block diagonal matrix
  friend SparseRealMatrix CreateBlockOffDiagonalMatrix(const SparseRealMatrix& matrix1, const SparseRealMatrix& matrix2, double coefficient1, double coefficient2);
  
  // extract a submatrix 
  //
  // nbrRow = number of rows for the submatrix
  // nbrColumn = number of columns for the submatrix
  // rowFlags = array that indicated if a row index is part of the submtrix
  // columnFlags = array that indicated if a column index is part of the submtrix
  // return value = extracted matrix
  SparseRealMatrix ExtractMatrix(int nbrRow, int nbrColumn, bool* rowFlags, bool* columnFlags);
  
  // extract a submatrix 
  //
  // nbrRow = number of rows for the submatrix
  // nbrColumn = number of columns for the submatrix and onto which index it has to be mapped (negative if it should not be kept)
  // rowKeptIndices = array that lists the row indices that have to be kept
  // columnKeptIndices = array that lists the column indices that have to be kept
  // return value = extracted matrix
  SparseRealMatrix ExtractMatrix(int nbrRow, int nbrColumn, int* rowKeptIndices, int* columnKeptIndices);

  // read matrix from a file 
  //
  // file = reference  on the input file stream
  // return value = true if no error occurs
  bool ReadMatrix (ifstream& file);

  // write matrix in a file 
  //
  // file = reference on the output file stream
  // return value = true if no error occurs
  virtual bool WriteMatrix (ofstream& file);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const SparseRealMatrix& P);

  // output the matrix in a sparse display (column formatted output)
  //
  // str = reference on output stream
  // error = numerical accuracy below which a matrix element is considered to be equal to zero (discarded fro sparse matrices)
  // return value = reference on output stream
  virtual ostream& PrintNonZero (ostream& str, double error = MACHINE_PRECISION);

  // output the matrix in a sparse display (column formatted output), using labels for the row and column indices
  //
  // str = reference on output stream
  // rowLabels = array of labels for the row indices
  // columnLabels = array of labels for the column indices
  // error = numerical accuracy below which a matrix element is considered to be equal to zero
  // return value = reference on output stream  
  virtual ostream& PrintNonZero (ostream& str, char** rowLabels, char** columnLabels, double error = MACHINE_PRECISION);

#ifdef USE_OUTPUT

 // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const SparseRealMatrix& P);

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

inline void SparseRealMatrix::GetMatrixElement(int i, int j, double& x) const
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

// lock the sparse matrix such that no additional element can be added (still SetMatrixElement/AddToMatrixElement can still be used to alter the existing matrix elements)
//

inline void SparseRealMatrix::LockMatrix()
{
  this->NbrMatrixElementPacketSize = 0l;
}

// get the matrix element
// i = position

inline double SparseRealMatrix::GetMatrixElement(int i)
{
  return this->MatrixElements[i];
}

// get the column index
// i = position

inline int SparseRealMatrix::GetColumnIndex(int i)
{
  return this->ColumnIndices[i];
}

// get the nbr of matrix elements
// i = position

inline long SparseRealMatrix::GetNbrMatrixElements()
{
  return this->NbrMatrixElements;
}



// find the position of a given column index in a row
//
// index = column index
// minPosition = minimum position of the row in the compressed row storage
// maxPosition = maximum position of the row in the compressed row storage
// return value = position of the column index (-1 if it does not exist)

inline long SparseRealMatrix::FindColumnIndexPosition(int index, long minPosition, long maxPosition) const
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
