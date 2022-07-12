////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of real tridiagoonal symmetric matrix                 //
//                                                                            //
//                        last modification : 16/01/2001                      //
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


#ifndef REALTRIDIAGONALSYMMETRICMATRIX_H
#define REALTRIDIAGONALSYMMETRICMATRIX_H


#include "config.h"
#include "Matrix/Matrix.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif
#include "Vector/RealVector.h"
#ifdef USE_POLYNOMIAL
#include "Polynomial/Polynomial.h"
#endif
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class Complex;
class ComplexVector;
class ComplexMatrix;
class ComplexTriDiagonalHermitianMatrix;
class RealUpperTriangularMatrix;
class RealLowerTriangularMatrix;


class RealTriDiagonalSymmetricMatrix : public Matrix
{

  friend class RealMatrix;
  friend class RealSymmetricMatrix;
  friend class RealBandDiagonalSymmetricMatrix;
  friend class ComplexMatrix;
  friend class HermitianMatrix;
  friend class ComplexTriDiagonalHermitianMatrix;
  friend class ComplexVector;
  friend class RealVector;
  friend ComplexTriDiagonalHermitianMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1,
						       const ComplexTriDiagonalHermitianMatrix& M2);
  friend ComplexTriDiagonalHermitianMatrix operator + (const ComplexTriDiagonalHermitianMatrix& M1,
						       const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexTriDiagonalHermitianMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1,
						       const ComplexTriDiagonalHermitianMatrix& M2);
  friend ComplexTriDiagonalHermitianMatrix operator - (const ComplexTriDiagonalHermitianMatrix& M1,
						       const RealTriDiagonalSymmetricMatrix& M2);
  friend RealMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator + (const RealMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend RealMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const RealMatrix& M2);
  friend RealMatrix operator - (const RealMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator + (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);
  friend ComplexMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const ComplexMatrix& M2);
  friend ComplexMatrix operator - (const ComplexMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

 protected:

  double* DiagonalElements;
  double* UpperDiagonalElements;
  GarbageFlag Flag;

  // dummy variable used if an element outside diagonal or upper(lower) diagonal is requested
  double Dummy;

 public:

  // default constructor
  //
  RealTriDiagonalSymmetricMatrix();

  // constructor for an empty matrix
  //
  // dimension = matrix dimension
  // zero = true if matrix has to be filled with zeros
  RealTriDiagonalSymmetricMatrix(int dimension, bool zero = false);

  // constructor from matrix elements (without duplicating datas)
  //
  // diagonal = pointer to diagonal element array
  // upperDiagonal = pointer to upper diagonal element arra
  // dimension = matrix dimension
  RealTriDiagonalSymmetricMatrix(double* diagonal, double* upperDiagonal, int dimension);

  // copy constructor (without duplicating datas)
  //
  // M = matrix to copy
  RealTriDiagonalSymmetricMatrix(const RealTriDiagonalSymmetricMatrix& M);

  // destructor
  //
  ~RealTriDiagonalSymmetricMatrix();

  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  RealTriDiagonalSymmetricMatrix& operator = (const RealTriDiagonalSymmetricMatrix& M);

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  Matrix* Clone ();

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

  // return refernce on real part of a given matrix element
  //
  // i = line position
  // j = column position
  // return value = reference on real part
  double& operator () (int i, int j);

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

  // copy matrix
  //
  // M = matrix to copy
  // return value = refence on current matrix
  RealTriDiagonalSymmetricMatrix& Copy (RealTriDiagonalSymmetricMatrix& M);

  // add two matrices
  //
  // M1 = first matrix
  // M2 = second matrix
  // return value = sum of the two matrices
  friend RealTriDiagonalSymmetricMatrix operator + (const RealTriDiagonalSymmetricMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

  // substract two matrices
  //
  // M1 = first matrix
  // M2 = matrix to substract to M1
  // return value = difference of the two matrices
  friend RealTriDiagonalSymmetricMatrix operator - (const RealTriDiagonalSymmetricMatrix& M1, const RealTriDiagonalSymmetricMatrix& M2);

  // multiply a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealTriDiagonalSymmetricMatrix operator * (const RealTriDiagonalSymmetricMatrix& M, double x);

  // multiply a matrix by a real number (left multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = product result
  friend RealTriDiagonalSymmetricMatrix operator * (double x, const RealTriDiagonalSymmetricMatrix& M);

  // divide a matrix by a real number (right multiplication)
  //
  // M = source matrix
  // x = real number to use
  // return value = division result
  friend RealTriDiagonalSymmetricMatrix operator / (const RealTriDiagonalSymmetricMatrix& M, double x);

  // add two matrices
  //
  // M = matrix to add to current matrix
  // return value = reference on current matrix
  RealTriDiagonalSymmetricMatrix& operator += (const RealTriDiagonalSymmetricMatrix& M);

  // substract two matrices
  //
  // M = matrix to substract to current matrix
  // return value = reference on current matrix
  RealTriDiagonalSymmetricMatrix& operator -= (const RealTriDiagonalSymmetricMatrix& M);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealTriDiagonalSymmetricMatrix& operator *= (double x);

  // divide a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  RealTriDiagonalSymmetricMatrix& operator /= (double x) ;

  // get a matrix element 
  // 
  // i = Row number
  // j = Column number
  // return value = matrix element M_(i,j)
  double GetElement(int i, int j);

  // access to i-th diagonal element
  // 
  // i = position 
  // return value = reference on i-th diagonal element
  double& DiagonalElement(int i);

  // access to i-th upper diagonal element
  // 
  // i = position 
  // return value = reference on i-th upper diagonal element 
  double& UpperDiagonalElement(int i);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  double Det ();

#ifdef USE_POLYNOMIAL
  // return matrix characteritic equation
  //
  // return value =  reference one polynomial corresponding to matrix characteritic equation  
  Polynomial& CharacteristicEquation();
#endif

  // Diagonalize RealTridiagonal Symmetric Matrix using QL algorithm with implicit shift
  // current matrix is replaced by its corresponding diagonalized matrix
  //
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& Diagonalize(int maxIter = 50);

  // Diagonalize RealTridiagonal Symmetric Matrix using QL algorithm with implicit shift, evaluating eigenvectors in a given base
  // current matrix is replaced by its corresponding diagonalized matrix
  //
  // Q = matrix initialized with corresponding base in which eigenvectors have to be calculated
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on current Matrix  
  RealTriDiagonalSymmetricMatrix& Diagonalize(ComplexMatrix& Q, int maxIter = 50);

  // Diagonalize RealTridiagonal Symmetric Matrix using QL algorithm with implicit shift, evaluating eigenvectors in a given base
  // current matrix is replaced by its corresponding diagonalized matrix
  //
  // Q = matrix initialized with corresponding base in which eigenvectors have to be calculated
  // maxIter = maximum number of iteration to fund an eigenvalue
  // return value = reference on current Matrix  
  RealTriDiagonalSymmetricMatrix& Diagonalize(RealMatrix& Q, int maxIter = 50);

  // find QR factorization
  //
  // Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
  // return value = upper triangular matrix corresponding to the QR factorization of the matrix
  RealUpperTriangularMatrix QRFactorization(RealMatrix& Q);

  // find QL factorization
  //
  // Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
  // return value = lower triangular matrix corresponding to the QL factorization of the matrix  
  RealLowerTriangularMatrix QLFactorization(RealMatrix& Q);

  // find QL factorization with shift (aka M - x 1) 
  //
  // Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
  // shift = shift value
  // return value = lower triangular matrix corresponding to the QL factorization of the matrix
  RealLowerTriangularMatrix QLFactorization(RealMatrix& Q, double shift);

  // find QL factorization and evaluate LQ (ie Qt H Q)
  //
  // Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
  // return value = Qt H Q
  RealTriDiagonalSymmetricMatrix QLConjugaison(RealMatrix& Q);

  // find QL factorization and evaluate LQ (ie Qt H Q), shifting initial matrix diagonal elements and shifting back after evaluating RQ
  //
  // Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
  // shift = shift value
  // return value = Qt H Q
  RealTriDiagonalSymmetricMatrix QLConjugaison(RealMatrix& Q, double shift);

  // find QR factorization and evaluate RQ (ie Qt H Q), shifting initial matrix diagonal elements and shifting back after evaluating RQ
  //
  // Q = matrix initialized with corresponding base in which unitary matrix of QR factorization has to be calculated
  // shift = shift value
  // return value = Qt H Q  
  RealTriDiagonalSymmetricMatrix ConjugateQR(RealMatrix& Q, double shift);
  
  // apply polynomial filter assumuing shifts are exacte shift (i.e. eigenvalues of the initial matrix)
  //
  // Q = unitary matrix encoding the transformation
  // shift = array of shift values
  // nbrShift = number of shifts
  // return value = filtered matrix stored as [[H 0], [0 D]] where D is a diagonal matrix with shifts as element and H a real tridiagonal matrix
  RealTriDiagonalSymmetricMatrix PolynomialFilterWithExactShift(RealMatrix& Q, double* shift, int nbrShift);

  // evaluate a normalized eigenvector for a given eigenvalue (supposing the eigenvalue is non-degenerate)
  //
  // eigenvalue = eigenvalue to use
  // eigenvector = vector where the eigenvector has to be stored
  // return value = reference on eigenvector  
  RealVector& Eigenvector(double eigenvalue, RealVector& eigenvector);

  // evaluate a normalized eigenvector for a given eigenvalue (supposing the eigenvalue is non-degenerate)
  //
  // eigenvalue = eigenvalue to use
  // eigenvector = vector where the eigenvector has to be stored
  // return value = reference on eigenvector
  ComplexVector& Eigenvector(double eigenvalue, ComplexVector& eigenvector);

  // Sort Matrix such that diagnonal elements are sort in increasing order (offdiagonal elements left unchanged)
  //
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& SortMatrixUpOrder();

  // Sort Matrix such that diagnonal elements are sort in increasing order (offdiagonal elements left unchanged) 
  // and apply corresponding transformation to column of a given real matrix 
  //
  // matrix = matrix on which transformation has to be applied
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& SortMatrixUpOrder(RealMatrix& matrix);

  // Sort Matrix such that diagnonal elements are sort in increasing order (offdiagonal elements left unchanged) 
  // and apply corresponding transformation to column of a given complex matrix 
  //
  // matrix = matrix on which transformation has to be applied
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& SortMatrixUpOrder(ComplexMatrix& matrix);

  // Sort Matrix such that diagnonal elements are sort in decreasing order (offdiagonal elements left unchanged)
  //
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& SortMatrixDownOrder();

  // Sort Matrix such that diagnonal elements are sort in decreasing order (offdiagonal elements left unchanged) 
  // and apply corresponding transformation to column of a given real matrix 
  //
  // matrix = matrix on which transformation has to be applied
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& SortMatrixDownOrder(RealMatrix& matrix);

  // Sort Matrix such that diagnonal elements are sort in increasing order and apply corresponding transformation to column
  // of a given complex matrix (offdiagonal elements left unchanged)
  //
  // Q = matrix on which transformation has to be applied
  // return value = reference on current Matrix
  RealTriDiagonalSymmetricMatrix& SortMatrix(ComplexMatrix& Q);

  // Output Stream overload
  //
  // Str = reference on output stream
  // P = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, const RealTriDiagonalSymmetricMatrix& P);

#ifdef USE_OUTPUT

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // P = matrix to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const RealTriDiagonalSymmetricMatrix& P);

#endif

};

#endif
