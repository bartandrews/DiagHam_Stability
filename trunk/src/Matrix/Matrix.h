////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          base class of for matrix                          //
//                                                                            //
//                        last modification : 05/01/2001                      //
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


#ifndef MATRIX_H
#define MATRIX_H


#include "config.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ofstream;
using std::ifstream;

#ifdef __MPI__
#include <mpi.h>
#endif


class Complex;
class RealMatrix;
class BlockDiagonalMatrix;
#ifdef USE_HILBERT_SPACE
class SubspaceSpaceConverter;
#endif


class Matrix
{

  friend class RealVector;
  friend class ComplexVector;

 private:
  
  double Dummy;

 protected:

  int NbrRow;
  int NbrColumn;

  int TrueNbrRow;
  int TrueNbrColumn;

  int MatrixType;

 public:

  enum MatrixProperties
  {
    RealElements = 0x0001,
    ComplexElements = 0x0002,
    LongRationalElements = 0x0004,
    IntegerElements = 0x0008,
    LongIntegerElements = 0x000c,
    Diagonal = 0x0010,
    TriDiagonal = 0x0020,
    Triangular = 0x0040,
    BlockDiagonal = 0x0080,
    Symmetric = 0x0100,
    Antisymmetric = 0x0200,
    Hermitian = 0x0400,
    AntiHermitian = 0x0800,    
    Block = 0x1000,
    Upper = 0x2000,
    Lower = 0x4000,
    BandDiagonal = 0x10000,
    Hessenberg = 0x20000,
    Sparse = 0x100000
  };

  // default constructor
  //
  Matrix();

  // virtual destructor
  //
  virtual ~Matrix();

  // return pointer on a clone matrix (without duplicating datas)
  //
  // retrun value = pointer on new matrix 
  virtual Matrix* Clone ();

#ifdef USE_HILBERT_SPACE
  // project matrix into a given subspace
  //
  // subspace = reference on subspace structure
  // return value = pointer to projected matrix
  virtual Matrix* Project (SubspaceSpaceConverter& subspace);  
#endif

  // get matrix type
  //
  // return value = matrix type 
  virtual int GetMatrixType () const;

  // get number of row
  //
  // return value = number of row
  virtual int GetNbrRow () const;

  // get number of column
  //
  // return value = number of column
  virtual int GetNbrColumn () const;

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  virtual void SetMatrixElement(int i, int j, double x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  virtual void SetMatrixElement(int i, int j, const Complex& x);

  // get a matrix element (real part if complex)
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  virtual void GetMatrixElement(int i, int j, double& x) const;

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  virtual void GetMatrixElement(int i, int j, Complex& x) const;

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
/*  template <class T>
  inline virtual void GetMatrixElement(int i, int j, T & x) const;*/


  // add a value to a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  virtual void AddToMatrixElement(int i, int j, double x);

  // add a value  a matrix element
  //
  // i = line position
  // j = column position
  // x = value to add to matrix element
  virtual void AddToMatrixElement(int i, int j, const Complex& x);

  // Resize matrix
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  virtual void Resize (int nbrRow, int nbrColumn);

  // Resize matrix and set to zero all elements that have been added
  //
  // nbrRow = new number of rows
  // nbrColumn = new number of columns
  virtual void ResizeAndClean (int nbrRow, int nbrColumn);

  // put all matrix elements to zero
  //
  virtual void ClearMatrix ();

  // return reference on real part of a given matrix element
  //
  // i = line position
  // j = column position
  // return value = reference on real part
  virtual double& operator () (int i, int j);

  // conjugate matrix with an unitary real matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  virtual Matrix* Conjugate (RealMatrix& UnitaryM);

  // conjugate matrix with an unitary block diagonal matrix (Ut M U)
  //
  // UnitaryM = unitary matrix to use
  // return value = pointer to conjugated matrix
  virtual Matrix* Conjugate (BlockDiagonalMatrix& UnitaryM);

  // multiply a matrix by a real number
  //
  // x = real number to use
  // return value = reference on current matrix
  //  virtual Matrix& operator *= (double x);

  // evaluate matrix trace
  //
  // return value = matrix trace 
  virtual double Tr ();

  // evaluate matrix determinant
  //
  // return value = matrix determinant 
  virtual double Det ();

  // evaluate matrix rank
  //
  // accuracy = numerical accuracy used to define linearly dependence 
  // return value = rank
  virtual int Rank(double accuracy = MACHINE_PRECISION);

  // test if a matrix is diagonal
  //
  // accuracy = numerical accuracy used to define a zero 
  // return value = true if the matrix is diagonal
  virtual bool IsDiagonal(double accuracy = MACHINE_PRECISION);

  // test if a matrix is the identity matrix
  //
  // accuracy = numerical accuracy used to define a zero 
  // return value = true if the matrix is diagonal
  virtual bool IsIdentity(double accuracy = MACHINE_PRECISION);

  // test if a matrix is symmetric
  //
  // accuracy = numerical accuracy used to define a zero 
  // return value = true if the matrix is symmetric
  virtual bool IsSymmetric(double accuracy = MACHINE_PRECISION);

  // test if a matrix is hermitian
  //
  // accuracy = numerical accuracy used to define a zero 
  // return value = true if the matrix is diagonal
  virtual bool IsHermitian(double accuracy = MACHINE_PRECISION);

  // test if a matrix is real
  //
  // accuracy = numerical accuracy used to define a zero 
  // return value = true if the matrix is real
  virtual bool IsReal(double accuracy = MACHINE_PRECISION);

  // compute the number of non-zero matrix elements (zero having strictly zero square norm)
  //
  // return value = number of non-zero matrix elements
  virtual long ComputeNbrNonZeroMatrixElements();

  // Output Stream overload
  //
  // str = reference on output stream
  // matrix = matrix to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const Matrix& matrix);

  // write matrix in a file 
  //
  // file = reference on the output file stream
  // return value = true if no error occurs
  virtual bool WriteMatrix (ofstream& file);

  // write matrix in a file 
  //
  // fileName = name of the file where the matrix has to be stored
  // return value = true if no error occurs
  virtual bool WriteMatrix (char* fileName);

  // write matrix in a file in ascii mode
  //
  // fileName = name of the file where the matrix has to be stored
  // gnuplotFlag = if true, write the matrix in a format compatible with gnuplot
  // return value = true if no error occurs
  virtual bool WriteAsciiMatrix (char* fileName, bool gnuplotFlag = false);

  // write matrix in a file in ascii mode, storing only its non zero elements, 
  // first column being the row index, second being the column index, the third is the matrix element real part and the fourth column the matrix element imaginary part
  //
  // fileName = name of the file where the matrix has to be stored
  // error = threshold below which a matrix element is considered to be null
  // zeroBased = indices are written starting from zero (i.e. C convention)
  // return value = true if no error occurs
  virtual bool SparseWriteAsciiMatrix (char* fileName, double error = 0.0, bool zeroBased = true);

  // read matrix from a file 
  //
  // file = reference  on the input file stream
  // return value = true if no error occurs
  virtual bool ReadMatrix (ifstream& file);

  // read matrix from a file 
  //
  // fileName = name of the file where the matrix has to be read
  // return value = true if no error occurs
  virtual bool ReadMatrix (char* fileName);

  // output the matrix in a sparse display (column formatted output)
  //
  // str = reference on output stream
  // error = numerical accuracy below which a matrix element is considered to be equal to zero
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

#ifdef __MPI__

  // send a matrix to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current matrix
  virtual Matrix& SendMatrix(MPI::Intracomm& communicator, int id);

  // broadcast a matrix to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // return value = reference on the current matrix
  virtual Matrix& BroadcastMatrix(MPI::Intracomm& communicator,  int id);

  // broadcast part of matrix to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // firstComponent = index of the column (or row) component (useless if the method is not called by the MPI process which broadcasts the matrix)
  // nbrComponent = number of column (or row) (useless if the method is not called by the MPI process which broadcasts the matrix)
  // return value = reference on the current matrix
  virtual Matrix& BroadcastPartialMatrix(MPI::Intracomm& communicator, int id, int firstComponent = 0, int nbrComponent = 0);

  // receive a matrix from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current matrix
  virtual Matrix& ReceiveMatrix(MPI::Intracomm& communicator, int id);

  // add current matrix to the current matrix of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current matrix
  virtual Matrix& SumMatrix(MPI::Intracomm& communicator, int id);

  // reassemble matrix from a scattered one
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current matrix
  virtual Matrix& ReassembleMatrix(MPI::Intracomm& communicator, int id);

  // create a new matrix on each MPI node which is an exact clone of the broadcasted one
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new matrix 
  virtual Matrix* BroadcastClone(MPI::Intracomm& communicator, int id);

  // create a new matrix on given MPI node which is an exact clone of the sent one but with only part of the data
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // firstComponent = index of the first column (or row)
  // nbrComponent = number of column (or row) to send
  // return value = reference on the current matrix
  virtual Matrix& SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent);

  // create a new matrix on given MPI node which is an exact clone of the sent one but with only part of the data
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new matrix 
  virtual Matrix* ReceivePartialClone(MPI::Intracomm& communicator, int id);

  // create a new matrix on each MPI node with same size and same type but non-initialized components
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the matrix
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new matrix 
  virtual Matrix* BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag = false);

#endif

};

// get matrix type
//
// return value = matrix type 

inline int Matrix::GetMatrixType () const
{
  return this->MatrixType;
}

// get number of row
//
// return value = number of row

inline int Matrix::GetNbrRow () const
{
  return this->NbrRow;
}

// get number of column
//
// return value = number of column

inline int Matrix::GetNbrColumn () const
{
  return this->NbrColumn;
}

/*
// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element
template <class T>
inline virtual void Matrix::GetMatrixElement(int i, int j, T & x) const
{
  if( 

}
*/

#endif
