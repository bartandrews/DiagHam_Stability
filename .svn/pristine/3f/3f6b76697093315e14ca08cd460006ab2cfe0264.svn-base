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


#include "Matrix/Matrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#ifdef USE_HILBERT_SPACE
#include "HilbertSpace/SubspaceSpaceConverter.h"
#endif
#include "GeneralTools/Endian.h"


#include <fstream>

using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//

Matrix::Matrix()
{
  this->NbrRow = 0;
  this->NbrColumn = 0;
  this->TrueNbrRow = 0;
  this->TrueNbrColumn = 0;
  this->MatrixType = 0;  
  this->Dummy = 0.0;  
}

// destructor
//

Matrix::~Matrix()
{
}

// return pointer on a clone matrix (without duplicating datas)
//
// retrun value = pointer on new matrix 

Matrix* Matrix::Clone ()
{
  return new Matrix ();
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void Matrix::SetMatrixElement(int i, int j, double x)
{
  return;
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element

void Matrix::SetMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// add a value to a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void Matrix::AddToMatrixElement(int i, int j, double x)
{
  return;
}

// get a matrix element (real part if complex)
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

void Matrix::GetMatrixElement(int i, int j, double& x) const
{
  return;
}

// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

void Matrix::GetMatrixElement(int i, int j, Complex& x) const
{
  return;
}

// add a value  a matrix element
//
// i = line position
// j = column position
// x = value to add to matrix element

void Matrix::AddToMatrixElement(int i, int j, const Complex& x)
{
  return;
}

// Resize matrix
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void Matrix::Resize (int nbrRow, int nbrColumn)
{
  if (nbrRow > this->TrueNbrRow)
    this->TrueNbrRow = nbrRow;
  if (nbrColumn > this->TrueNbrColumn)
    this->TrueNbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
}

// Resize matrix and set to zero all elements that have been added
//
// nbrRow = new number of rows
// nbrColumn = new number of columns

void Matrix::ResizeAndClean (int nbrRow, int nbrColumn)
{
  if (nbrRow > this->TrueNbrRow)
    this->TrueNbrRow = nbrRow;
  if (nbrColumn > this->TrueNbrColumn)
    this->TrueNbrColumn = nbrColumn;
  this->NbrRow = nbrRow;
  this->NbrColumn = nbrColumn;
}

#ifdef USE_HILBERT_SPACE
// project matrix into a given subspace
//
// subspace = reference on subspace structure
// return value = pointer to projected matrix

Matrix* Matrix::Project (SubspaceSpaceConverter& subspace)
{
  return 0;
}
#endif

// return refernce on real part of a given matrix element
//
// i = line position
// j = column position
// return value = reference on real part

double& Matrix::operator () (int i, int j)
{
  return this->Dummy;
}

// conjugate matrix with an unitary real matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = pointer to conjugated matrix

Matrix* Matrix::Conjugate (RealMatrix& UnitaryM)
{
  switch (this->MatrixType)
    {
      case (Matrix::RealElements | Matrix::Symmetric):
	return ((RealSymmetricMatrix&) *this).Conjugate(UnitaryM);
	break;
      case (Matrix::RealElements | Matrix::Antisymmetric):
	return ((RealAntisymmetricMatrix*) this)->Conjugate(UnitaryM);
	break;
    }
  return 0;
}

// conjugate matrix with an unitary block diagonal matrix (Ut M U)
//
// UnitaryM = unitary matrix to use
// return value = pointer to conjugated matrix

Matrix* Matrix::Conjugate (BlockDiagonalMatrix& UnitaryM)
{
  return 0;
}

// evaluate matrix trace
//
// return value = matrix trace 

double Matrix::Tr ()
{
  return 0.0;
}

// evaluate matrix determinant
//
// return value = matrix determinant 

double Matrix::Det ()
{
  return 1.0;
}

// Output Stream overload
//
// str = reference on output stream
// matrix = matrix to print
// return value = reference on output stream

ostream& operator << (ostream& str, const Matrix& matrix)
{
  switch (matrix.MatrixType)
    {
      case (Matrix::ComplexElements):
	str << ((ComplexMatrix&) matrix);
	break;
      case (Matrix::RealElements):
	str << ((RealMatrix&) matrix);
	break;
      case (Matrix::RealElements | Matrix::Symmetric):
	str << ((RealSymmetricMatrix&) matrix);
	break;
      case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
	str << ((RealDiagonalMatrix&) matrix);
	break;
      case (Matrix::RealElements | Matrix::Symmetric | Matrix::TriDiagonal):
	str << ((RealTriDiagonalSymmetricMatrix&) matrix);
	break;
      case (Matrix::RealElements | Matrix::Antisymmetric):
	str << ((RealAntisymmetricMatrix&) matrix);
	break;
      case (Matrix::ComplexElements | Matrix::Hermitian):
	str << ((HermitianMatrix&) matrix);
	break;
    }
  return str;
}


// write matrix in a file 
//
// fileName = name of the file where the matrix has to be stored
// return value = true if no error occurs

bool Matrix::WriteMatrix (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      int TmpType = Matrix::RealElements;
      WriteLittleEndian(File, TmpType);
      WriteLittleEndian(File, this->NbrRow);
      WriteLittleEndian(File, this->NbrColumn);
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    this->GetMatrixElement(i, j, Tmp);
	    WriteLittleEndian(File, Tmp);
	  }
    }
  else
    {
      int TmpType = Matrix::ComplexElements;
      WriteLittleEndian(File, TmpType);
      WriteLittleEndian(File, this->NbrRow);
      WriteLittleEndian(File, this->NbrColumn);
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    this->GetMatrixElement(i, j, Tmp);
	    WriteLittleEndian(File, Tmp.Re);
	    WriteLittleEndian(File, Tmp.Im);
	  }
    }
  File.close();
  return true;
}

// write matrix in a file in ascii mode
//
// fileName = name of the file where the matrix has to be stored
// return value = true if no error occurs

bool Matrix::WriteAsciiMatrix (char* fileName)
{
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      int ReducedNbrColumn = this->NbrColumn - 1;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < ReducedNbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      cout << Tmp << " ";
	    }
	  this->GetMatrixElement(i, ReducedNbrColumn, Tmp);
	  cout << Tmp << endl;
	}      
    }
  else
    {
      Complex Tmp;
      int ReducedNbrColumn = this->NbrColumn - 1;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < ReducedNbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      cout << Tmp.Re << " " << Tmp.Im << " ";
	    }
	  this->GetMatrixElement(i, ReducedNbrColumn, Tmp);
	  cout << Tmp.Re << " " << Tmp.Im << endl;
	}      
    }
  File.close();
  return true;
}

// read matrix from a file 
//
// fileName = name of the file where the matrix has to be read
// return value = true if no error occurs

bool Matrix::ReadMatrix (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  int TmpType = Matrix::RealElements;
  File.read ((char*) &(TmpType), sizeof(int));
  if (((this->MatrixType & TmpType & Matrix::RealElements) == 0) && ((this->MatrixType & TmpType & Matrix::ComplexElements) == 0))
    {
      File.close();
      return false;
    }
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      int TmpNbrRow;
      int TmpNbrColumn;
      ReadLittleEndian(File, TmpNbrRow);
      ReadLittleEndian(File, TmpNbrColumn);
      this->Resize(TmpNbrRow, TmpNbrColumn);
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    ReadLittleEndian(File, Tmp);
	    this->SetMatrixElement(i, j, Tmp);
	  }
    }
  else
    {
      int TmpNbrRow;
      int TmpNbrColumn;
      ReadLittleEndian(File, TmpNbrRow);
      ReadLittleEndian(File, TmpNbrColumn);
      this->Resize(TmpNbrRow, TmpNbrColumn);
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    ReadLittleEndian(File, Tmp.Re);
	    ReadLittleEndian(File, Tmp.Im);
	    this->SetMatrixElement(i, j, Tmp);
	  }
    }
  File.close();
  return true;
}

