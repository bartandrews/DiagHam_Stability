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

// put all matrix elements to zero
//

void Matrix::ClearMatrix ()
{
  this->ResizeAndClean(this->NbrRow, this->NbrColumn);
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
// file = reference on the output file stream
// return value = true if no error occurs

bool Matrix::WriteMatrix (ofstream& file)
{
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      int TmpType = Matrix::RealElements;
      WriteLittleEndian(file, TmpType);
      WriteLittleEndian(file, this->NbrRow);
      WriteLittleEndian(file, this->NbrColumn);
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    this->GetMatrixElement(i, j, Tmp);
	    WriteLittleEndian(file, Tmp);
	  }
    }
  else
    {
      int TmpType = Matrix::ComplexElements;
      WriteLittleEndian(file, TmpType);
      WriteLittleEndian(file, this->NbrRow);
      WriteLittleEndian(file, this->NbrColumn);
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    this->GetMatrixElement(i, j, Tmp);
	    WriteLittleEndian(file, Tmp.Re);
	    WriteLittleEndian(file, Tmp.Im);
	  }
    }
  return true;
}

// write matrix in a file 
//
// fileName = name of the file where the matrix has to be stored
// return value = true if no error occurs

bool Matrix::WriteMatrix (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  this->WriteMatrix(File);
  File.close();
  return true;
}

// write matrix in a file in ascii mode
//
// fileName = name of the file where the matrix has to be stored
// gnuplotFlag = if true, write the matrix in a format compatible with gnuplot
// return value = true if no error occurs

bool Matrix::WriteAsciiMatrix (char* fileName, bool gnuplotFlag)
{
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      if (gnuplotFlag == false)
	{
	  int ReducedNbrColumn = this->NbrColumn - 1;
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < ReducedNbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  File << Tmp << " ";
		}
	      this->GetMatrixElement(i, ReducedNbrColumn, Tmp);
	      File << Tmp << endl;
	    }      
	}
      else
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < this->NbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  File << i << " " << j << " " << Tmp << endl;
		}
	      File << endl;
	    }
	}
    }
  else
    {
      Complex Tmp;
      if (gnuplotFlag == false)
	{
	  int ReducedNbrColumn = this->NbrColumn - 1;
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < ReducedNbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  File << Tmp.Re << " " << Tmp.Im << " ";
		}
	      this->GetMatrixElement(i, ReducedNbrColumn, Tmp);
	      File << Tmp.Re << " " << Tmp.Im << endl;
	    }      
	}
      else
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < this->NbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  File << i << " " << j << " " << Tmp.Re << " " << Tmp.Im << endl;
		}
	      File << endl;
	    }
	}
    }
  File.close();
  return true;
}

// read matrix from a file 
//
// file = reference  on the input file stream
// return value = true if no error occurs

bool Matrix::ReadMatrix (ifstream& file)
{
  int TmpType = Matrix::RealElements;
  file.read ((char*) &(TmpType), sizeof(int));
  if (((this->MatrixType & TmpType & Matrix::RealElements) == 0) && ((this->MatrixType & TmpType & Matrix::ComplexElements) == 0))
    {
      file.close();
      cout << "error, trying to load the wrong type of matrix ";
      if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
	{
	  cout << "(complex data for real matrix)" << endl;
	}
      else
	{
	  cout << "(real data for complex matrix)" << endl;
	}
      return false;
    }
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      int TmpNbrRow;
      int TmpNbrColumn;
      ReadLittleEndian(file, TmpNbrRow);
      ReadLittleEndian(file, TmpNbrColumn);
      this->Resize(TmpNbrRow, TmpNbrColumn);
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    ReadLittleEndian(file, Tmp);
	    this->SetMatrixElement(i, j, Tmp);
	  }
    }
  else
    {
      int TmpNbrRow;
      int TmpNbrColumn;
      ReadLittleEndian(file, TmpNbrRow);
      ReadLittleEndian(file, TmpNbrColumn);
      this->Resize(TmpNbrRow, TmpNbrColumn);
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	for (int j = 0; j < this->NbrColumn; ++j)
	  {
	    ReadLittleEndian(file, Tmp.Re);
	    ReadLittleEndian(file, Tmp.Im);
	    this->SetMatrixElement(i, j, Tmp);
	  }
    }
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
  this->ReadMatrix(File);
  File.close();
  return true;
}

// evaluate matrix rank
//
// accuracy = numerical accuracy used to define linearly dependence 
// return value = rank

int Matrix::Rank(double accuracy)
{
  cout << "warning : rank calculation is not implemented for this type of matrix" << endl;
  return -1;
}

// test if a matrix is diagonal
//
// accuracy = numerical accuracy used to define a zero 
// return value = true if the matrix is diagonal

bool Matrix::IsDiagonal(double accuracy)
{
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < i; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp) > accuracy)
		return false;
	    }      
	  for (int j = i + 1; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  else
    {
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < i; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (Norm(Tmp) > accuracy)
		return false;
	    }      
	  for (int j = i + 1; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (Norm(Tmp) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  return true;
}

// test if a matrix is the identity matrix
//
// accuracy = numerical accuracy used to define a zero 
// return value = true if the matrix is diagonal

bool Matrix::IsIdentity(double accuracy)
{
  if (this->IsDiagonal() == false)
    {
      return false;
    }
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  this->GetMatrixElement(i, i, Tmp);
	  if (fabs(1.0 - Tmp) > accuracy)
	    return false;
	}      
      return true;
    }
  else
    {
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  this->GetMatrixElement(i, i, Tmp);
	  if (Norm(1.0 - Tmp) > accuracy)
	    return false;
	}
      return true;
    }
  return true;
}


// test if a matrix is symmetric
//
// accuracy = numerical accuracy used to define a zero 
// return value = true if the matrix is symmetric

bool Matrix::IsSymmetric(double accuracy)
{
  if (this->NbrRow != this->NbrColumn)
    return false;
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp1;
      double Tmp2;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < i; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp1);
	      this->GetMatrixElement(j, i, Tmp2);
	      if (fabs(Tmp1 - Tmp2) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  else
    {
      Complex Tmp1;
      Complex Tmp2;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  this->GetMatrixElement(i, i, Tmp1);
	  if (fabs(Tmp1.Im) > accuracy)
	    return false;
	  for (int j = 0; j < i; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp1);
	      this->GetMatrixElement(j, i, Tmp2);
	      if (Norm(Tmp1 - Tmp2) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  return true;
}

// test if a matrix is hermitian
//
// accuracy = numerical accuracy used to define a zero 
// return value = true if the matrix is diagonal

bool Matrix::IsHermitian(double accuracy)
{
  if (this->NbrRow != this->NbrColumn)
    return false;
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp1;
      double Tmp2;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < i; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp1);
	      this->GetMatrixElement(j, i, Tmp2);
	      if (fabs(Tmp1 - Tmp2) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  else
    {
      Complex Tmp1;
      Complex Tmp2;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  this->GetMatrixElement(i, i, Tmp1);
	  if (fabs(Tmp1.Im) > accuracy)
	    return false;
	  for (int j = 0; j < i; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp1);
	      this->GetMatrixElement(j, i, Tmp2);
	      if (Norm(Tmp1 - Conj(Tmp2)) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  return true;
}

// test if a matrix is real
//
// accuracy = numerical accuracy used to define a zero 
// return value = true if the matrix is real

bool Matrix::IsReal(double accuracy)
{
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      return true;
    }
  else
    {
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp.Im) > accuracy)
		return false;
	    }      
	}
      return true;
    }
  return true;
}

// compute the number of non-zero matrix elements (zero having strictly zero square norm)
//
// return value = number of non-zero matrix elements

long Matrix::ComputeNbrNonZeroMatrixElements()
{
  long TmpCount = 0l;
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (Tmp != 0.0)
	        ++TmpCount;
	    }      
	}
    }
  else
    {
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	        ++TmpCount;
	    }      
	}
    }
  return TmpCount;
}

// write matrix in a file in ascii mode, storing only its non zero elements, 
// first column being the row index, second being the column index, the third is the matrix element real part and the fourth column the matrix element imaginary part
//
// fileName = name of the file where the matrix has to be stored
// error = threshold below which a matrix element is considered to be null
// zeroBased = indices are written starting from zero (i.e. C convention)
// return value = true if no error occurs

bool Matrix::SparseWriteAsciiMatrix (char* fileName, double error, bool zeroBased)
{
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      if (zeroBased == true)
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < this->NbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  if (fabs(Tmp) > error)
		    File << i << " " << j << " " << Tmp << endl;
		}
	    }      
	}
      else
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < this->NbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  if (fabs(Tmp) > error)
		    File << (i + 1) << " " << (j + 1) << " " << Tmp << endl;
		}
	    }      
	}
    }
  else
    {
      Complex Tmp;
      if (zeroBased == true)
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < this->NbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  if (Norm(Tmp) > error)
		    File << i << " " << j << " " << Tmp.Re << " " << Tmp.Im << endl;
		}
	    }      
	}
      else
	{
	  for (int i = 0; i < this->NbrRow; ++i)
	    {
	      for (int j = 0; j < this->NbrColumn; ++j)
		{
		  this->GetMatrixElement(i, j, Tmp);
		  if (Norm(Tmp) > error)
		    File << (i + 1) << " " << (j + 1) << " " << Tmp.Re << " " << Tmp.Im << endl;
		}
	    }      
	}
    }
  File.close();
  return true;
}

// output the matrix in a sparse display (column formatted output)
//
// str = reference on output stream
// error = numerical accuracy below which a matrix element is considered to be equal to zero
// return value = reference on output stream  

ostream& Matrix::PrintNonZero (ostream& str, double error)
{
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp) > error)
		str << i << " " << j << " " << Tmp << endl;
	    }      
	}
    }
  else
    {
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (Norm(Tmp) > error)
		str << i << " " << j << " " << Tmp << endl;
	    }
	}
     }
  return str;  
}

// output the matrix in a sparse display (column formatted output), using labels for the row and column indices
//
// str = reference on output stream
// error = numerical accuracy below which a matrix element is considered to be equal to zero
// return value = reference on output stream  

ostream& Matrix::PrintNonZero (ostream& str, char** rowLabels, char** columnLabels, double error)
{
  if ((this->MatrixType & Matrix::RealElements) == Matrix::RealElements)
    {
      double Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (fabs(Tmp) > error)
		str << rowLabels[i] << " " << columnLabels[j] << " " << Tmp << endl;
	    }      
	}
    }
  else
    {
      Complex Tmp;
      for (int i = 0; i < this->NbrRow; ++i)
	{
	  for (int j = 0; j < this->NbrColumn; ++j)
	    {
	      this->GetMatrixElement(i, j, Tmp);
	      if (Norm(Tmp) > error)
		str << rowLabels[i] << " " << columnLabels[j] << " " << Tmp << endl;
	    }
	}
     }
  return str;  
}

#ifdef __MPI__

// send a matrix to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& Matrix::SendMatrix(MPI::Intracomm& communicator, int id)
{
  return *this;
}

// broadcast a matrix to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// return value = reference on the current matrix

Matrix& Matrix::BroadcastMatrix(MPI::Intracomm& communicator,  int id)
{
  return *this;
}

// broadcast part of matrix to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// firstComponent = index of the column (or row) component (useless if the method is not called by the MPI process which broadcasts the matrix)
// nbrComponent = number of column (or row) (useless if the method is not called by the MPI process which broadcasts the matrix)
// return value = reference on the current matrix
Matrix& Matrix::BroadcastPartialMatrix(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  return *this;
}

// receive a matrix from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current matrix

Matrix& Matrix::ReceiveMatrix(MPI::Intracomm& communicator, int id)
{
  return *this;
}

// add current matrix to the current matrix of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& Matrix::SumMatrix(MPI::Intracomm& communicator, int id)
{
  return *this;
}
 
// reassemble matrix from a scattered one
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current matrix

Matrix& Matrix::ReassembleMatrix(MPI::Intracomm& communicator, int id)
{
  return *this;
}

// create a new matrix on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new matrix 

Matrix* Matrix::BroadcastClone(MPI::Intracomm& communicator, int id)
{
  int Type = this->MatrixType;
  if (id != communicator.Get_rank())
    {
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      switch (Type)
	{
	case (Matrix::RealElements):
	  return new RealMatrix(communicator, id);
	  break;
	case (Matrix::ComplexElements):
	  return new ComplexMatrix(communicator, id);
	  break;
	default:
	  cout << "Matrix::BroadcastClone, matrix type not supported" << endl;
	  return 0;
	}
    }
  return 0;
}

// create a new matrix on given MPI node which is an exact clone of the sent one but with only part of the data
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// firstComponent = index of the first column (or row)
// nbrComponent = number of column (or row) to send
// return value = reference on the current matrix

Matrix& Matrix::SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  cout << "Matrix::SendPartialClone is not supported" << endl;
  return *this;
}

// create a new matrix on given MPI node which is an exact clone of the sent one but with only part of the data
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new matrix 

Matrix* Matrix::ReceivePartialClone(MPI::Intracomm& communicator, int id)
{
  cout << "Matrix::ReceivePartialClone is not supported" << endl;
  return 0;
}

// create a new matrix on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the matrix
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new matrix 

Matrix* Matrix::BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag)
{
  int Type = this->MatrixType;
  communicator.Bcast(&Type, 1, MPI::INT, id);  
  if (id != communicator.Get_rank())
    {
      switch (Type)
	{
	case (Matrix::RealElements):
	  return new RealMatrix(communicator, id);
	  break;
	case (Matrix::ComplexElements):
	  return new ComplexMatrix(communicator, id);
	  break;
	default:
	  return 0;
	}
    }
  return 0;
}

#endif
