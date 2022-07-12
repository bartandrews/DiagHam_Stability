////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//                                                                            //
//                        last modification : 04/01/2001                      //
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


#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <fstream>
#include <iostream>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//

RealVector::RealVector()
{
  this->VectorType = Vector::RealDatas;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->Components = 0;
  this->VectorId = 0;
}

// constructor for an empty real vector (all coordinates set to zero)
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

RealVector::RealVector(int size, bool zeroFlag)
{
  this->VectorType = Vector::RealDatas;
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Components = new double [this->Dimension + 1]; 
  this->Flag.Initialize();
  this->VectorId = 0;
  if (zeroFlag == true)
    for (int i = 0; i < this->Dimension; i++)
      {
	this->Components[i] = 0.0;
      }
}

// constructor from an array of doubles
//
// array = array of doubles with real in even position and imaginary part in odd position
// size = Vector Dimension
 
RealVector::RealVector(double* array, int size)
{
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Components = array;
  this->Flag.Initialize();
  this->VectorId = 0;
}

// copy constructor
//
// vector = vector to copy
// DuplicateFlag = true if datas have to be duplicated

RealVector::RealVector(const RealVector& vector, bool duplicateFlag)
{
  this->VectorType = Vector::RealDatas;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  if (duplicateFlag == false)
    {
      this->Components = vector.Components;
      this->Flag = vector.Flag;
    }
  else
    {
      if (vector.Dimension > 0)
	{
	  this->Flag.Initialize();
	  this->Components = new double [this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    this->Components[i] = vector.Components[i];
	}
      else
	{
	  this->Components = 0;
	}
    }
}

// copy constructor from a complex vector (keep only real part and datas are duplicated)
//
// vector = vector to copy

RealVector::RealVector(const ComplexVector& vector)
{
  this->VectorType = Vector::RealDatas;
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->VectorId = 0;
  if (this->Dimension > 0)
    {
      this->Components = new double[this->Dimension + 1];
      for (int i = 0; i < this->Dimension; ++i)
	{
	  this->Components[i] = vector.Components[i].Re;
	}
    }
  else
    this->Components = 0;
  this->Flag.Initialize();
}

// copy constructor from a vector (duplicate datas if necessary)
//
// vector = vector to copy

RealVector::RealVector(const Vector& vector)
{
  this->VectorType = Vector::RealDatas;
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->VectorId = vector.VectorId;
  if (vector.VectorType == Vector::RealDatas)
    {
      this->VectorType = Vector::RealDatas;
      this->Components = ((RealVector&) vector).Components;
      this->Flag = ((RealVector&) vector).Flag;
    }
  else
    if (vector.VectorType == Vector::ComplexDatas)
      {
	if (this->Dimension > 0)
	  {
	    this->Components = new double[this->Dimension + 1];
	    for (int i = 0; i < this->Dimension; ++i)
	      {
		this->Components[i] = ((ComplexVector&) vector).Components[i].Re;
	      }
	  }
	else
	  this->Components = 0;
	this->Flag.Initialize();
      }
    else
      {
	this->Components = 0;
	this->Flag.Initialize();
      }
}

#ifdef __MPI__

// constructor from informations sent using MPI
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the vector
// broadcast = true if the vector is broadcasted

RealVector::RealVector(MPI::Intracomm& communicator, int id, bool broadcast)
{
  this->VectorType = Vector::RealDatas;
  int TmpArray[3];
  if (broadcast == true)
    communicator.Bcast(TmpArray, 3, MPI::INT, id);      
  else
    communicator.Recv(TmpArray, 3, MPI::INT, id, 1);   
  this->Dimension = TmpArray[0];
  this->VectorId = TmpArray[1];
  this->Components = new double [this->Dimension + 1];
  if (TmpArray[2] == 1)
    for (int i = 0; i <= this->Dimension; ++i) 
      this->Components[i] = 0.0;
  else
    if (TmpArray[2] == 2)
      {
	if (broadcast == true)
	  communicator.Bcast(this->Components, this->Dimension, MPI::DOUBLE, id);      
	else
	  communicator.Recv(this->Components, this->Dimension, MPI::DOUBLE, id, 1);   
      }
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
}

#endif

// destructor
//

RealVector::~RealVector ()
{
  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

RealVector& RealVector::operator = (const RealVector& vector)
{
  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Flag = vector.Flag;
  this->VectorId = vector.VectorId;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.Dimension;
  return *this;
}

// assignement from a complex vector (keep only real part and datas are duplicated)
//
// vector = vector to assign
// return value = reference on current vector

RealVector& RealVector::operator = (const ComplexVector& vector)
{
  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->Components = new double[this->Dimension + 1];
  this->VectorId = 0;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i] = vector.Components[i].Re;
    }
  this->Flag.Initialize();
  return *this;
}

// assignement from a vector (duplicate datas if necessary)
//
// vector = vector to assign
// return value = reference on current vector

RealVector& RealVector::operator = (const Vector& vector)
{
  if (vector.VectorType == Vector::RealDatas)
    {
      return ((*this) = (RealVector&) vector);
    }
  else
    if (vector.VectorType == Vector::ComplexDatas)
      {
	return ((*this) = (ComplexVector&) vector);
      }
  return *this;
}

// Resize vector
//
// dimension = new dimension

void RealVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  double* TmpVector = new double [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    TmpVector[i] = this->Components[i];
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->Components = TmpVector;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void RealVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  double* TmpVector = new double [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    TmpVector[i] = this->Components[i];
  for (int i = this->Dimension; i < dimension; i++)
    TmpVector[i] = 0.0;  
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->Components = TmpVector;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

RealVector& RealVector::Copy (RealVector& vector, double coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  this->Localize();
  vector.Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] = vector.Components[i] * coefficient;
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* RealVector::EmptyClone(bool zeroFlag)
{
  return new RealVector(this->Dimension, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors

Vector* RealVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  RealVector* TmpVectors = new RealVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = RealVector(this->Dimension, zeroFlag);
  return TmpVectors;
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& RealVector::ClearVector ()
{
  this->Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] = 0.0;  
  this->Delocalize(true);
  return *this;
}

// change sign of a vector
//
// return value = reference on current vector

RealVector& RealVector::operator - ()
{
  this->Localize();
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i] *= -1;
    }
  this->Delocalize(true);
  return *this;
}

// return a new vector with opposite sign form a given source vector
//
// V1 = source vector
// return value = new vector

RealVector operator - (RealVector& V1)
{
  if (V1.Dimension != 0)
    {
      V1.Localize();
      double* TmpComponents = new double [V1.Dimension + 1];
      for (int i = 0; i < V1.Dimension; i++)
	TmpComponents[i] = -V1.Components[i];
      V1.Delocalize();
      return RealVector(TmpComponents, V1.Dimension);
    }
  else
    return RealVector();
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

double operator * (RealVector& V1, RealVector& V2)
{
  V1.Localize();
  V2.Localize();
  double x = V1.Components[0] * V2.Components[0];
  for (int i = 1; i < V1.Dimension; i++)
    x += V1.Components[i] * V2.Components[i];
  V1.Delocalize();
  V2.Delocalize();
  return x;
}

// do part of the scalar product between two vectors in a given range of indices
//
// vRight = right vector of the scalar product
// firstComponent = index of the first component to consider
// nbrComponent = number of components to consider
// step = increment between to consecutive indices to consider
// return value = result of the partial scalar product

double RealVector::PartialScalarProduct (RealVector& vRight, int firstComponent, int nbrComponent, int step)
{
  this->Localize();
  vRight.Localize();
  double x = this->Components[firstComponent] * vRight.Components[firstComponent];
  nbrComponent *= step;
  nbrComponent += firstComponent;
  firstComponent += step;
  for (; firstComponent < nbrComponent; firstComponent += step)
    x += this->Components[firstComponent] * vRight.Components[firstComponent];
  this->Delocalize();
  vRight.Delocalize();
  return x;
}

// sum two vectors
//
// V1 = vector to add
// return value = reference on current vector

RealVector& RealVector::operator += (RealVector& V1)
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  this->Localize();
  V1.Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] += V1.Components[i];
  this->Delocalize(true);
  V1.Delocalize();
  return *this;
}

// sum two vectors
//
// vector = vector to add
// return value = reference on current vector

Vector& RealVector::operator += (Vector& vector)
{
  if (vector.VectorType == Vector::RealDatas)
    return (*this += ((RealVector&) vector));
  return *this;
}

// substract two vectors
//
// V1 = first vector
// return value = reference on current vector

RealVector& RealVector::operator -= (RealVector& V1)
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  this->Localize();
  V1.Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] -= V1.Components[i];
  this->Delocalize(true);
  V1.Delocalize();
  return *this;
}

// sum two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

RealVector operator + (RealVector& V1, RealVector& V2)
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      V1.Localize();
      V2.Localize();
      double* TmpComponents = new double [V1.Dimension + 1];
      for (int i = 0; i < V1.Dimension; i++)
	TmpComponents[i] = V1.Components[i] + V2.Components[i];
      V1.Delocalize();
      V2.Delocalize();
      return RealVector(TmpComponents, V1.Dimension);
    }
  else
    return RealVector();
}

// substract two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

RealVector operator - (RealVector& V1, RealVector& V2)
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      V1.Localize();
      V2.Localize();
      double* TmpComponents = new double [V1.Dimension + 1];
      for (int i = 0; i < V1.Dimension; i++)
	TmpComponents[i] = V1.Components[i] - V2.Components[i];
      V1.Delocalize();
      V2.Delocalize();
      return RealVector(TmpComponents, V1.Dimension);
    }
  else
    return RealVector();
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

RealVector& RealVector::AddLinearCombination (const double& x, RealVector& V)
{
  if ((V.Dimension != this->Dimension))
    return *this;
  this->Localize();
  V.Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] += V.Components[i] * x;
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

RealVector& RealVector::AddLinearCombination (double x, RealVector& V, int firstComponent, 
					      int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
  this->Localize();
  V.Localize();
  for (int i = firstComponent; i < LastComponent; ++i)
    this->Components[i] += V.Components[i] * x;
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

RealVector& RealVector::AddLinearCombination (double x1, RealVector& v1, double x2, 
					      RealVector& v2)
{
  if ((v1.Dimension != this->Dimension) || (v2.Dimension != this->Dimension))
    return *this;
  this->Localize();
  v1.Localize();
  v2.Localize();
  for (int i = 0; i < this->Dimension; ++i)
    this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
  this->Delocalize(true);
  v1.Delocalize();
  v2.Delocalize();
  return *this;
}

// add a linear combination of two vectors to a given vector, for a given range of indices
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

RealVector& RealVector::AddLinearCombination (double x1, RealVector& v1, double x2, 
					      RealVector& v2, int firstComponent, 
					      int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
  this->Localize();
  v1.Localize();
  v2.Localize();
  for (int i = firstComponent; i < LastComponent; i++)
    this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
  this->Delocalize(true);
  v1.Delocalize();
  v2.Delocalize();
  return *this;
}

// multiply a vector with a real number on the right hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

RealVector operator * (RealVector& V1, double d)
{
  if (V1.Dimension != 0)
    {
      V1.Localize();
      double* TmpComponents = new double [V1.Dimension + 1];
      for (int i = 0; i < V1.Dimension; i++)
	TmpComponents[i] = V1.Components[i] * d ;
      V1.Delocalize();
      return RealVector(TmpComponents, V1.Dimension);
    }
  else
    return RealVector();
}

// multiply a vector with a real number on the left hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

RealVector operator * (double d, RealVector& V1)
{
  if (V1.Dimension != 0)
    {
      V1.Localize();
      double* TmpComponents = new double [V1.Dimension + 1];
      for (int i = 0; i < V1.Dimension; i++)
	TmpComponents[i] = V1.Components[i] * d ;
      V1.Delocalize();
      return RealVector(TmpComponents, V1.Dimension);
    }
  else
    return RealVector();
}

// multiply a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

RealVector& RealVector::operator *= (double d)
{
  this->Localize();
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] *= d;
  this->Delocalize(true);
  return *this;
}

// divide a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

RealVector& RealVector::operator /= (double d)
{
  this->Localize();
  double tmp = 1.0 / d;
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] *= tmp;
  this->Delocalize(true);
  return *this;
}

// left multiply a vector with a real symmetric matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

RealVector& RealVector::operator *= (const RealTriDiagonalSymmetricMatrix&  M)
{
  if ((this->Dimension != M.NbrRow) || (this->Dimension == 0))
    return *this;
  this->Localize();
  int ReducedDim = this->Dimension - 1;
  double Tmp1 = this->Components[0];
  double Tmp2;
  this->Components[0] *= M.DiagonalElements[0];
  this->Components[0] += this->Components[1] * M.UpperDiagonalElements[0];
  for (int i = 1; i < ReducedDim; i++)
    {
      Tmp2 = this->Components[i];
      this->Components[i] *= M.DiagonalElements[i];
      this->Components[i] += this->Components[i + 1] * M.UpperDiagonalElements[i] 
	+ Tmp1 * M.UpperDiagonalElements[i - 1]; 
      Tmp1 = Tmp2;
    }
  this->Components[this->Dimension - 1] *= M.DiagonalElements[this->Dimension - 1];
  this->Components[this->Dimension - 1] += Tmp1 * M.UpperDiagonalElements[this->Dimension - 2];
  this->Delocalize(true);
  return *this;
}

// left multiply a vector with a real tridiagonal symmetric matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

RealVector& RealVector::operator *= (const RealSymmetricMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  double* tmp = new double [this->Dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      tmp[i] = M.DiagonalElements[i] * this->Components[i];
      int j = 0;
      int pos = i - 1;
      for (; j < i; j++)
	{
	  tmp[i] += M.OffDiagonalElements[pos] * this->Components[j];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      pos++;
      j++;
      for (; j < this->Dimension; j++)
	{
	  tmp[i] += M.OffDiagonalElements[pos++] * this->Components[j];
	}
    }
  if (this->Flag.Shared() == true)
    {
      for (int i = 0; i < this->Dimension; i++)
	this->Components[i] = tmp[i];
    }
  else
    {
      delete[] this->Components;
      this->Components = tmp;
    }
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
  this->Delocalize(true);
  return *this;
}

// left multiply a vector with a symmetric matrix and use to store result 
// in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealSymmetricMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i] = M.DiagonalElements[i] * V.Components[i];
      int pos = i - 1;
      int j = 0;
      for (; j < i; j++)
	{
	  this->Components[i] += M.OffDiagonalElements[pos] * V.Components[j];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      pos++;
      j++;
      for (; j < this->Dimension; j++)
	this->Components[i] += M.OffDiagonalElements[pos++] * V.Components[j];
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real symmetric matrix and store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    {
      return *this;
    }
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  int j;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  double x;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; i++)
    {
      x = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->Components[i] = x;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      x += M.DiagonalElements[i] * V.Components[i];
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos2] * V.Components[j];
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->Components[i] = x;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->Components[i] = x;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    {
      return *this;
    }
  this->Localize();
  V.Localize();
  int pos3 = destStart;
  int pos;
  int j;
  int Inc =  M.NbrRow - 2 + M.Increment;
  double x;
  double* VectorPosition;
  for (int i = 0; i < M.NbrRow; i++)
    {
      pos = i - 1;
      VectorPosition = &(V.Components[sourceStart]);
      j = 0;
      x = 0;
      for (; j < i; j++)
	{
	  x += M.OffDiagonalElements[pos] * (*VectorPosition);
	  pos += Inc - j;
	  VectorPosition += sourceStep;
	}
      x += M.DiagonalElements[i] * (*VectorPosition);
      VectorPosition += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  x += M.OffDiagonalElements[pos++] * (*VectorPosition);
	  VectorPosition += sourceStep;
	}
      this->Components[pos3] = x;
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a symmetric matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceNbrComponent, int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    {
      return *this;
    }
  this->Localize();
  V.Localize();
  int pos3 = destStart;
  int pos;
  int j;
  int Inc =  M.NbrRow - 2 + M.Increment;
  double x;
  double* VectorPosition;
  for (int i = 0; i < M.NbrRow; i++)
    {
      pos = i - 1;
      VectorPosition = &(V.Components[sourceStart]);
      j = 0;
      x = 0;
      for (; j < i; j++)
	{
	  x += M.OffDiagonalElements[pos] * (*VectorPosition);
	  pos += Inc - j;
	  VectorPosition += sourceStep;
	}
      x += M.DiagonalElements[i] * (*VectorPosition);
      VectorPosition += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  x += M.OffDiagonalElements[pos++] * (*VectorPosition);
	  VectorPosition += sourceStep;
	}
      this->Components[pos3] = x;
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealSymmetricMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn) || (V.Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  int pos3 = 0;
  int Inc = M.NbrRow + M.Increment - 2;
  for (int i = 0; i < M.NbrRow; i++)
    {
      int pos = i - 1;
      int pos2 = 0;
      int j = 0;
      for (; j < i; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos] * V.Components[pos2++];
	  pos +=  Inc - j;
	}
      this->Components[pos3] += M.DiagonalElements[i] * V.Components[pos2++];
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos++] * V.Components[pos2++];
	}
      pos3++;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real symmetric matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    {
      return *this;
    }
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  int j;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  double x;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; i++)
    {
      x = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->Components[i] += x;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      x += M.DiagonalElements[i] * V.Components[i];
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos2] * V.Components[j];
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->Components[i] += x;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->Components[i] += x;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int pos;
  int j;
  int pos3 = destStart;
  int Inc = M.NbrRow + M.Increment - 2;
  double x;
  double* VectorPosition;
  for (int i = 0; i < M.NbrRow; i++)
    {
      pos = i - 1;
      j = 0;
      VectorPosition = &(V.Components[sourceStart]);
      x = 0;
      for (; j < i; j++)
	{
	  x += M.OffDiagonalElements[pos] * (*VectorPosition);
	  pos +=  Inc - j;
	  VectorPosition += sourceStep;
	}
      x += M.DiagonalElements[i] * (*VectorPosition);
      VectorPosition += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  x += M.OffDiagonalElements[pos++] * (*VectorPosition);
	  VectorPosition += sourceStep;
	}
      this->Components[pos3] += x;
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealSymmetricMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceStep, int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int pos;
  int j;
  int pos3 = destStart;
  int Inc = M.NbrRow + M.Increment - 2;
  double x;
  double* VectorPosition;
  for (int i = 0; i < M.NbrRow; i++)
    {
      pos = i - 1;
      j = 0;
      VectorPosition = &(V.Components[sourceStart]);
      x = 0;
      for (; j < i; j++)
	{
	  x += M.OffDiagonalElements[pos] * (*VectorPosition);
	  pos +=  Inc - j;
	  VectorPosition += sourceStep;
	}
      x += M.DiagonalElements[i] * (*VectorPosition);
      VectorPosition += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  x += M.OffDiagonalElements[pos++] * (*VectorPosition);
	  VectorPosition += sourceStep;
	}
      this->Components[pos3] += x;
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a antisymmetric matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealAntisymmetricMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i] = 0.0;
      int pos = i - 1;
      int j = 0;
      for (; j < i; j++)
	{
	  this->Components[i] -= M.OffDiagonalElements[pos] * V.Components[j];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      pos++;
      j++;
      for (; j < this->Dimension; j++)
	this->Components[i] += M.OffDiagonalElements[pos++] * V.Components[j];
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real antisymmetric matrix and store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    {
      return *this;
    }
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  int j;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  double x;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; i++)
    {
      x = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->Components[i] = x;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      ++j;
      for (; j < Last; ++j)
	{
	  x -= M.OffDiagonalElements[Pos2] * V.Components[j];
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->Components[i] = x;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x -= M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->Components[i] = x;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a antisymmetric matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				  int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int pos3 = destStart;
  for (int i = 0; i < M.NbrRow; i++)
    {
      int pos = i - 1;
      int pos2 = sourceStart;
      int j = 0;
      this->Components[pos3] = 0;
      for (; j < i; j++)
	{
	  this->Components[pos3] -= M.OffDiagonalElements[pos] * V.Components[pos2];
	  pos +=  M.NbrRow - j - 2 + M.Increment;
	  pos2 += sourceStep;
	}
      pos2 += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos++] * V.Components[pos2];
	  pos2 += sourceStep;
	}
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a antisymmetric matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				  int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int pos3 = destStart;
  for (int i = 0; i < M.NbrRow; i++)
    {
      int pos = i - 1;
      int pos2 = sourceStart;
      int j = 0;
      this->Components[pos3] = 0;
      for (; j < i; j++)
	{
	  this->Components[pos3] -= M.OffDiagonalElements[pos] * V.Components[pos2];
	  pos +=  M.NbrRow - j - 2 + M.Increment;
	  pos2 += sourceStep;
	}
      pos2 += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos++] * V.Components[pos2];
	  pos2 += sourceStep;
	}
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with an antisymmetric matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow) || (V.Dimension != M.NbrColumn))
    return *this;
  this->Localize();
  V.Localize();
  int pos3 = 0;
  for (int i = 0; i < M.NbrRow; i++)
    {
      int pos = i - 1;
      int pos2 = 0;
      int j = 0;
      for (; j < i; j++)
	{
	  this->Components[pos3] -= M.OffDiagonalElements[pos] * V.Components[pos2++];
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      pos2++;
      pos++;
      j++;
      for (; j < this->Dimension; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos++] * V.Components[pos2++];
	}
      pos3++;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real antisymmetric matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    {
      return *this;
    }
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  int j;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  double x;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; i++)
    {
      x = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos] * V.Components[j];
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->Components[i] += x;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x -= M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.OffDiagonalElements[Pos2] * V.Components[j];
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->Components[i] += x;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x -= M.OffDiagonalElements[Pos] * V.Components[j];
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->Components[i] += x;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with an antisymmetric matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				     int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int pos3 = destStart;
  for (int i = 0; i < M.NbrRow; i++)
    {
      int pos = i - 1;
      int pos2 = sourceStart;
      int j = 0;
      for (; j < i; j++)
	{
	  this->Components[pos3] -= M.OffDiagonalElements[pos] * V.Components[pos2];
	  pos += M.NbrRow - j - 2 + M.Increment;
	  pos2 += sourceStep;
	}
       pos2 += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos++] * V.Components[pos2];
	  pos2 += sourceStep;
	}
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with an antisymmetric matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealAntisymmetricMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				     int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (this->Dimension < ((M.NbrColumn - 1) * sourceStep + sourceStart)) || 
      (V.Dimension < ((M.NbrRow - 1) * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int pos3 = destStart;
  for (int i = 0; i < M.NbrRow; i++)
    {
      int pos = i - 1;
      int pos2 = sourceStart;
      int j = 0;
      for (; j < i; j++)
	{
	  this->Components[pos3] -= M.OffDiagonalElements[pos] * V.Components[pos2];
	  pos += M.NbrRow - j - 2 + M.Increment;
	  pos2 += sourceStep;
	}
       pos2 += sourceStep;
      pos++;
      j++;
      for (; j < M.NbrRow; j++)
	{
	  this->Components[pos3] += M.OffDiagonalElements[pos++] * V.Components[pos2];
	  pos2 += sourceStep;
	}
      pos3 += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn))
    return *this;
  this->Localize();
  V.Localize();
  if (this->Dimension != M.NbrRow)
   this->Resize(M.NbrRow);
  for (int i = 0; i < this->Dimension; i ++)
    {
      this->Components[i] = M.Columns[0].Components[i] * V.Components[0];
      for (int j = 1; j < V.Dimension; j++)
	this->Components[i] += M.Columns[j].Components[i] * V.Components[j];
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  int j = sourceStart;
  for (int i = 0; i < this->Dimension; ++i)
    {
      j = sourceStart;
      this->Components[i] = M.Columns[j].Components[i] * V.Components[j];
      ++j;
      for (; j < Last; ++j)
	{
	  this->Components[i] += M.Columns[j].Components[i] * V.Components[j];
	}
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  for (int i = 0; i < M.NbrRow; i ++)
    {
      int SourcePos = sourceStart;
      this->Components[DestPos] += M.Columns[1].Components[i] * V.Components[SourcePos];
      SourcePos += sourceStep;
      for (int j = 1; j < M.NbrColumn; j++)
	{
	  this->Components[DestPos] += M.Columns[j].Components[i] * V.Components[SourcePos];
	  SourcePos += sourceStep;
	}
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				  int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int Last = sourceNbrComponent / sourceStep + 1;
  for (int i = 0; i < Last; i ++)
    {
      int SourcePos = sourceStart;
      this->Components[DestPos] += M.Columns[1].Components[i] * V.Components[SourcePos];
      SourcePos += sourceStep;
      for (int j = 1; j < Last; j++)
	{
	  this->Components[DestPos] += M.Columns[j].Components[i] * V.Components[SourcePos];
	  SourcePos += sourceStep;
	}
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  for (int i = 0; i < M.NbrRow; i ++)
    {
      int SourcePos = sourceStart;
      for (int j = 0; j < M.NbrColumn; j++)
	{
	  this->Components[DestPos] += M.Columns[j].Components[i] * V.Components[SourcePos];
	  SourcePos += sourceStep;
	}
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  for (int i = 0; i < this->Dimension; ++i)
    {
      for (int j = sourceStart; j < Last; ++j)
	{
	  this->Components[i] += M.Columns[j].Components[i] * V.Components[j];
	}
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				     int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int Last = sourceNbrComponent / sourceStep + 1;
  if (Last > M.NbrRow)
    Last = M.NbrRow;
  for (int i = 0; i < Last; ++i)
    {
      int SourcePos = sourceStart;
      for (int j = 0; j < Last; ++j)
	{
	  this->Components[DestPos] += M.Columns[j].Components[i] * V.Components[SourcePos];
	  SourcePos += sourceStep;
	}
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealDiagonalMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  this->Localize();
  V.Localize();
  if (V.Dimension != M.NbrRow)
    V.Resize(M.NbrRow);
  for (int i = 0; i < V.Dimension; ++i)
    {
      this->Components[i] = M.DiagonalElements[i] * V.Components[i];
     }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    return *this;
  int Last = sourceStart + sourceNbrComponent;
  for (int i = sourceStart; i < Last; ++i)
    {
      this->Components[i] =  M.DiagonalElements[i] * V.Components[i];
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  for (int i = 0; i < M.NbrColumn; ++i)
    {
      this->Components[DestPos] = M.DiagonalElements[i] * V.Components[SourcePos];
      SourcePos += sourceStep;
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real diagonal matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				  int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  int Last = sourceNbrComponent / sourceStep + 1;
  if (Last > M.NbrRow)
    Last = M.NbrRow;
  for (int i = 0; i < Last; ++i)
    {
      this->Components[DestPos] = M.DiagonalElements[i] * V.Components[SourcePos];
      SourcePos += sourceStep;
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
 return *this;
}

// left multiply a vector with a real diagonal matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealDiagonalMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  for (int i = 0; i < V.Dimension; ++i)
    {
      this->Components[i] += M.DiagonalElements[i] * V.Components[i];
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// do a partial left multication of a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) 
      || (this->Dimension != M.NbrRow))
    return *this;
  this->Localize();
  V.Localize();
  int Last = sourceStart + sourceNbrComponent;
  for (int i = sourceStart; i < Last; ++i)
    {
      this->Components[i] +=  M.DiagonalElements[i] * V.Components[i];
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real diagonal matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  for (int i = 0; i < M.NbrColumn; ++i)
    {
      this->Components[DestPos] += M.DiagonalElements[i] * V.Components[SourcePos];
      SourcePos += sourceStep;
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a real diagonal matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const RealDiagonalMatrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				     int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  int Last = sourceNbrComponent / sourceStep + 1;
  if (Last > M.NbrRow)
    Last = M.NbrRow;
  for (int i = 0; i < Last; ++i)
    {
      this->Components[DestPos] += M.DiagonalElements[i] * V.Components[SourcePos];
      SourcePos += sourceStep;
      DestPos += destStep;
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a block-diagonal matrix and use to store result in current vector 
// (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::Multiply (const BlockDiagonalMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  this->Localize();
  V.Localize();
  if (V.Dimension != M.NbrRow)
    V.Resize(M.NbrRow);
  int SourcePos = 0;
  int DestPos = 0;
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlock(M.Blocks);
  while ((TmpM = IterBlock()))
    {
      this->Multiply (**TmpM, V, SourcePos, 1, DestPos, 1);
      SourcePos += (*TmpM)->GetNbrColumn();
      DestPos += (*TmpM)->GetNbrRow();
    } 
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a block-diagonal matrix and use to store result in 
// current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlock(M.Blocks);
  while ((TmpM = IterBlock()))
    {
      this->Multiply (**TmpM, V, SourcePos, sourceStep, DestPos, destStep);
      SourcePos += (*TmpM)->GetNbrColumn() * sourceStep;
      DestPos += (*TmpM)->GetNbrRow() * destStep;
    } 
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a block-diagonal matrix and use to store result in 
// current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
				  int sourceStep, int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlock(M.Blocks);
  while ((TmpM = IterBlock()))
    {
      this->Multiply (**TmpM, V, SourcePos, sourceStep, DestPos, destStep);
      SourcePos += (*TmpM)->GetNbrColumn() * sourceStep;
      DestPos += (*TmpM)->GetNbrRow() * destStep;
    } 
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a block-diagonal matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const BlockDiagonalMatrix&  M, RealVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  this->Localize();
  V.Localize();
  if (V.Dimension != M.NbrRow)
    V.Resize(M.NbrRow);
  int SourcePos = 0;
  int DestPos = 0;
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlock(M.Blocks);
  while ((TmpM = IterBlock()))
    {
      this->AddMultiply (**TmpM, V, SourcePos, 1, DestPos, 1);
      SourcePos += (*TmpM)->GetNbrColumn();
      DestPos += (*TmpM)->GetNbrRow();
    } 
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a block-diagonal matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceStep, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlock(M.Blocks);
  while ((TmpM = IterBlock()))
    {
      this->AddMultiply (**TmpM, V, SourcePos, sourceStep, DestPos, destStep);
      SourcePos += (*TmpM)->GetNbrColumn() * sourceStep;
      DestPos += (*TmpM)->GetNbrRow() * destStep;
    } 
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a block-diagonal matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const BlockDiagonalMatrix&  M, RealVector& V, int sourceStart, 
				     int sourceStep, int sourceNbrComponent, int destStart, int destStep)
{
  if ((this->Dimension == 0) || (V.Dimension != (M.NbrColumn * sourceStep + sourceStart)) 
      || (this->Dimension != (M.NbrRow * destStep + destStart)))
    return *this;
  this->Localize();
  V.Localize();
  int DestPos = destStart;
  int SourcePos = sourceStart;
  Matrix** TmpM;
  ListIterator<Matrix*> IterBlock(M.Blocks);
  while ((TmpM = IterBlock()))
    {
      this->AddMultiply (**TmpM, V, SourcePos, sourceStep, DestPos, destStep);
      SourcePos += (*TmpM)->GetNbrColumn() * sourceStep;
      DestPos += (*TmpM)->GetNbrRow() * destStep;
    } 
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// left multiply a vector with a matrix and use to store result in current vector 
// (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::Multiply (const Matrix&  M, RealVector& V)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->Multiply((BlockDiagonalMatrix&) M, V);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->Multiply((RealMatrix&) M, V);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->Multiply((RealDiagonalMatrix&) M, V);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->Multiply((RealSymmetricMatrix&) M, V);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->Multiply((RealAntisymmetricMatrix&) M, V);
      break;
    default:
      return *this;
    }
}

// left multiply a vector with a matrix and add result to current vector 
// (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const Matrix&  M, RealVector& V)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->AddMultiply((BlockDiagonalMatrix&) M, V);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->AddMultiply((RealMatrix&) M, V);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->AddMultiply((RealDiagonalMatrix&) M, V);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->AddMultiply((RealSymmetricMatrix&) M, V);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->AddMultiply((RealAntisymmetricMatrix&) M, V);
      break;
    default:
      return *this;
    }
}

// do a partial left multication of a vector with a real matrix and store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->Multiply((BlockDiagonalMatrix&) M, V, sourceStart, sourceNbrComponent);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->Multiply((RealMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->Multiply((RealDiagonalMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->Multiply((RealSymmetricMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->Multiply((RealAntisymmetricMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    default:
      return *this;
    }
}

// left multiply a vector with a matrix and use to store result in 
// current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, int destStart, int destStep)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->Multiply((BlockDiagonalMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->Multiply((RealMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->Multiply((RealDiagonalMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->Multiply((RealSymmetricMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->Multiply((RealAntisymmetricMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    default:
      return *this;
    }
}

// left multiply a vector with a matrix and use to store result in 
// current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::Multiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				  int sourceNbrComponent, int destStart, int destStep)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->Multiply((BlockDiagonalMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->Multiply((RealMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->Multiply((RealDiagonalMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      {
	return this->Multiply((RealSymmetricMatrix&) M, V, sourceStart,  sourceNbrComponent, sourceStep,destStart, destStep);
      }
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->Multiply((RealAntisymmetricMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    default:
      return *this;
    }
}

// do a partial left multication of a vector with a real matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceNbrComponent)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->AddMultiply((BlockDiagonalMatrix&) M, V, sourceStart, sourceNbrComponent);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->AddMultiply((RealMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->AddMultiply((RealDiagonalMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->AddMultiply((RealSymmetricMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->AddMultiply((RealAntisymmetricMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    default:
      return *this;
    }
}

// left multiply a vector with a matrix and add result to 
// current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				     int destStart, int destStep)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->AddMultiply((BlockDiagonalMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->AddMultiply((RealMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->AddMultiply((RealDiagonalMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->AddMultiply((RealSymmetricMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->AddMultiply((RealAntisymmetricMatrix&) M, V, sourceStart, sourceStep, destStart, destStep);
      break;
    default:
      return *this;
    }
}

// left multiply a vector with a matrix and add result to 
// current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

RealVector& RealVector::AddMultiply (const Matrix&  M, RealVector& V, int sourceStart, int sourceStep, 
				     int sourceNbrComponent, int destStart, int destStep)
{
  if ((M.MatrixType & Matrix::RealElements) == 0)
    return *this;
  if ((M.MatrixType & Matrix::BlockDiagonal) != 0)
    {
      return this->AddMultiply((BlockDiagonalMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
    }
  switch (M.MatrixType)
    {
    case (Matrix::RealElements):
      return this->AddMultiply((RealMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    case (Matrix::RealElements | Matrix::Symmetric | Matrix::Diagonal):
      return this->AddMultiply((RealDiagonalMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    case (Matrix::Symmetric | Matrix::RealElements):
      return this->AddMultiply((RealSymmetricMatrix&) M, V, sourceStart, sourceNbrComponent, sourceStep, destStart, destStep);
      break;
    case (Matrix::Antisymmetric | Matrix::RealElements):
      return this->AddMultiply((RealAntisymmetricMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    default:
      return *this;
    }
}

// get vector norm
//
// return value = vector norm

double RealVector::Norm()
{
  this->Localize();
  double x = 0.0;
  if (this->Dimension != 0)
    {
      for (int i = 0; i < this->Dimension; ++i)
	x += this->Components[i] * this->Components[i];
    }
  this->Delocalize();
  return sqrt(x);
}
  
// get square of vector norm
//
// return value = square of vector norm

double RealVector::SqrNorm ()
{
  this->Localize();
  double x = 0.0;
  if (this->Dimension != 0)
    {
      for (int i = 0; i < this->Dimension; i++)
	x += this->Components[i] * this->Components[i];
    }
  this->Delocalize();
  return x;
}
  
// normalize vector
//
// return value = reference on current vector

RealVector& RealVector::Normalize()
{
  this->Localize();
  double tmp = this->Components[0] * this->Components[0];
  for (int i = 1; i < this->Dimension; i ++)
    tmp += this->Components[i] * this->Components[i];
  tmp = 1.0 / sqrt(tmp);
  for (int i = 0; i < this->Dimension;)
    {
      this->Components[i++] *= tmp;
    }
  this->Delocalize(true);
  return *this;
}
  
// orthonormalized a vector with respect to a set of orthonormalized vectors
//
// vectors = vector array corresponding to the set
// nbrVectors = number of vectors in the set
// return value = resulting vector norm (can be used to see if vector is can be decomposed on vector set)

/*double RealVector::Orthonormalized (RealVector* vectors, int nbrVectors)
{
  double* Factors = new double []
  for (int i = 0; i )
  }*/

// Extract a subvector from a given vector
//
// FirstCoordinate = Coordinate where extraction has to begin
// LastCoordinate = Coordinate where extraction has to stop (extract also include this last coordinate)
// Step = distance to the next coordinate (1 means to take the following)
// return value = return corresponding subvector

RealVector RealVector::Extract(int FirstCoordinate, int LastCoordinate, int Step)
{
  if (this->Dimension == 0)
    return RealVector();
  this->Localize();
  RealVector TmpV ((int) ((LastCoordinate - FirstCoordinate + 1) / Step));
  for (int i = 0; i < TmpV.Dimension; i++)
    {
      TmpV.Components[i] = this->Components[FirstCoordinate];
      FirstCoordinate += Step;
    }
  this->Delocalize();
  return TmpV ; 
}
  
// Merge a subvector into a given vector
//
// V = vector to merge
// firstCoordinate = Coordinate where merge has to begin
// step = distance to the next coordinate in the destination vector (1 means to take the following)
// return value = reference to the current Vector

RealVector& RealVector::Merge(RealVector& V, int firstCoordinate, int step)
{
  if ((this->Dimension == 0) || (V.Dimension == 0))
    return *this;
  int Max = firstCoordinate + (V.Dimension * step);
  if (Max > this->Dimension)
    return *this;
  this->Localize();
  V.Localize();
  for (int i = firstCoordinate; i < Max; i += step)
    {
      this->Components[i] = V.Components[i - firstCoordinate];
      i++;
      this->Components[i] = V.Components[i - firstCoordinate];      
    }
  this->Delocalize(true);
  V.Delocalize();
  return *this;
}

// reverse elements of the current vector (i.e. exchanging i <-> N - i)
//
// return value = reference to the current Vector

RealVector& RealVector::ReverseVector()
{
  if (this->Dimension != 0)
    {
      this->Localize();
      int Max = this->Dimension - 1;
      double Tmp;
      for (int i = 0; i < this->Dimension/2; ++i)
	{
	  Tmp = this->Components[i];
	  this->Components[i] = this->Components[Max - i];
	  this->Components[Max - i] = Tmp;
	}
      this->Delocalize(true);
    }
  return *this;
}

// Output Stream overload
//
// str = reference on output stream
// v = vector to print
// return value = reference on output stream

ostream& operator << (ostream& str, RealVector& v)
{
  v.Localize();
  for (int i = 0; i < v.Dimension; ++i)
    {
      str << v.Components[i] << endl;
    }
  v.Delocalize();
  return str;
}

// output file stream overload
//
// file = reference on output file stream
// vector = reference on vector to save
// return value = reference on output file stream

/*ofstream& operator << (ofstream& file, const RealVector& vector)
{
  file.write ((char*) &(vector.Dimension), sizeof(int));
  file.write ((char*) vector.Components, sizeof(double) * vector.Dimension);
  return file;
}*/

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool RealVector::WriteVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->Dimension);
  for (int i = 0; i < this->Dimension; ++i)
    WriteLittleEndian(File, this->Components[i]);
  File.close();
  this->Delocalize();
  return true;
}

// write vector in a file in ascii mode
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool RealVector::WriteAsciiVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  int ReducedDimension = this->Dimension - 1;
  for (int i = 0; i < ReducedDimension; ++i)
    File << this->Components[i] << " ";
  File << this->Components[ReducedDimension] << endl;  
  File.close();
  this->Delocalize();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool RealVector::ReadVector (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  unsigned ZeroPos, MaxPos;
  File.seekg (0, ios::beg);
  ZeroPos = File.tellg();
  File.seekg (0, ios::end);
  MaxPos = File.tellg ();

  unsigned Length = MaxPos-ZeroPos-sizeof(int);  
  File.seekg (0, ios::beg);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (Length/sizeof(double)>(unsigned)TmpDimension)
    {      
      cout << "Error reading real vector "<<fileName<<": estimated length "<<Length/sizeof(double)<<" vs dimension "<<TmpDimension<<endl;
      if ((unsigned)TmpDimension*2==Length/sizeof(double))
	cout << "This could be a complex vector!"<<endl;
      exit(1);
    }
   
  this->Resize(TmpDimension);
  for (int i = 0; i < this->Dimension; ++i)
    ReadLittleEndian(File, this->Components[i]);
  File.close();
  return true;
}

// input file stream overload
//
// file = reference on input file stream
// vector = reference on vector to save
// return value = reference on output file stream

ifstream& operator >> (ifstream& file, RealVector& vector)
{
  file.read ((char*) &(vector.Dimension), sizeof(int));
  if (vector.Dimension > 0)
    {
      vector.Resize(vector.Dimension);
      file.read ((char*) (vector.Components), sizeof(double) * vector.Dimension);
    }
  else
    {
      vector.Dimension = 0;
      vector.Resize(vector.Dimension);      
    }
  return file;
}

#ifdef __MPI__

// send a vector to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& RealVector::SendVector(MPI::Intracomm& communicator, int id)
{
  communicator.Send(&this->VectorType, 1, MPI::INT, id, 1);
  communicator.Send(&this->Dimension, 1, MPI::INT, id, 1); 
  int Acknowledge = 0;
  communicator.Recv(&Acknowledge, 1, MPI::INT, id, 1);
  if (Acknowledge != 0)
    return *this;
  communicator.Send(this->Components, this->Dimension, MPI::DOUBLE, id, 1); 
  return *this;
}

// broadcast a vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// return value = reference on the current vector

Vector& RealVector::BroadcastVector(MPI::Intracomm& communicator,  int id)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  if (this->VectorType != TmpVectorType)
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    return *this;
  if (TmpDimension != this->Dimension)
    {
      this->Resize(TmpDimension);      
    }
  communicator.Bcast(this->Components, this->Dimension, MPI::DOUBLE, id);
  return *this;
}

// broadcast part of vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
// nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
// return value = reference on the current vector

Vector& RealVector::BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  communicator.Bcast(&firstComponent, 1, MPI::INT, id);
  communicator.Bcast(&nbrComponent, 1, MPI::INT, id);
  if (this->VectorType != TmpVectorType)
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    return *this;
  if (TmpDimension != this->Dimension)
    {
      this->Resize(TmpDimension);      
    }
  communicator.Bcast(this->Components + firstComponent, nbrComponent, MPI::DOUBLE, id);
  return *this;
}

// receive a vector from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current vector

Vector& RealVector::ReceiveVector(MPI::Intracomm& communicator, int id)
{
  int TmpVectorType = 0;
  int TmpDimension = 0;
  communicator.Recv(&TmpVectorType, 1, MPI::INT, id, 1);
  communicator.Recv(&TmpDimension, 1, MPI::INT, id, 1); 
  if (TmpVectorType != this->VectorType)
    {
      TmpDimension = 1;
      communicator.Send(&TmpDimension, 1, MPI::INT, id, 1);
      return *this;
    }
  else
    {
      if (TmpDimension != this->Dimension)
	{
	  this->Resize(TmpDimension);      
	}
      TmpDimension = 0;
      communicator.Send(&TmpDimension, 1, MPI::INT, id, 1);
    }
  communicator.Recv(this->Components, this->Dimension, MPI::DOUBLE, id, 1); 
  return *this;
}

// add current vector to the current vector of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& RealVector::SumVector(MPI::Intracomm& communicator, int id)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  if ((this->VectorType != TmpVectorType) || (this->Dimension != TmpDimension))
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    {
      return *this;
    }
  double* TmpComponents = 0;
  if (id == communicator.Get_rank())
    {
      TmpComponents = new double [this->Dimension];
    }
  communicator.Reduce(this->Components, TmpComponents, this->Dimension, MPI::DOUBLE, MPI::SUM, id);
  if (id == communicator.Get_rank())
    {
      for (int i = 0; i < this->Dimension; ++i)
	this->Components[i] = TmpComponents[i];
      delete[] TmpComponents;
    }
  return *this;
}

// create a new vector on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* RealVector::BroadcastClone(MPI::Intracomm& communicator, int id)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->VectorType, 1, MPI::INT, id);
      int TmpArray[3];
      TmpArray[0] = this->Dimension;
      TmpArray[1] = this->VectorId;
      TmpArray[2] = 2;
      communicator.Bcast(TmpArray, 3, MPI::INT, id);      
      communicator.Bcast(this->Components, this->Dimension, MPI::DOUBLE, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new RealVector(communicator, id);
    }
  return 0;
}

// create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// firstComponent = index of the first component 
// nbrComponent = number of component to send
// return value = reference on the current vector

Vector& RealVector::SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  communicator.Send(&this->VectorType, 1, MPI::INT, id, 1);
  int TmpArray[5];
  TmpArray[0] = nbrComponent;
  TmpArray[1] = this->VectorId;
  TmpArray[2] = 2;
  TmpArray[3] = firstComponent;
  TmpArray[4] = this->Dimension;
  communicator.Send(TmpArray, 5, MPI::INT, id, 1); 
  communicator.Send(this->Components + firstComponent, nbrComponent, MPI::DOUBLE, id, 1); 
  return *this;
}

// create a new vector on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* RealVector::BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->VectorType, 1, MPI::INT, id);
      int TmpArray[3];
      TmpArray[0] = this->Dimension;
      TmpArray[1] = this->VectorId;
      TmpArray[2] = 0;
      if (zeroFlag == true)
	{
	  TmpArray[2] = 1;
	}
      communicator.Bcast(TmpArray, 3, MPI::INT, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new RealVector(communicator, id);
    }
  return 0;
}

#endif
