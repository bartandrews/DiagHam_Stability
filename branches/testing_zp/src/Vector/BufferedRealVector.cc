////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//    where elements are stored on demand and tramsmitted to a master vector  //
//                           when the buffer is full                          //
//                                                                            //
//                        last modification : 19/07/2008                      //
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


#include "Vector/PartialRealVector.h"
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


using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//

PartialRealVector::PartialRealVector()
{
  this->VectorType = Vector::RealDatas | Vector::BufferedData;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->Components = 0;
  this->VectorId = 0;
  this->IndexShift = 0;
  this->RealDimension = 0;
}

// constructor for an empty real vector (all coordinates set to zero)
//
// size = effective vector dimension (i.e. number of stored elements)
// realSize = size of the full vector
// indexShift = index of the first element which is actually stored in the current partial vector
// zeroFlag = true if all coordinates have to be set to zero

PartialRealVector::PartialRealVector(int size, int realSize, int indexShift, bool zeroFlag)
{
  this->VectorType = Vector::RealDatas | Vector::PartialData;
  this->Dimension = size;
  this->IndexShift = indexShift;
  this->RealDimension = realSize;
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

// copy constructor
//
// vector = vector to copy
// DuplicateFlag = true if datas have to be duplicated

PartialRealVector::PartialRealVector(const PartialRealVector& vector, bool duplicateFlag)
{
  this->VectorType = Vector::RealDatas | Vector::PartialData;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
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

// destructor
//

PartialRealVector::~PartialRealVector ()
{
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

PartialRealVector& PartialRealVector::operator = (const PartialRealVector& vector)
{
  //  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
  if (this->Dimension != 0)
    {
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->Components;
	}
    }
  this->Flag = vector.Flag;
  this->VectorId = vector.VectorId;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.Dimension;
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

PartialRealVector& PartialRealVector::Copy (PartialRealVector& vector, double coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  this->Localize();
  vector.Localize();
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
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

Vector* PartialRealVector::EmptyClone(bool zeroFlag)
{
  return new PartialRealVector(this->Dimension, this->Dimension, 0, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors

Vector* PartialRealVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  PartialRealVector* TmpVectors = new PartialRealVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = PartialRealVector(this->Dimension, this->Dimension, 0, zeroFlag);
  return TmpVectors;
}


