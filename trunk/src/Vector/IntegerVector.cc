////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class for n dimensional integer vector using long           //
//                                                                            //
//                        last modification : 03/01/2022                      //
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


#include "Vector/IntegerVector.h"
#include "GeneralTools/Endian.h"

#include <fstream>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor  
// 

IntegerVector::IntegerVector()
{
  this->VectorType = Vector::IntegerData;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->LargeDimension = 0l;
  this->LargeTrueDimension = 0l;
  this->VectorId = 0;
  this->Components = 0;
}

// default constructor for an empty vector
// 
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

IntegerVector::IntegerVector(int size, bool zeroFlag)
{
  this->VectorType = Vector::IntegerData;
  this->Dimension = size;
  this->TrueDimension = size;
  this->LargeDimension = (long) size;
  this->LargeTrueDimension = (long) size;
  this->VectorId = 0;
  this->Components = new long [this->LargeDimension];
  this->Flag.Initialize();
}

// constructor for an empty integer vector bigger than 2^31
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

IntegerVector::IntegerVector(long size, bool zeroFlag)
{
  this->VectorType = Vector::IntegerData;
  this->LargeDimension = size;
  this->LargeTrueDimension = this->LargeDimension;
#ifdef  __64_BITS__
  if (this->LargeDimension < (1l << 31))
    this->Dimension = (int) size;
  else
    {
      this->Dimension = -1;
      this->VectorType |= Vector::LargeData;
    }
#else
  this->Dimension = (int) size;
#endif
  this->TrueDimension = this->Dimension;
  this->Components = new long [this->LargeDimension];
  this->Flag.Initialize();
}

// constructor from an array of long integer
//
// array = array of long integer to become Components of vector
// size = Vector Dimension
 
IntegerVector::IntegerVector(long* array, long size)
{
  this->LargeDimension = size;
  this->LargeTrueDimension = this->LargeDimension;
  this->VectorType = Vector::IntegerData;
  #ifdef  __64_BITS__
  if (this->LargeDimension < (1l << 31))
    this->Dimension = (int) size;
  else
    {
      this->Dimension = -1;
      this->VectorType |= Vector::LargeData;
    }
#else
  this->Dimension = (int) size;
#endif
  this->TrueDimension = this->Dimension;
  this->Components = array;
  this->Flag.Initialize();
  this->VectorId = 0;
}

// copy constructor
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

IntegerVector::IntegerVector(const IntegerVector& vector, bool duplicateFlag)
{
  this->VectorType = vector.VectorType;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->LargeDimension = vector.LargeDimension;
  this->LargeTrueDimension = vector.LargeTrueDimension;
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
	  this->Components = new long [this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    {
	      this->Components[i] = vector.Components[i];
	    }
	}
      else
	{
	  if (this->LargeDimension > 0l)
	    {
	      this->Flag.Initialize();
	      this->Components = new  long[this->LargeTrueDimension + 1]; 
	      for (long i = 0; i < this->LargeDimension; i++)
		{
		  this->Components[i] = vector.Components[i];
		}
	    }
	  else
	    {
	      this->Components = 0;
	    }
	}
    }
}

// destructor
//

IntegerVector::~IntegerVector ()
{
  if (((this->Dimension != 0) || (this->LargeDimension != 0l)) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
}

// Resize vector
//
// dimension = new dimension

void IntegerVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      this->LargeDimension = dimension;
      return;
    }
  long* TmpComponents = new long [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      TmpComponents[i] = this->Components[i];
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->LargeDimension = dimension;
  this->LargeTrueDimension = dimension;
  this->Components = TmpComponents;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void IntegerVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      this->LargeDimension = dimension;
      return;
    }
  long* TmpComponents = new long [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      TmpComponents[i] = this->Components[i];
    }
  for (int i = this->Dimension; i < dimension; i++)
    {
      TmpComponents[i] = 0l;
    }
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->LargeDimension = dimension;
  this->LargeTrueDimension = dimension;
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Components = TmpComponents;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& IntegerVector::ClearVector ()
{
  for (long i = 0l; i < this->LargeTrueDimension; ++i)
    {
      this->Components[i] = 0l;
    }
  return *this;
}

// put select vector components to zero
// start = start index
// nbrComponent = number of components to set to zero
// return value = reference on current vector

Vector& IntegerVector::ClearVectorSegment (long start, long nbrComponent)
{
  nbrComponent += start;
  for (; start < nbrComponent; ++start)
    {
      this->Components[start] = 0l;
    }
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

IntegerVector& IntegerVector::Copy (IntegerVector& vector)
{
  if ((this->Dimension != vector.Dimension)||(this->LargeDimension != vector.LargeDimension))
    this->Resize(vector.LargeDimension);
  this->Localize();
  vector.Localize();
  for (long i = 0; i < this->LargeDimension; i++)
    this->Components[i] = vector.Components[i];
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

IntegerVector& IntegerVector::Copy (IntegerVector& vector, const long& coefficient)
{
  if ((this->Dimension != vector.Dimension)||(this->LargeDimension != vector.LargeDimension))
    this->Resize(vector.LargeDimension);
  this->Localize();
  vector.Localize();
  for (long i = 0; i < this->LargeDimension; i++)
    this->Components[i] = vector.Components[i] * coefficient;
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* IntegerVector::EmptyClone(bool zeroFlag)
{
  return new IntegerVector(this->LargeDimension, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors
  
Vector* IntegerVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  IntegerVector* TmpVectors = new IntegerVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = IntegerVector(this->LargeDimension, zeroFlag);
  return TmpVectors;
}

// test if the vector is a null vector
//
// return value = true if the vector is a null vector

bool IntegerVector::IsNullVector()
{
  long i = 0l;
  while ((i < this->LargeDimension) && (this->Components[i] == 0l))
    ++i;
  if (i == this->LargeDimension)
    return true;
  else
    return false;  
}

// sum two vectors
//
// vector = vector to add
// return value = reference on current vector

IntegerVector& IntegerVector::operator += (IntegerVector& vector)
{
  if (this->Dimension == -1)
    for (long i = 0l; i < this->LargeDimension; ++i)
      this->Components[i] += vector[i];
  else
    {
      for (int i = 0; i < this->Dimension; i++)
	this->Components[i] += vector[i];
    }
  return *this;
}

// substract two vectors
//
// vector = vector to substract
// return value = reference on current vector

IntegerVector& IntegerVector::operator -= (IntegerVector& vector)
{
  if (this->Dimension == -1)
    for (long i = 0; i < this->LargeDimension; ++i)
      this->Components[i] -= vector[i];
  else
    {
      for (int i = 0; i < this->Dimension; i++)
	this->Components[i] -= vector[i];
    }
  return *this;
}

// multiply a vector with an integer number on the right hand side
//
// d = integer to use
// return value = reference on current vector

IntegerVector& IntegerVector::operator *= (const long& d)
{
  if (this->Dimension == -1)
    for (long i = 0; i < this->LargeDimension; ++i)
      this->Components[i] *= d;
  else
    {
      for (int i = 0; i < this->Dimension; i++)
	this->Components[i] *= d;
    }
  return *this;
}

// divide a vector by an integer number on the right hand side
//
// d = integer to use
// return value = reference on current vector

IntegerVector& IntegerVector::operator /= (const long& d)
{
  if (this->Dimension == -1)
    for (long i = 0; i < this->LargeDimension; ++i)
      this->Components[i] /= d;
  else
    {
      for (int i = 0; i < this->Dimension; i++)
	this->Components[i] /= d;
    }
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

IntegerVector& IntegerVector::AddLinearCombination (const long& x, IntegerVector& V)
{
  if ((V.Dimension != this->Dimension)||(V.LargeDimension != this->LargeDimension))
    return *this;
  if (V.Dimension == -1)
    for (long i = 0; i < this->LargeDimension; i++)
      this->Components[i] += V.Components[i] * x;
  else
    for (int i = 0; i < this->Dimension; i++)
      this->Components[i] += V.Components[i] * x;
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

IntegerVector& IntegerVector::AddLinearCombination (const long& x, IntegerVector& V, int firstComponent, 
							      int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    this->Components[i] += V.Components[i] * x;
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

IntegerVector& IntegerVector::AddLinearCombination (const long& x1, IntegerVector& v1, const long& x2, 
							      IntegerVector& v2)
{
  if ((v1.Dimension != this->Dimension) || (v2.Dimension != this->Dimension) ||
      (v1.LargeDimension != this->LargeDimension) || (v2.LargeDimension != this->LargeDimension))
    return *this;
  if (this->Dimension==-1)
    for (long i = 0; i < this->LargeDimension; ++i)
      this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
  else
    for (int i = 0; i < this->Dimension; ++i)
      this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
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

IntegerVector& IntegerVector::AddLinearCombination (const long& x1, IntegerVector& v1, const long& x2, 
						    IntegerVector& v2, int firstComponent, 
						    int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; i++)
    this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
  return *this;
}

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool IntegerVector::WriteVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->Dimension);
  if (this->Dimension == -1)
    {
      WriteLittleEndian(File, this->LargeDimension);
      for (long i = 0; i < this->LargeDimension; ++i)
	{
	  WriteLittleEndian(File, this->Components[i]);
	}
    }
  else
    for (int i = 0; i < this->Dimension; ++i)
      {
	WriteLittleEndian(File, this->Components[i]);
      }
  File.close();
  this->Delocalize();
  return true;
}

// write vector in a file in ascii mode
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool IntegerVector::WriteAsciiVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  if (this->Dimension == -1)
    {
      long ReducedDimension = this->LargeDimension - 1;
      for (long i = 0; i < ReducedDimension; ++i)
	File << this->Components[i] << "  ";
      File << this->Components[ReducedDimension] << endl;  
    }
  else
    {
      int ReducedDimension = this->Dimension - 1;
      for (int i = 0; i < ReducedDimension; ++i)
	File << this->Components[i] << "  ";
      File << this->Components[ReducedDimension] << endl;  
    }
  File.close();
  this->Delocalize();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool IntegerVector::ReadVector (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (TmpDimension > 0)
    {
      this->Resize(TmpDimension);
      for (int i = 0; i < this->Dimension; ++i)
	{
	  ReadLittleEndian(File, this->Components[i]);
	}
    }
  else
    {
      long TmpLargeDimension;
      ReadLittleEndian(File, TmpLargeDimension);
      this->Resize(TmpLargeDimension);
      for (long i = 0; i < this->LargeDimension; ++i)
	{
	  ReadLittleEndian(File, this->Components[i]);
	}
    }
  File.close();
  return true;
}

// read vector from a file, only within a given range of indices
//
// fileName = name of the file where the vector has to be read
// minIndex = index of the first component to read (if negative, start from the end of vector)
// maxIndex = index of the last component to read (negative or zero is considered as the last component)
// return value = true if no error occurs

bool IntegerVector::ReadVector (const char* fileName, long minIndex, long maxIndex)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (TmpDimension > 0)
    {
      if ((maxIndex >= TmpDimension) || (maxIndex <= 0l))
 	maxIndex = TmpDimension - 1;
      if (minIndex < 0l)      
 	minIndex += TmpDimension;
       this->Resize(maxIndex - minIndex + 1l);
      for (int i = 0; i < minIndex; ++i)
	{
	  ReadLittleEndian(File, this->Components[0]);
	}
      for (int i = minIndex; i <= maxIndex; ++i)
	{
	  ReadLittleEndian(File, this->Components[i - minIndex]);
	}
    }
  else
    {
      long TmpLargeDimension;
      ReadLittleEndian(File, TmpLargeDimension);
      if ((maxIndex >= TmpLargeDimension) || (maxIndex <= 0l))
 	maxIndex = TmpLargeDimension - 1l;
       if (minIndex < 0l)      
 	minIndex += TmpLargeDimension;
       this->Resize(maxIndex - minIndex + 1l);
      for (long i = 0l; i < minIndex; ++i)
	{
	  ReadLittleEndian(File, this->Components[0]);
	}
      for (long i = minIndex; i <= maxIndex; ++i)
	{
	  ReadLittleEndian(File, this->Components[i - minIndex]);
	}
    }
  File.close();
  return true;
}

// read vector dimension from a file, without loading the full vector 
//
// fileName = name of the file where the vector has to be read
// return value = vector dimension

long IntegerVector::ReadVectorDimension (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  if (TmpDimension > 0)
    {
      File.close();
      return ((long) TmpDimension);
    }
  long TmpLargeDimension = 0l;
  ReadLittleEndian(File, TmpLargeDimension);
  File.close();
  return TmpLargeDimension;
}

// test if a vector can be read from a file (matching the right type), without loading the full vector 
//
// fileName = name of the file where the vector has to be read
// return value = true if the vector can be read

bool IntegerVector::ReadVectorTest (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  
  std::streampos ZeroPos, MaxPos;
  File.seekg (0, ios::beg);
  ZeroPos = File.tellg();
  File.seekg (0, ios::end);
  MaxPos = File.tellg ();

  long Length = (long) (MaxPos - ZeroPos) - sizeof(int);  
  File.seekg (0, ios::beg);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (TmpDimension > 0)
    {
      File.close();
      Length /= sizeof(long);
      if (Length == ((long) TmpDimension))
	return true;
      else
	return false;
    }
  else
    {
      long TmpLargeDimension;
      ReadLittleEndian(File, TmpLargeDimension);
      File.close();
      Length -= sizeof (long);
      Length /= sizeof(long);
      if (Length == TmpLargeDimension)
	return true;
      else
	return false;
    }
  File.close();
  return false;
}

// Output Stream overload
//
// str = reference on output stream
// v = vector to print
// return value = referenceint GetNumSites(){return this->NSites;} on output stream

ostream& operator << (ostream& str, IntegerVector& v)
{
  for (long i = 0; i < v.LargeDimension; ++i)
    str << v.Components[i] << endl;
  return str;
}
