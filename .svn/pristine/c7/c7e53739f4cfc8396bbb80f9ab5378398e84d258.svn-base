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


#include "Vector/LongIntegerVector.h"
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

LongIntegerVector::LongIntegerVector()
{
  this->VectorType = Vector::LongIntegerData;
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

LongIntegerVector::LongIntegerVector(int size, bool zeroFlag)
{
  this->VectorType = Vector::LongIntegerData;
  this->Dimension = size;
  this->TrueDimension = size;
  this->LargeDimension = (long) size;
  this->LargeTrueDimension = (long) size;
  this->VectorId = 0;
#ifdef __GMP__
  this->Components = new mpz_t [this->LargeDimension];
  for (long i = 0l; i < this->LargeDimension; ++i)
    {
      mpz_init(this->Components[i]);
    }
#else
  this->Components = new LONGLONG [this->LargeDimension];
  if (zeroFlag == true)
    {
      for (long i = 0l; i < this->LargeDimension; ++i)
	{
	  this->Components[i] = (LONGLONG) 0l;
	}
    }
#endif  
  this->Flag.Initialize();
}

// constructor for an empty integer vector bigger than 2^31
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

LongIntegerVector::LongIntegerVector(long size, bool zeroFlag)
{
  this->VectorType = Vector::LongIntegerData;
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
#ifdef __GMP__
  this->Components = new mpz_t [this->LargeDimension];
  for (long i = 0l; i < this->LargeDimension; ++i)
    {
      mpz_init(this->Components[i]);
    }
#else
  this->Components = new LONGLONG [this->LargeDimension];
  if (zeroFlag == true)
    {
      for (long i = 0l; i < this->LargeDimension; ++i)
	{
	  this->Components[i] = (LONGLONG) 0l;
	}
    }
#endif  
  this->Flag.Initialize();
}

// constructor from an array of long integer
//
// array = array of long integer to become Components of vector
// size = Vector Dimension
 
LongIntegerVector::LongIntegerVector(long* array, long size)
{
  this->LargeDimension = size;
  this->LargeTrueDimension = this->LargeDimension;
  this->VectorType = Vector::LongIntegerData;
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
#ifdef __GMP__
  this->Components = new mpz_t [this->LargeDimension];
  for (long i = 0l; i < this->LargeDimension; ++i)
    {
      mpz_init_set_si(this->Components[i], array[i]);
    }
#else
  this->Components = new LONGLONG [this->LargeDimension];
  for (long i = 0l; i < this->LargeDimension; ++i)
    {
      this->Components[i] = (LONGLONG) array[i];
    }
#endif  
  this->Flag.Initialize();
  this->VectorId = 0;
}

// copy constructor
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

LongIntegerVector::LongIntegerVector(const LongIntegerVector& vector, bool duplicateFlag)
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
#ifdef __GMP__
	  this->Components = new mpz_t [this->TrueDimension + 1]; 
#else
	  this->Components = new LONGLONG [this->TrueDimension + 1]; 
#endif
	  for (int i = 0; i < this->Dimension; i++)
	    {
#ifdef __GMP__
	      mpz_init_set(this->Components[i], vector.Components[i]);
#else
	      this->Components[i] = vector.Components[i];
#endif
	    }
	}
      else
	{
	  if (this->LargeDimension > 0l)
	    {
	      this->Flag.Initialize();
#ifdef __GMP__
	      this->Components = new mpz_t [this->LargeTrueDimension + 1]; 
#else
	      this->Components = new LONGLONG [this->LargeTrueDimension + 1]; 
#endif
	      for (long i = 0; i < this->LargeDimension; i++)
		{
#ifdef __GMP__
		  mpz_init_set(this->Components[i], vector.Components[i]);
#else
		  this->Components[i] = vector.Components[i];
#endif		  
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

LongIntegerVector::~LongIntegerVector ()
{
  if (((this->Dimension != 0) || (this->LargeDimension != 0l)) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
#ifdef __GMP__
      for (long i =0l; i < this->LargeDimension; ++i)
	{
	  mpz_clear(this->Components[i]);
	}
#endif
      delete[] this->Components;
    }
}

// Resize vector
//
// dimension = new dimension

void LongIntegerVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      this->LargeDimension = dimension;
      return;
    }
#ifdef __GMP__
  mpz_t* TmpComponents = new mpz_t [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      mpz_init_set(TmpComponents[i], this->Components[i]);
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->Dimension; i++)
	{
	  mpz_clear(this->Components[i]);
	}
      delete[] this->Components;
    }
#else
  LONGLONG* TmpComponents = new LONGLONG [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      TmpComponents[i] = this->Components[i];
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
#endif
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

void LongIntegerVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      this->LargeDimension = dimension;
      return;
    }
#ifdef __GMP__
  mpz_t* TmpComponents = new mpz_t [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      mpz_init_set(TmpComponents[i], this->Components[i]);
    }
  for (int i = this->Dimension; i < dimension; i++)
    {
      mpz_init(TmpComponents[i]);
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->Dimension; i++)
	{
	  mpz_clear(this->Components[i]);
	}
      delete[] this->Components;
    }
#else
  LONGLONG* TmpComponents = new LONGLONG [dimension + 1];
  for (int i = 0; i < this->Dimension; i++)
    {
      TmpComponents[i] = this->Components[i];
    }
  for (int i = this->Dimension; i < dimension; i++)
    {
      TmpComponents[i] = (LONGLONG) 0l;
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
#endif
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->LargeDimension = dimension;
  this->LargeTrueDimension = dimension;
  this->Components = TmpComponents;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& LongIntegerVector::ClearVector ()
{
  for (long i = 0l; i < this->LargeTrueDimension; ++i)
    {
#ifdef __GMP__
      mpz_set_si(this->Components[i], 0l);
#else
      this->Components[i] = 0l;
#endif
    }
  return *this;
}

// put select vector components to zero
// start = start index
// nbrComponent = number of components to set to zero
// return value = reference on current vector

Vector& LongIntegerVector::ClearVectorSegment (long start, long nbrComponent)
{
  nbrComponent += start;
  for (; start < nbrComponent; ++start)
    {
#ifdef __GMP__
      mpz_set_si(this->Components[start], 0l);
#else
     this->Components[start] = 0l;
#endif
    }
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::Copy (LongIntegerVector& vector)
{
  if ((this->Dimension != vector.Dimension)||(this->LargeDimension != vector.LargeDimension))
    this->Resize(vector.LargeDimension);
  this->Localize();
  vector.Localize();
  for (long i = 0; i < this->LargeDimension; i++)
    {
#ifdef __GMP__
      mpz_init_set(this->Components[i], vector.Components[i]);
#else      
     this->Components[i] = vector.Components[i];
#endif
    }
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::Copy (LongIntegerVector& vector, const long& coefficient)
{
  if ((this->Dimension != vector.Dimension)||(this->LargeDimension != vector.LargeDimension))
    this->Resize(vector.LargeDimension);
  this->Localize();
  vector.Localize();
#ifdef __GMP__
  for (long i = 0; i < this->LargeDimension; i++)
    {      
      mpz_mul_si(this->Components[i], vector.Components[i], coefficient);
    }
#else      
  for (long i = 0; i < this->LargeDimension; i++)
    {      
      this->Components[i] = vector.Components[i] * coefficient;
    }
#endif
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* LongIntegerVector::EmptyClone(bool zeroFlag)
{
  return new LongIntegerVector(this->LargeDimension, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors
  
Vector* LongIntegerVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  LongIntegerVector* TmpVectors = new LongIntegerVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = LongIntegerVector(this->LargeDimension, zeroFlag);
  return TmpVectors;
}

// test if the vector is a null vector
//
// return value = true if the vector is a null vector

bool LongIntegerVector::IsNullVector()
{
  long i = 0l;
#ifdef __GMP__
  while ((i < this->LargeDimension) && (mpz_sgn(this->Components[i]) == 0))
#else      
    while ((i < this->LargeDimension) && (this->Components[i] == ((LONGLONG) 0l)))
#endif
    {
      ++i;
    }
  if (i == this->LargeDimension)
    return true;
  else
    return false;  
}

// test if the current vector is proportional to another one
//
// vector = reference on the vector to compare to
// return value = true if the two vectors are proportional 

bool LongIntegerVector::IsProportional(LongIntegerVector& vector)
{
  long FirstNonZeroCoefficient = 0l;
#ifdef __GMP__
  while ((FirstNonZeroCoefficient < this->LargeDimension) && (mpz_sgn(this->Components[FirstNonZeroCoefficient]) == 0))
#else
    while ((FirstNonZeroCoefficient < this->LargeDimension) && (this->Components[FirstNonZeroCoefficient] == ((LONGLONG) 0l)))
#endif    
    {
      ++FirstNonZeroCoefficient;
    }
  if (FirstNonZeroCoefficient == this->LargeDimension)
    {
      return false;
    }
  for (long i = 0l; i < FirstNonZeroCoefficient; ++i)
    {
#ifdef __GMP__
      if (mpz_sgn(vector.Components[i]) != 0)
#else
	if (vector.Components[i] != ((LONGLONG) 0l))
#endif    
	{
	  return false;
	}
    }
#ifdef __GMP__
  if (mpz_sgn(vector.Components[FirstNonZeroCoefficient]) == 0)
#else
    if (vector.Components[FirstNonZeroCoefficient] == ((LONGLONG) 0l))
#endif    
      {
	return false;
      }
  
#ifdef __GMP__
  mpz_t* SmallerVector = this->Components;
  mpz_t* LargerVector = vector.Components;
  if (mpz_cmpabs(this->Components[FirstNonZeroCoefficient], vector.Components[FirstNonZeroCoefficient]) > 0)
    {
      LargerVector = this->Components;
      SmallerVector = vector.Components;	  
    }
  if (mpz_divisible_p(LargerVector[FirstNonZeroCoefficient], SmallerVector[FirstNonZeroCoefficient]) == 0)
    {
      return false;
    }
  mpz_t TmpFactor;
  mpz_t TmpCoefficient;
  mpz_init(TmpFactor);
  mpz_init(TmpCoefficient);
  mpz_divexact(TmpFactor, LargerVector[FirstNonZeroCoefficient], SmallerVector[FirstNonZeroCoefficient]);
  for (long i = FirstNonZeroCoefficient + 1l; i < this->LargeDimension; ++i)
    {	  
      mpz_mul(TmpCoefficient, SmallerVector[i], TmpFactor);
      if (mpz_cmp(LargerVector[i], TmpCoefficient) != 0)
	{
	  mpz_clear(TmpCoefficient);  
	  mpz_clear(TmpFactor);  
	  return false;
	}
    }
  mpz_clear(TmpCoefficient);  
  mpz_clear(TmpFactor);  
  return true;
#else
  LONGLONG* SmallerVector = this->Components;
  LONGLONG* LargerVector = vector.Components;
  if (abs(this->Components[FirstNonZeroCoefficient]) > abs(vector.Components[FirstNonZeroCoefficient]))
    {
      LargerVector = this->Components;
      SmallerVector = vector.Components;	  
    }
  if ((LargerVector[FirstNonZeroCoefficient] % SmallerVector[FirstNonZeroCoefficient]) != 0)
    {
      return false;
    }
  LONGLONG TmpFactor = LargerVector[FirstNonZeroCoefficient] / SmallerVector[FirstNonZeroCoefficient];
  for (long i = FirstNonZeroCoefficient + 1l; i < this->LargeDimension; ++i)
    {	  
      if (LargerVector[i] != (TmpFactor * SmallerVector[i]))
	{
	  return false;
	}
    }
  return true;
#endif
}

// sum two vectors
//
// vector = vector to add
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::operator += (LongIntegerVector& vector)
{
  for (long i = 0l; i < this->LargeDimension; ++i)
    {
#ifdef __GMP__
      mpz_add(this->Components[i], this->Components[i], vector[i]);
#else      
      this->Components[i] += vector[i];
#endif
    }
  return *this;
}

// substract two vectors
//
// vector = vector to substract
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::operator -= (LongIntegerVector& vector)
{
  for (long i = 0; i < this->LargeDimension; ++i)
    {
#ifdef __GMP__
      mpz_sub(this->Components[i], this->Components[i], vector[i]);
#else      
      this->Components[i] -= vector[i];
#endif      
    }
  return *this;
}

// multiply a vector with an integer number on the right hand side
//
// d = integer to use
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::operator *= (const long& d)
{
  for (long i = 0; i < this->LargeDimension; ++i)
    {
#ifdef __GMP__
      mpz_mul_si(this->Components[i], this->Components[i], d);
#else
      this->Components[i] *= d;
#endif           
    }
  return *this;
}

// divide a vector by an integer number on the right hand side
//
// d = integer to use
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::operator /= (const long& d)
{
#ifdef __GMP__
  mpz_t Tmp;
  mpz_init_set_si(Tmp, d);
#endif           
  for (long i = 0; i < this->LargeDimension; ++i)
    {
#ifdef __GMP__
      mpz_divexact (this->Components[i], this->Components[i], Tmp);
#else
      this->Components[i] /= d;
#endif           
    }
#ifdef __GMP__
  mpz_clear(Tmp);
#endif           
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::AddLinearCombination (const long& x, LongIntegerVector& V)
{
  if ((V.Dimension != this->Dimension)||(V.LargeDimension != this->LargeDimension))
    return *this;
#ifdef __GMP__
  if (x >= 0l)
    {
      for (long i = 0; i < this->LargeDimension; i++)
	{
	  mpz_addmul_ui(this->Components[i], V.Components[i], (unsigned long) x);
	}
    }
  else
    {
      unsigned long Tmp = (unsigned long) -x;
      for (long i = 0; i < this->LargeDimension; i++)
	{
	  mpz_submul_ui(this->Components[i], V.Components[i], Tmp);
	}      
    }
#else
  for (long i = 0; i < this->LargeDimension; i++)
    {
      this->Components[i] += V.Components[i] * (LONGLONG) x;
    }      
#endif           
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::AddLinearCombination (const long& x, LongIntegerVector& V, int firstComponent, 
							    int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
#ifdef __GMP__
  if (x >= 0l)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  mpz_addmul_ui(this->Components[i], V.Components[i], (unsigned long) x);
	}
    }
  else
    {
      unsigned long Tmp = (unsigned long) -x;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  mpz_submul_ui(this->Components[i], V.Components[i], Tmp);
	}      
    }
#else
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->Components[i] += V.Components[i] * (LONGLONG) x;
    }      
#endif           
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

LongIntegerVector& LongIntegerVector::AddLinearCombination (const long& x1, LongIntegerVector& v1, const long& x2, 
							    LongIntegerVector& v2)
{
  if ((v1.Dimension != this->Dimension) || (v2.Dimension != this->Dimension) ||
      (v1.LargeDimension != this->LargeDimension) || (v2.LargeDimension != this->LargeDimension))
    return *this;
#ifdef __GMP__
#else
  for (long i = 0; i < this->LargeDimension; ++i)
    {
      this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
    }
#endif
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

LongIntegerVector& LongIntegerVector::AddLinearCombination (const long& x1, LongIntegerVector& v1, const long& x2, 
							    LongIntegerVector& v2, int firstComponent, 
							    int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
#ifdef __GMP__
#else
  for (int i = firstComponent; i < LastComponent; i++)
    {
      this->Components[i] += v1.Components[i] * x1 + v2.Components[i] * x2;
    }
#endif
  return *this;
}

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

// bool LongIntegerVector::WriteVector (const char* fileName)
// {
//   this->Localize();
//   ofstream File;
//   File.open(fileName, ios::binary | ios::out);
//   WriteLittleEndian(File, this->Dimension);
//   if (this->Dimension == -1)
//     {
//       WriteLittleEndian(File, this->LargeDimension);
//       for (long i = 0; i < this->LargeDimension; ++i)
// 	{
// 	  WriteLittleEndian(File, this->Components[i]);
// 	}
//     }
//   else
//     for (int i = 0; i < this->Dimension; ++i)
//       {
// 	WriteLittleEndian(File, this->Components[i]);
//       }
//   File.close();
//   this->Delocalize();
//   return true;
// }

// write vector in a file in ascii mode
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool LongIntegerVector::WriteAsciiVector (const char* fileName)
{
  this->Localize();
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
#ifdef __GMP__
  long ReducedDimension = this->LargeDimension - 1l;
  for (long i = 0; i < ReducedDimension; ++i)
    File << this->Components[i] << "  ";
  File << this->Components[ReducedDimension] << endl;  
#else
  long ReducedDimension = this->LargeDimension - 1l;
  for (long i = 0; i < ReducedDimension; ++i)
    File << (long) this->Components[i] << "  ";
  File << (long) this->Components[ReducedDimension] << endl;  
#endif
  File.close();
  this->Delocalize();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

// bool LongIntegerVector::ReadVector (const char* fileName)
// {
//   ifstream File;
//   File.open(fileName, ios::binary | ios::in);
//   if (!File.is_open())
//     {
//       cout << "Cannot open the file: " << fileName << endl;
//       return false;
//     }
  
//   int TmpDimension;
//   ReadLittleEndian(File, TmpDimension);

//   if (TmpDimension > 0)
//     {
//       this->Resize(TmpDimension);
//       for (int i = 0; i < this->Dimension; ++i)
// 	{
// 	  ReadLittleEndian(File, this->Components[i]);
// 	}
//     }
//   else
//     {
//       long TmpLargeDimension;
//       ReadLittleEndian(File, TmpLargeDimension);
//       this->Resize(TmpLargeDimension);
//       for (long i = 0; i < this->LargeDimension; ++i)
// 	{
// 	  ReadLittleEndian(File, this->Components[i]);
// 	}
//     }
//   File.close();
//   return true;
// }

// read vector from a file, only within a given range of indices
//
// fileName = name of the file where the vector has to be read
// minIndex = index of the first component to read (if negative, start from the end of vector)
// maxIndex = index of the last component to read (negative or zero is considered as the last component)
// return value = true if no error occurs

// bool LongIntegerVector::ReadVector (const char* fileName, long minIndex, long maxIndex)
// {
//   ifstream File;
//   File.open(fileName, ios::binary | ios::in);
//   if (!File.is_open())
//     {
//       cout << "Cannot open the file: " << fileName << endl;
//       return false;
//     }
  
//   int TmpDimension;
//   ReadLittleEndian(File, TmpDimension);

//   if (TmpDimension > 0)
//     {
//       if ((maxIndex >= TmpDimension) || (maxIndex <= 0l))
//  	maxIndex = TmpDimension - 1;
//       if (minIndex < 0l)      
//  	minIndex += TmpDimension;
//        this->Resize(maxIndex - minIndex + 1l);
//       for (int i = 0; i < minIndex; ++i)
// 	{
// 	  ReadLittleEndian(File, this->Components[0]);
// 	}
//       for (int i = minIndex; i <= maxIndex; ++i)
// 	{
// 	  ReadLittleEndian(File, this->Components[i - minIndex]);
// 	}
//     }
//   else
//     {
//       long TmpLargeDimension;
//       ReadLittleEndian(File, TmpLargeDimension);
//       if ((maxIndex >= TmpLargeDimension) || (maxIndex <= 0l))
//  	maxIndex = TmpLargeDimension - 1l;
//        if (minIndex < 0l)      
//  	minIndex += TmpLargeDimension;
//        this->Resize(maxIndex - minIndex + 1l);
//       for (long i = 0l; i < minIndex; ++i)
// 	{
// 	  ReadLittleEndian(File, this->Components[0]);
// 	}
//       for (long i = minIndex; i <= maxIndex; ++i)
// 	{
// 	  ReadLittleEndian(File, this->Components[i - minIndex]);
// 	}
//     }
//   File.close();
//   return true;
// }

// read vector dimension from a file, without loading the full vector 
//
// fileName = name of the file where the vector has to be read
// return value = vector dimension

// long LongIntegerVector::ReadVectorDimension (const char* fileName)
// {
//   ifstream File;
//   File.open(fileName, ios::binary | ios::in);
//   if (!File.is_open())
//     {
//       cout << "Cannot open the file: " << fileName << endl;
//       return false;
//     }
//   int TmpDimension;
//   ReadLittleEndian(File, TmpDimension);
//   if (TmpDimension > 0)
//     {
//       File.close();
//       return ((long) TmpDimension);
//     }
//   long TmpLargeDimension = 0l;
//   ReadLittleEndian(File, TmpLargeDimension);
//   File.close();
//   return TmpLargeDimension;
// }

// test if a vector can be read from a file (matching the right type), without loading the full vector 
//
// fileName = name of the file where the vector has to be read
// return value = true if the vector can be read

// bool LongIntegerVector::ReadVectorTest (const char* fileName)
// {
//   ifstream File;
//   File.open(fileName, ios::binary | ios::in);
//   if (!File.is_open())
//     {
//       cout << "Cannot open the file: " << fileName << endl;
//       return false;
//     }
  
//   std::streampos ZeroPos, MaxPos;
//   File.seekg (0, ios::beg);
//   ZeroPos = File.tellg();
//   File.seekg (0, ios::end);
//   MaxPos = File.tellg ();

//   long Length = (long) (MaxPos - ZeroPos) - sizeof(int);  
//   File.seekg (0, ios::beg);
//   int TmpDimension;
//   ReadLittleEndian(File, TmpDimension);

//   if (TmpDimension > 0)
//     {
//       File.close();
//       Length /= sizeof(long);
//       if (Length == ((long) TmpDimension))
// 	return true;
//       else
// 	return false;
//     }
//   else
//     {
//       long TmpLargeDimension;
//       ReadLittleEndian(File, TmpLargeDimension);
//       File.close();
//       Length -= sizeof (long);
//       Length /= sizeof(long);
//       if (Length == TmpLargeDimension)
// 	return true;
//       else
// 	return false;
//     }
//   File.close();
//   return false;
// }

// Output Stream overload
//
// str = reference on output stream
// v = vector to print
// return value = referenceint GetNumSites(){return this->NSites;} on output stream

ostream& operator << (ostream& str, LongIntegerVector& v)
{
  for (long i = 0; i < v.LargeDimension; ++i)
    {
#ifdef __GMP__
      str << v.Components[i] << endl;
#else
      str << (long) v.Components[i] << endl;
#endif
    }
  return str;
}
