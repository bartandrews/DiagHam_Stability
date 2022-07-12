////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class for n dimensional rational vector using long long         //
//                                                                            //
//                        last modification : 17/11/2010                      //
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


#ifndef LONGRATIONALVECTOR_H
#define LONGRATIONALVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/LongRational.h"

#include <iostream>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::ostream;


class LongRationalVector : public Vector
{

  friend class LongRationalMatrix;

 protected:

  // array for the vector components
  LongRational* Components;

  // garbage flag to avoid data duplication
  GarbageFlag Flag;

 public:

  // default constructor  
  // 
  LongRationalVector();

  // default constructor for an empty vector
  // 
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  LongRationalVector(int size, bool zeroFlag = false);

  // constructor for an empty rational vector bigger than 2^31
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  LongRationalVector(long size, bool zeroFlag = false);

  // constructor from integer arrays
  //
  // numerators = numerator array 
  // denominators = denominator array
  // size = vector Dimension 
  LongRationalVector(long* numerators, long* denominators, long size);

  // constructor from an array of long rational
  //
  // array = array of long rational to become Components of vector
  // size = Vector Dimension
  LongRationalVector(LongRational* array, long size);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  LongRationalVector(const LongRationalVector& vector, bool duplicateFlag = false);

  // destructor
  //
  ~LongRationalVector ();

  // Resize vector
  //
  // dimension = new dimension
  virtual void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  virtual void ResizeAndClean (int dimension);

  // put all vector components to zero
  //
  // return value = reference on current vector
  virtual Vector& ClearVector ();

  // put select vector components to zero
  // start = start index
  // nbrComponent = number of components to set to zero
  // return value = reference on current vector
  virtual Vector& ClearVectorSegment (long start, long nbrComponent);

  // copy a vector into another
  //
  // vector = vector to copy
  // return value = reference on current vector
  LongRationalVector& Copy (LongRationalVector& vector);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  LongRationalVector& Copy (LongRationalVector& vector, const LongRational& coefficient);

  // create a new vector with same size and same type but non-initialized components
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* EmptyClone(bool zeroFlag = false);

  // create an array of new vectors with same size and same type but non-initialized components
  //
  // nbrVectors = number of vectors to sreate
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to the array of new vectors
  virtual Vector* EmptyCloneArray(int nbrVectors, bool zeroFlag = false);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  virtual LongRational& operator [] (int i);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  virtual LongRational& operator [] (long i);

  // get the numerator of a given component
  // 
  // index = component index
  // return value = denominator
  long Num(int index);

  // get the numerator of a given component
  // 
  // index = component index
  // return value = denominator
  long Num(long index);

  // get the denominator of a given component
  // 
  // index = component index
  // return value = denominator
  long Den(int index);

  // get the denominator of a given component
  // 
  // index = component index
  // return value = denominator
  long Den(long index);

  // test if the vector is a null vector
  //
  // return value = true if the vector is a null vector
  bool IsNullVector();

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  LongRationalVector& operator += (LongRationalVector& vector);

  // substract two vectors
  //
  // vector = vector to substract
  // return value = reference on current vector
  LongRationalVector& operator -= (LongRationalVector& vector);

  // multiply a vector with a rational number on the right hand side
  //
  // d = rational to use
  // return value = reference on current vector
  LongRationalVector& operator *= (const LongRational& d);

  // divide a vector by a rational number on the right hand side
  //
  // d = rational to use
  // return value = reference on current vector
  LongRationalVector& operator /= (const LongRational& d);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  LongRationalVector& AddLinearCombination (const LongRational& x, LongRationalVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  LongRationalVector& AddLinearCombination (const LongRational& x, LongRationalVector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  LongRationalVector& AddLinearCombination (const LongRational& x1, LongRationalVector& v1, const LongRational& x2, LongRationalVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  LongRationalVector& AddLinearCombination (const LongRational& x1, LongRationalVector& v1, const LongRational& x2, 
					    LongRationalVector& v2, int firstComponent, int nbrComponent);

  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (const char* fileName);

  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (const char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (const char* fileName);

  // read vector from a file, only within a given range of indices
  //
  // fileName = name of the file where the vector has to be read
  // minIndex = index of the first component to read
  // maxIndex = index of the last component to read
  // return value = true if no error occurs
  bool ReadVector (const char* fileName, long minIndex, long maxIndex);

  // read vector dimension from a file, without loading the full vector 
  //
  // fileName = name of the file where the vector has to be read
  // return value = vector dimension
  long ReadVectorDimension (const char* fileName);

  // test if a vector can be read from a file (matching the right type), without loading the full vector 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if the vector can be read
  bool ReadVectorTest (const char* fileName);

  // Output stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, LongRationalVector& v);

  // print a single component
  //
  // str = reference on output stream
  // index = index of the component to print
  // return value  = reference on output stream
  virtual ostream& PrintComponent(ostream& str, long index);

};
 

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline LongRational& LongRationalVector::operator [] (int i)
{
  return this->Components[i];
}

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline LongRational& LongRationalVector::operator [] (long i)
{
  return this->Components[i];
}

// get the numerator of a given component
// 
// index = component index
// return value = denominator

inline long LongRationalVector::Num(int index)
{
  return this->Components[index].Num();
}

// get the numerator of a given component
// 
// index = component index
// return value = denominator

inline long LongRationalVector::Num(long index)
{
  return  this->Components[index].Num();
}

// get the denominator of a given component
// 
// index = component index
// return value = denominator

inline long LongRationalVector::Den(int index)
{
  return  this->Components[index].Den();
}

// get the denominator of a given component
// 
// index = component index
// return value = denominator

inline long LongRationalVector::Den(long index)
{
  return this->Components[index].Den();
}

// print a single component
//
// str = reference on output stream
// index = index of the component to print
// return value  = reference on output stream

inline ostream& LongRationalVector::PrintComponent(ostream& str, long index)
{
  str << this->Components[index];
  return str;
}

#endif
