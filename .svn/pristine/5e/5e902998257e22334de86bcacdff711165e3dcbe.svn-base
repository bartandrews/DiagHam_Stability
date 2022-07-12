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


#ifndef INTEGERVECTOR_H
#define INTEGERVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::ostream;
using std::cout;
using std::endl;


class IntegerVector : public Vector
{

  friend class IntegerMatrix;

 protected:

  // array for the vector components
  long* Components;

  // garbage flag to avoid data duplication
  GarbageFlag Flag;

 public:

  // default constructor  
  // 
  IntegerVector();

  // default constructor for an empty vector
  // 
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  IntegerVector(int size, bool zeroFlag = false);

  // constructor for an empty integer vector bigger than 2^31
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  IntegerVector(long size, bool zeroFlag = false);

  // constructor from an array of long integer
  //
  // array = array of long integer to become Components of vector
  // size = Vector Dimension
  IntegerVector(long* array, long size);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  IntegerVector(const IntegerVector& vector, bool duplicateFlag = false);

  // destructor
  //
  ~IntegerVector ();

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
  IntegerVector& Copy (IntegerVector& vector);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  IntegerVector& Copy (IntegerVector& vector, const long& coefficient);

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
  virtual long& operator [] (int i);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  virtual long& operator [] (long i);

  // test if the vector is a null vector
  //
  // return value = true if the vector is a null vector
  bool IsNullVector();

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  IntegerVector& operator += (IntegerVector& vector);

  // substract two vectors
  //
  // vector = vector to substract
  // return value = reference on current vector
  IntegerVector& operator -= (IntegerVector& vector);

  // multiply a vector with an integer number on the right hand side
  //
  // d = integer to use
  // return value = reference on current vector
  IntegerVector& operator *= (const long& d);

  // divide a vector by an integer number on the right hand side
  //
  // d = integer to use
  // return value = reference on current vector
  IntegerVector& operator /= (const long& d);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  IntegerVector& AddLinearCombination (const long& x, IntegerVector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  IntegerVector& AddLinearCombination (const long& x, IntegerVector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  IntegerVector& AddLinearCombination (const long& x1, IntegerVector& v1, const long& x2, IntegerVector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = second vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  IntegerVector& AddLinearCombination (const long& x1, IntegerVector& v1, const long& x2, 
				       IntegerVector& v2, int firstComponent, int nbrComponent);

  // compute the scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = scalar product
  friend long operator * (const IntegerVector& V1, const IntegerVector& V2);

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
  friend ostream& operator << (ostream& str, IntegerVector& v);

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

inline long& IntegerVector::operator [] (int i)
{
  return this->Components[i];
}

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline long& IntegerVector::operator [] (long i)
{
  return this->Components[i];
}

// print a single component
//
// str = reference on output stream
// index = index of the component to print
// return value  = reference on output stream

inline ostream& IntegerVector::PrintComponent(ostream& str, long index)
{
  str << this->Components[index];
  return str;
}

// compute the scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = scalar product

inline long operator * (const IntegerVector& V1, const IntegerVector& V2)
{
  long Tmp = 0l;
  for (long i = 0l; i < V1.LargeDimension; ++i)
    {
      Tmp += V1.Components[i] * V2.Components[i];
    }
  return Tmp;
}

#endif
