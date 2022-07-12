////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                    base class for n dimensional vector                     //
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


#ifndef VECTOR_H
#define VECTOR_H


#include "config.h"

#include <iostream>


using std::ostream;


class RealVector;
class ComplexVector;


class Vector
{

 private:

  friend class RealVector;
  friend class ComplexVector;

 protected:
  
  int Dimension;
  int TrueDimension;

  int VectorType;

 public:

  enum Type
    {
      RealDatas = 0x01,
      ComplexDatas = 0x02
    };

  // virtual destructor
  //
  virtual ~Vector ();

  // get vector norm
  //
  // return value = vector norm
  virtual double Norm();
  
  // get square of vector norm
  //
  // return value = square of vector norm
  virtual double SqrNorm ();
  
  // Get Vector dimension
  //
  // return value = vector dimension
  int GetVectorDimension();
  
  // Get Vector type
  //
  // return value = flag indicating vector type
  int GetVectorType();
  
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

  // create a new vector with same size and same type but non-initialized components
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* EmptyClone(bool zeroFlag = false);

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  virtual Vector& operator += (const Vector& vector);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  Vector& AddLinearCombination (double x, const Vector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  Vector& AddLinearCombination (double x, const Vector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  Vector& AddLinearCombination (double x1, const Vector& v1, double x2, const Vector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  Vector& AddLinearCombination (double x1, const Vector& v1, double x2, 
				const Vector& v2, int firstComponent, int nbrComponent);

  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, const Vector& v);

};
 
// Get Vector dimension
//
// return value = vector dimension

inline int Vector::GetVectorDimension()
{
  return this->Dimension;
}
  
// Get Vector type
//
// return value = flag indicating vector type

inline int Vector::GetVectorType()
{
  return this->VectorType;
}
  


#endif

