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


#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


using std::endl;


// virtual destructor
//

Vector::~Vector ()
{
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* Vector::EmptyClone(bool zeroFlag)
{
  return 0;
}

// get vector norm
//
// return value = vector norm

double Vector::Norm()
{
  return 0.0;
}
  
// get square of vector norm
//
// return value = square of vector norm

double Vector::SqrNorm ()
{
  return 0.0;
}
  
// Resize vector
//
// dimension = new dimension

void Vector::Resize (int dimension)
{
  return;
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void Vector::ResizeAndClean (int dimension)
{
  return;
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& Vector::ClearVector ()
{
  return *this;
}

// sum two vectors
//
// vector = vector to add
// return value = reference on current vector

Vector& Vector::operator += (const Vector& vector)
{
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

Vector& Vector::AddLinearCombination (double x, const Vector& V)
{
  switch ((V.VectorType & this->VectorType))
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x, (RealVector&) V);
      break;
    case (Vector::ComplexDatas):
      return *this;
      break;
    default:
      return *this;
    }
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

Vector& Vector::AddLinearCombination (double x, const Vector& V, int firstComponent, int nbrComponent)
{
  switch ((V.VectorType & this->VectorType))
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x, (RealVector&) V, firstComponent, nbrComponent);
      break;
    case (Vector::ComplexDatas):
      return *this;
      break;
    default:
      return *this;
    }
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

Vector& Vector::AddLinearCombination (double x1, const Vector& v1, double x2, const Vector& v2)
{
  switch (((v1.VectorType & this->VectorType) & v2.VectorType))
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x1, (RealVector&) v1, x1, (RealVector&) v2);
      break;
    case (Vector::ComplexDatas):
      return *this;
      break;
    default:
      return *this;
    }
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

Vector& Vector::AddLinearCombination (double x1, const Vector& v1, double x2, 
				      const Vector& v2, int firstComponent, int nbrComponent)
{
  switch (((v1.VectorType & this->VectorType) & v2.VectorType))
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x1, (RealVector&) v1, x1, (RealVector&) v2, firstComponent, nbrComponent);
      break;
    case (Vector::ComplexDatas):
      return *this;
      break;
    default:
      return *this;
    }
  return *this;
}

// Output Stream overload
//
// str = reference on output stream
// v = vector to print
// return value = reference on output stream

ostream& operator << (ostream& str, const Vector& v)
{
  switch (v.VectorType)
    {
      case (Vector::RealDatas):
	str << (RealVector&) v;
	break;
    case (Vector::ComplexDatas):
      str << (ComplexVector&) v;
      break;
    default:
      str << "unknown vector type " << v.VectorType << endl; 
    }
  return str;
}

