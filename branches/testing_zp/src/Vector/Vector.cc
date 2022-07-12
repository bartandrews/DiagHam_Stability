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
#include "Vector/PartialRealVector.h"
#include "Vector/PartialComplexVector.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::cout;
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

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors

Vector* Vector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
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

Vector& Vector::operator += (Vector& vector)
{
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

Vector& Vector::AddLinearCombination (double x, Vector& V)
{
  switch ((V.VectorType & this->VectorType) & Vector::DataTypeMask)
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x, (RealVector&) V);
      break;
    case (Vector::ComplexDatas):
      return ((ComplexVector&) (*this)).AddLinearCombination (x, (ComplexVector&) V);
      break;
    default:
      cout <<"Vector type not recognized!"<<endl;
      return *this;
    }
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

Vector& Vector::AddLinearCombination (double x, Vector& V, int firstComponent, int nbrComponent)
{
  switch ((V.VectorType & this->VectorType)  & Vector::DataTypeMask)
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x, (RealVector&) V, firstComponent, nbrComponent);
      break;
    case (Vector::ComplexDatas):
      return ((ComplexVector&) (*this)).AddLinearCombination (x, (ComplexVector&) V, firstComponent, nbrComponent);      break;
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

Vector& Vector::AddLinearCombination (double x1, Vector& v1, double x2, Vector& v2)
{
  switch (((v1.VectorType & this->VectorType) & v2.VectorType)  & Vector::DataTypeMask)
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x1, (RealVector&) v1, x1, (RealVector&) v2);
      break;
    case (Vector::ComplexDatas):
      return ((ComplexVector&) (*this)).AddLinearCombination (x1, (ComplexVector&) v1, x1, (ComplexVector&) v2);      
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

Vector& Vector::AddLinearCombination (double x1, Vector& v1, double x2, 
				      Vector& v2, int firstComponent, int nbrComponent)
{
  switch (((v1.VectorType & this->VectorType) & v2.VectorType)  & Vector::DataTypeMask)
    {
    case (Vector::RealDatas):
      return ((RealVector&) (*this)).AddLinearCombination (x1, (RealVector&) v1, x1, (RealVector&) v2, firstComponent, nbrComponent);
      break;
    case (Vector::ComplexDatas):
      return ((ComplexVector&) (*this)).AddLinearCombination (x1, (ComplexVector&) v1, x1, (ComplexVector&) v2, firstComponent, nbrComponent);
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

ostream& operator << (ostream& str, Vector& v)
{
  switch (v.VectorType & Vector::DataTypeMask)
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

#ifdef __MPI__

// send a vector to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& Vector::SendVector(MPI::Intracomm& communicator, int id)
{
  return *this;
}

// broadcast a vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// return value = reference on the current vector

Vector& Vector::BroadcastVector(MPI::Intracomm& communicator,  int id)
{
  return *this;
}

// broadcast part of vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
// nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
// return value = reference on the current vector

Vector& Vector::BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  return *this;
}

// receive a vector from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current vector

Vector& Vector::ReceiveVector(MPI::Intracomm& communicator, int id)
{
  return *this;
}

// add current vector to the current vector of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& Vector::SumVector(MPI::Intracomm& communicator, int id)
{
  return *this;
}

// create a new vector on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* Vector::BroadcastClone(MPI::Intracomm& communicator, int id)
{
  int Type = this->VectorType;
  if (id != communicator.Get_rank())
    {
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      switch (Type & Vector::DataTypeMask)
	{
	case (Vector::RealDatas):
	  return new RealVector(communicator, id);
	  break;
	case (Vector::ComplexDatas):
	  return new ComplexVector(communicator, id);
	  break;
	default:
	  return 0;
	}
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

Vector& Vector::SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  switch (this->VectorType & Vector::DataTypeMask)
    {
    case (Vector::RealDatas):
      return ((RealVector*) this)->SendPartialClone(communicator, id, firstComponent, nbrComponent);
      break;
    case (Vector::ComplexDatas):
      return ((ComplexVector*) this)->SendPartialClone(communicator, id, firstComponent, nbrComponent);
      break;
    default:
      return *this;
    }
  return *this;
}

// create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* Vector::ReceivePartialClone(MPI::Intracomm& communicator, int id)
{
  int Type = this->VectorType;
  if (id != communicator.Get_rank())
    {
      communicator.Recv(&Type, 1, MPI::INT, id, 1);  
      switch (Type & Vector::DataTypeMask)
	{
	case (Vector::RealDatas):
	  return new PartialRealVector(communicator, id, false);
	  break;
	case (Vector::ComplexDatas):
	  return new PartialComplexVector(communicator, id, false);
	  break;
	default:
	  return 0;
	}
    }
  return 0;
}

// create a new vector on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* Vector::BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag)
{
  int Type = this->VectorType;
  communicator.Bcast(&Type, 1, MPI::INT, id);  
  if (id != communicator.Get_rank())
    {
      switch (Type & Vector::DataTypeMask)
	{
	case (Vector::RealDatas):
	  return new RealVector(communicator, id);
	  break;
	case (Vector::ComplexDatas):
	  return new ComplexVector(communicator, id);
	  break;
	default:
	  return 0;
	}
    }
  return 0;
}

#endif

