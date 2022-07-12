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

#ifdef __MPI__
#include <mpi.h>
#endif

#include <iostream>

using std::ostream;

class RealVector;
class ComplexVector;
class RealPtrVector;

class Vector
{

 private:

  friend class RealVector;
  friend class RealPtrVector;
  friend class ComplexVector;
  friend class DelocalizedRealVector;

 protected:

  // dimension as it appears to the user
  int Dimension;
  // dimension of the allocated array
  int TrueDimension;

  // dimension as it appears to the user for vector bigger than 2^31
  long LargeDimension;
  // dimension of the allocated array for vector bigger than 2^31
  long LargeTrueDimension;

  // flag indacting vector type
  int VectorType;

  // vector id
  int VectorId;


 public:

  enum Type
    {
      RealDatas = 0x01,
      ComplexDatas = 0x02,      
      RealPtrDatas = 0x04,
      RationalData = 0x08,
      LongRationalData = 0x1008,
      IntegerData = 0x10000,
      LongIntegerData = 0x11000,
      DataTypeMask = 0x0f,
      NonLocalDatas = 0x10,
      DistributedDatas = 0x20,
      PartialData = 0x100,
      LargeData = 0x200
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
  int GetVectorDimension() const;
  
  // Get Vector dimension for vector bigger than 2^31
  //
  // return value = vector dimension
  long GetLargeVectorDimension() const;
  
  // check if the vector is bigger than 2^31
  //
  // return value = true if the vector is bigger than 2^31
  bool IsLargeVector() const;
  
  // Get Vector type
  //
  // return value = flag indicating vector type
  int GetVectorType() const;
  
  // get vector id
  //
  // return value = vector id
  virtual int GetVectorId() const;

  // set vector id
  //
  // id = vector new id
  virtual void SetVectorId(int id);

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

  // sum two vectors
  //
  // vector = vector to add
  // return value = reference on current vector
  virtual Vector& operator += (Vector& vector);

  // add a linear combination to a given vector
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  Vector& AddLinearCombination (double x, Vector& V);

  // add a linear combination to a given vector, for a given range of indices
  //
  // x = multiplicative coefficient
  // V = vector to add
  // return value = reference on current vector
  Vector& AddLinearCombination (double x, Vector& V, int firstComponent, int nbrComponent);

  // add a linear combination of two vectors to a given vector
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // return value = reference on current vector
  Vector& AddLinearCombination (double x1, Vector& v1, double x2, Vector& v2);

  // add a linear combination of two vectors to a given vector, for a given range of indices
  //
  // x1 = multiplicative coefficient of first vector
  // v1 = first vector to add
  // x2 = multiplicative coefficient of first vector
  // v2 = first vector to add
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on current vector
  Vector& AddLinearCombination (double x1, Vector& v1, double x2, 
				Vector& v2, int firstComponent, int nbrComponent);

  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, Vector& v);

  // print a single component
  //
  // str = reference on output stream
  // index = index of the component to print
  // return value  = reference on output stream
  virtual ostream& PrintComponent(ostream& str, long index);

  // output the vector in a sparse display
  //
  // str = reference on output stream
  // error = numerical accuracy below which a vector component is considered to be equal to zero
  // return value = reference on output stream  
  virtual ostream& PrintNonZero(ostream& str, double error = MACHINE_PRECISION);

  // output the vector in a sparse display, using labels for the component indices
  //
  // str = reference on output stream
  // componentLabels = array of labels for the component indices
  // error = numerical accuracy below which a vector component is considered to be equal to zero
  // return value = reference on output stream  
  virtual ostream& PrintNonZero(ostream& str, char** componentLabels, double error = MACHINE_PRECISION);

  // localize the current vector to the current process
  // 
  virtual void Localize() const;

  // delocalize the current vector from the current process
  // 
  // transfertFlag = indicates if the current vector datas have to sent to the vector real location
  virtual void Delocalize(bool transfertFlag = false) const;

#ifdef __MPI__

  // send a vector to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current vector
  virtual Vector& SendVector(MPI::Intracomm& communicator, int id);

  // broadcast a vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // return value = reference on the current vector
  virtual Vector& BroadcastVector(MPI::Intracomm& communicator,  int id);

  // broadcast part of vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
  // nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
  // return value = reference on the current vector
  virtual Vector& BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent = 0, int nbrComponent = 0);

  // receive a vector from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current vector
  virtual Vector& ReceiveVector(MPI::Intracomm& communicator, int id);

  // add current vector to the current vector of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  virtual Vector& SumVector(MPI::Intracomm& communicator, int id);

  // reassemble vector from a scattered one
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  virtual Vector& ReassembleVector(MPI::Intracomm& communicator, int id);

  // create a new vector on each MPI node which is an exact clone of the broadcasted one
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* BroadcastClone(MPI::Intracomm& communicator, int id);

  // create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // firstComponent = index of the first component 
  // nbrComponent = number of component to send
  // return value = reference on the current vector
  virtual Vector& SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent);

  // create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* ReceivePartialClone(MPI::Intracomm& communicator, int id);

  // scatter this vector across all MPI nodes with the given load balancing information
  // 
  // communicator = reference on the communicator to use
  // mininumIndices = lowest index for each thread
  // maximumIndices = largest index for each thread
  // id = id of the process to send the vector
  // return value = reference on the current vector
  Vector& ScatterPartialClones(MPI::Intracomm& communicator, long *minimumIndices, long *maximumIndices, int id);

  // create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
  // using efficient implementation with Scatterv
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which scatters the vector
  // return value = pointer to new vector 
  Vector* ReceiveScatteredClone(MPI::Intracomm& communicator, int id);


  // create a new vector on each MPI node with same size and same type but non-initialized components
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  virtual Vector* BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag = false);

#endif

};
 
// Get Vector dimension
//
// return value = vector dimension

inline int Vector::GetVectorDimension() const
{
  return this->Dimension;
}
  
// Get Vector dimension for vector bigger than 2^31
//
// return value = vector dimension

inline long Vector::GetLargeVectorDimension() const
{
  return this->LargeDimension;
}
  
// check if the vector is bigger than 2^31
//
// return value = true if the vector is bigger than 2^31

inline bool Vector::IsLargeVector() const
{
  if ((this->VectorType & Vector::LargeData) == 0)
    return false;
  else
    return true;
}

// Get Vector type
//
// return value = flag indicating vector type

inline int Vector::GetVectorType() const
{
  return this->VectorType;
}
  
// get vector id
//
// return value = vector id

inline int Vector::GetVectorId() const
{
  return this->VectorId;
}

// set vector id
//
// id = vector new id

inline void Vector::SetVectorId(int id)
{
  this->VectorId = id;
}


// localize the current vector on the current process
// 
inline void Vector::Localize() const
{
}

// delocalize the current vector on the current process
// 
// transfertFlag = indicates if the current vector datas have to sent to the vector real location

inline void Vector::Delocalize(bool transfertFlag) const
{
}



#endif
