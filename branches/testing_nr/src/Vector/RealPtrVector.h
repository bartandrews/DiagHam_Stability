////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//   class for n dimensional vector each entry pointing to a real number       //
//                                                                            //
//                        last modification : 10/11/2005                      //
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


#ifndef REALPTRVECTOR_H
#define REALPTRVECTOR_H


#include "config.h"
#include "Vector/Vector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ifstream;
using std::ofstream;


class RealVector;

class RealPtrVector : public Vector
{

  friend class RealVector;
  friend class Vector;

  friend double operator * (RealPtrVector& V1, RealVector& V2);
  friend double operator * (RealVector& V1, RealPtrVector& V2);

 protected:
  
  double **Components;
  GarbageFlag Flag;

 public:

  // default constructor
  //
  RealPtrVector();

  // constructor for an empty real vector 
  //
  // size = Vector Dimension 
  // zeroFlag = true if all coordinates have to be set to zero
  RealPtrVector(int size, bool zeroFlag = false);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  RealPtrVector(const RealPtrVector& vector, bool duplicateFlag = false);

  // copy constructor from a vector (duplicate datas if necessary)
  //
  // vector = vector to copy
  RealPtrVector(const Vector& vector);

#ifdef __MPI__
  // constructor from informations sent using MPI
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts or sends the vector
  // broadcast = true if the vector is broadcasted
  RealPtrVector(MPI::Intracomm& communicator, int id, bool broadcast = true);
#endif

  // destructor
  //
  virtual ~RealPtrVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  RealPtrVector& operator = (const RealPtrVector& vector);

  // assignement from a vector (duplicate datas if necessary)
  //
  // vector = vector to assign
  // return value = reference on current vector
  RealPtrVector& operator = (const Vector& vector);

  // Resize vector
  //
  // dimension = new dimension
  void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (int dimension);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  RealPtrVector& Copy (RealPtrVector& vector);



  // scalar product between two vectors
  //
  // V1 = first vector
  // V2 = second vector
  // return value = result of scalar product
  friend double operator * (RealPtrVector& V1, RealPtrVector& V2);


  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  double*& operator [] (int i);

  // Output Stream overload
  //
  // str = reference on output stream
  // v = vector to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, RealPtrVector& v);


#ifdef __MPI__

  // send a vector to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& SendVector(MPI::Intracomm& communicator, int id);

  // broadcast a vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // return value = reference on the current vector
  Vector& BroadcastVector(MPI::Intracomm& communicator,  int id);

  // broadcast part of vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
  // nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
  // return value = reference on the current vector
  Vector& BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent = 0, int nbrComponent = 0);

  // receive a vector from a MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the source MPI process
  // return value = reference on the current vector
  Vector& ReceiveVector(MPI::Intracomm& communicator, int id);

  // add current vector to the current vector of a given MPI process
  // 
  // communicator = reference on the communicator to use 
  // id = id of the destination MPI process
  // return value = reference on the current vector
  Vector& SumVector(MPI::Intracomm& communicator, int id);

  // create a new vector on each MPI node which is an exact clone of the broadcasted one
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* BroadcastClone(MPI::Intracomm& communicator, int id);

  // create a new vector on each MPI node with same size and same type but non-initialized components
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag = false);

#endif

};

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double*& RealPtrVector::operator [] (int i)
{
  return this->Components[i];
}
 

#endif
