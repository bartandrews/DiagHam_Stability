////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//             where only part of the elements are actually stored            //
//                                                                            //
//                        last modification : 23/12/2007                      //
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


#ifndef PARTIALREALVECTOR_H
#define PARTIALREALVECTOR_H


#include "config.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;


class PartialRealVector : public RealVector
{

 protected:
  
  // index of the first element which is actually stored in the current partial vector
  int IndexShift;
  // size of the full vector
  int RealDimension;

 public:

  // default constructor
  //
  PartialRealVector();

  // constructor for an empty real vector 
  //
  // size = effective vector dimension (i.e. number of stored elements)
  // realSize = size of the full vector
  // indexShift = index of the first element which is actually stored in the current partial vector
  // zeroFlag = true if all coordinates have to be set to zero
  PartialRealVector(int size, int realSize, int indexShift, bool zeroFlag = false);

  // constructor from an array of doubles
  //
  // array = array of doubles with real in even position and imaginary part in odd position
  // size = effective vector dimension (i.e. number of stored elements)
  // realSize = size of the full vector
  // indexShift = index of the first element which is actually stored in the current partial vector
  PartialRealVector(double* array, int size, int realSize, int indexShift);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  PartialRealVector(const PartialRealVector& vector, bool duplicateFlag = false);

#ifdef __MPI__

  // constructor from informations sent using MPI
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts or sends the vector
  // broadcast = true if the vector is broadcasted
  PartialRealVector(MPI::Intracomm& communicator, int id, bool broadcast = true);
#endif

  // destructor
  //
  ~PartialRealVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  PartialRealVector& operator = (const PartialRealVector& vector);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  PartialRealVector& Copy (PartialRealVector& vector, double coefficient = 1.0);

  // create a new vector with same size and same type but without duplicating datas
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* EmptyClone(bool zeroFlag = false);

  // create an array of new vectors with same size and same type but non-initialized components
  //
  // nbrVectors = number of vectors to sreate
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to the array of new vectors
  Vector* EmptyCloneArray(int nbrVectors, bool zeroFlag = false);

  // Get Vector dimension
  //
  // return value = vector dimension
  virtual int GetVectorDimension();

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  double& operator [] (int i);

  // write vector in a file 
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteVector (char* fileName);

  // write vector in a file in ascii mode
  //
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  bool WriteAsciiVector (char* fileName);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (char* fileName);

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

// Get Vector dimension
//
// return value = vector dimension

inline int PartialRealVector::GetVectorDimension()
{
  return this->RealDimension;
}
  
// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& PartialRealVector::operator [] (int i)
{
  return this->Components[i - this->IndexShift];
}
 

#endif

