////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//    where elements are stored on demand and tramsmitted to a master vector  //
//                           when the buffer is full                          //
//                                                                            //
//                        last modification : 19/07/2008                      //
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


#ifndef BUFFEREDREALVECTOR_H
#define BUFFEREDREALVECTOR_H


#include "config.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"

#ifdef __SMP__
#include <pthread.h>
#endif
#include <iostream>
#include <fstream>


using std::ostream;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;


class BufferedRealVector : public RealVector
{

 protected:
  
  // buffer size
  int BufferSize;

  // buffer for indices for current operations
  int* CurrentIndexBuffer;
  // buffer for indices that has to be transmitted
  int* TransmittedIndexBuffer;
  // buffer for components for current operations
  int* CurrentComponentBuffer;
  // buffer for components that has to be transmitted
  int* TransmittedComponentBuffer;

  // mutex to lock transmitted buffers
#ifdef __SMP__
  pthread_mutex_t* BufferMutex;
#endif
  // current position if Current*Buffer
  int CurrentPosition;

  // size of the full vector
  int RealDimension;


 public:

  // default constructor
  //
  BufferedRealVector();

  // constructor for an empty real vector 
  //
  // realSize = size of the full vector
  // bufferSize = memory in bytes that can be allocated for buffers
  // zeroFlag = true if all coordinates have to be set to zero
  BufferedRealVector(int realSize, long bufferSize, bool zeroFlag = false);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  PartialRealVector(const BufferedRealVector& vector, bool duplicateFlag = false);

  // destructor
  //
  ~BufferedRealVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  BufferedRealVector& operator = (const BufferedRealVector& vector);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  BufferedRealVector& Copy (BufferedRealVector& vector, double coefficient = 1.0);

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

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  double& operator [] (int i);

 protected:

  // resize buffers
  //
  // size = new number of elements per buffer
  void ResizeBuffer();

  // transmit buffrered component to the master vector
  //
  void TransmitBuffer();

};

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& BufferedRealVector::operator [] (int i)
{
  if (this->CurrentPosition < this->BufferSize)
    {
      this->CurrentIndexBuffer[this->CurrentPosition] = i;      
      return this->CurrentComponentBuffer[this->CurrentPosition++];
    }
  else
    {
#ifdef __SMP__
      pthread_mutex_lock(this->BufferMutex);
#endif
      int* TmpIndexBuffer = this->CurrentIndexBuffer;
      this->CurrentIndexBuffer = this->TransmittedIndexBuffer;
      this->TransmittedIndexBuffer = TmpIndexBuffer;
      double* TmpComponentBuffer = this->CurrentComponentBuffer;
      this->CurrentComponentBuffer = this->TransmittedComponentBuffer;
      this->TransmittedComponentBuffer = TmpComponentBuffer;
#ifdef __SMP__
      pthread_mutex_unlock(this->BufferMutex);
#endif
      this->TransmitBuffer();
      this->CurrentPosition = 0;
      this->CurrentIndexBuffer[this->CurrentPosition] = i;      
      return this->CurrentComponentBuffer[this->CurrentPosition++];      
    }
}
 

#endif

