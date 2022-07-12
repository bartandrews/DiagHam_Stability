////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of Abstract Architecture                      //
//                                                                            //
//                        last modification : 10/04/2002                      //
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


#ifndef ABSTRACTARCHITECTURE_H
#define ABSTRACTARCHITECTURE_H


#include "config.h"


class Vector;
class RealVector;
class ComplexVector;

class AbstractArchitecture
{

 protected:

  // dimension of the Hilbert space on which the architecture has to work
  long HilbertSpaceDimension;

  // architecture ID
  int ArchitectureID;

 public:

  enum ArchitectureIDs
    {
      MonoProcessor = 0x01,
      SMP = 0x02,
      SimpleMPI = 0x30,
      WithCommunicator = 0x10,
      MixedMPISMP = 0x32
    };

  
  // destructor
  //
  virtual ~AbstractArchitecture();

  // get ID of the architecture
  //
  // return value = architecture ID
  int GetArchitectureID();
  
  // get the amount of memory available for the local architecture
  //
  // return value = amount of memory in byte (negative if the information is not available)
  virtual long GetLocalMemory();
  
  // get typical range of indices on which the local architecture acts
  //
  // minIndex = reference on the minimum index on which the local architecture can act
  // maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
  //            architecture doesn't support this feature)
  virtual void GetTypicalRange (long& minIndex, long& maxIndex);
  
  // get a new real vector with memory alloaction depending on the architecture
  //
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual RealVector* GetNewRealVector ();
  
  // get a new real vector with memory alloaction depending on the architecture
  //
  // dimension = dimension of the requested vector
  // zeroFlag = true if all vector entries has to be set to zero
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual RealVector* GetNewRealVector (long dimension, bool zeroFlag = false);
  
  // get a new complex vector with memory alloaction depending on the architecture
  //
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual ComplexVector* GetNewComplexVector ();
  
  // get a new complex vector with memory alloaction depending on the architecture
  //
  // dimension = dimension of the requested vector
  // zeroFlag = true if all vector entries has to be set to zero
  // return value = pointer to the requested vector (zero if an error occurs)
  virtual ComplexVector* GetNewComplexVector (long dimension, bool zeroFlag = false);
  
  // set dimension of the Hilbert space on which the architecture has to work
  // 
  // dimension = dimension of the Hilbert space
  virtual void SetDimension (long dimension);

  // request a given amount of memory for an array of Type element
  //
  // architecture = reference on the architecture to which memory will be asked
  // pointer = reference on pointer that will be used for the memory allocation
  // size = number of elements of type Type
  // return value = reference on the pointer
  template <class Type>
  Type*& New (Type*& pointer,  unsigned long size);
  
  // delete an array that has been requested by the New function
  //
  // architecture = reference on the architecture to which memory has been asked
  // pointer = reference on pointer that of the memory allocation
  template <class Type>
  void Delete (Type*& pointer);
  
  // get a temporary file name
  //
  // return value = string corresponding to a temporary file name
  virtual char* GetTemporaryFileName();
  
 protected:

  // indicate an allocation of memory to the architecture
  //
  // pointer = pointer to the memory zone which will be allocated
  // memory = amount of requested memory in bytes
  virtual void AllocateMemory (void* pointer, unsigned long memory);

  // indicate an deallocation of memory to the architecture
  //
  // pointer = pointer to the memory zone which will be free
  virtual void DeallocateMemory (void* pointer);

};


// request a given amount of memory for an array of Type element
//
// pointer = reference on pointer that will be used for the memory allocation
// size = number of elements of type Type
// return value = reference on the pointer

template <class Type>
inline Type*& AbstractArchitecture::New (Type*& pointer, unsigned long size)
{
  pointer = new Type [size];
  this->AllocateMemory(pointer, size * sizeof(Type));
  return pointer;
}
  
// delete an array that has been requested by the New function
//
// architecture = reference on the architecture to which memory has been asked
// pointer = reference on pointer that of the memory allocation

template <class Type>
inline void AbstractArchitecture::Delete (Type*& pointer)
{
  delete[] pointer;
  this->DeallocateMemory(pointer);
}

// indicate an allocation of memory to the architecture
//
// pointer = pointer to the memory zone which will be allocated
// memory = amount of requested memory in bytes

inline void AbstractArchitecture::AllocateMemory (void* pointer, unsigned long memory)
{
}

// indicate an deallocation of memory to the architecture
//
// pointer = pointer to the memory zone which will be free

inline void AbstractArchitecture::DeallocateMemory (void* pointer)
{
}

// get ID of the architecture
//
// return value = architecture ID

inline int AbstractArchitecture::GetArchitectureID()
{
  return this->ArchitectureID;
}

#endif
