////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                         class author: Gunnar Moeller                       //
//                                                                            //
//                 class for an array which has unique entries                //
//                                                                            //
//                        last modification : 13/02/2008                      //
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

#ifndef REALUNIQUEARRAY_H
#define REALUNIQUEARRAY_H

#include "config.h"

#include "GeneralTools/GarbageFlag.h"

#include <iostream>

#include <pthread.h>


using std::ofstream;
using std::ifstream;

class RealUniqueArray
{
 protected:
  
  // array with elements
  double *Elements;
  // size of array
  unsigned InternalSize;
  // number of elements stored
  unsigned NbrElements;

  // mutex to lock write access to array Elements
#ifdef __SMP__
  pthread_mutex_t* BufferMutex;
#endif

  // garbage flag
  GarbageFlag Flag;

 public:
  // standard constructor
  // internalSize = minimum table size to allocate
  RealUniqueArray(unsigned internalSize=100);

  // copy constructor
  RealUniqueArray(RealUniqueArray &array, bool duplicateFlag=false);

  // destructor
  ~RealUniqueArray();

  // Insert element
  // element = new element to be inserted
  // returns : index of this element  
  unsigned InsertElement(double element);

  // search entry
  // value = value to be searched for
  // returns : index of the element, or -1 if not found
  unsigned SearchElement(double value);

  // get number of elements
  unsigned GetNbrElements(){ return NbrElements;}

  // empty all elements
  // disallocate = flag indicating whether all memory should be unallocated
  // internalSize = minimum table size to allocate (only used if disallocating)
  void Empty(bool disallocate = false, unsigned internalSize = 100);

  // Access an element
  double& operator [] (unsigned i);

  // Write to file
  // file = open stream to write to
  void WriteArray(ofstream &file);

  // Read from file
  // file = open stream to read from
  void ReadArray(ifstream &file);
};


// return array's i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& RealUniqueArray::operator [] (unsigned i)
{
  return this->Elements[i];
}

#endif
