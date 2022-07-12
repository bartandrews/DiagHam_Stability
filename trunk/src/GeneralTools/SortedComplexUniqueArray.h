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

#ifndef SORTEDCOMPLEXUNIQUEARRAY_H
#define SORTEDCOMPLEXUNIQUEARRAY_H

#include "config.h"

#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"

#ifdef __MPI__
#include <mpi.h>
#endif

#include <iostream>

using std::ostream;
using std::ofstream;
using std::ifstream;

class SortedComplexUniqueArray
{
 public:
  // define the type of indices - int is best for MPI routines
  typedef int ElementIndexType;

 protected:
  // array with elements
  Complex *Elements;
  // tolerance for taking elements to be the same (squared)
  double ToleranceSqr;
  // size of array
  ElementIndexType InternalSize;
  // number of elements stored
  ElementIndexType NbrElements;

  // garbage flag
  GarbageFlag Flag;

  // flag indicating how many entries have been sorted
  ElementIndexType Sorted;

  // Universal ID, for some checking of MPI communications
  static const int UniversalID = 31415;

  // flag indicating whether to keep elements sorted
  bool KeepSorted;

  // flag indicating whether current order needs to be preserved
  bool KeepOrder;

 public:
  // standard constructor
  SortedComplexUniqueArray(double tolerance = MACHINE_PRECISION, ElementIndexType internalSize=128, bool keepSorted=true);

  // copy constructor
  SortedComplexUniqueArray(SortedComplexUniqueArray &array, bool duplicateFlag = false);

#ifdef __MPI__

  // constructor array from informations sent using MPI
  //
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts or sends the array
  // broadcast = true if the vector is broadcasted  
  SortedComplexUniqueArray(MPI::Intracomm& communicator, int id, bool broadcast = false);

#endif

  // destructor
  ~SortedComplexUniqueArray();

  // Insert element
  // element = new element to be inserted
  // returns : index of this element  
  ElementIndexType InsertElement(const Complex &element);

  // search entry using a binary search
  // value = value to be searched for
  // @param[out] index : index of the element, or -1 if not found
  // return : true if element was found, false otherwise.
  bool SearchElement(const Complex &value, ElementIndexType &index);

  // search entry, performing a linear search
  // value = value to be searched for
  // @param[out] index : index of the element, or -1 if not found
  // @param enhanceTolerance : increase tolerance by this factor
  // return : true if element was found, false otherwise.
  bool CarefulSearchElement(const Complex &value, ElementIndexType &index, double enhanceTolerance=1.0);

  // search entry closest to a given value
  // value = value to be searched for
  // @param[out] index : index of a nearby element, or the element itself, if found.
  // return : true if the exact element was found, false otherwise.
  bool NearbyEntry(const Complex &value, ElementIndexType &index);

  // get number of elements
  ElementIndexType GetNbrElements(){ return NbrElements;}

  // fix the order of elements that are already in the array
  void FixOrder(){ this->KeepOrder = true;}

  // empty all elements
  // disallocate = flag indicating whether all memory should be unallocated
  // internalSize = minimum table size to allocate (only used if disallocating)
  void Empty(bool disallocate = false, ElementIndexType internalSize = 100);

  // Access an element
  Complex& operator [] (ElementIndexType i);

  // Sort the entries
  void SortEntries();
  
  // Test if the array is sorted
  bool IsSorted();
  
  // Merge data with another UniqueArray
  void MergeArray(SortedComplexUniqueArray &a);

  // Test all entries
  // search for entries of the given array
  // for the search on (this->), make sure their indices match the search result
  // result: true if all entries are found, false otherwise
  bool TestAllEntries() {return this->TestAllEntries(*this);}
  bool TestAllEntries(SortedComplexUniqueArray &a);

  // Resize array to new internal size
  //
  // dimension = new dimension
  void IncreaseInternalSize (ElementIndexType size);

  // Write to file
  // file = open stream to write to
  bool WriteArray(const char*filename);
  bool WriteArray(ofstream &file);

  // Read from file
  // file = open stream to read from
  bool ReadArray(const char*filename);
  bool ReadArray(ifstream &file);
   
  // output stream overload
  friend ostream& operator << (ostream& Str, const SortedComplexUniqueArray &A);

#ifdef __MPI__

  // create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current array

  void SendClone(MPI::Intracomm& communicator, int id);

  // send entries to a given MPI process
  // 
  // communicator = reference on the communicator to use
  // id = id of the destination MPI process
  // return value = reference on the current vector

  void SendArray(MPI::Intracomm& communicator, int id);

  // broadcast a vector to all MPI processes associated to the same communicator
  // 
  // communicator = reference on the communicator to use 
  // id = id of the MPI process which broadcasts the vector
  // return value = true if operation was successful

  bool BroadcastArray(MPI::Intracomm& communicator,  int id);
    
  // merge all data on master node and broadcast to clones
  // 
  // communicator = reference on the communicator to use 
  // return value = reference on the current vector
  
  bool MergeAcrossNodes(MPI::Intracomm& communicator);

#endif // end MPI interface

};

// Test object
void TestClassSortedComplexUniqueArray(SortedComplexUniqueArray::ElementIndexType samples=2048, bool keepSorted=true);


// return array's i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline Complex& SortedComplexUniqueArray::operator [] (ElementIndexType i)
{
  return this->Elements[i];
}

#endif
