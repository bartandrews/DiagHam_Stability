////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                         class author: Gunnar Möller                        //
//                                                                            //
// class for an array which has unique entries for single-threaded insertion  //
//                                                                            //
//                        last modification : 25/12/2016                      //
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


#include "SortedComplexUniqueArray.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::endl;
using std::max;

// flag for testing
#define TESTING_SCUA

// standard constructor
SortedComplexUniqueArray::SortedComplexUniqueArray(double tolerance, ElementIndexType internalSize, bool keepSorted)
{
  if (internalSize<=0)
    internalSize = 4;
  this->InternalSize=internalSize;
  this->ToleranceSqr=tolerance*tolerance;
  this->NbrElements=0;
  if (internalSize>0)
    {
      this->Elements=new Complex[internalSize];
      this->Flag.Initialize();
    }
  this->Sorted = 0;
  this->KeepSorted = keepSorted;
  this->KeepOrder = false;
}

SortedComplexUniqueArray::SortedComplexUniqueArray(SortedComplexUniqueArray &array, bool duplicateFlag)
{
  this->InternalSize=array.InternalSize;
  this->ToleranceSqr=array.ToleranceSqr;
  this->NbrElements=array.NbrElements;
  if (duplicateFlag == false)
    {
      this->Elements = array.Elements;
      this->Flag = array.Flag;
    }
  else
    {
      if (this->InternalSize>0)
	{
	  this->Elements=new Complex[InternalSize];
	  for (ElementIndexType i=0; i<NbrElements; ++i)
	    this->Elements[i]=array.Elements[i];
	  this->Flag.Initialize();
	}
    }
  this->Sorted = array.Sorted;
  this->KeepSorted = array.KeepSorted;
  this->KeepOrder = array.KeepOrder;
}


#ifdef __MPI__

// constructor array from informations sent using MPI
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the array
// broadcast = true if the vector is broadcasted

SortedComplexUniqueArray::SortedComplexUniqueArray(MPI::Intracomm& communicator, int id, bool broadcast)
{
  int TmpArray[6];
  if (broadcast == true)
    communicator.Bcast(TmpArray, 6, MPI::INT, id);      
  else
    communicator.Recv(TmpArray, 6, MPI::INT, id, 1);
  
  int TmpDimension = TmpArray[0];
  this->NbrElements = (ElementIndexType) TmpDimension;
  if (TmpArray[1] != this->UniversalID)
    {
      cout << "Unexpected ID in SortedComplexUniqueArray::ComplexUniqueArray - aborting" <<endl;
      exit(1);
    }
  this->Elements = new Complex [this->NbrElements];
  this->InternalSize = this->NbrElements;
  if (TmpArray[2] == 1)
    for (int i = 0; i < TmpDimension; ++i) 
      this->Elements[i] = 0.0;
  else
    if (TmpArray[2] == 2)
      {
	if (broadcast == true)
	  communicator.Bcast(this->Elements, 2*TmpDimension, MPI::DOUBLE, id);      
	else
	  communicator.Recv(this->Elements, 2*TmpDimension, MPI::DOUBLE, id, 1);   
      }
  if (broadcast == true)
    communicator.Bcast(&this->ToleranceSqr, 1, MPI::DOUBLE, id);      
  else
    communicator.Recv(&this->ToleranceSqr, 1, MPI::DOUBLE, id, 1);
  this->Sorted = (ElementIndexType) TmpArray[3];
  this->KeepSorted = (bool) TmpArray[4];
  this->KeepOrder = (bool) TmpArray[5];
  this->Flag.Initialize();
}

#endif


// destructor
SortedComplexUniqueArray::~SortedComplexUniqueArray()
{
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
}

// Insert element
// element = new element to be inserted
// returns : index of this element  
SortedComplexUniqueArray::ElementIndexType SortedComplexUniqueArray::InsertElement(const Complex& element)
{
  ElementIndexType index;
  if (this->SearchElement(element, index))
    {
      //cout << "Found element"<<element<<" as index "<<index<<endl;
      return index;
    }
  // element not found
  if (NbrElements < InternalSize)
    {
      this->Elements[NbrElements]=element;
      ++NbrElements;      
    }
  else
    {
      if (this->InternalSize < (std::numeric_limits<ElementIndexType>::max() / 2))
	this->InternalSize*=2;
      else
	{
	  if (this->InternalSize == std::numeric_limits<ElementIndexType>::max())
	    {
	      cout << "Array overflow in SortedComplexUniqueArray: cannot store more entries"<<endl;
	      exit(1);
	    }
	  else
	    this->InternalSize = std::numeric_limits<ElementIndexType>::max();
	}
      Complex *newElements= new Complex[InternalSize];
      for (ElementIndexType i=0; i<NbrElements; ++i)
	newElements[i]=Elements[i];
      newElements[NbrElements]=element;
      ++NbrElements;
      Complex *tmpElements=this->Elements;
      this->Elements=newElements;
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	delete [] tmpElements;
    }
  index=NbrElements-1;
  if ((!this->KeepOrder) && this->KeepSorted && this->NbrElements > this->Sorted+12)
    {
      this->SortEntries();
      if (!this->SearchElement(element, index))
	{
	  cout << "Error: did not find the element that was just inserted"<<endl;
	  exit(1);
	}
    }
#ifdef TESTING_SCUA
  ElementIndexType test;
  if (!this->SearchElement(element, test))
    {
      cout << "Error: did not find the element that was just inserted ("<<element<<")"<<endl;
    }
  if (test != index)
    {
      cout << "Error: inconsistent indices for the element that was just inserted ("<<element<<")"<<endl;
    }
#endif

  return index;
}

// search entry
// value = value to be searched for
// @param[out] index : index of the element, if found.
// return : true if element was found, false otherwise.
bool SortedComplexUniqueArray::SearchElement(const Complex &value, ElementIndexType &index)
{
  ElementIndexType start=0;
  if (this->Sorted>3)
    {
      ElementIndexType PosMax = this->Sorted - 1;
      ElementIndexType PosMin = 0;
      ElementIndexType PosMid = (PosMin + PosMax) >> 1;
      Complex CurrentState = this->Elements[PosMid];
      // cout << "Searching "<<value<<"...";
      while ((PosMin != PosMid) && (SqrNorm(CurrentState - value) >= this->ToleranceSqr))
	{
	  if (CurrentState > value)
	    {
	      PosMax = PosMid;
	    }
	  else
	    {
	      PosMin = PosMid;
	    } 
	  PosMid = (PosMin + PosMax) >> 1;
	  CurrentState = this->Elements[PosMid];
	  // cout << "PosMid="<<PosMid<<", CurrentState="<<CurrentState<<", PosMin="<<PosMin<<", PosMax="<<PosMax<< endl;
	}
      if (SqrNorm(CurrentState - value) < this->ToleranceSqr)
	{
	  index = PosMid;
	  // cout << "Found "<<this->Elements[index]<<endl;
	  return true;
	}
      else
	{
	  index = PosMax; // ?? 
	  if (SqrNorm(this->Elements[PosMax] - value) < this->ToleranceSqr)
	    {
	      // cout << "Found "<<this->Elements[PosMax]<<endl;
	      return true;
	    }
	  else 
	    {
	      // cout << "Not found."<<endl;
	      start = this->Sorted-1;
	    }
	}
    }
  // cout << "Sequential search from "<<start<<", this->NbrElements="<<this->NbrElements<<endl;
  for (ElementIndexType i=start; i<this->NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-value)<this->ToleranceSqr)
	{
	  index = i;
	  return true;
	}
    }
  // element not found
  // cout << "element not found"<<endl;
  index = 0;
  return false;
}


// search entry closest to a given value
// value = value to be searched for
// @param[out] index : index of a nearby element, or the element itself, if found.
// return : true if the exact element was found, false otherwise.
bool SortedComplexUniqueArray::NearbyEntry(const Complex &value, ElementIndexType &index)
{
  ElementIndexType start=0;
  double distance = 1e300;
  if (this->Sorted>3)
    {
      ElementIndexType PosMax = this->Sorted - 1;
      ElementIndexType PosMin = 0;
      ElementIndexType PosMid = (PosMin + PosMax) >> 1;
      Complex CurrentState = this->Elements[PosMid];
      // cout << "Searching "<<value<<"...";
      while ((PosMin != PosMid) && (SqrNorm(CurrentState - value) >= this->ToleranceSqr))
	{
	  if (CurrentState > value)
	    {
	      PosMax = PosMid;
	    }
	  else
	    {
	      PosMin = PosMid;
	    } 
	  PosMid = (PosMin + PosMax) >> 1;
	  CurrentState = this->Elements[PosMid];
	  //cout << "PosMid="<<PosMid<<", CurrentState="<<CurrentState<<", PosMin="<<PosMin<<", PosMax="<<PosMax<< endl;
	}
      if (SqrNorm(CurrentState - value) < this->ToleranceSqr)
	{
	  index = PosMid;
	  // cout << "Found "<<this->Elements[index]<<endl;
	  return true;
	}
      else
	{
	  index = PosMax; // ?? 
	  if (SqrNorm(this->Elements[PosMax] - value) < this->ToleranceSqr)
	    {
	      // cout << "Found "<<this->Elements[PosMax]<<endl;
	      return true;
	    }
	  else 
	    {
	      // cout << "Not found."<<endl;
	      start = this->Sorted-1;
	      if ((distance = SqrNorm(this->Elements[PosMax] - value)) > SqrNorm(CurrentState - value))
		{
		  index = PosMid;
		  distance = SqrNorm(CurrentState - value);
		}
	      else
		index = PosMax;
	    }
	}
    }
  for (ElementIndexType i=start; i<this->NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-value)<distance)
	{
	  index = i;
	  if (SqrNorm(Elements[i]-value)<this->ToleranceSqr)
	    return true;
	}
    }
  return false;
}


// search entry, performing a linear search
// value = value to be searched for
// @param[out] index : index of the element, or -1 if not found
// return : true if element was found, false otherwise.
bool SortedComplexUniqueArray::CarefulSearchElement(const Complex &value, ElementIndexType &index, double enhanceTolerance)
{
  double myTolerance  =  enhanceTolerance*enhanceTolerance*this->ToleranceSqr;
  for (ElementIndexType i=0; i<this->NbrElements; ++i)
    {
      if (SqrNorm(Elements[i]-value)<this->ToleranceSqr)
	{
	  index = i;
	  return true;
	}
    }
  // element not found
  index = 0;
  return false;
}


// empty all elements
// disallocate = flag indicating whether all memory should be unallocated
// internalSize = minimum table size to allocate (only used if disallocating)
void SortedComplexUniqueArray::Empty(bool disallocate, ElementIndexType internalSize)
{
  if (disallocate)
    {
      if ((this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete [] Elements;
	}
      this->InternalSize=internalSize;
      this->Elements = new Complex[InternalSize];
    }
  this->NbrElements=0;
}

// apply a simple sort algorithm to the existing entries of the table
void SortedComplexUniqueArray::SortEntries()
{
  if (this->Sorted==this->NbrElements) return;
  ElementIndexType inc = std::floor(NbrElements/2.0 + 0.5);
  // if (this->Sorted>inc) inc=this->Sorted-1;
  Complex tmpC;
  while (inc > 0)
    {
      for (ElementIndexType i = inc; i< NbrElements; ++i)
	{
	  tmpC = this->Elements[i];
	  ElementIndexType j = i;
	  while ((j>=inc) && (this->Elements[j-inc] > tmpC) )
	    {
	      this->Elements[j] = this->Elements[j - inc];
	      j = j - inc;
	    }
	  this->Elements[j] = tmpC;
	}
      inc = std::floor(inc / 2.2 + 0.5);
    }
  this->Sorted=this->NbrElements;

#ifdef TESTING_SCUA
  if (!this->IsSorted())
    {
      cout << "Error: sort algorithm left unsorted array."<<endl;
      exit(1);
    }
#endif //TESTING_SCUA
}


// Test if the array is sorted
bool SortedComplexUniqueArray::IsSorted()
{
  if (this->NbrElements<2) return true;
  for (ElementIndexType i=0; i<this->NbrElements-1; ++i)
    if (this->Elements[i]>=this->Elements[i+1])
      {
	cout << "Unexpected order at position "<<i<<": (this->Elements["<<i<<"] = "<<this->Elements[i]<<" < this->Elements["<<i+1<<"] = "<<this->Elements[i+1]<<", this->NbrElements="<<this->NbrElements<<endl; 
	return false;
      }
  return true;
}


// Merge data with another UniqueArray
void SortedComplexUniqueArray::MergeArray(SortedComplexUniqueArray &a)
{
  a.SortEntries();
  this->SortEntries();

  // cout << "this ="<<*this<<endl;
  // cout << "a ="<<a << endl;
  // cout << "Merging arrays with "<<this->NbrElements<<" and "<<a.NbrElements<<" entries "<<endl;
  long newPos = a.NbrElements+this->NbrElements;
  if (newPos > std::numeric_limits<ElementIndexType>::max())
    {
      cout << "Error merged array size exceeds maximum"<< endl;
      exit(1);
    }
  ElementIndexType tmpInternalSize = (ElementIndexType)newPos;
  Complex *newElements = new Complex[newPos];
  newPos=0;
  ElementIndexType myPos = 0;
  for (ElementIndexType theirPos=0; theirPos<a.NbrElements; ++theirPos)
    {
      // definite insert elements in this-> that are smaller than the next element of a and not within tolerance
      while (myPos < this->NbrElements && this->Elements[myPos] < a.Elements[theirPos] && SqrNorm(this->Elements[myPos] - a.Elements[theirPos]) >= this->ToleranceSqr)
	{
	  // cout << "Insert this->Elements["<<myPos<<"] at 1 with this->Elements["<<myPos<<"] < a.Elements["<<theirPos<<"]"<<endl;
	  newElements[newPos++] = this->Elements[myPos++];
	}
      // also insert any elements that are equal or approximately equal
      while (myPos < this->NbrElements && SqrNorm(this->Elements[myPos] - a.Elements[theirPos]) < this->ToleranceSqr)
	{
	  //  cout << "Insert this->Elements["<<myPos<<"] at 2"<<endl;
	  newElements[newPos++] = this->Elements[myPos++];
	}
      // cout << "myPos = "<<myPos << " this->NbrElements = "<<this->NbrElements << " this->Elements[myPos]="<<this->Elements[myPos]<<", a.Elements[theirPos] ="<<a.Elements[theirPos]<<endl;
      if (newPos > 0)
	{
	  if (SqrNorm(a.Elements[theirPos] - newElements[newPos-1]) < this->ToleranceSqr)
	    {
	      // next element on a is identical within accuracy to the last entry on the new array
	      // cout << "Omitting identical element a.Elements["<<theirPos<<"] at 1"<<endl;
	      // do nothing - counter is increased in loop.
	    }
	  else
	    {
	      //cout << "Insert a.Elements["<<theirPos<<"] at 1"<<endl;
	      newElements[newPos++] = a.Elements[theirPos];
	    }
	}
      else
	{
	  //cout << "Insert a.Elements["<<theirPos<<"] at 2"<<endl;
	  newElements[newPos++] = a.Elements[theirPos];
	}
    }
  
  while (myPos < this->NbrElements) //  && this->Elements[myPos] < a.Elements[theirPos] && SqrNorm(this->Elements[myPos] - a.Elements[theirPos]) >= this->ToleranceSqr)
    {
      // cout << "Insert this->Elements["<<myPos<<"] at 3 with this->Elements["<<myPos<<"]"<<endl;
      newElements[newPos++] = this->Elements[myPos++];
    }
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
  this->Flag = GarbageFlag();
  this->Flag.Initialize(); // start newly allocated array
  this->InternalSize = tmpInternalSize;
  this->Elements = newElements;
  this->NbrElements = newPos;
  // cout << "Unique entries retained: "<<this->NbrElements<<endl;
  this->Sorted=NbrElements;
}

// Test all entries
// search for entries and make sure their indices match the search result
// result: true if all entries are found, false otherwise
bool SortedComplexUniqueArray::TestAllEntries(SortedComplexUniqueArray &a)
{
  ElementIndexType index;
  bool success=true;
  for (ElementIndexType i=0; i<a.GetNbrElements(); ++i)
    {
      if (! this->SearchElement(a.Elements[i], index))
	{
	  cout << "Element " << i << " not found during check"<<endl;
	  success=false;
	}
      if (&a == this)
	if (index != i)
	  {
	    cout << "Discrepancy in search for index "<< index <<" on element " << i << " during self-check."<<endl;
	    success=false;
	  }
    }
  return success;
}


// Resize array to new internal size
//
// dimension = new dimension

void SortedComplexUniqueArray::IncreaseInternalSize (ElementIndexType size)
{
  if (size <= this->InternalSize) // do nothing if size is sufficient
    {
      return;
    }
  Complex* TmpElements = new Complex [size];
  for (int i = 0; i < this->NbrElements; i++)
    TmpElements[i] = this->Elements[i];
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Elements;
    }
  this->InternalSize = size;
  this->Elements = TmpElements;
  this->Flag = GarbageFlag();
  this->Flag.Initialize();
}


// Write to file
// filename = file to open and write to
bool SortedComplexUniqueArray::WriteArray(const char*filename)
{
  ofstream file;
  if (filename!=NULL)
    {
      file.open(filename, std::ios::out | std::ios::binary );
      if (file.is_open())
	return WriteArray(file);
      cout << "Error writing file "<<filename<<endl;
      return false;
    }
  return false;
}


// write to file
// file = open stream to write to
bool SortedComplexUniqueArray::WriteArray(ofstream &file)
{
  WriteLittleEndian(file, this->NbrElements);
  WriteLittleEndian(file, this->ToleranceSqr);
  WriteLittleEndian(file, this->Sorted);
  WriteLittleEndian(file, this->KeepSorted);
  WriteLittleEndian(file, this->KeepOrder);
  for (ElementIndexType i = 0; i < this->NbrElements; ++i)
    WriteLittleEndian(file, this->Elements[i]);  
  return true;
}

// Read from file
// filename = file to open and read from
bool SortedComplexUniqueArray::ReadArray(const char*filename)
{
  ifstream file;
  if (filename!=NULL)
    {
      file.open(filename, std::ios::in | std::ios::binary );
      if (file.is_open())
	return ReadArray(file);
      cout << "Error reading file "<<filename<<endl;
      return false;
    }
  return false;
}


// Read from file
// file = open stream to read from
bool SortedComplexUniqueArray::ReadArray(ifstream &file)
{
  if ( (this->InternalSize!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] Elements;
    }
  ElementIndexType TmpDimension;
  ReadLittleEndian(file, TmpDimension);
  ReadLittleEndian(file, this->ToleranceSqr);
  ReadLittleEndian(file, this->Sorted);
  ReadLittleEndian(file, this->KeepSorted);
  ReadLittleEndian(file, this->KeepOrder);
  this->InternalSize=TmpDimension;
  this->NbrElements=TmpDimension;
  this->Elements=new Complex[TmpDimension];
  for (ElementIndexType i = 0; i < this->NbrElements; ++i)
    ReadLittleEndian(file, this->Elements[i]);
  this->Flag = GarbageFlag();
  this->Flag.Initialize(); // start newly allocated array
  return true;
}



// Output Stream overload
//

ostream& operator << (ostream& Str, const SortedComplexUniqueArray& A)
{
  for (SortedComplexUniqueArray::ElementIndexType i = 0; i < A.NbrElements; ++i)
    Str << A.Elements[i]<<endl;
  return Str;
}


#ifdef __MPI__

// create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current array

void SortedComplexUniqueArray::SendClone(MPI::Intracomm& communicator, int id)
{
  if (2*this->NbrElements > std::numeric_limits<int>::max())
    {
      cout << "Error: cannot send unique arrays larger than max(int)"<<endl;
    }
  int TmpArray[6];
  TmpArray[0] = (int)this->NbrElements;
  TmpArray[1] = this->UniversalID; // an ad-hoc number to be checked
  TmpArray[2] = 2;
  TmpArray[3] = (int)this->Sorted;
  TmpArray[4] = (int)this->KeepSorted;
  TmpArray[5] = (int)this->KeepOrder;
  communicator.Send(TmpArray, 6, MPI::INT, id, 1); 
  communicator.Send(this->Elements, 2*NbrElements, MPI::DOUBLE, id, 1); 
  communicator.Send(&this->ToleranceSqr, 1, MPI::DOUBLE, id, 1); 
}

// send entries to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current vector

void SortedComplexUniqueArray::SendArray(MPI::Intracomm& communicator, int id)
{
  int TmpArray[2] = {(int) this->NbrElements, (int) this->Sorted};
  communicator.Send(&TmpArray, 2, MPI::INT, id, 1); 
  // int Acknowledge = 0;
  // communicator.Recv(&Acknowledge, 1, MPI::INT, id, 1);
  // if (Acknowledge != 0)
  //   return;
  communicator.Send(this->Elements, 2*this->NbrElements, MPI::DOUBLE, id, 1); 
}

// broadcast the entries of the array on node "id" to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the array
// return value = true if operation was successful

bool SortedComplexUniqueArray::BroadcastArray(MPI::Intracomm& communicator,  int id)
{
  if (2*this->NbrElements > std::numeric_limits<int>::max())
    {
      cout << "Error: cannot broadcast unique arrays larger than max(int)"<<endl;
      return false;
    }
  int TmpArray[2] = {(int) this->NbrElements, (int) this->Sorted};
  communicator.Bcast(&TmpArray, 2, MPI::INT, id);
  if (TmpArray[0] > this->InternalSize)
    {
      this->IncreaseInternalSize(TmpArray[0]);
    }
  this->NbrElements = (ElementIndexType) TmpArray[0];
  this->Sorted = (ElementIndexType) TmpArray[1];
  communicator.Bcast(this->Elements, 2*this->NbrElements, MPI::DOUBLE, id);
  return true;
}

// merge all data on master node and broadcast to clones
// 
// communicator = reference on the communicator to use 
// return = true if successfully merged

bool SortedComplexUniqueArray::MergeAcrossNodes(MPI::Intracomm& communicator)
{
  if (communicator.Get_size()==1) return true;
  if (2*this->NbrElements > std::numeric_limits<int>::max())
    {
      cout << "Error: cannot merge unique arrays larger than max(int)"<<endl;
    }
  if (this->KeepOrder)
    {
      cout << "Attention: merge algorithm will scramble order of array entries"<<endl;
    }
  int TmpNbrElements = (int) this->NbrElements;
  
#ifdef TESTING_SCUA
  SortedComplexUniqueArray OldArray(*this, true);
#endif

  int Acknowledge = 0;
  if (communicator.Get_rank() != 0)
    this->SendClone(communicator, 0);
  else
    {
      int NbrMPINodes = communicator.Get_size();
      //cout << "Master="<<*this<<"done Master"<<endl;
      for (int id = 1; id < NbrMPINodes; ++id)
	{
	  SortedComplexUniqueArray TmpArray(communicator, id);
	  //cout << "TmpArray["<<id<<"]="<<TmpArray<<"done Slave "<<id<<endl;
	  this->MergeArray(TmpArray);
	  if (this->NbrElements > std::numeric_limits<int>::max())
	    {
	      cout << "Error: cannot merge unique arrays larger than max(int)"<<endl;
	      return false;
	    }
#ifdef TESTING_SCUA
	  this->TestAllEntries(TmpArray);
#endif
	}
    }
  bool Rst = this->BroadcastArray(communicator, 0);
#ifdef TESTING_SCUA
  this->TestAllEntries(OldArray);
#endif
  return Rst;
}

#endif


// Test object
void TestClassSortedComplexUniqueArray(SortedComplexUniqueArray::ElementIndexType samples, bool keepSorted)
{
  double randNorm=1.0/(double)RAND_MAX;

  double precision = 1e-13;
  SortedComplexUniqueArray a1(precision, samples>>1, keepSorted);
  SortedComplexUniqueArray a2(precision, samples>>1, keepSorted);
  
  // insert half the elements as independent numbers
  for (int i=0; i<samples; ++i)
    {
      a1.InsertElement( Complex(randNorm*(double)std::rand(),randNorm*(double)std::rand()) );
      a2.InsertElement( Complex(randNorm*(double)std::rand(),randNorm*(double)std::rand()) );
      // cout << "a1.GetNbrElements()="<<a1.GetNbrElements()<<endl;
    }
  // count identical entries
  SortedComplexUniqueArray::ElementIndexType common=0, index;
  for (int i=0; i<samples; ++i)
    {
      if (a2.SearchElement(a1[i],index))
	{
	  ++common;
	  cout << "Random common entry at index "<<index<<endl;
	}
    }

  // insert other half as same number
  Complex TmpC;
  for (int i=0; i<samples; ++i)
    {
      TmpC.Re=randNorm*(double)std::rand();
      TmpC.Im=randNorm*(double)std::rand();
      a1.InsertElement(TmpC);
      a2.InsertElement(TmpC);
    }
  
  SortedComplexUniqueArray a3(a1, true);
  SortedComplexUniqueArray a4(a2, true);
  
  cout << "Copied array 3 has "<< a3.GetNbrElements() <<" entries"<<endl;
  cout << "Copied array 4 has "<< a4.GetNbrElements() <<" entries"<<endl;

  a3.MergeArray(a2);

  cout << "Merged array has "<< a3.GetNbrElements() <<" entries"<<endl;
  
  if (a3.GetNbrElements() < 3*samples - common)
    {
      cout << "Unexpected number of entries in a3: "<< a3.GetNbrElements() << " vs " << 3*samples - common << " expected "<<endl;
    }

  for (int i=0; i<2*samples; ++i)
    {
      a4.InsertElement(a1[i]);
    }

  cout << "Manually merged array has "<< a4.GetNbrElements() <<" entries"<<endl;
  
  if (a4.GetNbrElements() < 3*samples - common)
    {
      cout << "Unexpected number of entries in a4: "<< a4.GetNbrElements() << " vs " << 3*samples - common << " expected "<<endl;
    }

  // check all entries are present:
  for (SortedComplexUniqueArray::ElementIndexType i=0; i<a4.GetNbrElements(); ++i)
    {
      if (! a3.SearchElement(a4[i], index))
	{
	  cout << "Element " << i << " of array 4, "<<a4[i] << ", missing from array 3"<<endl;
	}
    }

  for (SortedComplexUniqueArray::ElementIndexType i=0; i<a3.GetNbrElements(); ++i)
    {
      if (! a4.SearchElement(a3[i], index))
	{
	  cout << "Element "<<i<<" of array 3, "<<a3[i] << ", missing from array 4"<<endl;
	}
    }

  // check all entries are found:
  for (SortedComplexUniqueArray::ElementIndexType i=0; i<a3.GetNbrElements(); ++i)
    {
      if (! a3.SearchElement(a3[i], index))
	{
	  cout << "Element " << i << " not found by search on array 3"<<endl;
	}
      if (index != i)
	cout << "Discrepancy in search for index "<< index <<" on element " << i << " (array 3)."<<endl;
    }

  for (SortedComplexUniqueArray::ElementIndexType i=0; i<a4.GetNbrElements(); ++i)
    {
      if (! a4.SearchElement(a4[i], index))
	{
	  cout << "Element " << i << " not found by search on array 4"<<endl;
	}
      if (index != i)
	cout << "Discrepancy in search for index "<< index <<" on element " << i << " (array 4)."<<endl;
    }

    cout << "Tests completed."<<endl;

}


