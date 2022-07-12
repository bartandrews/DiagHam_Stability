////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2008 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//     implements an array of small unsigned integers in packed storage       //
//                                                                            //
//                        last modification : 16/04/2008                      //
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

#include "SmallIntegerArray.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <iostream>
using std::cout;
using std::endl;

// default constructor
//
SmallIntegerArray::SmallIntegerArray(int largestInteger)
{
  this->NbrEntries=0;
  this->NbrWords=0;  
  this->NbrBitsPerEntry=getHighestBit(largestInteger);
}

// constructor for an empty array
// nbrEntries = length of array, set to be >=1
// largestInteger = largest integer to hold in each field
//
SmallIntegerArray::SmallIntegerArray(int nbrEntries, int largestInteger)
{
  this->NbrBitsPerEntry=getHighestBit(largestInteger);
  this->NbrEntries=( nbrEntries>0 ? nbrEntries : 1);
  int EntriesPerWord=32/NbrBitsPerEntry;
  this->NbrWords=NbrEntries/EntriesPerWord;
  if (NbrEntries%EntriesPerWord!=0) ++this->NbrWords;
  this->InternalArray=new unsigned [this->NbrWords];
}

// constructor from given content as integers
// nbrEntries = length of array
// largestInteger = largest integer to hold in each field
// allEntries = array with entries to be stored
// 
SmallIntegerArray::SmallIntegerArray(int nbrEntries, int largestInteger, unsigned *allEntries)
{
  this->NbrBitsPerEntry=getHighestBit(largestInteger);
  int EntriesPerWord=32/NbrBitsPerEntry;
  this->NbrEntries=(nbrEntries>0 ? nbrEntries : 1);
  this->NbrWords=NbrEntries/EntriesPerWord;
  if (NbrEntries%EntriesPerWord!=0) ++this->NbrWords;
  this->InternalArray=new unsigned [this->NbrWords];
  this->SetElements(allEntries);
}

// copy constructor (full duplication)
SmallIntegerArray::SmallIntegerArray( const SmallIntegerArray &array)
{
  this->NbrWords=array.NbrWords;
  this->InternalArray=new unsigned [this->NbrWords];
  for (int i=0; i<NbrWords; ++i)
    this->InternalArray[i]=array.InternalArray[i];
  this->NbrEntries=array.NbrEntries;
  this->NbrBitsPerEntry=array.NbrBitsPerEntry;
}

//destructor
SmallIntegerArray::~SmallIntegerArray()
{
  if (NbrEntries>0)
    delete [] this->InternalArray;
}


// augment an array by an additional element
// array = initial part of array
// toAppend = new entry
SmallIntegerArray::SmallIntegerArray( const SmallIntegerArray &array, unsigned toAppend)
{
  this->NbrEntries = array.NbrEntries+1;
  this->NbrBitsPerEntry=array.NbrBitsPerEntry;
  int EntriesPerWord=32/array.NbrBitsPerEntry;
  this->NbrWords=NbrEntries/EntriesPerWord;
  if (NbrEntries%EntriesPerWord!=0) ++this->NbrWords;
  this->InternalArray=new unsigned [this->NbrWords];
  for (int i=0; i<array.NbrWords; ++i)
    this->InternalArray[i]=array.InternalArray[i];
  this->SetElement(NbrEntries-1,toAppend);
}

// assignment operator
//
// array = array to assign
// return value = reference on current array
SmallIntegerArray& SmallIntegerArray::operator = (const SmallIntegerArray& array)
{
  if (this->NbrEntries>0)
    delete [] this->InternalArray;
  this->NbrWords=array.NbrWords;
  this->InternalArray=new unsigned [this->NbrWords];
  for (int i=0; i<NbrWords; ++i)
    this->InternalArray[i]=array.InternalArray[i];
  this->NbrEntries=array.NbrEntries;
  this->NbrBitsPerEntry=array.NbrBitsPerEntry;
  return *this;
}

// access function: reading
void SmallIntegerArray::GetElement(int i, unsigned &value)
{
  int EntriesPerWord=32/NbrBitsPerEntry;
  int Word = i/EntriesPerWord;
  int Entry = i%EntriesPerWord;
  unsigned Mask = ( ( ( (0x1u << ((Entry+1)*NbrBitsPerEntry-1)) - 1) | (0x1u << ((Entry+1)*NbrBitsPerEntry-1)) ) ^ ((0x1u << ((Entry*NbrBitsPerEntry))) -1));
  value =(InternalArray[Word] & Mask) >> (Entry*NbrBitsPerEntry);
}

unsigned SmallIntegerArray::GetElement(int i)
{
  int EntriesPerWord=32/NbrBitsPerEntry;
  int Word = i/EntriesPerWord;
  int Entry = i%EntriesPerWord;
  unsigned Mask = ( ( ( (0x1u << ((Entry+1)*NbrBitsPerEntry-1)) - 1) | (0x1u << ((Entry+1)*NbrBitsPerEntry-1)) ) ^ ((0x1u << ((Entry*NbrBitsPerEntry))) -1));
  return (InternalArray[Word] & Mask) >> (Entry*NbrBitsPerEntry);
}

// access function: writing
void SmallIntegerArray::SetElement(int i, unsigned value)
{
  int EntriesPerWord=32/NbrBitsPerEntry;
  int Word = i/EntriesPerWord;
  int Entry = i%EntriesPerWord;
  unsigned Mask = ( ( ( (0x1u << ((Entry+1)*NbrBitsPerEntry-1)) - 1) | (0x1u << ((Entry+1)*NbrBitsPerEntry-1)) ) ^ ((0x1u << ((Entry*NbrBitsPerEntry))) -1));
  InternalArray[Word] &= ~Mask;
  InternalArray[Word] |= (Mask&(value<<Entry*NbrBitsPerEntry));
}


// fast access for all elements
// values = array of unsigned integers, supposed to be sufficiently large for answer
void SmallIntegerArray::GetElements(unsigned *values)
{
  int EntriesPerWord=32/NbrBitsPerEntry;  
  int TotalCount=0;
  for (int w=0; w<NbrWords; ++w)
    {
      int Mask = ((0x1u << NbrBitsPerEntry) - 1);
      int Shift=0;
      for (int i=0; (i<EntriesPerWord)&&(TotalCount<NbrEntries); ++i, ++TotalCount)
      {
	values[TotalCount]=(InternalArray[w] & Mask)>>Shift;
	Mask <<= NbrBitsPerEntry;
	Shift += NbrBitsPerEntry;
      }
    }
}

// fast access for all elements
// values = array of unsigned integers, supposed to be sufficiently large for answer
void SmallIntegerArray::SetElements(unsigned *values)
{
  int EntriesPerWord=32/NbrBitsPerEntry;  
  int TotalCount=0;  
  for (int w=0; w<NbrWords; ++w)
    {
      InternalArray[w]=0x0u;
      int Mask = ((0x1u << NbrBitsPerEntry) -1);
      int Shift=0;
      for (int i=0; (i<EntriesPerWord)&&(TotalCount<NbrEntries); ++i, ++TotalCount)
      {
	InternalArray[w] |= (values[TotalCount]&Mask)<<Shift;
	Shift += NbrBitsPerEntry;
      }
    }
}

bool operator == (const SmallIntegerArray& a1,const SmallIntegerArray& a2)
{
  bool result=(a1.NbrWords==a2.NbrWords);
  int i=0;
  while ((result)&&(i<a1.NbrWords))
    {
      result = (a1.InternalArray[i]==a2.InternalArray[i]);
      ++i;
    }
  return result;
}

bool operator != (const SmallIntegerArray& a1,const SmallIntegerArray& a2)
{
  bool result=(a1.NbrWords!=a2.NbrWords);
  int i=0;
  while ((result==false)&&(i<a1.NbrWords))
    {
      result = (a1.InternalArray[i]!=a2.InternalArray[i]);
      ++i;
    }
  return result;
}

bool operator < (const SmallIntegerArray& a1,const SmallIntegerArray& a2)
{
  if (a1.NbrWords<a2.NbrWords) return true;
  else if (a1.NbrWords<a2.NbrWords)
    return false;
  else // a1.NbrWords==a2.NbrWords
    {
      for (int w=a1.NbrWords-1;w>=0; --w)
	if (a1.InternalArray[w]<a2.InternalArray[w]) return true;
	else if (a1.InternalArray[w]>a2.InternalArray[w]) return false;
      return false; // equality!
    }  
}

bool operator > (const SmallIntegerArray& a1,const SmallIntegerArray& a2)
{
  if (a1.NbrWords>a2.NbrWords) return true;
  else if (a1.NbrWords>a2.NbrWords)
    return false;
  else // a1.NbrWords==a2.NbrWords
    {
      for (int w=a1.NbrWords-1;w>=0; --w)
	if (a1.InternalArray[w]>a2.InternalArray[w]) return true;
	else if (a1.InternalArray[w]<a2.InternalArray[w]) return false;
      return false; // equality!
    }
}

bool operator <= (const SmallIntegerArray& a1,const SmallIntegerArray& a2)
{
  return (!(a1>a2));
}

bool operator >= (const SmallIntegerArray& a1,const SmallIntegerArray& a2)
{
  return (!(a1<a2));
}

ostream& operator << (ostream & Str, SmallIntegerArray& a)
{  
  if (a.NbrEntries>0)
    {
      Str << a.GetElement(0);
      for (int i=1; i<a.NbrEntries; ++i)
	Str << " " << a.GetElement(i);
    }
  else
    Str << "--";
  return Str;
}
