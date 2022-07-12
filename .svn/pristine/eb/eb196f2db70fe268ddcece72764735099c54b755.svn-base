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


#ifndef SMALLINTEGERARRAY_H
#define SMALLINTEGERARRAY_H

#include <iostream>
using std::ostream;

class SmallIntegerArray
{
 private:
  // internal storage
  unsigned long int *InternalArray;
  // number of words used
  int NbrWords;
  // number of entries
  int NbrEntries;
  // number of bits per entry
  int NbrBitsPerEntry;

 public:

  // default constructor
  // largestInteger = largest integer to hold in each field
  SmallIntegerArray(int largestInteger=1);
  
  // constructor for an empty array
  // nbrEntries = length of array
  // largestInteger = largest integer to hold in each field
  //
  SmallIntegerArray(int nbrEntries, int largestInteger);

  // constructor from given content as integers
  // nbrEntries = length of array
  // largestInteger = largest integer to hold in each field
  // allEntries = array with entries to be stored
  //
  SmallIntegerArray(int nbrEntries, int largestInteger, unsigned *allEntries);

  // copy constructor (full duplication)
  //
  SmallIntegerArray( const SmallIntegerArray &array);

  // augment an array by an additional element
  // array = initial part of array
  // toAppend = new entry
  //
  SmallIntegerArray( const SmallIntegerArray &array, unsigned toAppend);

  //destructor
  //
  ~SmallIntegerArray();

  

  // assignment operator
  //
  // array = array to assign
  // return value = reference on current array
  //
  SmallIntegerArray& operator = (const SmallIntegerArray& array);


  // access functions: reading
  //
  void GetElement(int i, unsigned &value) const;
  unsigned GetElement(int i) const;

  // access function: writing
  //
  void SetElement(int i, unsigned value);

  // fast access for all elements
  // values = array of unsigned integers, supposed to be sufficiently large for answer
  //
  void GetElements(unsigned *values);

  // fast access for all elements
  // values = array of unsigned integers, supposed to be sufficiently large for answer
  //
  void SetElements(unsigned *values);

  // request length of array
  //
  int GetNbrElements() { return NbrEntries;}

  // logical operations
  //
  friend bool operator == (const SmallIntegerArray& a1, const SmallIntegerArray& a2);
  friend bool operator != (const SmallIntegerArray& a1, const SmallIntegerArray& a2);
  friend bool operator < (const SmallIntegerArray& a1,const SmallIntegerArray& a2);
  friend bool operator > (const SmallIntegerArray& a1,const SmallIntegerArray& a2);
  friend bool operator <= (const SmallIntegerArray& a1,const SmallIntegerArray& a2);
  friend bool operator >= (const SmallIntegerArray& a1,const SmallIntegerArray& a2);

  // output stream overload
  //
  friend ostream& operator << (ostream & Str, const SmallIntegerArray& a);
  
  
};

#endif
