////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2008 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//     implements an array of bits coded on unsigned integers                 //
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


#ifndef BITMARKERS
#define BITMARKERS

#include "config.h"
#include <cstdlib>

class BitMarkers
{
 protected:

  // length of bitstring
  unsigned Length;

  // number of words
  unsigned NbrWords;

  // array of bit 'markers'
  unsigned *MarkerArray;
  
 public:
  // standard constructor - create empty array
  BitMarkers();
  
  // constructor
  BitMarkers(unsigned length);

  // copy constructor
  BitMarkers(const BitMarkers &markers);

  // assignment operator
  BitMarkers& operator = (const BitMarkers &markers);

  // destructor
  ~BitMarkers();

  // accessor methods
  
  // mark bit with given index
  // bitIndex = bit affected
  void MarkBit(unsigned bitIndex);


  // request if bit is marked / unmarked
  // bitIndex = bit affected
  bool IsMarked(unsigned bitIndex);
  bool IsNotMarked(unsigned bitIndex);

  // query if all bits are set
  bool HaveMarkedAll();

  // reset all bits to zero
  void ResetMarked();
};


// mark bit with given index
// bitIndex = bit affected
inline void BitMarkers::MarkBit(unsigned bitIndex)
{
  unsigned int * pt = this->MarkerArray+(bitIndex>>5);
  (*pt) |= (0x1u <<(bitIndex&0x1fu));
}

// request if bit is not marked
// bitIndex = bit affected
inline bool BitMarkers::IsNotMarked(unsigned bitIndex)
{
  unsigned int * pt = this->MarkerArray+(bitIndex>>5);
  if ( (*pt) & (0x1u <<(bitIndex&0x1fu)))
    return false;
  else return true;
}

// request if bit is marked
// bitIndex = bit affected
inline bool BitMarkers::IsMarked(unsigned bitIndex)
{
  unsigned int * pt = this->MarkerArray+(bitIndex>>5);
  if ( (*pt) & (0x1u <<(bitIndex&0x1fu)))
    return true;
  else return false;
}


// query if all bits are set
inline bool BitMarkers::HaveMarkedAll()
{
  unsigned Test32 = 0xffffffff;
  for (unsigned w=0; w<NbrWords-1;++w)
    if (this->MarkerArray[w] != Test32) return false;
  if ((Length&0x1fu)==0)
    {
      if (this->MarkerArray[NbrWords-1] != Test32) return false;
      else return true;
    }
  else
    {
      Test32 = (1u<<(Length&0x1fu))-1;
      if (this->MarkerArray[NbrWords-1] !=Test32) return false;
      else return true;
    }
}

// reset all bits to zero
inline void BitMarkers::ResetMarked()
{
  for (unsigned w=0; w<NbrWords;++w) 
    this->MarkerArray[w]=0u;
}

#endif
