////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2006 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//          various unsigned long tools related to QHE Hilbert spaces         //
//                                                                            //
//                        last modification : 23/01/2006                      //
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


#include "config.h"
#include "UnsignedIntegerTools.h"
#include <iostream>
#include <cstdlib>


// return next bigger word with same number of bits set.
//
// i = word from which next bigger word has to be built
// return value = next bigger word

unsigned long nextone (unsigned long i)
{
  unsigned long bit=1ul;
  int count=-1;
  if (i==0ul) return 0ul;
  /* find first one bit */
  while (!(bit&i)) bit<<=1;
  /* find next zero bit */
  while (bit&i)
    {
      count++;
      bit<<=1;
    }
  if (!bit) {printf("overflow in nextone");return (~0ul);}
  
  i &= (~(bit-1));           /* clear lower bits */
  
  i |= bit | ((1u<<count)-1); /* put them in new places */
  
  return i;
}


// return next smaller word with same number of bits set.
//
// i = word from which next smaller word has to be built
// return value = next smaller word

unsigned long lastone (unsigned long i)
{
#ifdef  __64_BITS__
  int exitOn=64;
#else
  int exitOn=32;
#endif

  unsigned long  bit = 0x1ul;
  int count=0, spacing=-1;
  if (i==0x0ul) return 0ul;
  
  /* count set bits on the left */
  while (bit&i) { ++count; bit<<=1;}
  exitOn-=count;


  /* find next non-zero bit */
  if (!(i&(~(bit-1)))) return 0x0ul; // if non existant, exit
  while (!(bit&i))
    {
      bit<<=1;
      ++spacing;
    }
  
  i &= ~((bit<<1)-1);           /* clear lower bits */

  i |= (bit>>1) | ( ( (0x1l<<(count+spacing))-1) ^ ( (0x1ul<<(spacing))-1) ); /* put them in new places */
  
  return i;
}


// count bits in word
//
// i = word to test
// return value = number of bits

int bitcount (unsigned long i)
{
  int result=0;
  while(i) {
    result++;
    i &= i-1;
  }
  return result;
}

// get number of biggest bit set, numbering from 1 as lsb to 64/32, returns 0 if argument is 0

int leftmostBit(unsigned long i)
{
  if (i==0x0ul) return 0;
  else {
#ifdef  __64_BITS__
    int maxBit=64;
#else
    int maxBit=32;
#endif
    unsigned long bit=0x1ul << (maxBit-1);
    while (!(i&bit))
      {
	bit >>= 1;
	--maxBit;
      }
    return maxBit;
  }
}
	
