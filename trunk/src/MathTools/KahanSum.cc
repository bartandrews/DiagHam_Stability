////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                         class Author: Gunnar Möller                        //
//                                                                            //
//              class implementing Kahan summation (single sum)               //
//                                                                            //
//                        last modification : 04/02/2017                      //
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
#include "KahanSum.h"

#include <iostream>
using std::ostream;


// default constructor 
//
KahamSum<T>::KahamSum(): Overflow(0.0), Sum(0.0)
{}

// constructor from a double
//
// x = value to assign to the sum
KahamSum::KahamSum(double x): Overflow(0.0), Sum(x)
{
}

// assignment
//
// sum = sum to assign (including overflow)
// return value = reference on current sum
KahamSum& KahamSum::operator = (KahamSum& sum)
{
  this->Overflow=sum.Overflow;
  this->Sum=sum.Sum;
  return *this;
}

// assignment
//
// sum = sum to assign (setting overflow to zero)
// return value = reference on current sum
KahamSum& KahamSum::operator = (double x)
{
  this->Overflow=0.0;
  this->Sum=x;
  return *this;
}

// add a term, using Kahan summation
//
// x = term to add
// return value = reference on current sum
KahamSum& KahamSum::operator += (double x)
{
  this->Overflow-=x;  // overflow now contains (previous overflow - input) = -NewEntry
  x=this->Sum;  // copy overlaps to input array
  this->Sum-=this->Overflow;  // subtract overflow from final sum
  this->Overflow+=(this->Sum-x);  // adjust remaining overflow
  return *this;
}

// subtract a term, using Kahan summation
//
// x = term to subtract
// return value = reference on current sum
KahamSum& KahamSum::operator -= (double x)
{
  this->Overflow+=x;  // overflow now contains (previous overflow - input) = -NewEntry
  x=this->Sum;  // copy overlaps to input array
  this->Sum-=this->Overflow;  // subtract overflow from final sum
  this->Overflow+=(this->Sum-x);  // adjust remaining overflow
  return *this;
}

// multiply by a double
//
// x = integer to use
// return value = reference on current sum
KahamSum& KahamSum::operator *= (double x)
{
  this->Sum*=x; 
  this->Overflow*=x;
  return *this;
}

// output stream, printing internal representation of the sum
// str = stream to write to
// sum = Kahan sum object
// return = str
ostream& operator << (ostream& str, KahamSum &sum)
{
  str << sum.Sum << "(";
  if (sum.Overflow>0) str << "+";
  str << sum.Overflow << ")";
  return str;
}
