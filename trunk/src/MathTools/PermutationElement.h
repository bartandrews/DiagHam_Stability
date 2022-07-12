////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2008 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//  allows the calculation of permutations of groups and their multiplicity   //
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


#ifndef SIMPLEPERMUTATIONS_H
#define SIMPLEPERMUTATIONS_H

#include "GeneralTools/SmallIntegerArray.h"


class PermutationElement
{
public:

  SmallIntegerArray Value;
  unsigned long Multiplicity;
  
  PermutationElement();

  PermutationElement(SmallIntegerArray value, unsigned long multiplicity);
  
  PermutationElement(const PermutationElement &per);

  ~PermutationElement();

  // assignment operator
  //
  // per = permutation to assign
  // return value = reference on current permutation
  PermutationElement& operator = (const PermutationElement& array);

  // multiply on multiplicity
  PermutationElement& operator *= (const unsigned long mult);

  friend bool operator == (const PermutationElement& a1, const PermutationElement& a2);
  friend bool operator < (const PermutationElement& a1,const PermutationElement& a2);
  friend bool operator > (const PermutationElement& a1,const PermutationElement& a2);

  friend ostream& operator << ( ostream &Str, PermutationElement PE);
  
};

#endif
