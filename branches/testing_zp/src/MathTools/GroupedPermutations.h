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


#ifndef GROUPEDPERMUTATIONS_H
#define GROUPEDPERMUTATIONS_H

#include "GeneralTools/SmallIntegerArray.h"
#include "GeneralTools/OrderedList.h"
#include "GeneralTools/UnsignedIntegerTools.h"


class PermutationElement;

class GroupedPermutations
{
 private:
  // number of groups  
  int NbrGroups;

  // number of elements per groups  
  int ElementsPerGroup;

  // total number of elements
  int NbrElements;

  // number of bits required to represent NbrElements
  int NbrBitsForElements;

  // number of bits required to represent groups
  int NbrBitsPerGroup;

  // flag indicating whether order of groups matters
  bool OrderedGroups;

  // number of different permutations
  int NbrPermutations;
  
  // permutations
  SmallIntegerArray *Permutations;

  // multiplicities
  unsigned long *Multiplicities;

  // temporary space for interface with SmallIntegerArray
  unsigned *MyArray;

  // map of groups (temporary, unless having ordered group)
  unsigned *MapOfGroups;
  unsigned *InverseMapOfGroups;

  // count of group members (temporary)
  unsigned *CountOfGroups;

  // internal storage whilst growing list of elements:
  OrderedList<PermutationElement> PermutationList;
  

 public:
  // default constructor
  // nbrGroups = number of groups
  // elementsPerGroup = number of elements per group
  // orderedGroups = flag indicating whether order of groups matters
  GroupedPermutations(int nbrGroups, unsigned elementsPerGroup, bool orderedGroups=false);

  // destructor
  ~GroupedPermutations();

  // get Array of permutations
  SmallIntegerArray* GetPermutations() {return this->Permutations;}

  // get Array of multiplicities
  unsigned long* GetMultiplicities(){return this->Multiplicities;}

  // get length of arrays
  unsigned GetNbrPermutations() {return this->NbrPermutations; }

 private:

  // central recursive function that generates all different permutations
  void CentralRecursion(SmallIntegerArray remainingElements, SmallIntegerArray permutation, unsigned long multiplicity);

  // translate internal form of permutations to a canonic expression
  SmallIntegerArray GetPermutationString(SmallIntegerArray &permutation);

  // get an initial string without permutations
  SmallIntegerArray GetInitialString();

  
  
};

#endif
