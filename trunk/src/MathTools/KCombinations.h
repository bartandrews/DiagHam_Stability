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


#ifndef KCOMBINATIONS_H
#define KCOMBINATIONS_H

#include "config.h"

#include <iostream>
using std::ostream;

class KCombinations
{
 private:
  // number of elements  
  int NbrElements;

  // number of choices
  int NbrChoices;

  // number of different permutations
  int NbrCombinations;
  
  // permutations
  int **Combinations;

 public:
  // default constructor
  // nbrChoices = number of elements to choose
  // nbrElements = total number of elements
  // orderedGroups = flag indicating whether order of groups matters
  KCombinations(int nbrChoices, int nbrElements);

  // destructor
  ~KCombinations();

  // get Array of permutations
  int** GetCombinations() {return this->Combinations;}

  // get a particukar combination
  int* GetCombination(int i) {return this->Combinations[i];}

  // get length of arrays
  int GetNbrCombinations() {return this->NbrCombinations; }

  // output stream overload
  friend ostream& operator << (ostream & Str, KCombinations& a);

 private:

  // central recursive function that generates all different permutations
  // currentPos = current position to be chosen (start from highest, stop after placing zero)
  // minVal = minimum value for current position
  // pos = place where to start storing things
  // return = total number of quantum numbers
  int GenerateCombinations(int currentPos, int minVal, int pos);

};

#endif
