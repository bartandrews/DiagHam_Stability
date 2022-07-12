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


class SimplePermutations
{
 private:
  // total number of elements
  int NbrElements;

  // number of different permutations
  long NbrPermutations;
  
  // permutations
  unsigned **Permutations;

  // space used internally, only
  unsigned *MyPermutation;

  // current index to be filled
  long CurrentPermutation;

  // cleanup flag
  bool CleanUp;
  
 public:
  // default constructor
  // nbrElements = number of elements
  SimplePermutations(int nbrElements);

  // destructor
  ~SimplePermutations();

  // get Array of permutations
  inline unsigned* GetPermutation(int nbrPermutation);

  // get number of permutations
  long GetNbrPermutations(){return NbrPermutations;}

  // check out Permutations to an external method
  // if this method is called, the caller is responsible for cleaning up
  unsigned** CheckOutPermutations();

 private:

  // central recursion for generating permutations
  void Permute(const int start, const int n, int sign);

  // rotate permutation left
  void RotateLeft(const int start, const int n, int &sign);

  // swap two members
  void Swap(const int i, const int j, int &sign);

};

// get Array of permutations
unsigned* SimplePermutations::GetPermutation(int nbrPermutation)
{
  return this->Permutations[nbrPermutation];
}


#endif
