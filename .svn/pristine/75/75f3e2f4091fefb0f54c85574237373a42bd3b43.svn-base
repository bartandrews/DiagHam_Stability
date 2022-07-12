////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//                       class of abstract spin collection                    //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#ifndef GENERICSUNSPINCOLLECTION_H
#define GENERICSUNSPINCOLLECTION_H


#include "config.h"
#include "AbstractSUNSpinCollection.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>
using std::cout;
using std::endl;


// highest level of SU(N) possible in this implementation, MAXLEVEL needs to be set to 2^BITS
#define MAXLEVEL 8

// nbr of bits coding a single spins, accord with MAXLEVEL
#define BITS 3

// Mask for coding bits in lowest position
#define MASK 0x7ul



class GenericSUNSpinCollection : public AbstractSUNSpinCollection
{
 protected:
  
  // Level of SU(N) implementation
  int LevelN;

  // Number of spins
  int NbrSpins;

  // Quantum numbers associated with the operators of the Cartan algebra
  int *CartanQuantumNumbers;

  // array describing each state
  unsigned long* StateDescription;

  // number of spins indexed:
  int LookUpTableDepth;

  // Shift of description required to extract coding segment
  int LookUpTableShift;
  
  // array with entries of the LookUpTable
  int* LookUpTable;

  // garbage flag
  GarbageFlag Flag;

 public:

  // create Hilbert-space for a collection of spins
  // levelN = level of SU(N) symmetry
  // nbrSpins = number of spins in collection
  // cartanQuantumNumbers = quantum numbers of the first N-1 operators of the Cartan-algebra
  // memory = amount of memory granted for precalculations
  GenericSUNSpinCollection (int levelN, int nbrSpins, int *cartanQuantumNumbers, unsigned long memory);

  // copy constructor
  GenericSUNSpinCollection (GenericSUNSpinCollection &collection);

  // destructor
  //
  virtual ~GenericSUNSpinCollection ();

  // assignment operator
  GenericSUNSpinCollection& operator = (GenericSUNSpinCollection &collection);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  
  AbstractHilbertSpace* Clone();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // permutation operator of two spins
  // index = index of state to perform on
  // s1 = index of spin 1
  // s2 = index of spin 2
  // return = index of final state
  int SpinPermutation(int index, int s1, int s2);

  // cyclic permutation of three spins
  // index = index of state to perform on
  // s1 = index of spin 1
  // s2 = index of spin 2
  // s3 = index of spin 3
  // return = index of final state
  int CyclicSpinPermutation(int index, int s1, int s2, int s3);

  // cyclic permutation of k spins
  // index = index of state to perform on
  // numS = number of spins to permute
  // si = indices of spins
  // return = index of final state
  int CyclicSpinPermutation(int index, int numS, int *si);
  
  // get diagonal terms for an S*S interaction (counting instances for connections with same prefactor)
  // index = number of state to be considered
  int S2DiagonalElements(int index);

  // get off-diagonal terms for an S*S interaction (counting instances for connections with same prefactor)
  // index = number of state to be considered
  // spin = number of spin associated with first spin operator
  // targetIndices = states connected to (size: N-1)
  //
  void S2OffDiagonalElements(int index, int spin, int *targetIndices);  
  
  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state);

  // evaluate the Hilbert space dDimension
  // nbrSpins = number of spins
  // levelN = SU(N) level
  // cartanQuantumNumbers = full set of eigenvalues of the Cartan operators
  //
  long EvaluateHilbertSpaceDimension(int nbrSpins, int levelN, int *cartanQuantumNumbers);

  // generate all states corresponding to the constraints
  //
  // nbrSpin = number of spin to be considered next
  // levelN = highest remaining SU(N)-level
  // pos = position in StateDescription array where to store states
  // t_i = number of spins remaining with <T_i> = t_i
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrSpin, int levelN, long pos,
		      int t0, int t1, int t2, int t3, int t4, int t5, int t6, int t7);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
  // print but do not check any issues
  ostream& PrintStateOnly (ostream& Str, int state);
    
 private:

  // generate key from coding bits
  // codingBits = indexed part of a StateDescription
  // 
  inline unsigned GenerateKey(unsigned long codingBits);


};

// generate key from coding bits
// codingBits = indexed part of a StateDescription
// 
unsigned GenericSUNSpinCollection::GenerateKey(unsigned long codingBits)
{
  unsigned Key=(codingBits>>((this->LookUpTableDepth-1)*BITS));
  for (int k=this->LookUpTableDepth-2; k>0; --k)
    {
      Key*=this->LevelN;
      Key+=(codingBits>>(k*BITS))&MASK;      
    }
  return Key;
}

#endif


