////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin                    //
//                  including the Haldane squeezing technique                 //
//                  and the  Lz<->-Lz and Sz<->-Sz symmetries                 //
//                                                                            //
//                        last modification : 26/06/2009                      //
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


#ifndef FERMIONONSPHEREWITHSPINHALDANELZSZSYMMETRY_H
#define FERMIONONSPHEREWITHSPINHALDANELZSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>
#include <map>
#include <string>


class FermionOnSphereWithSpinHaldaneLzSzSymmetry :  public FermionOnSphereWithSpinLzSzSymmetry
{

 protected:

  // array of root partitions describing the squeezed basis
  unsigned long* RootPartitions;
  // number of root partitions
  int NbrRootPartitions;

  // three temporary arrays used during Hilbert space generation
  unsigned long* TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;
  
  // array storing the number of permutation used to texture a state
  int* NbrPermutations;
  // array storing the permutations used to texture a state
  unsigned long** Permutations;

 public:

  // default constructor
  //
  FermionOnSphereWithSpinHaldaneLzSzSymmetry();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = twice the maximum Lz value reached by a fermion
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  // rootPartitions = array of root partitions describing the squeezed basis
  // nbrRootPartitions = number of root partitions
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinHaldaneLzSzSymmetry (int nbrFermions, int lzMax,
					      bool minusSzParity, bool minusLzParity,
					      int** rootPartitions, int nbrRootPartitions, 
					      unsigned long memory = 10000000);
					      
  // textureless constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = twice the maximum Lz value reached by a fermion
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  // texturelessRootPartition = root partition describing the squeezed basis, spin texture has to be added on top of it   
  // nbrRootPartitions = number of root partitions
  // texturelessFlag = flag to indicate textureless squeezing.
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinHaldaneLzSzSymmetry (int nbrFermions, int lzMax,
					      bool minusSzParity, bool minusLzParity,
					      int** texturelessRootPartition, int nbrRootPartitions, bool texturelessFlag,					      
					      unsigned long memory = 10000000);					      

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinHaldaneLzSzSymmetry(const FermionOnSphereWithSpinHaldaneLzSzSymmetry& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinHaldaneLzSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinHaldaneLzSzSymmetry& operator = (const FermionOnSphereWithSpinHaldaneLzSzSymmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a given state from symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector  
  virtual RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis);

  // convert a given state from the usual n-body basis to the symmetric basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  virtual RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);


  // generate all squeezed states from a root partition
  // 
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateSqueezedStates(int lzMax, unsigned long referenceState, long pos, long& memory);
  
   // generate all squeezed states from a textureless root partition
  // 
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateSqueezedTexturelessStates(int lzMax, unsigned long referenceState, long pos, BosonOnSphereShort *bosonSpace, long& memory);
  
  // compute the projection of the product of a monomial in the two lowest LL and the halperin 110 state
  //
  // slater = array where the monomial representation of the slater determinant for half the number of particles is stored
  // monomial = array where the monomial representation is stored
  // sortingMap = map in which the generated states and their coefficient will be stored
  // nbrPermutations = number of different permutations
  // permutations1 = array where are stored the permutations of the spin up
  // permutations2 = array where are stored the permutations of the spin down
  // initialCoef = inital coefficient in front of the monomial
  virtual void MonomialsTimesPolarizedSlaterProjection(unsigned long * slater, unsigned long * monomial, map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2, double initialCoef);
  
  // evaluate all permutations requested to sapply spin texture
  //
  void EvaluatePermutations();
  
  // convert a given coefficient of state from the usual n-body basis to the symmetric basis
  //
  // coefficient = coefficient of configuration in usual n-body basis
  // stateDescription = configuration in usual n-body basis
  // return value = one if coefficient is to be used afterwards and zero otherwise
  virtual  int ConvertToSymmetricNbodyBasis(double& coefficient, unsigned long& stateDescription);

};



#endif


