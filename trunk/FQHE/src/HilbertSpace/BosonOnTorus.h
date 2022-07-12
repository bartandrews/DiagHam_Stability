////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of bosons on torus                         //
//                                                                            //
//                        last modification : 03/09/2002                      //
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


#ifndef BOSONONTORUS_H
#define BOSONONTORUS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


class BosonOnTorus :  public ParticleOnTorus
{

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // maximum momentum value reached by a boson
  int MaxMomentum;
  // number of Lz values in a state
  int NbrLzValue;

  // index of the momentum orbit
  int MomentumConstraint;
  // index of the momentum orbit
  bool MomentumConstraintFlag;
  int GCDMaxMomentum;

  // array describing each state
  int** StateDescription;
  // array giving maximum Lz value reached for a boson in a given state
  int* StateMaxMomentum;

  // multiplicative factors used during key construction
  int* KeyMultiplicationTable;
  // keys associated to each state
  int* Keys;
  // indicate position of the first state with a given number of boson having a given maximum Lz value
  int* MomentumMaxPosition;

  // indicates how many different states are store for each sector (a sector is given by its lzmax and the number of bosons that are at lzmax)
  int* KeyInvertSectorSize;
  // 
  int** KeyInvertTable;
  //
  int** KeyInvertTableNbrIndices;
  //
  int*** KeyInvertIndices;

  // temporary state used when applying operators
  int* TemporaryState;

  // pointer to the target space when an index is require after applying basic operation
  BosonOnTorus* TargetSpace;


 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  BosonOnTorus (int nbrBosons, int maxMomentum);

  // constructor with a constraint of the total momentum of states
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // momentumConstraint = index of the momentum orbit
  BosonOnTorus (int nbrBosons, int maxMomentum, int momentumConstraint);

  // constructor from full datas (with no constraint on the total momentum)
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // hilbertSpaceDimension = Hilbert space dimension
  // stateDescription = array describing each state
  // stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state
  BosonOnTorus (int nbrBosons, int maxMomentum, int hilbertSpaceDimension, 
		int** stateDescription, int* stateMaxMomentum);

  // constructor from full datas
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // momentumConstraint = index of the momentum orbit
  // hilbertSpaceDimension = Hilbert space dimension
  // stateDescription = array describing each state
  // stateMaxMomentum = array giving maximum Lz value reached for a boson in a given state
  BosonOnTorus (int nbrBosons, int maxMomentum, int momentumConstraint, int hilbertSpaceDimension, 
		int** stateDescription, int* stateMaxMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorus(const BosonOnTorus& bosons);

  // destructor
  //
  ~BosonOnTorus ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorus& operator = (const BosonOnTorus& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnTorus* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // return number of particles of the target space
  //
  // return value = number of particles of the target space
  virtual int GetTargetNbrParticles();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value of a given state
  //
  // index = state index
  // return value = state momentum
  int GetMomentumValue(int index);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

  // return matrix representation of the annihilation operator a_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& A (int i, Matrix& M);

  // return matrix representation ofthw creation operator a^+_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& Ad (int i, Matrix& M);

  // apply an annihilation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // coefficient = will be multiplied by the prefactor of the bosonic ladder operator
  // return value = index in the target space
  virtual int A (int i, int m, double &coefficient);

  // apply a creation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // coefficient = will be multiplied by the prefactor of the bosonic ladder operator
  // return value = index in the target space
  virtual int Ad (int i, int m, double &coefficient);

  // get a pointer to the target space
  // return value = target space
  //  ParticleOnTorus* GetTargetSpace()
  //  { return (ParticleOnTorus*) this->TargetSpace;}

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState);

 private:

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  int FindStateIndex(int* stateDescription, int lzmax);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate look-up table associated to current Hilbert space
  // 
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = key associated to the state
  int GenerateKey(int* stateDescription, int lzmax);
    
  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentMaxMomentum = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentMaxMomentum = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos, int currentMomentum);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorus::GetParticleStatistic()
{
  return ParticleOnTorus::BosonicStatistic;
}

#endif


