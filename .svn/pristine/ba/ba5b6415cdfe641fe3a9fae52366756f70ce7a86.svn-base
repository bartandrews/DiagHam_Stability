////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//                                                                            //
//                        last modification : 10/09/2003                      //
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


#ifndef FERMIONONTORUSWITHMAGNETICTRANSLATIONS_H
#define FERMIONONTORUSWITHMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"

using std::cout;
using std::endl;
using std::dec;
using std::hex;


class FermionOnTorusWithMagneticTranslations :  public ParticleOnTorusWithMagneticTranslations
{

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;

  // maximum momentum value reached by a fermion
  int MaxMomentum;
  // number of momentum values in a state (= MaxMomentum +1)
  int NbrMomentum;
  // GCD of MaxMomentum and NbrFermions (momemta are defined modulo MomentumModulo)
  int MomentumModulo;
  // momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  int XMomentum;
  // momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  int YMomentum;

  // value that has to be substracted to the momentum for each translation of the canonical form research
  int MomentumIncrement;
  // shift that has to be done on a state for each translation of the canonical form research
  int StateShift;
  // complementary shift (with respect to MaxMomentum) to StateShift
  int ComplementaryStateShift;
  // mask corresponding to StateShift
  unsigned long MomentumMask;

  // array describing each state 
  unsigned long* StateDescription;
  // array giving maximum momentum value reached for a fermion in a given state
  int* StateMaxMomentum;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given maxMomentum sector
  int LookUpTableMemorySize;
  // shift used in each maxMomentum sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used maxMomentum value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;
  // a table containing parity of the sum of 1 bits for all integer ranging from 0 to 2^MaximumSignLookUp - 1 (1 if odd)
  int* NbrParticleLookUpTable;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // sign due to state reordering when applying translation operator 
  unsigned long* ReorderingSign;
  // array of unsigned long where each bit describes sign associated to each translation of the orbit representant (0 for +, 1 for -) with respect to N-body ordering convention
//  int* StateSignature;

  // array containing for each state the sign due to fermion reordering when translating state (1 bit to 0 if sign is negative)
//  unsigned long* TranslationSign;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions 
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  FermionOnTorusWithMagneticTranslations (int nbrFermions, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithMagneticTranslations(const FermionOnTorusWithMagneticTranslations& fermions);

  // destructor
  //
  ~FermionOnTorusWithMagneticTranslations ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithMagneticTranslations& operator = (const FermionOnTorusWithMagneticTranslations& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value in the y direction of a given state
  //
  // index = state index
  // return value = state momentum in the y direction
  int GetYMomentumValue(int index);

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

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

  // apply a^+_m a_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation operator
  // return value =  resulting multiplicative factor 
  double AdA (int index, int m);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 private:

  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // yMomentum = state momentum value in the y direction
  // return value = canonical form of a state description
  unsigned long FindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint
  unsigned long FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

 // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  int FindNumberXTranslation(unsigned long stateDescription);

  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
  // return value = true if the state satisfy the x momentum constraint
  bool TestXMomentumConstraint(unsigned long stateDescription, int maxMomentum);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate look-up table associated to sign calculations
  // 
  void GenerateSignLookUpTable();

  // generate all states corresponding to the constraints
  // 
  // return value = hilbert space dimension
  int GenerateStates();

  // generate all states corresponding to the constraints (without taking into the canonical form) 
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentYMomentum = current value of the momentum in the y direction
  // return value = position from which new states have to be stored
  int RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentUMomentum);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnTorusWithMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithMagneticTranslations::FermionicStatistic;
}


#endif


