////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                     class of bosons on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 08/09/2009                      //
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


#ifndef BOSONONSPHERETWOLANDAULEVELS_H
#define BOSONONSPHERETWOLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "GeneralTools/ArrayTools.h"

#include <iostream>

#include <map>

using std::map;

class FermionOnSphereTwoLandauLevels;

class BosonOnSphereTwoLandauLevels :  public ParticleOnSphereWithSpin
{
  
  friend class FermionOnSphereTwoLandauLevels;
	friend class FermionOnSphereWithSpin;
  
 protected:

  // number of fermions
  int NbrBosons;
  // number of fermions plus 1
  int IncNbrBosons;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a boson
  int LzMax;
  // number of bosons with spin up / down
  int NbrBosonsUp;
  int NbrBosonsDown;
  // number of Lz values in a stat
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;

  // array describing each state
  unsigned long* StateDescription;
  
  // array giving maximum Lz value reached for a boson in a given state
  int* StateLzMax;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state
  unsigned long* TemporaryState;
  // Lz maximum value associated to temporary state
  int TemporaryStateLzMax;
  
  // temporary state used when applying ProdA operator
  unsigned long* ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdATemporaryStateLzMax;

  // maximum Lz value reached by a boson with a spin up
  int LzMaxUp;
  // maximum Lz value reached by a boson with a spin down
  int LzMaxDown;
  // shift to apply on the spin up part
  int LzShiftUp;
  // shift to apply on the spin down part
  int LzShiftDown;
  // sum of LzShiftUp and LzShiftDown
  int LzTotalShift;
  // numer of bit occupied by the down part
  int UpStateShift;

  // holds instance of class to calculate the Clebsch Gordan coefficients
  ClebschGordanCoefficients*	CGCoefficients;
  
  // array to store the precalculated PseudoPotentials
  double *PseudoPotentials; 

 public:

  // default constructor
  //
  BosonOnSphereTwoLandauLevels();

  // basic constructor with contraint on the number of particles per Landau level
  // 
  // nbrBosonsUp = number of bosons in N=1 LL
  // nbrBosonsDown = number of bosons in N=0 LL
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a boson with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a boson with a spin down
  // memory = amount of memory granted for precalculations
  BosonOnSphereTwoLandauLevels (int nbrBosonsUp, int nbrBosonsDown, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);

  // basic constructor with no contraint on the number of particles per spin component
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a boson with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a boson with a spin down
  // memory = amount of memory granted for precalculations
  BosonOnSphereTwoLandauLevels (int nbrBosons, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereTwoLandauLevels(const BosonOnSphereTwoLandauLevels& bosons);

  // destructor
  //
  ~BosonOnSphereTwoLandauLevels ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereTwoLandauLevels& operator = (const BosonOnSphereTwoLandauLevels& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
   
  // project out any configurations that have particles on levels other than lll
  //
  // inputVector = vector to apply the projection to
  // outputVector = projected vector
  // finalSpace = reference to space of output vector
  void  ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, BosonOnSphereShort * finalSpace);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);
  
  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, BosonOnSphere& upStateSpace, RealVector& downState, BosonOnSphere& downStateSpace);
  
  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(2) index for the first annihilation operator
  // sigma2 = SU(2) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AdAd (int index, int n1, int n2);
  
  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AuAu (int index, int n1, int n2);
  
  // apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AuAd (int index, int n1, int n2);
  
  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(2) index for the first creation operator
  // sigma2 = SU(2) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAdd (int m1, int m2, double& coefficient);
  
  
  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, double& coefficient);
  
  // compute the number of particles in each Landau level
  //
  // state = ID of the state to handle
  // LLOccupationConfiguration = array where the decomposition will be store
  void  LandauLevelOccupationNumber(int state, int* lLOccupationConfiguration);
  
  // Calculate normalisation needed for this config
  // 
  // index = index of state.
  // return value = index used for bosonic representation
  double GetConfigNorm(long index);
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = binary representation of state to print
  // return value = reference on current output stream 
  ostream& PrintStateBinary (ostream& Str, unsigned long state);
	
  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem  
  ///virtual RealMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient = false);
  
  // compute the projection of the product of a bosonic state and the halperin 110 state
  //
  // bosonState = real vector where the bosonic state is stored
  // outputVector = real vector where the result has to be stored
  // fermionSpace = pointer to the fermionic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void BosonicStateTimePolarizedSlaters(RealVector& bosonState, RealVector& outputVector, FermionOnSphere * fermionSpace , FermionOnSphereWithSpin* finalSpace, int firstComponent,int nbrComponent, unsigned long** slaterPermutations, double *slaterSigns, int nbrSlaterPermutations);

  // compute the projection of the product of a bosonic state and the halperin 110 state
  //
  // bosonState = real vector where the bosonic state is stored
  // outputVector = real vector where the result has to be stored
  // fermionSpace = pointer to the fermionic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void BosonicStateTimePolarizedSlatersLzSymmetry(RealVector& bosonState, RealVector& outputVector, FermionOnSphere * fermionSpace , FermionOnSphereWithSpin* finalSpace, int firstComponent,int nbrComponent, unsigned long** slaterPermutations, double *slaterSigns, int nbrSlaterPermutations);

  
  // compute the projection of the product of a bosonic state and the halperin 110 state
  //
  // bosonState = real vector where the bosonic state is stored
  // outputVector = real vector where the result has to be stored
  // fermionSpace = pointer to the fermionic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void BosonicStateTimePolarizedSlatersLzSzSymmetry(RealVector& bosonState, RealVector& outputVector, FermionOnSphere * fermionSpace , FermionOnSphereWithSpin* finalSpace, int firstComponent,int nbrComponent, unsigned long** slaterPermutations, double *slaterSigns, int nbrSlaterPermutations);
	
  // compute the projection of the product of a bosonic state and the halperin 110 state
  //
  // bosonState = real vector where the bosonic state is stored
  // outputVector = 	real vector where the result has to be stored
  // fermionSpace = pointer to the fermionic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void BosonicStateTimePolarizedSlaters(RealVector& bosonState, RealVector& outputVector, FermionOnSphere * fermionSpaceUp ,FermionOnSphere * fermionSpaceDown,  FermionOnSphereWithSpin* finalSpace, int indexUp , int indexDown, int firstComponent,int nbrComponent);

  // compute  the lz symmetric state with given index
  //
  // index = index of configuration
  // return = lz symmstric configuration
  inline int GetSymmetricStateIndex (int index);
  
  // remove all zeros from the vector and remove corresponding state information
  //
  // initialState = reference to the vector in question  
  // return = new size of the space
  int RemoveZeros(RealVector& initialState, bool lzSym = false);
  
 protected:

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where the initial up bosonic state is stored
  // initialStateLzMax = reference on the initial down bosonic state maximum Lz value
  // return value = corresponding fermionic state
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  void FermionToBoson(unsigned long initialState, int initialStateLzMax,unsigned long*& finalState, int& finalStateLzMax);
  
  // convert a fermionic state to its monomial representation
  //
  // index = index of the fermionic state
  // finalState = reference on the array where the monomial representation has to be stored
  void GetMonomial(long index, unsigned long*& finalState);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial  bosonic state
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long* initialState, int initialStateLzMax, unsigned long*& finalState);
  
  // convert a fermionic state to its monomial Landau representation
  //
  // index = index of the fermionic state
  // finalState = reference on the array where the monomial representation has to be stored
  void GetMonomialLandau(long index, unsigned long*& finalState);
  
  // convert a bosonic state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  unsigned long ConvertFromMonomial(unsigned long* initialState);
  
  // works out the maximum possible totallz that is left
  //
  // nbrBosons = the number of bosons that are left
  // pos = the index of the position we are on in filling with bosons where 0 is the largest Lz on the second Landau level
  // return value = maximum possible totallz that can be result
  long MaxLzLeft(int nbrBosons, int pos);
  
  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription);
  
  // calculate the pseudo potentials denoted V^S_J in the literature
  //
  // S double the max angular momentum
  // J the relative angular momentum
  // return value = The pseudo potential V^S_J
  double CalculatePseudoPotential(int S, int J);
  
  // calculate the number of ways of choosing c elements from n
  // 
  // n the number of options
  // c the number of choices
  // return value = number of ways of choosing c elements from n
  unsigned long NChooseC(int n , int c );
  
  // evaluate Hilbert space dimension without constraint on the number of particles per level
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMaxUp = momentum maximum value for a boson
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax, int totalLz, long pos);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  //virtual void GenerateLookUpTable(unsigned long memory);
  
  // Get the Lz value in second LL from the index given. index starts off at 0 on largest Lz of second LL (Up)
  // 
  // index = index used in boson representation
  // return value = equivalent Lz value with appropriate shift
  int GetLzFromIndexU(int index);
  
  // Get the Lz value in LLL from the index given. index starts off at 0 on largest Lz of second LL (Up)
  // 
  // index = index used in boson representation
  // return value = equivalent Lz value with appropriate shift
  int GetLzFromIndexD(int index);
  
  // Get the index value from the Lz for the second LL (labelled with Up)
  // 
  // lz = shifted Lz value
  // return value = index used for bosonic representation
  int GetIndexFromLzU(int lz);
  
  // Get the index value from the Lz for the LLL (labelled with Down)
  // 
  // lz = shifted Lz value
  // return value = index used for bosonic representation
  int GetIndexFromLzD(int lz);
	
  // generate the different states that appear in the product of a slater in the lowest Landau level and a Slater determinant in the two Landau levels
  //
  // sortingMap = map in which the generated states and their coefficient will be stored
  // slater = array where the Slater determinant is stored in its monomial representation
  // state = array where the obtained state is stored in its monomial representation
  // slaterSpace = pointer to the Hilbert Space which the Slater determinant belongs to
  // index = index of the particle being examinate
  // coef = coefficient of the state being generate
  virtual void GeneratesDifferentState(map <unsigned long, double> & sortingMap, unsigned long * slater, unsigned long * state, FermionOnSphereTwoLandauLevels * slaterSpace, int index, double coef);
  
};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereTwoLandauLevels::GetParticleStatistic()
{
  return ParticleOnSphere::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialStateUp = reference on the array where the initial up bosonic state is stored
// initialStateLzMaxUp = reference on the initial up bosonic state maximum Lz value
// initialStateDown = reference on the array where the initial down bosonic state is stored
// initialStateLzMaxDown = reference on the initial down bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereTwoLandauLevels::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
{
  unsigned long TmpState = 0x0ul;
  unsigned long bosons, Mask;
  int bosons_placed = 0;

  for (int i = 0; i <= initialStateLzMax; ++i)
    {
      if ( initialState[i] > 0 ) 
	{
	  bosons = (0x1ul << initialState[i]) - 1;
#ifdef  __64_BITS__  	  
	  Mask = (bosons << (63 - (initialState[i] - 1 ) - (i + bosons_placed) )); //this places bosons starting from the MSB.
#else	  
	  Mask = (bosons << (31 - (initialState[i] - 1 ) - (i + bosons_placed) )); //this places bosons starting from the MSB.
#endif	  
	  bosons_placed += initialState[i];
	  TmpState |= Mask;
	}
    }	 
  return TmpState;
}

// convert a fermionic state into its bosonic  counterpart
// 
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnSphereTwoLandauLevels::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
{
  int pos = 0;
  finalState[pos] = 0;
  unsigned long mask;
  bool last = false;
#ifdef  __64_BITS__   
  int bit = 63;
  while ( bit >= 0 && pos < this->NbrLzValue && bit >= (63 - initialStateLzMax) )
#else    
  int bit = 31;
  while ( bit >= 0 && pos < this->NbrLzValue && bit >= (31 - initialStateLzMax) )     
#endif    
    {
      mask = 0x1ul << bit; 
      if ( (mask & initialState) > 0 ) 
        {
	  finalState[pos]++;
	  last = true;
	}	
      else 
        {
	  if ( last ) last = false;
	  pos++;
	  finalState[pos] = 0;
	}
       bit--;
    }
  if ( pos < this->NbrLzValue ) 
    {
      finalStateLzMax = pos;
    }
  else
    {
      finalStateLzMax = this->NbrLzValue-1;
    } 
}

// convert a fermionic state to its monomial representation
//
// index = index of the fermionic state
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereTwoLandauLevels::GetMonomial(long index, unsigned long*& finalState)
{
  int Index = this->NbrBosons-1; //start filling from last element as will have lowest value.
  unsigned long InitialState = this->StateDescription[index];
  int Pos  = 0;
  int StateLzMax = this->StateLzMax[index]; //number of bits away from most significant bits.
  int TmpLz = 0; // using labels as in bosonic representation where 0 corresponds to largest lz on SLL. 
  while (Pos <= StateLzMax)
    {
#ifdef  __64_BITS__      
      while ((Pos <= StateLzMax) && (((InitialState >> (63 - Pos)) & 0x1ul) != 0x0ul))
#else	
      while ((Pos <= StateLzMax) && (((InitialState >> (31 - Pos)) & 0x1ul) != 0x0ul))	
#endif	
	{
	  finalState[Index--] = (unsigned long) TmpLz;
	  ++Pos;
	}
#ifdef  __64_BITS__	
      while ((Pos <= StateLzMax) && (((InitialState >> (63 - Pos)) & 0x1ul) == 0x0ul))
#else	
      while ((Pos <= StateLzMax) && (((InitialState >> (31 - Pos)) & 0x1ul) == 0x0ul))	
#endif	
	{
	  ++TmpLz;
	  ++Pos;
	}
    }
}

// convert a fermionic state to its monomial representation
//
// index = index of the fermionic state
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereTwoLandauLevels::GetMonomialLandau(long index, unsigned long*& finalState)
{
  int Index = this->NbrBosons-1; //start filling from last element as will have lowest value.
  unsigned long InitialState = this->StateDescription[index];
  int Pos  = 0;
  int StateLzMax = this->StateLzMax[index]; //number of bits away from most significant bits.
  int TmpLz = 0; // using labels as in bosonic representation where 0 corresponds to largest lz on SLL. 
  while (Pos <= StateLzMax)
    {
#ifdef  __64_BITS__      
      while ((Pos <= StateLzMax) && (((InitialState >> (63 - Pos)) & 0x1ul) != 0x0ul))
#else	
	while ((Pos <= StateLzMax) && (((InitialState >> (31 - Pos)) & 0x1ul) != 0x0ul))	
#endif	
	  {
	    finalState[Index--] = (unsigned long) TmpLz;
	    ++Pos;
	  }
#ifdef  __64_BITS__	
      while ((Pos <= StateLzMax) && (((InitialState >> (63 - Pos)) & 0x1ul) == 0x0ul))
#else	
	while ((Pos <= StateLzMax) && (((InitialState >> (31 - Pos)) & 0x1ul) == 0x0ul))	
#endif	
	  {
	    ++TmpLz;
	    ++Pos;
	  }
    }
  
  for (int i = 0 ; i < this->NbrBosons ; i++)
		{
		  if (finalState[i] < this->LzMaxUp + 1)
		    finalState[i] = (GetLzFromIndexU(finalState[i]) <<1) + 1;
		  else
		    finalState[i] = GetLzFromIndexD(finalState[i])<<1;
		}
  SortArrayDownOrdering(finalState,this->NbrBosons);
}

// convert a bosonic state to its monomial representation
//
// initialState = initial  bosonic state
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereTwoLandauLevels::ConvertToMonomial(unsigned long* initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int Index = 0;
  for (int i = initialStateLzMax; i >= 0; --i)
    for (unsigned long j = 0l; j < initialState[i]; ++j)
      finalState[Index++] = i;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereTwoLandauLevels::ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int Index = this->NbrBosons-1; //start filling from last element as will have lowest value.
  int Pos  = 0;
  int TmpLz = 0; // using labels as in bosonic representation where 0 corresponds to largest lz on SLL. 
  while (Pos <= initialStateLzMax)
    {
#ifdef  __64_BITS__        
      while ((Pos <= initialStateLzMax) && (((initialState >> (63 - Pos)) & 0x1ul) != 0x0ul))
#else
      while ((Pos <= initialStateLzMax) && (((initialState >> (31 - Pos)) & 0x1ul) != 0x0ul))
#endif
	{
	  finalState[Index--] = (unsigned long) TmpLz;
	  ++Pos;
	}
#ifdef  __64_BITS__	
      while ((Pos <= initialStateLzMax) && (((initialState >> (63 - Pos)) & 0x1ul) == 0x0ul))
#else	
      while ((Pos <= initialStateLzMax) && (((initialState >> (31 - Pos)) & 0x1ul) == 0x0ul))
#endif	
	{
	  ++TmpLz;
	  ++Pos;
	}
    }
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnSphereTwoLandauLevels::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long Tmp = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
#ifdef  __64_BITS__     
    Tmp |= 0x1ul << (63 - (initialState[this->NbrBosons-i-1]+i)); 
#else
    Tmp |= 0x1ul << (31 - (initialState[this->NbrBosons-i-1]+i)); 
#endif
  return Tmp;
}



// Get the Lz value in second LL from the index given. index starts off at 0 on largest Lz of second LL (Up)
// 
// index = index used in boson representation
// return value = equivalent Lz value with appropriate shift

inline int BosonOnSphereTwoLandauLevels::GetLzFromIndexU(int index)
{
  // for the second (Up) LL the indices go from 0 for highest Lz to this->LzMaxUp. does not check index is within expected range
  return this->LzMaxUp - index; 
}

// Get the Lz value in LLL from the index given. index starts off at 0 on largest Lz of second LL (Up)
// 
// index = index used in boson representation
// return value = equivalent Lz value with appropriate shift

inline int BosonOnSphereTwoLandauLevels::GetLzFromIndexD(int index)
{
  // for the LLL (Down) the indices go from this->LzMaxUp+1 for highest Lz to this->LzMaxUp+1 + this->LzMaxDown. does not check index is within expected range
  //return this->LzMaxDown + 1 - (index - (this->LzMaxUp+1)); 
  //return this->LzMaxDown + this->LzMaxUp + 2 - index;
  return (this->LzMaxUp << 1 ) - index;
}

// Get the index value from the Lz for the second LL (labelled with Up)
// 
// lz = shifted Lz value
// return value = index used for bosonic representation

inline int BosonOnSphereTwoLandauLevels::GetIndexFromLzU(int lz)
{
  return this->LzMaxUp - lz;
}
  
// Get the index value from the Lz for the LLL (labelled with Down)
// 
// lz = shifted Lz value
// return value = index used for bosonic representation

inline int BosonOnSphereTwoLandauLevels::GetIndexFromLzD(int lz)
{
  //return this->LzMaxUp + 1 + this->LzMaxDown + 1 - lz;
  return (this->LzMaxUp << 1 ) - lz;
}

// Calculate normalisation needed for this config
// 
// index = index of state.
// return value = index used for bosonic representation

inline double BosonOnSphereTwoLandauLevels::GetConfigNorm(long index)
{
  double Value = 1.0;
  this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  
  for ( int i = 0 ; i <= this->ProdATemporaryStateLzMax ; i++ ) 
    {
      if ( this->ProdATemporaryState[i] > 0 ) 
	{
	  Value *= (double)this->ProdATemporaryState[i];
	}
    }
  return sqrt(Value);
}

// compute  the lz symmetric state with given index
//
// index = index of configuration
// return = lz symmstric configuration

inline int BosonOnSphereTwoLandauLevels::GetSymmetricStateIndex (int index)
{
  int Index = this->NbrBosons-1; //start filling from last element as will have lowest value.
  unsigned long InitialState = this->StateDescription[index];
  int Pos  = 0;
  int StateLzMax = this->StateLzMax[index]; //number of bits away from most significant bits.
  int TmpLz = 0; // using labels as in bosonic representation where 0 corresponds to largest lz on SLL. 
  while (Pos <= StateLzMax)
    {
#ifdef  __64_BITS__      
      while ((Pos <= StateLzMax) && (((InitialState >> (63 - Pos)) & 0x1ul) != 0x0ul))
#else	
	while ((Pos <= StateLzMax) && (((InitialState >> (31 - Pos)) & 0x1ul) != 0x0ul))	
#endif	
	  {
	    this->TemporaryState[Index--] = (unsigned long) TmpLz;
	    ++Pos;
	  }
#ifdef  __64_BITS__	
      while ((Pos <= StateLzMax) && (((InitialState >> (63 - Pos)) & 0x1ul) == 0x0ul))
#else	
	while ((Pos <= StateLzMax) && (((InitialState >> (31 - Pos)) & 0x1ul) == 0x0ul))	
#endif	
	  {
	    ++TmpLz;
	    ++Pos;
	  }
    }
  
  for (int i = 0 ; i < this->NbrBosons ; i++)
    {      
      if (this->TemporaryState[i] < this->LzMaxUp + 1)
	this->TemporaryState[i] = GetIndexFromLzU(this->LzMax - GetLzFromIndexU(this->TemporaryState[i]));
      else
	this->TemporaryState[i] = GetIndexFromLzD(this->LzMax -GetLzFromIndexD(this->TemporaryState[i]));
    }
  
  SortArrayDownOrdering(this->TemporaryState,this->NbrBosons);
	
// 	for (int i = 0 ; i < this->NbrBosons ; i++)
// 		{
// 			cout <<this->TemporaryState[i]<<" ";
// 		}
// 		cout <<endl;
// 	this->TemporaryStateLzMax = this->NbrLzValue - 1;
//   while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
//     --this->TemporaryStateLzMax;
  return this->FindStateIndex(this->ConvertFromMonomial(this->TemporaryState));	
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereTwoLandauLevels::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  cout << "warning : AsigmaAsigma not defined in BosonOnSphereTwoLandauLevels" << endl;
  return 0.0;
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereTwoLandauLevels::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  cout << "warning : AdsigmaAdsigma not defined in BosonOnSphereTwoLandauLevels" << endl;
  return this->HilbertSpaceDimension;
}


#endif
