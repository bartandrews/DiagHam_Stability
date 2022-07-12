////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 19/05/2009                      //
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


#ifndef FERMIONONSPHERETWOLANDAULEVELS_H
#define FERMIONONSPHERETWOLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"

#include <iostream>
#include <map>

using std::map;

class BosonOnSphereTwoLandauLevels;

class FermionOnSphereTwoLandauLevels :  public FermionOnSphereWithSpin
{

  friend class BosonOnSphereTwoLandauLevels;
 
 protected:
  
  // maximum Lz value reached by a fermion with a spin up
  int LzMaxUp;
  // maximum Lz value reached by a fermion with a spin down
  int LzMaxDown;
  // shift to apply on the spin up part
  int LzShiftUp;
  // shift to apply on the spin down part
  int LzShiftDown;
  // sum of LzShiftUp and LzShiftDown
  int LzTotalShift;
  
 public:
  
  // default constructor
  //
  FermionOnSphereTwoLandauLevels();
  
  // basic constructor with contraint on the number of particles per Landau level component
  // 
  // nbrFermionsUp = number of fermions in level N=1
  // nbrFermionsDown = number of fermions in level N=0
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
  // memory = amount of memory granted for precalculations
  FermionOnSphereTwoLandauLevels (int nbrFermionsUp, int nbrFermionsDown, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);

  // basic constructor with no contraint on the number of particles per spin component
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
  // memory = amount of memory granted for precalculations
  FermionOnSphereTwoLandauLevels (int nbrFermions, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);
  
  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereTwoLandauLevels(const FermionOnSphereTwoLandauLevels& fermions);
  
  // destructor
  //
  ~FermionOnSphereTwoLandauLevels ();
  
  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereTwoLandauLevels& operator = (const FermionOnSphereTwoLandauLevels& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

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

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient);

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace);

  // compute the projection of the product of a Slater determinant and a monomial 
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable,int nbrVariable, map <unsigned long, double> & sortingMap, FermionOnSphere* finalSpace);
  
  // compute the projection of the product of a Slater determinant and a monomial 
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlaterProjection(unsigned long* slater,unsigned long* monomial,unsigned long* variable,int nbrVariable,  map <unsigned long, LongRational> & sortingMap, FermionOnSphere* finalSpace);
  
  // compute the projection of the product of a bosonic state and a fermionic state
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed	
  // reverseFluxFlag = true if it a reverse flux attachment
  virtual void BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, int firstComponent, int nbrComponent,  bool reverseFluxFlag = false);
  
  // compute the projection of the product of a bosonic state and a fermionic state
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed	
  virtual void BosonicStateTimeFermionicState(LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, int firstComponent,int nbrComponent);
  
  // compute the projection of the product of a bosonic state and a fermionic state using the lz->-lz symmetry
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed  
  // reverseFluxFlag = true if it a reverse flux attachment
  virtual void BosonicStateTimeFermionicStateSymmetric(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector,BosonOnSphereShort* bosonSpace, FermionOnSphere* finalSpace, int firstComponent, int nbrComponent, bool reverseFluxFlag = false);
  
  // compute the projection of the product of a bosonic state and a fermionic state using the lz->-lz symmetry
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed  
  virtual void BosonicStateTimeFermionicStateSymmetric(LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, BosonOnSphereShort* bosonSpace, FermionOnSphere* finalSpace, int firstComponent, int nbrComponent);
  
  // compute the product of a bosonic state and a fermionic state belonging in two Landau levels
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed	
  virtual void BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, BosonOnSphereShort* bosonSpace, FermionOnSphereTwoLandauLevels* finalSpace, int firstComponent, int nbrComponent);

  // compute the product of a bosonic state and a fermionic state belonging in two Landau levels
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed	
  virtual void BosonicStateTimeFermionicState(LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, BosonOnSphereShort* bosonSpace, FermionOnSphereTwoLandauLevels* finalSpace, int firstComponent, int nbrComponent);
  
  // compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  // reverseFluxFlag = true if it a reverse flux attachment
  virtual void LLLFermionicStateTimeFermionicState(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent, bool reverseFluxFlag = false);
	
  // compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  // reverseFluxFlag = true if it a reverse flux attachment
  virtual void LLLFermionicStateTimeFermionicState(LongRationalVector& lllFermionState, LongRationalVector& fermionState, LongRationalVector& outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent, bool reverseFluxFlag = false);
  
  // compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels using lz->-lz symmetry
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  // reverseFluxFlag = true if it a reverse flux attachment
  virtual void LLLFermionicStateTimeFermionicStateSymmetric(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent, bool reverseFluxFlag = false);
  
  // compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels using lz->-lz symmetry
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  // reverseFluxFlag = true if it a reverse flux attachment
  virtual void LLLFermionicStateTimeFermionicStateSymmetric(LongRationalVector& lllFermionState, LongRationalVector& fermionState, LongRationalVector& outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent, bool reverseFluxFlag = false);
	
  // compute the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void LLLFermionicStateTimeFermionicState(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, 
						   FermionOnSphere* lllFermionSpace, BosonOnSphereTwoLandauLevels * finalSpace, int firstComponent, int nbrComponent);
	  
  // compute the number of particles in each Landau level
  //
  // state = ID of the state to handle
  // lLOccupationConfiguration = array where the decomposition will be store
  virtual void LandauLevelOccupationNumber(int state, int* lLOccupationConfiguration);
  
  // project out any configurations that have particles on levels other than lll
  //
  // inputVector = vector to apply the projection to
  // outputVector = projected vector
  // finalSpace = reference to space of output vector space
  void  ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, FermionOnSphere * finalSpace);

  // compute the projection of the product of a fermionic state in the lowest Landau level and a fermionic state in the two lowest Landau levels
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void FermionicStateTimeFermionicState(RealVector& fermionState1, RealVector& fermionState2, RealVector& outputVector, FermionOnSphereTwoLandauLevels*  fermionSpace2 , 
						BosonOnSphereShort* finalSpace, int firstComponent,int nbrComponent);


protected:
  
  // evaluate Hilbert space dimension
  //
  // nbrFermionsUp = number of fermions with spin up
  // nbrFermionsDown = number of fermions with spin down
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz);
  
  // evaluate Hilbert space dimension without constraint on the number of particles per level
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermionsUp = number of fermions with spin up
  // nbrFermionsDown = number of fermions with spin down
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz, long pos);
  
  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMaxUp = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateFullStates(int nbrFermions, int lzMax, int totalLz, long pos);

  // convert a fermionic state to its monomial representation
  //
  // state = initial fermionic state in its fermionic representation
  // monomial = reference on the array where the monomial representation has to be stored  
  virtual void ConvertToMonomial(unsigned long state , unsigned long*& monomial);
  
  // convert a fermionic state from its monomial representation
  //
  // monomial = the array where the monomial representation is stored
  // return value = fermionic state in its fermionic representation  
  virtual  unsigned long ConvertFromMonomial(unsigned long* monomial);
	
  // convert a fermionic state to its monomial representation
  //
  // index = place in StateDescription of the fermionic state
  // monomial = reference on the array where the monomial representation has to be stored  
  virtual void ConvertToMonomial(int index, unsigned long*& monomial);
		
  // convert a fermionic state to its monomial variable representation
  //
  // state = initial fermionic state in its fermionic representation
  // slater = reference on the array where the monomial representation has to be stored
  // nbrVariable = number of fermions in the second Landau level
  // variable = reference on the array where the indice of fermions in the second Landau level has to be stored  
  virtual void ConvertToMonomialVariable(unsigned long state, unsigned long*& slater,int& nbrVariable, unsigned long *& variable);
	
  // generate the different states that appear in the product of a monomial and a Slater determinant in the two Landau levels
  //
  // sortingMap = map in which the generated states and their coefficient will be stored
  // slater = array where the Slater determinant is stored in its monomial representation
  // state = array where the obtained state is stored in its monomial representation
  // slaterSpace = pointer to the Hilbert Space which the Slater determinant belongs to
  // index = index of the particle being examinate
  // coef = coefficient of the state being generate	
  virtual void GeneratesDifferentState( map <unsigned long, double> & sortingMap ,unsigned long* slater,unsigned long* state,FermionOnSphereTwoLandauLevels * slaterSpace, int index, double coef);
							 
  // generate the different states that appear in the product of a monomial and a Slater determinant in the two Landau levels
  //
  // sortingMap = map in which the generated states and their coefficient will be stored
  // slater = array where the Slater determinant is stored in its monomial representation
  // state = array where the obtained state is stored in its monomial representation
  // slaterSpace = pointer to the Hilbert Space which the Slater determinant belongs to
  // index = index of the particle being examinate
  // coef = coefficient of the state being generate	
  virtual void GeneratesDifferentState( map <unsigned long, LongRational> & sortingMap , unsigned long* slater, unsigned long* state,FermionOnSphereTwoLandauLevels * slaterSpace, int index, LongRational coef);
  
  // compute the product of a monomial and a Slater determinant belonging in two Landau levels
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlater(unsigned long* slater,unsigned long* monomial, map <unsigned long, double> & sortingMap,FermionOnSphereTwoLandauLevels * finalSpace);
  
  // compute the product of a monomial and a Slater determinant belonging in two Landau levels
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  // return value = number of different obtained states
  virtual void MonomialsTimesSlater(unsigned long* slater,unsigned long* monomial, map <unsigned long, LongRational> & sortingMap ,FermionOnSphereTwoLandauLevels * finalSpace);
  
  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // finalStates = array where the obtained states are stored in their fermionic representation
  // sortingMap = map in which the generated states and their coefficient will be stored
  virtual void SlaterTimesSlaterProjection(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, map <unsigned long, double> & sortingMap, BosonOnSphereShort* finalSpace);
	
  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels with reverse flux attachment
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  void SlaterTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, map <unsigned long, double> & sortingMap, BosonOnSphereShort* finalSpace);

  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void SlaterTimesSlaterProjection(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, map <unsigned long, LongRational> & sortingMap, BosonOnSphereShort* finalSpace);
  
  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels with reverse flux attachment
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void SlaterTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, map <unsigned long, LongRational> & sortingMap, BosonOnSphereShort* finalSpace);
  
  // convert to standard fermionic representation
  //
  // monomial = reference to state stored in monomial format
  // landaulevel = the relevant landau level
  // return value = fermionic representation
  unsigned long ConvertFromPowerLandauRepresentation(unsigned long* monomial,unsigned long* landauLevel);
	
  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void SlaterTimesSlater(unsigned long* slater, unsigned long* lllslater, map <unsigned long, double> & sortingMap, BosonOnSphereTwoLandauLevels* finalSpace);

  // find state index. not using lookup table at the moment
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);
	
  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);

  // compute the product and the projection of a Slater determinant and a monomial with reverse flux attachment
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // landau =  array where the landau level of fermions is stored
  // sortingMap = map in which the generated states and their coefficient will be stored
  // binomialsCoefficient = binomials coefficient needed in the computation
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* monomial, map <unsigned long, double> & sortingMap, BinomialCoefficients& binomialsCoefficient, FermionOnSphere* finalSpace);


	// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in three Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

	virtual void SecondLandauLevelSlaterTimesSlaterProjection(unsigned long* slater1, unsigned long* slater2, map <unsigned long, double> & sortingMap, BosonOnSphereShort* finalSpace);
	

};

// convert a fermionic state to its monomial representation
//
// index = place in StateDescription of the fermionic state
// monomial = reference on the array where the monomial representation has to be stored

inline void FermionOnSphereTwoLandauLevels::ConvertToMonomial(int index , unsigned long*& monomial)
{
  this->ConvertToMonomial(this->StateDescription[index],monomial);
}

// convert a fermionic state to its monomial representation
//
// state = initial fermionic state in its fermionic representation
// monomial = reference on the array where the monomial representation has to be stored

inline void FermionOnSphereTwoLandauLevels::ConvertToMonomial(unsigned long state , unsigned long*& monomial)
{
  int Pos=0;
  for(int i= 2 * this->LzMax + 1; i > 0; --i)
    {
      if (((state >> i) & 0x1ul) != 0x0ul)
	{
	  monomial[Pos] = i;
	  ++Pos;
	}
    }
}

// convert a fermionic state from its monomial representation
//
// monomial = the array where the monomial representation is stored
// return value = fermionic state in its fermionic representation

inline unsigned long FermionOnSphereTwoLandauLevels::ConvertFromMonomial(unsigned long* monomial)
{
  unsigned long TmpState = 0x0ul;
  for (int i = 0; i < this->NbrFermions; i++)
    {
      TmpState |= 0x1ul << monomial[i];
    }
  return TmpState;
}

// convert a fermionic state to its monomial variable representation
//
// state = initial fermionic state in its fermionic representation
// slater = reference on the array where the monomial representation has to be stored
// nbrVariable = number of fermions in the second Landau level
// variable = reference on the array where the indice of fermions in the second Landau level has to be stored

inline void FermionOnSphereTwoLandauLevels::ConvertToMonomialVariable(unsigned long state, unsigned long*& slater,int & nbrVariable,unsigned long *& variable)
{
  unsigned long Tmp;
  int Pos=0;
  for(int i= this->LzMax; i >= 0;i--)
    {
      Tmp=  ((state >> (i << 1)) & 0x3ul);
      switch (Tmp)
	{
	case 0x1l: 
	  {
	    slater[Pos]=i;
	    Pos++;
	    break;
	  }
	case 0x2l:
	  {
	    slater[Pos]=i;
	    variable[nbrVariable]=Pos;
	    nbrVariable++;
	    Pos++;
	    break;
	  }
	case 0x3l:
	  {
	    slater[Pos]=i;
	    variable[nbrVariable]=Pos;
	    nbrVariable++;
	    Pos++;
	    slater[Pos]=i;
	    Pos++;
	    break;
	  }
	  break;
	default : 
	  {
	    break;
	  }
	}
    }
}

inline unsigned long FermionOnSphereTwoLandauLevels::ConvertFromPowerLandauRepresentation(unsigned long* monomial,unsigned long* landauLevel)
{
  unsigned long TmpState = 0x0ul;
  for (int i = 0; i < this->NbrFermions; i++)
    {
      TmpState |= 0x1ul << (((monomial[i] + (this->LzShiftDown * (1-landauLevel[i]))) << 1) + landauLevel[i]);
    }
  return TmpState;
}

#endif


