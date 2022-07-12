////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including four               //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 14/10/2010                      //
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


#ifndef FERMIONONSPHEREFOURLANDAULEVELS_H
#define FERMIONONSPHEREFOURLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>
#include <map>

using std::map;

class FermionOnSphere;
class BosonOnSphereShort;


class FermionOnSphereFourLandauLevels : public FermionOnSphereWithSU4Spin
{

 protected:
  
  
 public:
  
  // default constructor
  // 
  FermionOnSphereFourLandauLevels ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereFourLandauLevels (int nbrFermions, int totalLz, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereFourLandauLevels(const FermionOnSphereFourLandauLevels& fermions);

  // destructor
  //
  virtual ~FermionOnSphereFourLandauLevels ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereFourLandauLevels& operator = (const FermionOnSphereFourLandauLevels& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

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
  virtual void BosonicStateTimeFermionicState( RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, BosonOnSphereShort* bosonSpace, FermionOnSphere * finalSpace, int firstComponent, int nbrComponent , bool reverseFluxFlag = false);
  
  // compute the product and the projection of a Slater determinant and a monomial 
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable, int nbrVariable, map <unsigned long , double > & sortingMap, FermionOnSphere* finalSpace);
	
  // compute the product and the projection of a Slater determinant and a monomial with reverse flux attachment
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // landau =  array where the landau level of fermions is stored
  // sortingMap = map in which the generated states and their coefficient will be stored
  // binomialsCoefficient = binomials coefficient needed in the computation
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* monomial, unsigned long* variable, map <unsigned long, double > & sortingMap, BinomialCoefficients& binomialsCoefficient, FermionOnSphere* finalSpace);
  
	
  // compute the projection of the product of a bosonic state and a fermionic state
  //
  // bosonState = real vector where the bosonic state is stored
  // fermionState = real vector where the fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // bosonSpace = pointer to the bosonic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void BosonicStateTimeFermionicState(LongRationalVector& bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, BosonOnSphereShort* bosonSpace , FermionOnSphere * finalSpace, int firstComponent, int nbrComponent);
  
  // compute the product and the projection of a Slater determinant and a monomial 
  // 
  // slater = array where the slater is stored in its monomial representation
  // monomial = array where the monomial is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable, int nbrVariable, map <unsigned long, LongRational > & sortingMap, FermionOnSphere* finalSpace);
	

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);
	
  // compute the number of particles in each Landau level
  //
  // state = ID of the state to handle
  // lLOccupationConfiguration = array where the decomposition will be store
  virtual void LandauLevelOccupationNumber(int state, int * lLOccupationConfiguration);
  
  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in three Landau levels
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void SlaterTimesSlaterProjection(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, map <unsigned long , double> & sortingMap, BosonOnSphereShort* finalSpace);
  
  // compute the projection of the product of a bosonic state and a fermionic state
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void LLLFermionicStateTimeFermionicState( RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, FermionOnSphere* lllFermionSpace,BosonOnSphereShort* finalSpace, int firstComponent,int nbrComponent);
  
  // compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in three Landau levels
  //
  // slater = array where the slater determinant in the two landau levels is stored in its monomial representation
  // lllslater = array where the slater determinant in the LLL is stored in its monomial representation
  // variable = reference on the array where the indice of fermions in the second Landau level is stored
  // nbrVariable = number of fermions in the second Landau level
  // sortingMap = map in which the generated states and their coefficient will be stored
  // finalSpace = pointer to the final HilbertSpace
  virtual void SlaterTimesSlaterProjection(unsigned long* slater,unsigned long* lllslater,unsigned long * variable,int nbrVariable, map <unsigned long, LongRational > & sortingMap, BosonOnSphereShort* finalSpace);
  
  // compute the projection of the product of a bosonic state and a fermionic state
  //
  // lllFermionState = real vector where the lowest Landau level fermionic state is stored
  // fermionState = real vector where the two Landau level fermionic state is stored
  // outputVector = real vector where the result has to be stored
  // lllFermionSpace = pointer to the lowest Landau level Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed
  virtual void LLLFermionicStateTimeFermionicState(LongRationalVector& lllFermionState, LongRationalVector& fermionState, LongRationalVector & outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent);
  
 protected:
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // nbrFluxQuanta = number of flux quanta
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // nbrFluxQuanta = number of flux quanta
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz, long pos);

  // convert a fermionic state to its monomial variable representation
  //
  // state = initial fermionic state in its fermionic representation
  // slater = reference on the array where the monomial representation has to be stored
  // nbrVariable = number of fermions in the second Landau level
  // variable = reference on the array where the indice of fermions in the second Landau and the third level has to be stored
  virtual void ConvertToMonomialVariable(unsigned long state, unsigned long*& slater, int& nbrVariable, unsigned long*& variable);
	
  // convert a fermionic state to its monomial variable representation
  //
  // state = initial fermionic state in its fermionic representation
  // slater = reference on the array where the monomial representation has to be stored
  // landau = array where landau level of each fermions has to be stored
  virtual void ConvertToMonomialLandau(unsigned long state, unsigned long*& landau, unsigned long*& variable);

};

// convert a fermionique state to its monomial variable representation
//
// state = initial fermionic state in its fermionic representation
// slater = reference on the array where the monomial representation has to be stored
// nbrVariable = number of fermions in the second Landau level
// variable = reference on the array where the indice of fermions in the second Landau and the third level has to be stored

inline void FermionOnSphereFourLandauLevels::ConvertToMonomialVariable(unsigned long state, unsigned long*& slater, int& nbrVariable, unsigned long*& variable)
{
  int TmpPos = 0;
  for(int i = this->LzMax; i >= 0; i--)
    {
      switch ((state >> (i << 2)) & 0xful)
	{
	case 0x1ul: 
	  {
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0x2ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0x3ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0x4ul: 
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable]= 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0x5ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0x6ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0x7ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0x8ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0x9ul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0xaul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0xbul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0xcul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0xdul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	case 0xeul:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	  }
	  break;
	case 0xful:
	  {
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 2;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos + 1;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    variable[nbrVariable] = 4 * TmpPos;
	    nbrVariable++;
	    TmpPos++;
	    slater[TmpPos] = i;
	    TmpPos++;
	  }
	  break;
	default : 
	  {
	    break;
	  }
	}
    }
}

// convert a fermionic state to its monomial variable representation
//
// state = initial fermionic state in its fermionic representation
// slater = reference on the array where the monomial representation has to be stored
// landau = array where landau level of each fermions has to be stored

inline void FermionOnSphereFourLandauLevels::ConvertToMonomialLandau(unsigned long state, unsigned long*& slater, unsigned long*& landau)
{
  int TmpPos = 0;
  for(int i = this->LzMax; i >= 0; i--)
    {
      switch ((state >> (i << 2)) & 0xful)
	{
	case 0x1ul: 
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0x2ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	  }
	  break;
	case 0x3ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0x4ul: 
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	  }
	  break;
	case 0x5ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0x6ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	  }
	  break;
	case 0x7ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0x8ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	  }
	  break;
	case 0x9ul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0xaul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	  }
	  break;
	case 0xbul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0xcul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	  }
	  break;
	case 0xdul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	case 0xeul:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	  }
	  break;
	case 0xful:
	  {
	    slater[TmpPos] = i;
	    landau[TmpPos] = 3;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 2;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 1;
	    TmpPos++;
	    slater[TmpPos] = i;
	    landau[TmpPos] = 0;
	    TmpPos++;
	  }
	  break;
	default : 
	  {
	    break;
	  }
	}
    }
}

#endif


