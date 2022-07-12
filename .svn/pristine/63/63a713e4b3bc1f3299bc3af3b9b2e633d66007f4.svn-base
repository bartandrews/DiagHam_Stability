////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere including three               //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 15/12/2011                      //
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


#ifndef BOSONONSPHERETHREELANDAULEVELS_H
#define BOSONONSPHERETHREELANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"
#include "HilbertSpace/BosonOnSphereShort.h"

#include <iostream>
#include <map>

using std::map;


class BosonOnSphereThreeLandauLevels : public BosonOnSphereWithSU3Spin
{

 protected:

  
 public:
  
  // default constructor
  // 
  BosonOnSphereThreeLandauLevels ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // memory = amount of memory granted for precalculations
  BosonOnSphereThreeLandauLevels (int nbrBosons, int totalLz, int lzMax, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereThreeLandauLevels(const BosonOnSphereThreeLandauLevels& bosons);

  // destructor
  //
  virtual ~BosonOnSphereThreeLandauLevels ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereThreeLandauLevels& operator = (const BosonOnSphereThreeLandauLevels& bosons);
  
  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:
  
  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int nbrFluxQuanta, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta
  // lzMax = momentum maximum value for a boson in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int nbrFluxQuanta, int lzMax, int totalLz, long pos);

  // convert a bosonic state to its monomial variable representation
  //
  // state1 = initial bosonic state in its fermionic representation for the first Landau level
  // state2 = initial bosonic state in its fermionic representation for the second Landau level
  // state3 = initial bosonic state in its fermionic representation for the third Landau level
  // slater = reference on the array where the monomial representation has to be stored
  // nbrVariable = number of bosons in the second Landau level
  // variable = reference on the array where the indice of bosons in the second Landau and the third level has to be stored
  virtual void ConvertToMonomialVariable(unsigned long state1, unsigned long state2, unsigned long state3, 
					 unsigned long*& slater, int& nbrVariable, unsigned long*& variable);
	
  // convert a bosonic state to its monomial variable representation
  //
  // state1 = initial bosonic state in its fermionic representation for the first Landau level
  // state2 = initial bosonic state in its fermionic representation for the second Landau level
  // state3 = initial bosonic state in its fermionic representation for the third Landau level
  // slater = reference on the array where the monomial representation has to be stored
  // landau = array where landau level of each bosons has to be stored
  virtual void ConvertToMonomialLandau(unsigned long state1, unsigned long state2, unsigned long state3,
				       unsigned long*& slater, unsigned long*& landau);

};

// convert a bosonique state to its monomial variable representation
//
// state1 = initial bosonic state in its fermionic representation for the first Landau level
// state2 = initial bosonic state in its fermionic representation for the second Landau level
// state3 = initial bosonic state in its fermionic representation for the third Landau level
// slater = reference on the array where the monomial representation has to be stored
// nbrVariable = number of bosons in the second Landau level
// variable = reference on the array where the indice of bosons in the second Landau and the third level has to be stored

inline void BosonOnSphereThreeLandauLevels::ConvertToMonomialVariable(unsigned long state1, unsigned long state2, unsigned long state3,
								      unsigned long*& slater, int& nbrVariable, unsigned long*& variable)
{
/*   unsigned long Tmp; */
/*   int TmpPos=0; */
/*   for(int i=this->LzMax;i>=0;i--) */
/*     { */
/*       Tmp=((state >> (3*i)) & ((unsigned long) 0x7)); */
/*       switch (Tmp) */
/* 	{ */
/* 	case 0x1l:  */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x2l: */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x3l: */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos]=i; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x4l:  */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos+0x01ul; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x5l: */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos+0x01ul; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos]=i; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x6l: */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos+0x01ul; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x7l: */
/* 	  { */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos+0x01ul; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos]=i; */
/* 	    variable[nbrVariable]=2*TmpPos; */
/* 	    nbrVariable++; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos]=i; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	  break; */
/* 	default :  */
/* 	  { */
/* 	    break; */
/* 	  } */
/* 	} */
/*     } */
}

// convert a bosonic state to its monomial variable representation
//
// state1 = initial bosonic state in its fermionic representation for the first Landau level
// state2 = initial bosonic state in its fermionic representation for the second Landau level
// state3 = initial bosonic state in its fermionic representation for the third Landau level
// slater = reference on the array where the monomial representation has to be stored
// landau = array where landau level of each bosons has to be stored

inline void BosonOnSphereThreeLandauLevels::ConvertToMonomialLandau(unsigned long state1, unsigned long state2, unsigned long state3,
								    unsigned long*& slater, unsigned long*& landau)
{
/*   unsigned long Tmp; */
/*   int TmpPos = 0; */
/*   for(int i = this->LzMax; i >= 0; i--) */
/*     { */
/*       Tmp = ((state >> (3*i)) & ((unsigned long) 0x7)); */
/*       switch (Tmp) */
/* 	{ */
/* 	case 0x1l:  */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 0; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x2l: */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 1; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x3l: */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 1; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 0; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x4l:  */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 2; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x5l: */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 2; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 0; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x6l: */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 2; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 1; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	case 0x7l: */
/* 	  { */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 2; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 1; */
/* 	    TmpPos++; */
/* 	    slater[TmpPos] = i; */
/* 	    landau[TmpPos] = 0; */
/* 	    TmpPos++; */
/* 	    break; */
/* 	  } */
/* 	  break; */
/* 	default :  */
/* 	  { */
/* 	    break; */
/* 	  } */
/* 	} */
/*     } */
}

#endif

