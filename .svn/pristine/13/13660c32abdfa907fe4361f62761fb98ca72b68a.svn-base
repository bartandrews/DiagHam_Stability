////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              class of Fermions                             //
//                                                                            //
//                        last modification : 11/05/2001                      //
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


#ifndef FERMIONS_H
#define FERMIONS_H


#include "config.h"
#include "GeneralTools/List.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


#include <iostream>


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class Fermions :  public AbstractHilbertSpace
{

 protected:

  int NbrSite;
  int NbrFermion;
  int Sz;

  bool FixedQuantumNumberFlag;
  
  int* LookUpTable;
  unsigned long LookUpTableMask;
  int LookUpTableSize;
  int LookUpPosition;

  unsigned long* FermionDescription;
  unsigned long* Parity;

 public:


  // default constructor
  //
  Fermions ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz nor on
  // total number of fermions 
  //
  // nbrSite = number of site
  // memorySize = memory size in bytes allowed for look-up table
  Fermions (int nbrSite, int memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz and a
  // given number of fermions
  //
  // nbrSite = number of site
  // n = number of fermions
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Fermions (int nbrSite, int n, int sz, int memorySize) ;

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on fermions to copy
  Fermions (const Fermions& fermions);

  // destructor
  //
  ~Fermions ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on fermions to copy
  // return value = reference on current fermions
  Fermions& operator = (const Fermions& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  int TotalSz (int index);

  // get number of particles for a given state
  //
  // index = index of the state to test
  // return value = number of particles
  int GetNumberParticle (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  int FindStateIndex(unsigned long state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 private:

  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // fermionDescription = array describing states
  // parity = array describing parity
  // nbrSite = number of spin sites
  // n = number of fermions
  // sz = twice the value of total Sz component
  // fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
  // lookUpTable = look-up table
  // lookUpTableSize = look-Up table size
  // lookUpTablePosition = last position described by the look-Up table
  // lookUpTableMask = look-Up table mask  
  Fermions (int hilbertSpaceDimension, unsigned long* fermionDescription, unsigned long* parity,
	    int nbrSite, int n, int sz, bool fixedQuantumNumberFlag, int* lookUpTable, 
	    int lookUpTableSize, int lookUpPosition, unsigned long lookUpTableMask);
  
  // generate all states
  //
  // statePosition = position for the new states
  // sitePosition = site where fermion has to be changed
  // currentStateDescription = description of current state
  // return value = number of generated states
  int GenerateStates(int statePosition, int sitePosition, unsigned long currentStateDescription);

  // generate all states corresponding to a given total Sz and a total number of fermions
  //
  // statePosition = position for the new states
  // sitePosition = site where fermion has to be changed
  // currentStateDescription = description of current state
  // currentN = current number of fermions
  // currentSz = total Sz value of current state
  // return value = number of generated states
  int GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription, 
		     int currentN, int currentSz);

  // Evaluate Hilbert space dimension
  // 
  // nbrSite = number of sites
  // nUp = number of fermion with spin up
  // nDown = number of fermion with spin down
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension (int nbrSite, int nUp, int nDown);

  // evaluate parity of a given state
  //
  // state = description of the state to evaluate
  // return value = parity
  unsigned long EvaluateParity(unsigned long state);

};

#endif


