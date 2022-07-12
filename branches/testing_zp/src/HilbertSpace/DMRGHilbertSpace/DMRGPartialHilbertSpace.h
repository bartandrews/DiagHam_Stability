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


#ifndef DMRGPARTIALHILBERTSPACE_H
#define DMRGPARTIALHILBERTSPACE_H


#include "config.h"
#include "GeneralTools/List.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "HilbertSpace/SpaceDecomposition.h"
#include "GeneralTools/GarbageFlag.h"


class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class DMRGPartialHilbertSpace :  public AbstractHilbertSpace
{

 protected:

  List<AbstractQuantumNumber*> ListQuantumNumber;

  RealDiagonalMatrix Energies;

  SpaceDecomposition Decomposition;

 public:

  // default constructor
  //
  DMRGPartialHilbertSpace ();
  
  // constructor from datas 
  //
  // energies = matrix conatining energies
  // quantumNumbers = list of quantum numbers
  // decomposition = space decompostion with respect to the subspaces
  DMRGPartialHilbertSpace (RealDiagonalMatrix&  energies, List<AbstractQuantumNumber*> quantumNumbers, 
			   const SpaceDecomposition& decomposition);
  
  // copy constructor (without duplicating datas)
  //
  // space = reference on DMRG partial Hilbert space to copy
  DMRGPartialHilbertSpace (const DMRGPartialHilbertSpace& space);

  // destructor
  //
  ~DMRGPartialHilbertSpace ();

  // assignement (without duplicating datas)
  //
  // space = reference on DMRG partial Hilbert space to copy
  // return value = reference on current fermions
  DMRGPartialHilbertSpace& operator = (const DMRGPartialHilbertSpace& space);

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

  // get energy of a given state
  //
  // index = index of the state to test
  // return value = energy
  double GetEnergy (int index);

  // get description of a given subspace
  //
  // index = index of the subspace to describe
  // return value = reference on subspace description
  SubspaceSpaceConverter& GetSubspaceDescription (int index);
  
  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

};

#endif


