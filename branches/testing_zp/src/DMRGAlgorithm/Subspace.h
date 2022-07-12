////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of decribing subspace                       //
//                                                                            //
//                        last modification : 26/04/2001                      //
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


#ifndef SUBSPACE_H
#define SUBSPACE_H


#include "config.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


class Subspace
{

 protected:

  int SubspaceDimension;

  AbstractQuantumNumber* QuantumNumber;
  SubspaceSpaceConverter Converter;

  RealMatrix VectorsInSubspace; 

  RealDiagonalMatrix Hamiltonian;

 public: 

  // default constructor
  //
  Subspace ();

  // constructor for a not yet defined subspace coresponding to a given quantum number
  // 
  // Q = pointer to quantum number
  Subspace (AbstractQuantumNumber* Q);

  // copy constructor
  //
  // subspace = subspace to copy
  Subspace (const Subspace& subspace);

  // destructor
  //
  ~Subspace();

  // assignment 
  //
  // subspace = subspace to assign
  // return value = reference to current subspace
  Subspace& operator = (const Subspace& subspace);

  // get subspace dimension
  //
  // return value = subspace dimension
  int GetSubspaceDimension ();

  // set subspace dimension
  // 
  // dimension = new subspace dimension
  //  void SetSubspaceDimension (int dimension);

  // get reference to associated subspace-space converter
  //
  // rerturn value = reference to converter
  SubspaceSpaceConverter& GetConverter ();

  // get reference to diagonalized hamiltonian associated to the current subspace
  //
  // return value = reference to hamiltonian
  RealDiagonalMatrix& GetHamiltonian ();
  
  // get reference to Transformation Matrix from eigenvector base to canonical subspace base 
  //
  // return value = reference to coresponding matrix
  RealMatrix& GetMatrix ();

};

#endif
