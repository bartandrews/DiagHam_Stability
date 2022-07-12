////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//                       class of Hamiltonian H=S*S on lattice                //
//                                                                            //
//                        last modification : 31/07/2009                      //
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

#ifndef SUNSPINONLATTICEQUADRATICHAMILTONIAN_H
#define SUNSPINONLATTICEQUADRATICHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/GenericSUNSpinCollection.h"
#include "Hamiltonian/AbstractSUNSpinOnLatticeHamiltonian.h"
#include "Tools/LatticeConnections.h"


class SUNSpinOnLatticeQuadraticHamiltonian : public AbstractSUNSpinOnLatticeHamiltonian
{
 protected:
  // definitions of the underlying lattice
  LatticeConnections *Lattice;
  
  // prefactor of cyclic permutation terms
  double CyclicPrefactor;
 public:
  // constructor
  // space = Hilbert space for problem
  // lattice = class providing size and geometry / connections on lattice
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // cyclic = prefactor of optional cyclic permutations around plaquettes
  //
  SUNSpinOnLatticeQuadraticHamiltonian(GenericSUNSpinCollection *space, LatticeConnections *lattice, 
				       AbstractArchitecture* architecture, long memory = -1,
				       char* precalculationFileName = 0, double cyclic = 0.0);

  // destructor
  //
  ~SUNSpinOnLatticeQuadraticHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();

  
 protected:

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionTerms();

  

};

#endif
