////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                         to particles on a torus with                       //
//                                                                            //
//                        last modification : 27/06/2003                      //
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


#ifndef ABSTRACTQHEONTORUSHAMILTONIAN_H
#define ABSTRACTQHEONTORUSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorus.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEOnTorusHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{

  friend class QHEParticlePrecalculationOperation;
  
 protected:
  
  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

  // number of Lz values in a state
  int NbrLzValue;
  
  // flag to indicate that only some terms in the interaction should be kept
  bool FilterInteractionFlag;
  // Indices of the annihilation (or creation) terms that have to be kept. Indices are stored a linearized index \sum_{j=1}^{nbr_nbody} i_j * N_\phi^j 
  int* FilterInteractionIndices;
  // number of kept indices 
  int NbrFilterInteractionIndices;

 public:

  // default constructor
  //
  AbstractQHEOnTorusHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnTorusHamiltonian() = 0;

 protected:

  // get all the indices that should appear in the annihilation/creation operators
  //
  virtual void GetIndices();

  // find all the indices that have to be kept when filtering the hamiltionian
  //
  // filterInteractionFile = name of the file that describe which terms in the interaction should be kept
  virtual void FindFilteredIndices(char* filterInteractionFile);

};


#endif
