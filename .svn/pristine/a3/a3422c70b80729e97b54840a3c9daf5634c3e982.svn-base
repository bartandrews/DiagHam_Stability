////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2017 Nicolas Regnault                  //
//                                                                            //
//                         class Author: Gunnar Möller                        //
//                                                                            //
//              class providing interface for interactions on lattice         //
//                                                                            //
//                        last modification : 16/01/2017                      //
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


#ifndef ABSTRACTTIGHTBINDINGINTERACTION_H
#define ABSTRACTTIGHTBINDINGINTERACTION_H


#include "config.h"

#include "AbstractTightBindingModel.h"

/* #include <iostream> */
/* using std::ostream; */


class AbstractTightBindingInteraction
{

 public:

  // type definition for functional dependency of interactions 
  // as a generic function pointer without additional parameters
  typedef double (*InteractionType)(double);

  // default destructor
  virtual ~AbstractTightBindingInteraction(){}

  // accessor routine to get the magnitude of the interaction
  // obtain Amplitude for a given separation
  // s = sublattice of the initial site
  // dR1 = displacement final site in unit cells along 1-direction
  // dR2 = displacement final site in unit cells along 2-direction
  // s2 = sublattice index of the final site
  virtual double GetAmplitude(int s, int dR1, int dR2, int s2) = 0;

 protected:

  // calculate the periodized version of the interaction over images in adjacent unit cells
  // s = sublattice of the initial site
  // dR1 = displacement final site in unit cells along 1-direction
  // dR2 = displacement final site in unit cells along 2-direction
  // s2 = sublattice index of the final site
  // nbrImages = maximum number of images to sum over
  virtual double SumToConvergence(int s, int dR1, int dR2, int s2, int nbrImages=300);

};

#endif  // ABSTRACTTIGHTBINDINGINTERACTION_H
