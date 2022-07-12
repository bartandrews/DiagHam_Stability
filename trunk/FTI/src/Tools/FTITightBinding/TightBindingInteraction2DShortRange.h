////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                           DiagHam  version 0.01                            //
//                                                                            //
//                  Copyright (C) 2001-2017 Nicolas Regnault                  //
//                                                                            //
//                       class Author: Gunnar Möller                          //
//                                                                            //
//              class providing periodic interactions on a 2D lattice         //
//                                                                            //
//                      last modification : 16/01/2017                        //
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


#ifndef TIGHTBINDINGINTERACTION2DSHORTRANGE_H
#define TIGHTBINDINGINTERACTION2DSHORTRANGE_H


#include "config.h"
#include "TightBindingInteraction2DShortRange.h"
#include "Abstract2DTightBindingModel.h"
#include "AbstractTightBindingInteraction.h"

/* #include <iostream> */
/* using std::ostream; */


class TightBindingInteraction2DShortRange : public AbstractTightBindingInteraction
{
 protected:
  Abstract2DTightBindingModel *TightBindingModel; // pointer to a tight-binding model
  double TwistAngle; 	  // angle between unit vectors (radians)
  int NbrCells1;          // no. of sites in sim. cell in 1 direction
  int NbrCells2;          // no. of sites in sim. cell in 2 direction
  int NbrImages;             // number of images used for summation
    
  RealVector BasisVector1;         // first basis vector
  RealVector BasisVector2;         // second basis vector
        
  RealMatrix U;           // vector storing the potential terms (1st index: sublattice, 2nd index: linearized index of 2nd site)


 public:

  // constructor
  TightBindingInteraction2DShortRange(Abstract2DTightBindingModel *tightBinding, InteractionType interaction, int nbrImages);
  
  // default destructor
  virtual ~TightBindingInteraction2DShortRange();

  // accessor routine to get the magnitude of the interaction
  // obtain Amplitude for a given separation
  // s = sublattice of the initial site
  // dR1 = displacement final site in unit cells along 1-direction
  // dR2 = displacement final site in unit cells along 2-direction
  // s2 = sublattice index of the final site
  virtual double GetAmplitude(int s, int dR1, int dR2, int s2);

 protected:

  // calculate the periodized version of the interaction over images in adjacent unit cells
  // s = sublattice of the initial site
  // dR1 = displacement final site in unit cells along 1-direction
  // dR2 = displacement final site in unit cells along 2-direction
  // s2 = sublattice index of the final site
  // nbrImages = maximum number of images to sum over
  double SumToConvergence(int s, int dR1, int dR2, int s2, int nbrImages=300);

};

#endif  // TIGHTBINDINGINTERACTION2DSHORTRANGE_H
