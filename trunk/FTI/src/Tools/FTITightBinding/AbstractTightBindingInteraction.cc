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


#include "config.h"
#include <iostream>

#include "AbstractTightBindingInteraction.h"



// calculate the periodized version of the interaction over images in adjacent unit cells
// s = sublattice of the initial site
// dR1 = displacement final site in unit cells along 1-direction
// dR2 = displacement final site in unit cells along 2-direction
// s2 = sublattice index of the final site
// nbrImages = maximum number of images to sum over
double AbstractTightBindingInteraction::SumToConvergence(int s, int dR1, int dR2, int s2, int nbrImages)
{
  std::cout << "Attention, dummy method AbstractTightBindingInteraction::SumToConvergence called"<<std::endl;
}
