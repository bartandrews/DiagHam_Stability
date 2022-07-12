////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                          class of bosons on square lattice                 //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 16/09/2011                      //
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


#ifndef PARTICLEONSQUARELATTICEWANNIERINTERFACE_H
#define PARTICLEONSQUARELATTICEWANNIERINTERFACE_H

#include "config.h"

//#include "HilbertSpace/ParticleOnTorus.h"

class ParticleOnTorus;

class ParticleOnSquareLatticeWannierInterface
{

 protected:

 public:

  // virtual destructor
  virtual ~ParticleOnSquareLatticeWannierInterface();

  // get total momentum in the linearized momentum index (moduly N_phi)
  // index = state to operate on
  virtual int GetLinearizedMomentum(int index) = 0;

  // get total momentum in the linearized momentum index (moduly N_phi)
  // state = index of state to operate on
  virtual int ProjectToTorus(ParticleOnTorus *torusSpace, int state) = 0;

  
};

#endif
