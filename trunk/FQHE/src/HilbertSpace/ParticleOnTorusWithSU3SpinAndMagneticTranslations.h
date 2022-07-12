////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of particle with SU3 spin on a torus                   //
//                taking into account magnetic translations                   //
//                                                                            //
//                        last modification : 18/06/2012                      //
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


#ifndef PARTICLEONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONS_H
#define PARTICLEONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU3Spin.h"

#include <iostream>


using std::ostream;


class Matrix;


class ParticleOnTorusWithSU3SpinAndMagneticTranslations :  public ParticleOnSphereWithSU3Spin
{

 public:

  // virtual destructor
  //
  virtual ~ParticleOnTorusWithSU3SpinAndMagneticTranslations ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
  
};


#endif
