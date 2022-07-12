////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//             class of hilbert space of two periodic 3d particles             //
//                                                                            //
//                        last modification : 18/10/2004                      //
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


#ifndef PERIODICTHREEDTWOPARTICLES_H
#define PERIODICTHREEDTWOPARTICLES_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/PeriodicThreeDOneParticle.h"
#include "HilbertSpace/ThreeDTwoParticles.h"


class PeriodicThreeDTwoParticles : public ThreeDTwoParticles
{

 protected:

 public:

  // default constructor
  //
  PeriodicThreeDTwoParticles ();

  // constructor
  //
  // firstParticle = pointer to the first abstract particle
  // secondParticle = pointer to the second abstract particle
  PeriodicThreeDTwoParticles (PeriodicThreeDOneParticle* firstParticle, PeriodicThreeDOneParticle* secondParticle);

  // copy constructor
  //
  // space = reference on Hilbert space to copy
  PeriodicThreeDTwoParticles (const PeriodicThreeDTwoParticles& space);

  // destructor
  //
  virtual ~PeriodicThreeDTwoParticles();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignement
  //
  // space = reference on Hilbert space to assign
  // return value = reference on current Hilbert space
  PeriodicThreeDTwoParticles& operator = (const PeriodicThreeDTwoParticles& space);
  
  // get the Hilbert space description of the first particle
  //
  // return = pointer to the Hilbert space
  ThreeDOneParticle* GetFirstParticleSpace ();

  // get the Hilbert space description of the second particle
  //
  // return = pointer to the Hilbert space
  ThreeDOneParticle* GetSecondParticleSpace ();
  
};

// get the Hilbert space description of the first particle
//
// return = pointer to the Hilbert space

inline ThreeDOneParticle* PeriodicThreeDTwoParticles::GetFirstParticleSpace ()
{
  PeriodicThreeDOneParticle* firstParticle = (PeriodicThreeDOneParticle*) this->FirstParticle->Clone ();
  return firstParticle;
}

// get the Hilbert space description of the second particle
//
// return = pointer to the Hilbert space

inline ThreeDOneParticle* PeriodicThreeDTwoParticles::GetSecondParticleSpace ()
{
  PeriodicThreeDOneParticle* secondParticle = (PeriodicThreeDOneParticle*) this->SecondParticle->Clone ();
  return secondParticle;
}

#endif
