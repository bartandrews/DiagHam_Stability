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


#include "config.h"
#include "HilbertSpace/PeriodicThreeDTwoParticles.h"


using std::cout;
using std::endl;


// default constructor
//

PeriodicThreeDTwoParticles::PeriodicThreeDTwoParticles ()
{
}

// constructor from two abstract particles
//
// firstParticle = pointer to the first abstract particle
// secondParticle = pointer to the second abstract particle

PeriodicThreeDTwoParticles::PeriodicThreeDTwoParticles (PeriodicThreeDOneParticle* firstParticle, PeriodicThreeDOneParticle* secondParticle)
{
  this->FirstParticle = (PeriodicThreeDOneParticle*) firstParticle;
  this->SecondParticle = (PeriodicThreeDOneParticle*) secondParticle;
  this->HilbertSpaceDimension = this->FirstParticle->GetHilbertSpaceDimension () * this->SecondParticle->GetHilbertSpaceDimension ();
}

// copy constructor
//
// space = reference on Hilbert space to copy

PeriodicThreeDTwoParticles::PeriodicThreeDTwoParticles (const PeriodicThreeDTwoParticles& space)
{
  this->FirstParticle = space.FirstParticle;
  this->SecondParticle = space.SecondParticle;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

PeriodicThreeDTwoParticles::~PeriodicThreeDTwoParticles ()
{
  delete (PeriodicThreeDOneParticle*) this->FirstParticle;
  delete (PeriodicThreeDOneParticle*) this->SecondParticle;
  this->FirstParticle = NULL; this->SecondParticle = NULL;
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PeriodicThreeDTwoParticles& PeriodicThreeDTwoParticles::operator = (const PeriodicThreeDTwoParticles& space)
{
  this->FirstParticle = space.FirstParticle;
  this->SecondParticle = space.SecondParticle;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PeriodicThreeDTwoParticles::Clone ()
{
  return new PeriodicThreeDTwoParticles (*this);
}
