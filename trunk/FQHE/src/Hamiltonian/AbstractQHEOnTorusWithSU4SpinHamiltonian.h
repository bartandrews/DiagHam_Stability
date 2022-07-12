////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract quantum Hall hamiltonian associated          //
//                   to particles on a torus with SU(6) spin                  //
//                                                                            //
//                        last modification : 26/06/2012                      //
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


#ifndef ABSTRACTQHEONTORUSWITHSU4SPINHAMILTONIAN_H
#define ABSTRACTQHEONTORUSWITHSU4SPINHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSU4SpinHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;

class AbstractQHEOnTorusWithSU4SpinHamiltonian : public AbstractQHEOnSphereWithSU4SpinHamiltonian
{
  
 protected:
  
  
  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;
  
		
 public:

  // destructor
  //
  virtual ~AbstractQHEOnTorusWithSU4SpinHamiltonian() = 0;
  
  
};

#endif
