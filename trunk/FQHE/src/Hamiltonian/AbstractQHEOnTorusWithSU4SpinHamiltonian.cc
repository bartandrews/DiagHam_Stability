////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract quantum Hall hamiltonian associated          //
//                   to particles on a torus with SU(4) spin                  //
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


#include "config.h"
#include "Hamiltonian/AbstractQHEOnTorusWithSU4SpinHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"


#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

AbstractQHEOnTorusWithSU4SpinHamiltonian::~AbstractQHEOnTorusWithSU4SpinHamiltonian()
{
}
