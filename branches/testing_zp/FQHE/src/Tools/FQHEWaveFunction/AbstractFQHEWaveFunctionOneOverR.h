////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of abstract wave function for FQHE that evalute both         //
//              the wave function and the wave function time the Coulomb      //
//                                interaction term                            //
//                                                                            //
//                        last modification : 23/10/2006                      //
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


#ifndef ABSTRACTFQHEWAVEFUNCTIONONEOVERR_H
#define ABSTRACTFQHEWAVEFUNCTIONONEOVERR_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class AbstractFQHEWaveFunctionOneOverR: public Abstract1DComplexFunction
{

 public:

  // destructor
  //
  virtual ~AbstractFQHEWaveFunctionOneOverR();

  // evaluate the norm to the square of the wave function at a given point time the coulomb term (assume the coordinates are those provides by the previous operator() method call)
  //
  // return value = corresponding numerical value
  virtual double CoulombContribution() = 0;

};

#endif
