////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of random number generator based on the stdlib rand function      //
//                                                                            //
//                        last modification : 15/09/2004                      //
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
#include "MathTools/RandomNumber/RanluxRandomNumberGenerator.h"


// constructor
//
// seed = optional seed definition (0 if none)

RanluxRandomNumberGenerator::RanluxRandomNumberGenerator(unsigned int& seed)
{
  if (seed != 0)
    srand (seed);
  this->NbrGeneratedNumbers = 0ul;
}

// copy constructor
//
// generator = generator to copy

RanluxRandomNumberGenerator::RanluxRandomNumberGenerator(const RanluxRandomNumberGenerator& generator)
{
  this->NbrGeneratedNumbers = generator.NbrGeneratedNumbers;
}

// destructor
//

RanluxRandomNumberGenerator::~RanluxRandomNumberGenerator()
{
}

// clone random number generator 
//
// return value = clone of the random number generator

AbstractRandomNumberGenerator* RanluxRandomNumberGenerator::Clone ()
{
  return new RanluxRandomNumberGenerator(*this);
}

