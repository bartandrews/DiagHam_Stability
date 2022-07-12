////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract random number generator                //
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
#include "MathTools/RandomNumber/AbstractRandomNumberGenerator.h"

#include <ctime>
#include <cmath>

// virtual destructor
//

AbstractRandomNumberGenerator::~AbstractRandomNumberGenerator()
{
}


// set seed of the random number generator to system time
//

void AbstractRandomNumberGenerator::UseTimeSeed()
{
  this->SetSeed(time(NULL));
}


// get real random number with gaussian distribution (uses multiple calls to generator)
//
// return value = random number

double AbstractRandomNumberGenerator::GetGaussianRandomNumber ()
{
  double fac, rsq, v1, v2;
  if (this->iset == 0)
    {
     
mark: v1=2.*GetRealRandomNumber()-1.0;      
      v2=2.*GetRealRandomNumber()-1.0;
      rsq=v1*v1+v2*v2;
      if (rsq >=1.0 || rsq == 0.0 ) goto mark;
      fac = sqrt(-2.*log(rsq)/rsq);
      this->gset=v1*fac;
      this->iset=1;
      return (v2*fac);
    }
  else 
    {
      this->iset=0;
      return this->gset;
    }
}
