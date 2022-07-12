////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of the equivalent of clebsch gordan coefficients          //
//          for addition of two angular monenta for Landau states in the LLL  //
//                       and on the thin annulus geometry                     //
//                                                                            //
//                        last modification : 07/11/2021                      //
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


#ifndef CLEBSCHGORDANTHINANNULUSCOEFFICIENTS_H
#define CLEBSCHGORDANTHINANNULUSCOEFFICIENTS_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/ClebschGordanDiskCoefficients.h"

#include <iostream>


using std::ostream;


class ClebschGordanThinAnnulusCoefficients : public ClebschGordanDiskCoefficients
{

 public:

  // default constructor
  //
  ClebschGordanThinAnnulusCoefficients();

  // constructor 
  //
  // mmax = maximum angular momentum for a single particle
  ClebschGordanThinAnnulusCoefficients(int mmax);

  // copy constructor (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to copy
  ClebschGordanThinAnnulusCoefficients (const ClebschGordanThinAnnulusCoefficients& coefficients);

  // destructor
  //
  ~ClebschGordanThinAnnulusCoefficients ();

  // assignment (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to assign
  // return value = reference on current Clebsch Gordan coefficients
  ClebschGordanThinAnnulusCoefficients& operator = (const ClebschGordanThinAnnulusCoefficients& coefficients);


 protected:
  
  // evaluate all Clebsch Gordan coefficients using Schulten Gordon recursion algorithm
  //
  void EvaluateClebschGordanThinAnnulusCoefficients();

};

#endif


