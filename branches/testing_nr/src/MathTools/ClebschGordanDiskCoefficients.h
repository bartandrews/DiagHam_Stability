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
//                              and on the disk geometry                      //
//                                                                            //
//                        last modification : 05/07/2008                      //
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


#ifndef CLEBSCHGORDANDISKCOEFFICIENTS_H
#define CLEBSCHGORDANDISKCOEFFICIENTS_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ClebschGordanDiskCoefficients
{

 private:

  // maximum angular momentum for a single particle
  int  MaximumMomentum;

  // Clebsch Gordan coefficient array accessed as Coefficients[j1][j2][j3]
  double*** Coefficients;
  // garbage flag associated to coefficient array
  GarbageFlag Flag;


  // position associated to the projection of first angular momentum for current iterator
  int M1;
  // position associated to the projection of second angular momentum for current iterator
  int M2;
  // resulting angular momentum for current iterator  
  int J;
  // position associated to the resulting angular momentum  for current iterator  
  int CurrentPosition;

 public:

  // default constructor
  //
  ClebschGordanDiskCoefficients();

  // constructor 
  //
  // mmax = maximum angular momentum for a single particle
  ClebschGordanDiskCoefficients(int mmax);

  // copy constructor (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to copy
  ClebschGordanDiskCoefficients (const ClebschGordanDiskCoefficients& coefficients);

  // destructor
  //
  ~ClebschGordanDiskCoefficients ();

  // assignment (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to assign
  // return value = reference on current Clebsch Gordan coefficients
  ClebschGordanDiskCoefficients& operator = (const ClebschGordanDiskCoefficients& coefficients);

  // get a particular coefficient (without testing if m1, m2 and j are valid)
  //
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum 
  // j = resulting angular momentum
  // return value = corresponding Clebsch Gordan coefficient
  double GetCoefficient (int m1, int m2, int j);

  // print a particular coefficient (without testing if m1, m2 and j are valid)
  //
  // str = reference on output stream
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum 
  // j = resulting angular momentum
  // return value = reference on output stream
  ostream& PrintCoefficient (ostream& str, int m1, int m2, int j);

 private:
  
  // evaluate all Clebsch Gordan coefficients using Schulten Gordon recursion algorithm
  //
  void EvaluateClebschGordanDiskCoefficients();

};

#endif


