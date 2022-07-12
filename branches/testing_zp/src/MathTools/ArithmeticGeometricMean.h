////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of clebsch gordan coefficients                   //
//                                                                            //
//                        last modification : 19/06/2002                      //
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


#ifndef ARITHMETICGEOMETRICMEAN_H
#define ARITHMETICGEOMETRICMEAN_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ArithmeticGeometricMean
{

 private:

  // first series of numbers a_n
  double *A;

  // second series of numbers a_n
  double *B;

  // third series of numbers a_n
  double *C;

  // number of values calculated
  int Length;

  // precision required
  double Precision;

  // garbage flag
  GarbageFlag Flag;

  // do one step in the recursion
  // a = on entry a_n on return: a_n+1
  // b = on entry b_n ... 
  // c = on entry c_n ...
  void Recurse(double &a, double &b, double &c);
	       
    
 public:

  // default constructor
  //
  ArithmeticGeometricMean();

  // constructor 
  //
  // a0, b0, c0 = initial triple of numbers
  // precision = required precision
  ArithmeticGeometricMean(double a0, double b0, double c0, double precision=1e-15);

  // copy constructor (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to copy
  ArithmeticGeometricMean (const ArithmeticGeometricMean& agm);

  // destructor
  //
  ~ArithmeticGeometricMean ();

  // get number of non-zero values in series A_n, B_n, C_n
  int GetLength() { return Length;}
  
  // get number A_n
  double GetA (int n) {return A[n];}

  // get number B_n
  double GetB (int n) {return B[n];}

  // get number C_n
  double GetC (int n) {return C[n];}

};


#endif
