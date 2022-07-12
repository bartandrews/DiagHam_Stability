////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                a set of functions usefull for integer algebra              //
//                                                                            //
//                        last modification : 11/09/2003                      //
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


#ifndef INTEGERALGEBRATOOLS_H
#define INTEGERALGEBRATOOLS_H


// find greatest common divider
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD
int FindGCD(int m, int n);

// find greatest common divider (recurisive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD
int RecursiveFindGCD(int m, int n);

// get all binomial coefficients up to a given number of element
//
// n = maximum number of elements
// return value = array containing the binomial coefficients (first index (i) corresponding to the number of elements, the second index going from 0 to i)
long** GetBinomialCoefficients (int n); 

// get all dimensions of irreducible representations of the symmetric group
//
// n = maximum number of elements
// return value = array containing binomial coeffcient (first index corresponds to the number of elements, the second is the number of indices,
//                if the number of indices is zero, the dimension is then equal to one)
long** GetIrreducibleRepresentationDimensionSymmetricGroup (int n);

#endif
