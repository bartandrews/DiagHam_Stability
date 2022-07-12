////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//       class implementing Symmetrizer for complex function //
//                                                                            //
//                        last modification : 06/12/2007                      //
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


#ifndef SYMMETRIZEDCOMPLEXFUNCTION_H
#define SYMMETRIZEDCOMPLEXFUNCTION_H

#include "config.h"
#include "MathTools/Complex.h"
#include "Vector/RealVector.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"

class SymmetrizedComplexFunction : public Abstract1DComplexFunction
{
 public:
  // default constructor
  SymmetrizedComplexFunction();

  // construct Symmetric form of given function
  // originalFunction = function to be AS
  // nbrCoordinates = number of coordinates to be permuted
  // dimCoordinate = number of elements per coordinate
  SymmetrizedComplexFunction(Abstract1DComplexFunction *originalFunction, int nbrCoordinates, int dimCoordinate);

  // copy constructor
  SymmetrizedComplexFunction(SymmetrizedComplexFunction &F);

  // destructor
  ~SymmetrizedComplexFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();
  
  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

    
 protected:
  // size of input
  int NbrCoordinates;
  // number of dimensions per coordinate
  int CoordinateDimension;
  
  Abstract1DComplexFunction *OriginalFunction;
  // table for permutations
  int *Permutations;

  double FactorialPrefactor;
  
  // Complex to store sum of Symmetrizer
  Complex SSum;

  // counter for testing
  int NbrTerms;
  
  RealVector InternalCoordinates;
  RealVector *GivenCoordinates;

 private:
  // client function to generate permutations (from http://www.bearcave.com/random_hacks/permute.html)
  void AddSTerm(int sign);
  void Swap(const int i, const int j, int &sign);
  void RotateLeft(const int start, const int n, int &sign);
  void Permute(const int start, const int n, int sign);
};

#endif // SYMMETRIZEDCOMPLEXFUNCTION_H
