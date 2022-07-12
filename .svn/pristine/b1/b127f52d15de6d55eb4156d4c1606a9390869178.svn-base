////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of abstract 1D complex function                   //
//                                                                            //
//                        last modification : 01/09/2004                      //
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


#ifndef ABSTRACT1DCOMPLEXFUNCTION_H
#define ABSTRACT1DCOMPLEXFUNCTION_H


#include "config.h"
#include "MathTools/Complex.h"


class RealVector;


class Abstract1DComplexFunction
{

 public:

  enum ComplexFunctionProperties
    {
      Basic = 0x1u,
      OnSphere = 0x2u,
      Trial = 0x4u,
      Normalizable = 0x10000
    };
    

  // virtual destructor
  //
  virtual ~Abstract1DComplexFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone () = 0;

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x) = 0;

  // change the normalization of the funtion by a multiplicative factor
  // factor = factor to be multiplied
  virtual void Renormalize(double factor);

  // get function properties, and possible extensions of interface 
  // 
  virtual unsigned GetProperties() {return Abstract1DComplexFunction::Basic;}

  // set wavefunction to one for a given set of particle coordinates
  virtual void AdaptNorm(RealVector& x);
  
  // utility function to set the right dynamic interval for Monte-Carlo
  virtual void AdaptAverageMCNorm(int thermalize = 500 , int average = 1000);

  
};

#endif
