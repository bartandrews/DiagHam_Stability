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


#ifndef ABSTRACT1DCOMPLEXTRIALFUNCTIONONSPHERE_H
#define ABSTRACT1DCOMPLEXTRIALFUNCTIONONSPHERE_H


#include "config.h"
#include "MathTools/Complex.h"
#include "Abstract1DComplexTrialFunction.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


class RealVector;


class Abstract1DComplexTrialFunctionOnSphere : public Abstract1DComplexTrialFunction
{
 protected:
  double *TrialParameters;
  int NbrParameters;
  
 public:

  // virtual destructor
  //
  virtual ~Abstract1DComplexTrialFunctionOnSphere();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone () = 0;

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x) = 0;

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv) = 0;

  // get a value of the wavefunction for the last set of coordinates, but with different variational parameters
  // parameters =  alternative set of parameters
  virtual Complex GetForOtherParameters( double *parameters) = 0;

  // do many evaluations of the function, storing the result in the vector results given in the call
  // result = vector of leading dimension of the array coefficients for returning values
  // x = positions to evaluate the wavefuntion in
  // format for passing parameters as [nbrSet][nbrParameter],
  virtual void GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients) = 0;

  // access internal values of parameters
  virtual double *GetTrialParameters(){return this->TrialParameters;}

  // get number of parameters
  virtual int GetNbrParameters() {return this->NbrParameters;}
  
  // set new values of the trial coefficients (keeping the initial number of parameters)
  virtual void SetTrialParameters(double * coefficients) = 0;


  
};

#endif
