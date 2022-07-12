////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of Laughlin wave function on disk with one Jain quasielectron      //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 13/11/2008                      //
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


#ifndef FQHEDISKLAUGHLINONEQUASIELECTRONJAINWAVEFUNCTION_H
#define FQHEDISKLAUGHLINONEQUASIELECTRONJAINWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class FQHEDiskLaughlinOneJainQuasielectronWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // inverse value of the filling factor
  int InvFillingFactor;

  // invert of the maximum x value
  double InvScale;

  // quasielectron position 
  Complex ZElectron;

  // gaussian weight associated to the quasielectron
  double ElectronWeight;
  Complex* TmpElectronWeight;

  // temporary array where terms (and square of) of the Jastrow factor will be stored 
  Complex** TmpJastrow;
  Complex** TmpSqrJastrow;


 public:

  // constructor
  //
  // nbrParticles = number of particles
  // zElectron = quasielectron position
  // scale = typical sytem size
  // invFillingFactor = inverse value of the filling factor
  FQHEDiskLaughlinOneJainQuasielectronWaveFunction(int nbrParticles, Complex zElectron, double scale, int invFillingFactor);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHEDiskLaughlinOneJainQuasielectronWaveFunction(const FQHEDiskLaughlinOneJainQuasielectronWaveFunction& function);

  // destructor
  //
   ~FQHEDiskLaughlinOneJainQuasielectronWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

};

#endif
