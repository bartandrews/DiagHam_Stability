////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Halperin (m1, m2, n) wave function on disk             //
//              that evalute both the wave function and the wave              //
//                  function time the Coulomb interaction term                //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 23/10/2006                      //
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


#ifndef HALPERINONDISKWAVEFUNCTIONONEOVERR_H
#define HALPERINONDISKWAVEFUNCTIONONEOVERR_H


#include "config.h"
#include "Tools/FQHEWaveFunction/AbstractFQHEWaveFunctionOneOverR.h"


class HalperinOnDiskWaveFunctionOneOverR: public AbstractFQHEWaveFunctionOneOverR
{

 protected:

  // number of particles with spin up
  int NbrSpinUpParticles;
  // number of particles with spin down
  int NbrSpinDownParticles;
  // total number of particles
  int TotalNbrParticles;

  // m1 index ( (m1, m2, n) Halperin wave function)
  int M1Index;
  // m2 index ( (m1, m2, n) Halperin wave function)
  int M2Index;
  // n index ( (m1, m2, n) Halperin wave function)
  int NIndex;

  // number of Jastrow factors
  int NbrFactors;

  // array where the Jastrow factors that have been evaluated during the previous operator() method call have been stored
  Complex* JastrowFactors;
  
  // temporaray array used during 1/r |Psi(...)|^2 calculations
  double* TemporaryFactors;

 public:

  // constructor for the (m1, m2, n) Halperin wave function
  //
  // nbrSpinUpParticles = number of particles with spin up
  // nbrSpinDownParticles = number of particles with spin down
  // m1Index = m1 index
  // m2Index = m2 index
  // nIndex = n index
  HalperinOnDiskWaveFunctionOneOverR(int nbrSpinUpParticles, int nbrSpinDownParticles, int m1Index, int m2Index, int nIndex);

  // copy constructor
  //
  // function = reference on the wave function to copy
  HalperinOnDiskWaveFunctionOneOverR(const HalperinOnDiskWaveFunctionOneOverR& function);

  // destructor
  //
   ~HalperinOnDiskWaveFunctionOneOverR();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // evaluate the norm to the square of the wave function at a given point time the coulomb term (assume the coordinates are those provides by the previous operator() method call)
  //
  // return value = corresponding numerical value
  double CoulombContribution();

};

#endif
