////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of S^- spin operator on the torus with magnetic translations    //
//                                                                            //
//                        last modification : 03/02/2015                      //
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


#ifndef PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSSMINUSOPERATOR_H
#define PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONSSMINUSOPERATOR_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"

#include <iostream>


class MathematicaOutput;


class ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the particles
  ParticleOnTorusWithSpinAndMagneticTranslations* Particle;
  
  // maximum momentum value reached by a particle in the state
  int MaxMomentum;
  // momentum value in the x direction (modulo GCD of NbParticles and MaxMomentum)
  int XMomentum;

  //array containing all the phase factors that are needed when computing matrix elements
  Complex* ExponentialFactors;

 public:
  
  // constructor from default data
  //
  // particle = hilbert space associated to the particles
  // maxMomentum = number of flux quanta
  // xMomentum= momentum along the x direction
  ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator(ParticleOnTorusWithSpinAndMagneticTranslations* particle, int maxMomentum, int xMomentum);

  // copy constructor
  //
  // oper = operator to copy
  ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator(ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator& oper);

  // destructor
  //
  ~ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator();
  
  // clone operator without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractOperator* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which operator acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where operator acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);
  
};

#endif
