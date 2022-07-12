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


#include "config.h"
#include "Operator/ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default data
//
// particle = hilbert space associated to the particles
// maxMomentum = number of flux quanta
// xMomentum= momentum along the x direction

ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator(ParticleOnTorusWithSpinAndMagneticTranslations* particle, int maxMomentum, int xMomentum)
{
  this->Particle= (ParticleOnTorusWithSpinAndMagneticTranslations*) (particle->Clone());
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->ExponentialFactors = new Complex[this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator(ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator& oper)
{
  this->Particle = (ParticleOnTorusWithSpinAndMagneticTranslations*) (oper.Particle->Clone());
  this->MaxMomentum = oper.MaxMomentum;
  this->XMomentum = oper.XMomentum;
  this->ExponentialFactors = new Complex[this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
}

// destructor
//

ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::~ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator()
{
  delete this->Particle;
  delete[] this->ExponentialFactors;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::Clone ()
{
  return new ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusWithSpinAndMagneticTranslationsSMinusOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int Index = 0;
  double Coefficient = 0.0;
  int TargetDim = this->Particle->GetTargetHilbertSpaceDimension();
  int NbrTranslations;
  for (int i = firstComponent; i < Last; ++i)
    {
      Complex& Tmp = vSource[i];
      for (int j = 0; j < this->MaxMomentum; ++j)
	{
	  Index = this->Particle->AddAu(i, j, Coefficient, NbrTranslations);
	  if (Index < TargetDim)
	    {
	      vDestination[Index] += Tmp * this->ExponentialFactors[NbrTranslations] * Coefficient;		  
	    }
	}
    }
  return vDestination;
}
  


