////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on torus Kx momentum operator            //
//                                                                            //
//                        last modification : 09/02/2012                      //
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
#include "Operator/ParticleOnTorusKxOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default data
//
// particle = hilbert space associated to the particles

ParticleOnTorusKxOperator::ParticleOnTorusKxOperator(ParticleOnTorus* particle)
{
  this->Particle= (ParticleOnTorus*) (particle->Clone());
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnTorusKxOperator::ParticleOnTorusKxOperator(ParticleOnTorusKxOperator& oper)
{
  this->Particle = (ParticleOnTorus*) (oper.Particle->Clone());
}

// destructor
//

ParticleOnTorusKxOperator::~ParticleOnTorusKxOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusKxOperator::Clone ()
{
  return new ParticleOnTorusKxOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusKxOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorus*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusKxOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusKxOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnTorusKxOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  double Element = 0.0;
  int Last = firstComponent + nbrComponent;
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int i = firstComponent; i < Last; ++i)
	{
	  Element += V1[this->Particle->ApplyXMagneticTranslation(i)] * V2[i];
	}
    }
  else
    {
      for (int i = firstComponent; i < Last; ++i)
	{
	  Element += V1[this->Particle->ApplyXMagneticTranslation(i)] * V2[i];
	}
    }
  return Complex(Element);
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnTorusKxOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  Complex Element = 0.0;
  int Last = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Last; ++i)
    {
      Element += Conj(V1[this->Particle->ApplyXMagneticTranslation(i)]) * V2[i];
    }
  return Element;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnTorusKxOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      for (int i = firstComponent; i < Last; ++i)
	{
	  vDestination[this->Particle->ApplyXMagneticTranslation(i)] += vSource[i];
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      for (long i = firstComponent; i < Last; ++i)
	{
	  vDestination[this->Particle->ApplyXMagneticTranslation(i)] += vSource[i];
	}
    }
  return vDestination;
}
  


// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusKxOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      double TmpSign;
      if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
	{
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      int TmpIndex = this->Particle->ApplyXMagneticTranslation(i, TmpSign);
	      vDestination[TmpIndex] += TmpSign * vSource[i];
	    }
	}
      else
	{
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      long TmpIndex = this->Particle->ApplyXMagneticTranslation(i, TmpSign);
	      vDestination[TmpIndex] += TmpSign * vSource[i];
	    }
	}
    }
  else
    {
      if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
	{
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      vDestination[this->Particle->ApplyXMagneticTranslation(i)] += vSource[i];
	    }
	}
      else
	{
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      vDestination[this->Particle->ApplyXMagneticTranslation(i)] += vSource[i];
	    }
	}
    }
  return vDestination;
}
  


