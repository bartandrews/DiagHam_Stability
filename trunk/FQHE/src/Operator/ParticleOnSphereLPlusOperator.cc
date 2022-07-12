////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of particle on sphere lowering momentum L operator           //
//                                                                            //
//                        last modification : 06/06/2006                      //
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
#include "Operator/ParticleOnSphereLPlusOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles with totalLz  Lz momentum, target space has to be fixed to hilbert space totalLz - 1 Lz momentum
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a fermion

ParticleOnSphereLPlusOperator::ParticleOnSphereLPlusOperator(ParticleOnSphere* particle, int totalLz, int lzMax)
{
  this->Particle = (ParticleOnSphere*) (particle->Clone());
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->Coefficients = RealVector(this->LzMax + 1);
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->Coefficients[i] = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax + 2))))));
    }
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnSphereLPlusOperator::ParticleOnSphereLPlusOperator(const ParticleOnSphereLPlusOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->TotalLz = oper.TotalLz;
  this->LzMax = oper.LzMax;
  this->Coefficients = oper.Coefficients;
}

// destructor
//

ParticleOnSphereLPlusOperator::~ParticleOnSphereLPlusOperator()
{
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereLPlusOperator::Clone ()
{
  return new ParticleOnSphereLPlusOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereLPlusOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereLPlusOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereLPlusOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnSphereLPlusOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = firstComponent + nbrComponent;
  int TargetDim = this->Particle->GetTargetHilbertSpaceDimension();
  double Element = 0.0;
  int Index = 0;
  double Coefficient = 0.0;
  double Tmp;
  for (int i = firstComponent; i < Dim; ++i)
    {
      Tmp = V2[i];
      for (int j = 0; j < this->LzMax; ++j)
	{
	  Index = this->Particle->AdA(i, j + 1, j, Coefficient);
	  if (Index < TargetDim)
	    {
	      Element += V1[Index] * Tmp * Coefficient * this->Coefficients[j];		  
	    }
	}
    }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereLPlusOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int Index = 0;
  double Coefficient = 0.0;
  int TargetDim = this->Particle->GetTargetHilbertSpaceDimension();
  double Tmp;
  for (int i = firstComponent; i < Last; ++i)
    {
      Tmp = vSource[i];
      for (int j = 0; j < this->LzMax; ++j)
	{
	  Index = this->Particle->AdA(i, j + 1, j, Coefficient);
	  if (Index < TargetDim)
	    {
	      vDestination[Index] += Tmp * Coefficient * this->Coefficients[j];		  
	    }
	}
    }
  return vDestination;
}
  
