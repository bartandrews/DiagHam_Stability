////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of particle on sphere sqare of total momentum operator         //
//                                                                            //
//                        last modification : 06/03/2005                      //
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
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


#include <iostream>


using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// lzMax = maximum Lz value reached by a fermion
// factor = multiplicative factor in front of the L^2 operator
// memory = amount of memory (in bytes) that can be used for precalculations (none if memory < 0)

ParticleOnSphereSquareTotalMomentumOperator::ParticleOnSphereSquareTotalMomentumOperator(ParticleOnSphere* particle, int lzMax, double factor, long memory)
{
  this->Particle = (ParticleOnSphere*) (particle->Clone());
  this->LzMax = lzMax;
  RealMatrix Coefficients (this->LzMax + 1, this->LzMax + 1);
  this->TotalLz = 0;
  for (int k = 0; k <= this->LzMax; ++k)
    this->TotalLz += ((2 * k) - this->LzMax) * (int) this->Particle->AdA(0, k);
  this->Shift = 0.25 * factor * ((double) (this->TotalLz * this->TotalLz));
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax + 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	{
     	  Coefficients(i, j) = TmpCoefficient;
	}
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax - 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	{
     	  Coefficients(j, i) *= 0.5 * TmpCoefficient;
	}
    }

  double Coefficient = 2.0;
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Coefficient = -2.0;
  this->TwoBodyCoefficients = new double [this->LzMax * this->LzMax];
  double* TmpTwoBodyCoefficients = this->TwoBodyCoefficients;
  for (int j = 0; j < this->LzMax; ++j)
    for (int k = 1; k <= this->LzMax; ++k)
      {
	(*TmpTwoBodyCoefficients) = factor * Coefficient * Coefficients(j, k);
	++TmpTwoBodyCoefficients;
      }

  this->OneBodyCoefficients = new double [this->LzMax + 1];
  this->OneBodyCoefficients[0] = factor * Coefficients(0, 1);
  for (int k = 1; k < this->LzMax; ++k)
    this->OneBodyCoefficients[k] = factor * (Coefficients(k, k + 1) + Coefficients(k - 1, k));
  this->OneBodyCoefficients[this->LzMax] = factor * Coefficients(this->LzMax - 1, this->LzMax);
  
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnSphereSquareTotalMomentumOperator::ParticleOnSphereSquareTotalMomentumOperator(const ParticleOnSphereSquareTotalMomentumOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->LzMax = oper.LzMax;
  this->Shift = oper.Shift;
  int Lim = this->LzMax * this->LzMax;
  this->TwoBodyCoefficients = new double [Lim];
  this->OneBodyCoefficients = new double [this->LzMax + 1];
  for (int k = 0; k < Lim; ++k)
    this->TwoBodyCoefficients[k] = oper.TwoBodyCoefficients[k];
  for (int k = 0; k <= this->LzMax; ++k)
    this->OneBodyCoefficients[k] = oper.OneBodyCoefficients[k];  
}

// destructor
//

ParticleOnSphereSquareTotalMomentumOperator::~ParticleOnSphereSquareTotalMomentumOperator()
{
  delete this->Particle;
  delete[] this->TwoBodyCoefficients;
  delete[] this->OneBodyCoefficients;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereSquareTotalMomentumOperator::Clone ()
{
  return new ParticleOnSphereSquareTotalMomentumOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereSquareTotalMomentumOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereSquareTotalMomentumOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereSquareTotalMomentumOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnSphereSquareTotalMomentumOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  double Element = 0.0;
  int Index = 0;
  double Coefficient = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      double* TmpTwoBodyCoefficients = this->TwoBodyCoefficients;
      for (int j = 0; j < this->LzMax; ++j)
	for (int  k = 1; k <= this->LzMax; ++k)
	  {
	    Index = this->Particle->AdAdAASafe(i, k - 1, j + 1, k, j, Coefficient);
	    if (Index != this->Particle->GetHilbertSpaceDimension())
	      {
		Element += V1[Index] * V2[i] * Coefficient * (*TmpTwoBodyCoefficients);		  
	      }
	    ++TmpTwoBodyCoefficients;
	  }
      Coefficient = 0.0;
      for (int k = 0; k <= this->LzMax; ++k)
	Coefficient += this->OneBodyCoefficients[k] * this->Particle->AdA(i, k);
      Element += V1[i] * V2[i] * (Coefficient + this->Shift); 
    }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and add result to another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereSquareTotalMomentumOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Index = 0;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      double* TmpTwoBodyCoefficients = this->TwoBodyCoefficients;
      for (int j = 0; j < this->LzMax; ++j)
	for (int k = 1; k <= this->LzMax; ++k)
	  {
	    Index = this->Particle->AdAdAA(i, k - 1, j + 1, k, j, Coefficient);
	    if (Index != this->Particle->GetHilbertSpaceDimension())
	      {
		vDestination[Index] += vSource[i] * Coefficient * (*TmpTwoBodyCoefficients);		  
	      }
	    ++TmpTwoBodyCoefficients;
	  }
      Coefficient = 0.0;
      for (int k = 0; k <= this->LzMax; ++k)
	Coefficient += this->OneBodyCoefficients[k] * this->Particle->AdA(i, k);
      vDestination[i] += vSource[i] * (Coefficient + this->Shift);
    }
  return vDestination;
}
  


