////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of particle on sphere square of total isospin operator         //
//                                                                            //
//                        last modification : 01/01/2007                      //
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
#include "Operator/ParticleOnSphereSquareTotalIsospinOperator.h"
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
// nbrParticles = number of particles
// totalIz = twice the total isospin projection (i.e Iz) value

ParticleOnSphereSquareTotalIsospinOperator::ParticleOnSphereSquareTotalIsospinOperator(ParticleOnSphereWithSU4Spin* particle, int lzMax, int nbrParticles, int totalIz)
{
  this->Particle = (ParticleOnSphereWithSU4Spin*) (particle->Clone());
  this->LzMax = lzMax;
  this->TotalIz = totalIz;
  this->NbrParticles = nbrParticles;
  this->Shift = (0.5 * ((double) this->NbrParticles)) + (0.25 * ((double) (this->TotalIz * this->TotalIz)));
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnSphereSquareTotalIsospinOperator::ParticleOnSphereSquareTotalIsospinOperator(const ParticleOnSphereSquareTotalIsospinOperator& oper)
{
  this->Particle = (ParticleOnSphereWithSU4Spin*) (oper.Particle->Clone());
  this->LzMax = oper.LzMax;
  this->Shift = oper.Shift;
  this->NbrParticles = oper.NbrParticles;
  this->TotalIz = oper.TotalIz;
}

// destructor
//

ParticleOnSphereSquareTotalIsospinOperator::~ParticleOnSphereSquareTotalIsospinOperator()
{
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereSquareTotalIsospinOperator::Clone ()
{
  return new ParticleOnSphereSquareTotalIsospinOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereSquareTotalIsospinOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSU4Spin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereSquareTotalIsospinOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereSquareTotalIsospinOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSquareTotalIsospinOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Element = 0.0;
  int Index = 0;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  for (int i = 0; i < Dim; ++i)
    {
      for (int j = 0; j <= this->LzMax; ++j)
	for (int  k = 0; k <= this->LzMax; ++k)
	  {
	    Coefficient = this->Particle->AupAum(i, k, j);
	    if (Coefficient != 0.0)
	      {
		Index = this->Particle->AdupAdum(j, k, Coefficient2);
		if (Index != Dim)
		  Element += V1[Index] * V2[i] * Coefficient * Coefficient2;		  
	      }
	    Coefficient = this->Particle->AdpAdm(i, k, j);
	    if (Coefficient != 0.0)
	      {
		Index = this->Particle->AddpAddm(j, k, Coefficient2);
		if (Index != Dim)
		  Element += V1[Index] * V2[i] * Coefficient * Coefficient2;		  
	      }
	    Coefficient = this->Particle->AumAdp(i, k, j);
	    if (Coefficient != 0.0)
	      {
		Index = this->Particle->AdupAddm(k, j, Coefficient2);
		if (Index != Dim)
		  Element -= V1[Index] * V2[i] * Coefficient * Coefficient2;		  
	      }
	    Coefficient = this->Particle->AupAdm(i, k, j);
	    if (Coefficient != 0.0)
	      {
		Index = this->Particle->AdumAddp(k, j, Coefficient2);
		if (Index != Dim)
		  Element -= V1[Index] * V2[i] * Coefficient * Coefficient2;		
	      }
	  }
    }
  Element += this->Shift; 
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSquareTotalIsospinOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator for a given range of indices 
// and add result to another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereSquareTotalIsospinOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  return vDestination;
}
  


