////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of particle on lattice density-density operator          //
//                                                                            //
//                        last modification : 10/12/2002                      //
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
#include "Operator/ParticleOnLatticeDensityDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/BosonOnLattice.h"

#include <iostream>


using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// siteIndex1 = index of first site to be considered
// siteIndex2 = index of second site to be considered
ParticleOnLatticeDensityDensityOperator::ParticleOnLatticeDensityDensityOperator(ParticleOnLattice* particle, int siteIndex1, int siteIndex2)
{
  this->Particle = particle;
  this->CreationIndex1 = siteIndex1;
  this->CreationIndex2 = siteIndex2;
  this->AnnihilationIndex1 = siteIndex2;
  this->AnnihilationIndex2 = siteIndex1;
}



// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex1 = index of the leftmost creation operator
// creationIndex2 = index of the rightmost creation operator
// annihilationIndex1 = index of the leftmost annihilation operator
// annihilationIndex2 = index of the rightmost annihilation operator

ParticleOnLatticeDensityDensityOperator::ParticleOnLatticeDensityDensityOperator(ParticleOnLattice* particle,
										 int creationIndex1, int creationIndex2,
										 int annihilationIndex1, int annihilationIndex2)
{
  this->Particle = particle;
  this->CreationIndex1 = creationIndex1;
  this->CreationIndex2 = creationIndex2;
  this->AnnihilationIndex1 = annihilationIndex1;
  this->AnnihilationIndex2 = annihilationIndex2;
}

ParticleOnLatticeDensityDensityOperator::ParticleOnLatticeDensityDensityOperator(const ParticleOnLatticeDensityDensityOperator& oper)
{
  this->Particle = oper.Particle;
  this->CreationIndex1 = oper.CreationIndex1;
  this->CreationIndex2 = oper.CreationIndex2;
  this->AnnihilationIndex1 = oper.AnnihilationIndex1;
  this->AnnihilationIndex2 = oper.AnnihilationIndex2;
}


// destructor
//

ParticleOnLatticeDensityDensityOperator::~ParticleOnLatticeDensityDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeDensityDensityOperator::Clone ()
{
  return new ParticleOnLatticeDensityDensityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeDensityDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeDensityDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeDensityDensityOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnLatticeDensityDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  double Coefficient = 0.0;
  double Element = 0.0;
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      Index = TmpParticle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      if (Index != FullDim)
	Element += V1[Index] * V2[i] * Coefficient;      
    }
  delete TmpParticle;
  return Complex(Element);
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnLatticeDensityDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      Index = TmpParticle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      if (Index != FullDim)
	Element += Conj(V1[Index]) * V2[i] * Coefficient;  
    }
  delete TmpParticle;
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
ComplexVector& ParticleOnLatticeDensityDensityOperator::LowLevelAddMultiply(ComplexVector& vSource,
									    ComplexVector& vDestination, 
									    int firstComponent, int nbrComponent)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  int Last = firstComponent + nbrComponent;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = TmpParticle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      if (Index != Dim)
	vDestination[Index] += vSource[i] * Coefficient;
    }
  delete TmpParticle;
  return vDestination;
}
  

