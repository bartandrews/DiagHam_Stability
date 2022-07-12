////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 09/04/2008                      //
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
#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include <iostream>

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator
//
ParticleOnLatticeOneBodyOperator::ParticleOnLatticeOneBodyOperator(ParticleOnLattice* particle,
								   int creationIndex, int annihilationIndex)
{
  this->Particle = (ParticleOnLattice*) (particle->Clone());
  this->CreationIndex = creationIndex;
  this->AnnihilationIndex = annihilationIndex;
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeOneBodyOperator::ParticleOnLatticeOneBodyOperator(const ParticleOnLatticeOneBodyOperator& oper)
{
  this->Particle = (ParticleOnLattice*) (oper.Particle->Clone());
  this->CreationIndex = oper.CreationIndex;
  this->AnnihilationIndex = oper.AnnihilationIndex;
}

// destructor
//

ParticleOnLatticeOneBodyOperator::~ParticleOnLatticeOneBodyOperator()
{
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeOneBodyOperator::Clone ()
{
  return new ParticleOnLatticeOneBodyOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeOneBodyOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeOneBodyOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeOneBodyOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}


// change indices of creation / annihilation operators
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator
void ParticleOnLatticeOneBodyOperator::SetCreationAnnihilationIndex (int creationIndex, int annihilationIndex)
{
  this->CreationIndex=creationIndex;
  this->AnnihilationIndex=annihilationIndex;
}

  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnLatticeOneBodyOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particle->AdA(i, this->CreationIndex, this->AnnihilationIndex, Coefficient);
      if ((Index<Dim)&&(Coefficient!=0.0))
	Element += V1[Index] * V2[i] * Coefficient;
    }
  return Complex(Element);
}


Complex ParticleOnLatticeOneBodyOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particle->AdA(i, this->CreationIndex, this->AnnihilationIndex, Coefficient);
      if ((Index<Dim)&&(Coefficient!=0.0))
	Element += (Conj(V1[Index]) * V2[i] * Coefficient);
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

ComplexVector& ParticleOnLatticeOneBodyOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = this->Particle->AdA(i, this->CreationIndex, this->AnnihilationIndex, Coefficient);
      if ((Index<Dim)&&(Coefficient!=0.0))
	{
	  vDestination[Index].Re = vSource[i].Re * Coefficient;
	  vDestination[Index].Im = vSource[i].Im * Coefficient;
	}
    }
  return vDestination;
}


