////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on sphere n-body operator               //
//                                                                            //
//                        last modification : 13/12/2004                      //
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
#include "Operator/ParticleOnSphereNBodyOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndices = array containing the indices of the creation operators from left to right (a+_0.....a+_(n-1))
// annihilationIndices = array containing the indices of the annihilation operators from left to right (a_0.....a_(n-1))
// nbrNBody = number of creation (or annihilation) operators

ParticleOnSphereNBodyOperator::ParticleOnSphereNBodyOperator(ParticleOnSphere* particle, int* creationIndices,
							     int* annihilationIndices, int nbrNBody)
{
  this->Particle = (ParticleOnSphere*) (particle->Clone());
  this->NbrNBody = nbrNBody;
  this->CreationIndices = new int [this->NbrNBody];  
  this->AnnihilationIndices = new int [this->NbrNBody];
  for (int i = 0; i < this->NbrNBody; ++i)
    {
      this->CreationIndices[i] = creationIndices[i];
      this->AnnihilationIndices[i] = annihilationIndices[i];
    }
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnSphereNBodyOperator::ParticleOnSphereNBodyOperator(const ParticleOnSphereNBodyOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->NbrNBody = oper.NbrNBody;
  this->CreationIndices = new int [this->NbrNBody];  
  this->AnnihilationIndices = new int [this->NbrNBody];
  for (int i = 0; i < this->NbrNBody; ++i)
    {
      this->CreationIndices[i] = oper.CreationIndices[i];
      this->AnnihilationIndices[i] = oper.AnnihilationIndices[i];
    }
}

// destructor
//

ParticleOnSphereNBodyOperator::~ParticleOnSphereNBodyOperator()
{
  delete this->Particle;
  delete[] this->CreationIndices;
  delete[] this->AnnihilationIndices;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereNBodyOperator::Clone ()
{
  return new ParticleOnSphereNBodyOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereNBodyOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereNBodyOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereNBodyOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereNBodyOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particle->ProdAdProdA(i, this->CreationIndices, this->AnnihilationIndices, this->NbrNBody, Coefficient);
      Element += V1[Index] * V2[i] * Coefficient;
    }
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereNBodyOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereNBodyOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							    int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = this->Particle->ProdAdProdA(i, this->CreationIndices, this->AnnihilationIndices, this->NbrNBody, Coefficient);
      vDestination[Index] = vSource[i] * Coefficient;
    }
  return vDestination;
}


