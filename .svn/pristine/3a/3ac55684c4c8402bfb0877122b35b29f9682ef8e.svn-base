////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of particle on lattice translation operator            //
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
#include "Operator/ParticleOnLatticeTranslationOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// rx = x-component of the desired translation (generally: multiples of first lattice vector)
// ry = y-component of the desired translation (generally: multiples of second lattice vector)
ParticleOnLatticeTranslationOperator::ParticleOnLatticeTranslationOperator(ParticleOnLattice* particle, int rx, int ry)
{
  this->Particle = particle;
  this->Rx=rx;
  this->Ry=ry;  
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeTranslationOperator::ParticleOnLatticeTranslationOperator(const ParticleOnLatticeTranslationOperator& oper)
{
  this->Particle = oper.Particle;
  this->Rx=oper.Rx;
  this->Ry=oper.Ry;
}

// destructor
//

ParticleOnLatticeTranslationOperator::~ParticleOnLatticeTranslationOperator()
{
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeTranslationOperator::Clone ()
{
  return new ParticleOnLatticeTranslationOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeTranslationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeTranslationOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeTranslationOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}

// set components of translation vector
//
// rx = x-component of the desired translation
// ry = y-component of the desired translation
void ParticleOnLatticeTranslationOperator::SetTranslationComponents(int rx, int ry)
{
  this->Rx=rx;
  this->Ry=ry;  
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnLatticeTranslationOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = firstComponent + nbrComponent;
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  Complex TranslationPhase;  
  Complex Element = 0.0;
  int Index;
  for (int i = firstComponent; i < Dim; ++i)
    {
      Index = TmpParticle->TranslateState(i, this->Rx, this->Ry, TranslationPhase);      
      Element += V1[Index] * V2[i] * TranslationPhase;
    }
  delete TmpParticle;
  return Element;
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
  
Complex ParticleOnLatticeTranslationOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  Complex TranslationPhase;
  Complex Element = 0.0;
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      Index = TmpParticle->TranslateState(i, this->Rx, this->Ry, TranslationPhase);
      Element += Conj(V1[Index]) * V2[i] * TranslationPhase;
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

ComplexVector& ParticleOnLatticeTranslationOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Last = (int) (firstComponent + nbrComponent);
  int Index;
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  Complex TranslationPhase;
  //std::cout << "R=("<<Rx<<", "<<Ry<<")"<<std::endl<<vSource<<std::endl<<"Norm:"<<vSource.Norm()<<std::endl;
  //std::cout.precision(8);
  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Index = TmpParticle->TranslateState(i, this->Rx, this->Ry, TranslationPhase);      
      //std::cout << "Translated("<<i<<")="<<TranslationPhase<<"* ["<<Index<<"]"<<std::endl;
      
      //std::cout << "Source("<<i<<")="<<Norm(vSource[i])<<", Target("<<Index<<")="<<Norm(vSource[Index])<<" d="<<(Arg(vSourc e[i]/vSource[Index]))/M_PI<<"pi, TranslationPhase="<<Arg(TranslationPhase)/M_PI<<"pi"<<std::endl;
      
      //std::cout << "Ignoring translation Phase!"<<std::endl; TranslationPhase=1.0;
      // permuted i and Index on the following lines!
      vDestination.Re(i) += vSource[Index].Re * TranslationPhase.Re - vSource[Index].Im * TranslationPhase.Im;
      vDestination.Im(i) += vSource[Index].Re * TranslationPhase.Im + vSource[Index].Im * TranslationPhase.Re;
    }
  //std::cout << vDestination<<std::endl<<"Norm:"<<vDestination.Norm()<<std::endl;
  delete TmpParticle;
  return vDestination;
}


