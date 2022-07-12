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
#include "Operator/ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include <iostream>

using std::cout;
using std::endl;
  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator
//
ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator(FermionOnLatticeRealSpaceAnd1DTranslation * particle,  int creationIndex, int annihilationIndex)
{
  this->Particle =  (FermionOnLatticeRealSpaceAnd1DTranslation*) (particle->Clone()); 
  this->CreationIndex = creationIndex;
  this->AnnihilationIndex = annihilationIndex;
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator (const ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator & oper)
{
  this->Particle =  (FermionOnLatticeRealSpaceAnd1DTranslation*) (oper.Particle->Clone()); 
  this->CreationIndex = oper.CreationIndex;
  this->AnnihilationIndex = oper.AnnihilationIndex;
}

// destructor
//

ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::~ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator()
{
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::Clone ()
{
  return new ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator (*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (FermionOnLatticeRealSpaceAnd1DTranslation*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace * ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}


// change indices of creation / annihilation operators
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator

void ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::SetCreationAnnihilationIndex (int creationIndex, int annihilationIndex)
{
  this->CreationIndex=creationIndex;
  this->AnnihilationIndex=annihilationIndex;
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
  
Complex ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  int NbrTranslationX=0;
  double Coefficient = 0.0;
  Complex Element = 0.0;
  FermionOnLatticeRealSpaceAnd1DTranslation*  TmpParticle = (FermionOnLatticeRealSpaceAnd1DTranslation*) this->Particle->Clone();

  int AnnihilationOrbital;
  int AnnihilationXPosition;
  TmpParticle->GetLinearizedIndex(this->AnnihilationIndex, AnnihilationXPosition, AnnihilationOrbital);

  int CreationOrbitalIndex;
  int CreationXPosition;
  TmpParticle->GetLinearizedIndex(this->CreationIndex, CreationXPosition,CreationOrbitalIndex);

  int TmpNbrUnitCells = TmpParticle->GetMaxXMomentum();
  int* ShiftedCreationIndices = new int [TmpNbrUnitCells];
  int* ShiftedAnnihilationIndices = new int [TmpNbrUnitCells];
  Complex* RightPhases = new Complex [TmpNbrUnitCells];
  int Tmp = 0;  
  for (int i = 0; i < TmpParticle->GetMaxXMomentum(); ++i)
    {
      ShiftedCreationIndices[i] = TmpParticle->GetLinearizedIndexSafe(CreationXPosition - i, CreationOrbitalIndex);
      ShiftedAnnihilationIndices[i] = TmpParticle->GetLinearizedIndexSafe(AnnihilationXPosition - i, AnnihilationOrbital);
      RightPhases[i] = Phase (2.0 * M_PI * ((((double) (i * TmpParticle->GetKxMomentum())) / ((double) TmpParticle->GetMaxXMomentum()))));
    }

  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
    for (int j = 0; j < TmpNbrUnitCells; ++j)
      {
      Index = TmpParticle->AdA(i, ShiftedCreationIndices[j], ShiftedAnnihilationIndices[j], Coefficient,NbrTranslationX);
      if (Index < FullDim)
      {
        NbrTranslationX -= j;
        if (NbrTranslationX < 0)
           NbrTranslationX += TmpParticle->GetMaxXMomentum();
	Element += Conj(V1[Index]) * V2[i] * ((Complex) Coefficient) * RightPhases[j] * RightPhases[NbrTranslationX] / ((double) TmpNbrUnitCells);
      }
 	}
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

ComplexVector& ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  int NbrTranslationX = 0;
  double Coefficient = 0.0;
  Complex Element = 0.0;
  FermionOnLatticeRealSpaceAnd1DTranslation*  TmpParticle = (FermionOnLatticeRealSpaceAnd1DTranslation*) this->Particle->Clone();

  int AnnihilationOrbital;
  int AnnihilationXPosition;
  TmpParticle->GetLinearizedIndex(this->AnnihilationIndex, AnnihilationXPosition, AnnihilationOrbital);

  int CreationOrbitalIndex;
  int CreationXPosition;
  TmpParticle->GetLinearizedIndex(this->CreationIndex, CreationXPosition,CreationOrbitalIndex);

  int TmpNbrUnitCells = TmpParticle->GetMaxXMomentum();
  int* ShiftedCreationIndices = new int [TmpNbrUnitCells];
  int* ShiftedAnnihilationIndices = new int [TmpNbrUnitCells];
  Complex* RightPhases = new Complex [TmpNbrUnitCells];
  int Tmp = 0;  
  for (int i = 0; i < TmpParticle->GetMaxXMomentum(); ++i)
    {
      ShiftedCreationIndices[i] = TmpParticle->GetLinearizedIndexSafe(CreationXPosition - i, CreationOrbitalIndex);
      ShiftedAnnihilationIndices[i] = TmpParticle->GetLinearizedIndexSafe(AnnihilationXPosition - i, AnnihilationOrbital);
      RightPhases[i] = Phase (2.0 * M_PI * ((((double) (i * TmpParticle->GetKxMomentum())) / ((double) TmpParticle->GetMaxXMomentum()))));
    }

  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
    for (int j = 0; j < TmpNbrUnitCells; ++j)
      {
      Index = TmpParticle->AdA(i, ShiftedCreationIndices[j], ShiftedAnnihilationIndices[j], Coefficient,NbrTranslationX);
      if (Index < FullDim)
      {
        NbrTranslationX -= j;
        if (NbrTranslationX < 0)
           NbrTranslationX += TmpParticle->GetMaxXMomentum();
        vDestination[Index] += vSource[i] * ((Complex) Coefficient) * RightPhases[j] * RightPhases[NbrTranslationX] / ((double) TmpNbrUnitCells);;      
      }
 	}
}
  delete TmpParticle;
  return vDestination;
}   
