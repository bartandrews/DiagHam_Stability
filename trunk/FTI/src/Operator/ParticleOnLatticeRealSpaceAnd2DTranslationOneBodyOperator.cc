////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                           Class author: Cecile Repellin                    //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 10/04/2018                      //
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
#include "Operator/ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator.h"
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
ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator(ParticleOnSphere* particle,  int creationIndex, int annihilationIndex)
{
  this->Particle =  (ParticleOnSphere*) (particle->Clone()); 
  this->OperatorIndexDagger = creationIndex;
  this->OperatorIndex = annihilationIndex;
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator (const ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator & oper)
{
  this->Particle =  (ParticleOnSphere*) (oper.Particle->Clone()); 
  this->OperatorIndexDagger = oper.OperatorIndexDagger;
  this->OperatorIndex = oper.OperatorIndex;
}

// destructor
//

ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::~ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator()
{
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::Clone ()
{
  return new ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator (*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace * ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}


// change indices of creation / annihilation operators
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator

void ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::SetCreationAnnihilationIndex (int creationIndex, int annihilationIndex)
{
  this->OperatorIndexDagger=creationIndex;
  this->OperatorIndex=annihilationIndex;
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
  
Complex ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  
  double Coefficient = 0.0;
  Complex Element = 0.0;
  FermionOnLatticeRealSpaceAnd2DTranslation*  TmpParticle = (FermionOnLatticeRealSpaceAnd2DTranslation*) this->Particle->Clone();
  
  int CreationOrbitalIndex;
  int CreationXPosition;
  int CreationYPosition;
  TmpParticle->GetLinearizedIndex(this->OperatorIndexDagger, CreationXPosition, CreationYPosition, CreationOrbitalIndex);
  
  int AnnihilationOrbitalIndex;
  int AnnihilationXPosition;
  int AnnihilationYPosition;
  TmpParticle->GetLinearizedIndex(this->OperatorIndexDagger, AnnihilationXPosition, AnnihilationYPosition, AnnihilationOrbitalIndex);

  int TmpNbrUnitCells = TmpParticle->GetMaxXMomentum() * TmpParticle->GetMaxYMomentum();
  int* ShiftedCreationIndices = new int [TmpNbrUnitCells];
  int* ShiftedAnnihilationIndices = new int [TmpNbrUnitCells];
  Complex* RightPhases = new Complex [TmpNbrUnitCells];
  Complex** LeftPhases = new Complex* [TmpParticle->GetMaxXMomentum()];
  int Tmp = 0;  
  for (int i = 0; i < TmpParticle->GetMaxXMomentum(); ++i)
    {
      LeftPhases[i] = new Complex[TmpParticle->GetMaxYMomentum()];
      for (int j = 0; j < TmpParticle->GetMaxYMomentum(); ++j)
	{  
	  ShiftedCreationIndices[Tmp] = TmpParticle->GetLinearizedIndexSafe(CreationXPosition - i, CreationYPosition - j, CreationOrbitalIndex);
	  ShiftedAnnihilationIndices[Tmp] = TmpParticle->GetLinearizedIndexSafe(AnnihilationXPosition - i, AnnihilationYPosition - j, AnnihilationOrbitalIndex);
	  RightPhases[Tmp] = Phase (2.0 * M_PI * ((((double) (i * TmpParticle->GetKxMomentum())) / ((double) TmpParticle->GetMaxXMomentum()))
						  + (((double) (j * TmpParticle->GetKyMomentum())) / ((double) TmpParticle->GetMaxYMomentum())))) / ((double) TmpNbrUnitCells);
	  LeftPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * TmpParticle->GetKxMomentum())) / ((double) TmpParticle->GetMaxXMomentum())
						   + ((double) (j * TmpParticle->GetKyMomentum())) / ((double) TmpParticle->GetMaxYMomentum())));
	  ++Tmp;
	}
    }

  int NbrTranslationX;
  int NbrTranslationY;  
  
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
    for (int j = 0; j < TmpNbrUnitCells; ++j)
      {
	Index = TmpParticle->AdA(i, ShiftedCreationIndices[j], ShiftedAnnihilationIndices[j], Coefficient,NbrTranslationX, NbrTranslationY);
	if (Index < FullDim)
	{
	  NbrTranslationX -= j / TmpParticle->GetMaxYMomentum();
	  if (NbrTranslationX < 0)
	      NbrTranslationX += TmpParticle->GetMaxXMomentum();
	  NbrTranslationY -= j % TmpParticle->GetMaxYMomentum();
	  if (NbrTranslationY < 0)
	      NbrTranslationY += TmpParticle->GetMaxYMomentum();
	  Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;  
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

ComplexVector& ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  
  double Coefficient = 0.0;
  Complex Element = 0.0;
  FermionOnLatticeRealSpaceAnd2DTranslation*  TmpParticle = (FermionOnLatticeRealSpaceAnd2DTranslation*) this->Particle->Clone();

  int AnnihilationOrbitalIndex;
  int AnnihilationXPosition;
  int AnnihilationYPosition;
  TmpParticle->GetLinearizedIndex(this->OperatorIndex, AnnihilationXPosition, AnnihilationYPosition, AnnihilationOrbitalIndex);

  int CreationOrbitalIndex;
  int CreationXPosition;
  int CreationYPosition;
  TmpParticle->GetLinearizedIndex(this->OperatorIndexDagger, CreationXPosition, CreationYPosition, CreationOrbitalIndex);

  int TmpNbrUnitCells = TmpParticle->GetMaxXMomentum();
  int* ShiftedCreationIndices = new int [TmpNbrUnitCells];
  int* ShiftedAnnihilationIndices = new int [TmpNbrUnitCells];
  Complex* RightPhases = new Complex [TmpNbrUnitCells];
  Complex** LeftPhases = new Complex* [TmpParticle->GetMaxXMomentum()];
  int Tmp = 0;  
  for (int i = 0; i < TmpParticle->GetMaxXMomentum(); ++i)
    {
      LeftPhases[i] = new Complex[TmpParticle->GetMaxYMomentum()];
      for (int j = 0; j < TmpParticle->GetMaxYMomentum(); ++j)
	{
	  ShiftedCreationIndices[Tmp] = TmpParticle->GetLinearizedIndexSafe(CreationXPosition - i, CreationYPosition - j, CreationOrbitalIndex);
	  ShiftedAnnihilationIndices[Tmp] = TmpParticle->GetLinearizedIndexSafe(AnnihilationXPosition - i, AnnihilationYPosition - j, AnnihilationOrbitalIndex);
	  RightPhases[Tmp] = Phase (2.0 * M_PI * (((double) (i * TmpParticle->GetKxMomentum())) / ((double) TmpParticle->GetMaxXMomentum())
						   + ((double) (j * TmpParticle->GetKyMomentum())) / ((double) TmpParticle->GetMaxYMomentum()))) / ((double) TmpNbrUnitCells);
	  LeftPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * TmpParticle->GetKxMomentum())) / ((double) TmpParticle->GetMaxXMomentum())
						   + ((double) (j * TmpParticle->GetKyMomentum())) / ((double) TmpParticle->GetMaxYMomentum())));
	  ++Tmp;
	}
    }

  int NbrTranslationX;
  int NbrTranslationY;  
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
    for (int j = 0; j < TmpNbrUnitCells; ++j)
      {
      Index = TmpParticle->AdA(i, ShiftedCreationIndices[j], ShiftedAnnihilationIndices[j], Coefficient, NbrTranslationX, NbrTranslationY);
      if (Index != Dim)
	{
	  NbrTranslationX -= j / TmpParticle->GetMaxYMomentum();
	  if (NbrTranslationX < 0)
	    NbrTranslationX += TmpParticle->GetMaxXMomentum();
	  NbrTranslationY -= j % TmpParticle->GetMaxYMomentum();
	  if (NbrTranslationY < 0)
	    NbrTranslationY += TmpParticle->GetMaxYMomentum();
	  vDestination[Index] += vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * ((Complex) Coefficient);
	}
      }
    }
  delete TmpParticle;
  return vDestination;
}   
