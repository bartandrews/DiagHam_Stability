////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class density-density operator for particle with spin          //
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
#include "Operator/ParticleOnSphereWithSpinDensityDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for up, 1 for plus)
// creationMomentumIndex2 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for up, 1 for plus)
// annihilationMomentumIndex1 = momentum index of the leftmost annihilation operator (from 0 to 2S)
// annihilationSymmetryIndex1 = symmetry index of the leftmost annihilation operator (0 for up, 1 for plus)
// annihilationMomentumIndex2 = momentum index of the rightmost annihilation operator(from 0 to 2S)
// annihilationSymmetryIndex2 = symmetry index of the rightmost annihilation operator (0 for up, 1 for plus)

ParticleOnSphereWithSpinDensityDensityOperator::ParticleOnSphereWithSpinDensityDensityOperator(ParticleOnSphereWithSpin* particle, 
											       int creationMomentumIndex1, int creationSymmetryIndex1,
											       int creationMomentumIndex2, int creationSymmetryIndex2,
											       int annihilationMomentumIndex1, int annihilationSymmetryIndex1,
											       int annihilationMomentumIndex2, int annihilationSymmetryIndex2)
{
  this->Particle= particle;
  this->SignFactor = 1.0;
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->CreationSymmetryIndex1 = creationSymmetryIndex1;
  this->CreationMomentumIndex2 = creationMomentumIndex2;
  this->CreationSymmetryIndex2 = creationSymmetryIndex2;
  if (this->CreationSymmetryIndex1 > this->CreationSymmetryIndex2)
    {
      int Tmp = this->CreationSymmetryIndex1;
      this->CreationSymmetryIndex1 = this->CreationSymmetryIndex2;
      this->CreationSymmetryIndex2 = Tmp;
      Tmp = this->CreationMomentumIndex1;
      this->CreationMomentumIndex1 = this->CreationMomentumIndex2;
      this->CreationMomentumIndex2 = Tmp;   
      if (this->Particle->GetParticleStatistic() == AbstractQHEParticle::FermionicStatistic)
	this->SignFactor *= 1.0;
    }
  this->AnnihilationMomentumIndex1 = annihilationMomentumIndex1;
  this->AnnihilationSymmetryIndex1 = annihilationSymmetryIndex1;
  this->AnnihilationMomentumIndex2 = annihilationMomentumIndex2;
  this->AnnihilationSymmetryIndex2 = annihilationSymmetryIndex2;
  if (this->AnnihilationSymmetryIndex1 > this->AnnihilationSymmetryIndex2)
    {
      int Tmp = this->AnnihilationSymmetryIndex1;
      this->AnnihilationSymmetryIndex1 = this->AnnihilationSymmetryIndex2;
      this->AnnihilationSymmetryIndex2 = Tmp;
      Tmp = this->AnnihilationMomentumIndex1;
      this->AnnihilationMomentumIndex1 = this->AnnihilationMomentumIndex2;
      this->AnnihilationMomentumIndex2 = Tmp;   
      if (this->Particle->GetParticleStatistic() == AbstractQHEParticle::FermionicStatistic)
	this->SignFactor *= 1.0;
    }
}

// destructor
//

ParticleOnSphereWithSpinDensityDensityOperator::~ParticleOnSphereWithSpinDensityDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinDensityDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinDensityDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinDensityDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinDensityDensityOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnSphereWithSpinDensityDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  double Element = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 12) | (this->CreationSymmetryIndex2 << 8) | (this->AnnihilationSymmetryIndex1 << 4) | this->AnnihilationSymmetryIndex2;
  switch (SymmetryIndex)
    {
    case 0x0000 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    Coefficient = this->Particle->AuAu(i, this->AnnihilationMomentumIndex1, this->AnnihilationMomentumIndex2);
	    if (Coefficient != 0.0)
	      {
		int Index = this->Particle->AduAdu(this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient2);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * Coefficient * Coefficient2;      
	      }
	  }
	break;
      }
    case 0x0101 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    Coefficient = this->Particle->AuAd(i, this->AnnihilationMomentumIndex1, this->AnnihilationMomentumIndex2);
	    if (Coefficient != 0.0)
	      {
		int Index = this->Particle->AduAdd(this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient2);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * Coefficient * Coefficient2;      
	      }
	  }
	break;
      }
    case 0x1111 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    Coefficient = this->Particle->AdAd(i, this->AnnihilationMomentumIndex1, this->AnnihilationMomentumIndex2);
	    if (Coefficient != 0.0)
	      {
		int Index = this->Particle->AddAdd(this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient2);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * Coefficient * Coefficient2;      
	      }
	  }
	break;
      }
    }
  return Complex(this->SignFactor * Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinDensityDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 12) | (this->CreationSymmetryIndex2 << 8) | (this->AnnihilationSymmetryIndex1 << 4) | this->AnnihilationSymmetryIndex2;
  switch (SymmetryIndex)
    {
    case 0x0000 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    Coefficient = this->Particle->AuAu(i, this->AnnihilationMomentumIndex1, this->AnnihilationMomentumIndex2);
	    if (Coefficient != 0.0)
	      {
		int Index = this->Particle->AduAdu(this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient2);
		if (Index != Dim)
		  vDestination[Index] += this->SignFactor * vSource[i] * Coefficient * Coefficient2;
	      }
	  }
	break;
      }
    case 0x0101 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    Coefficient = this->Particle->AuAd(i, this->AnnihilationMomentumIndex1, this->AnnihilationMomentumIndex2);
	    if (Coefficient != 0.0)
	      {
		int Index = this->Particle->AduAdd(this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient2);
		if (Index != Dim)
		  vDestination[Index] += this->SignFactor * vSource[i] * Coefficient * Coefficient2;
	      }
	  }
	break;
      }
    case 0x1111 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    Coefficient = this->Particle->AdAd(i, this->AnnihilationMomentumIndex1, this->AnnihilationMomentumIndex2);
	    if (Coefficient != 0.0)
	      {
		int Index = this->Particle->AddAdd(this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient2);
		if (Index != Dim)
		  vDestination[Index] += this->SignFactor * vSource[i] * Coefficient * Coefficient2;
	      }
	  }
	break;
      }
    }
  return vDestination;
}

