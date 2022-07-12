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
#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationMomentumIndex = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex = symmetry index of the leftmost creation operator (0 for up, 1 for plus)
// annihilationMomentumIndex = momentum index of the leftmost annihilation operator (from 0 to 2S)
// annihilationSymmetryIndex = symmetry index of the leftmost annihilation operator (0 for up, 1 for plus)

ParticleOnSphereWithSpinDensityOperator::ParticleOnSphereWithSpinDensityOperator(ParticleOnSphereWithSpin* particle, 
										 int creationMomentumIndex, int creationSymmetryIndex,
										 int annihilationMomentumIndex, int annihilationSymmetryIndex)
{
  this->Particle= particle;
  this->CreationMomentumIndex = creationMomentumIndex;
  this->CreationSymmetryIndex = creationSymmetryIndex;
  this->AnnihilationMomentumIndex = annihilationMomentumIndex;
  this->AnnihilationSymmetryIndex = annihilationSymmetryIndex;
}

// destructor
//

ParticleOnSphereWithSpinDensityOperator::~ParticleOnSphereWithSpinDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinDensityOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex << 1) | this->AnnihilationSymmetryIndex;
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = 0; i < Dim; ++i)
	  {
	    int Index = this->Particle->AddAd(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
    case 3 :
      {
	for (int i = 0; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAu(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
    case 1 :
      {
	for (int i = 0; i < Dim; ++i)
	  {
	    int Index = this->Particle->AddAu(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
    case 2 :
      {
	for (int i = 0; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAd(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
   }
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinDensityOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
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

RealVector& ParticleOnSphereWithSpinDensityOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex << 1) | this->AnnihilationSymmetryIndex;
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AddAd(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] = vSource[i] * Coefficient;      
	  }
      }
      break;
    case 3 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAu(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] = vSource[i] * Coefficient;      
	  }
	break;
      }
    case 1 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AddAu(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] = vSource[i] * Coefficient;      
	  }
 	break;
     }
    case 2 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAd(i, this->CreationMomentumIndex, this->AnnihilationMomentumIndex, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] = vSource[i] * Coefficient;      
	  }
	break;
      }
    }
  return vDestination;
}

