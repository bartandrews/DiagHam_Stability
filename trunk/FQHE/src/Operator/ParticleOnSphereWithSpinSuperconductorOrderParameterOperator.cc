////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class superconductor order parameter operator              // 
//                           for particle with spin                           //
//                                                                            //
//                        last modification : 03/07/2014                      //
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
#include "Operator/ParticleOnSphereWithSpinSuperconductorOrderParameterOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::ParticleOnSphereWithSpinSuperconductorOrderParameterOperator()
{
  this->Particle = 0;
  this->CreationMomentumIndex1 = 0;
  this->CreationSymmetryIndex1 = 0;
  this->CreationMomentumIndex2 = 0;
  this->CreationSymmetryIndex2 = 0;
  this->CombinationFlag = false;
  this->CreationSymmetryIndex1SecondTerm = 0;
  this->CreationSymmetryIndex2SecondTerm = 0;
  this->CombinationSign = 0.0;  
}

// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for up, 1 for plus)
// creationMomentumIndex2 = momentum index of the rightmost creation operator (from 0 to 2S)
// creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for up, 1 for plus)

ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(ParticleOnSphereWithSpin* particle, 
															   int creationMomentumIndex1, int creationSymmetryIndex1,
															   int creationMomentumIndex2, int creationSymmetryIndex2)
{
  this->Particle = (ParticleOnSphereWithSpin*) (particle->Clone());
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->CreationSymmetryIndex1 = creationSymmetryIndex1;
  this->CreationMomentumIndex2 = creationMomentumIndex2;
  this->CreationSymmetryIndex2 = creationSymmetryIndex2;
  this->CombinationFlag = false;
  this->CreationSymmetryIndex1SecondTerm = creationSymmetryIndex1;
  this->CreationSymmetryIndex2SecondTerm = creationSymmetryIndex2;
  this->CombinationSign = 0.0;  
}

// constructor for operator such as a+_{sigma_1,i_1} a+_{sigma_2,i_2} +/- a+_{sigma_3,i_1} a+_{sigma_4,i_2}
//
// particle = hilbert space associated to the particles with the small number of particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for down, 1 for up)
// creationMomentumIndex2 = momentum index of the rightmost creation operator (from 0 to 2S)
// creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for down, 1 for up)
// creationSymmetryIndex1SecondTerm = symmetry index of the rightmost creation operator for the second term
// creationSymmetryIndex2SecondTerm = symmetry index of the leftmost creation operator for the second term
// sign = sign in front of the second term

ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(ParticleOnSphereWithSpin* particle,  int creationMomentumIndex1, int creationSymmetryIndex1,
															   int creationMomentumIndex2, int creationSymmetryIndex2,
															   int creationSymmetryIndex1SecondTerm, int creationSymmetryIndex2SecondTerm, 
															   double sign)
{
  this->Particle = (ParticleOnSphereWithSpin*) (particle->Clone());
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->CreationSymmetryIndex1 = creationSymmetryIndex1;
  this->CreationMomentumIndex2 = creationMomentumIndex2;
  this->CreationSymmetryIndex2 = creationSymmetryIndex2;
  this->CombinationFlag = true;
  this->CreationSymmetryIndex1SecondTerm = creationSymmetryIndex1SecondTerm;
  this->CreationSymmetryIndex2SecondTerm = creationSymmetryIndex2SecondTerm;
  this->CombinationSign = sign;
}

// copy constructor
//

ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(ParticleOnSphereWithSpinSuperconductorOrderParameterOperator& oper)
{
  this->Particle = (ParticleOnSphereWithSpin*) (oper.Particle->Clone());
  this->CreationMomentumIndex1 = oper.CreationMomentumIndex1;
  this->CreationSymmetryIndex1 = oper.CreationSymmetryIndex1;
  this->CreationMomentumIndex2 = oper.CreationMomentumIndex2;
  this->CreationSymmetryIndex2 = oper.CreationSymmetryIndex2;
  this->CombinationFlag = oper.CombinationFlag;
  this->CombinationSign = oper.CombinationSign;  
  this->CreationSymmetryIndex1SecondTerm = oper.CreationSymmetryIndex1SecondTerm;
  this->CreationSymmetryIndex2SecondTerm = oper.CreationSymmetryIndex2SecondTerm;
  
}

// destructor
//

ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::~ParticleOnSphereWithSpinSuperconductorOrderParameterOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::Clone ()
{
  return new ParticleOnSphereWithSpinSuperconductorOrderParameterOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetTargetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 1) | this->CreationSymmetryIndex2;
  double Sign = 1.0;
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != FullDim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
    case 3 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != FullDim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
    case 1 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
	    if (Index != FullDim)
	      Element += V1[Index] * V2[i] * (Sign * Coefficient);      
	  }
      }
      break;
    case 2 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != FullDim)
	      Element += V1[Index] * V2[i] * Coefficient;      
	  }
      }
      break;
   }
  if (this->CombinationFlag == true)
    {
      SymmetryIndex = (this->CreationSymmetryIndex1SecondTerm << 1) | this->CreationSymmetryIndex2SecondTerm;
      switch (SymmetryIndex)
	{
	case 0 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * this->CombinationSign * Coefficient;      
	      }
	  }
	  break;
	case 3 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * this->CombinationSign * Coefficient;      
	      }
	  }
	  break;
	case 1 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * (Sign * this->CombinationSign * Coefficient);      
	      }
	  }
	  break;
	case 2 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != FullDim)
		  Element += V1[Index] * V2[i] * this->CombinationSign * Coefficient;      
	      }
	  }
	  break;
	}
    }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
											      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetTargetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 1) | this->CreationSymmetryIndex2;
  double Sign = 1.0;
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * Coefficient;      
	  }
      }
      break;
    case 3 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * Coefficient;      
	  }
	break;
      }
    case 1 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * (Sign * Coefficient);      
	  }
 	break;
     }
    case 2 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * Coefficient;      
	  }
	break;
      }
    }
  if (this->CombinationFlag == true)
    {
      SymmetryIndex = (this->CreationSymmetryIndex1SecondTerm << 1) | this->CreationSymmetryIndex2SecondTerm;
      switch (SymmetryIndex)
	{
	case 0 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * Coefficient;      
	      }
	  }
	  break;
	case 3 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * Coefficient;      
	      }
	    break;
	  }
	case 1 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * (Sign * Coefficient);      
	      }
	    break;
	  }
	case 2 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * Coefficient;      
	      }
	break;
	  }
	}
    }
  return vDestination;
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetTargetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 1) | this->CreationSymmetryIndex2;
  double Sign = 1.0;
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != FullDim)
	      Element += Conj(V1[Index]) * V2[i] * Coefficient;      
	  }
      }
      break;
    case 3 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != FullDim)
	      {
		Element += Conj(V1[Index]) * V2[i] * Coefficient;      
	      }
	  }
      }
      break;
    case 1 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
	    if (Index != FullDim)
	      Element += Conj(V1[Index]) * V2[i] * (Sign * Coefficient);      
	  }
      }
      break;
    case 2 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != FullDim)
	      Element += Conj(V1[Index]) * V2[i] * Coefficient;      
	  }
      }
      break;
   }
  if (this->CombinationFlag == true)
    {
      SymmetryIndex = (this->CreationSymmetryIndex1SecondTerm << 1) | this->CreationSymmetryIndex2SecondTerm; 
      switch (SymmetryIndex)
	{
	case 0 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != FullDim)
		  Element += Conj(V1[Index]) * V2[i] * (this->CombinationSign * Coefficient);      
	      }
	  }
	  break;
	case 3 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != FullDim)
		  {
		    Element += Conj(V1[Index]) * V2[i] * (this->CombinationSign * Coefficient);      
		  }
	      }
	  }
	  break;
	case 1 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
		if (Index != FullDim)
		  Element += Conj(V1[Index]) * V2[i] * (this->CombinationSign * Sign * Coefficient);      
	      }
	  }
	  break;
	case 2 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != FullDim)
		  Element += Conj(V1[Index]) * V2[i] * (this->CombinationSign * Coefficient);      
	      }
	  }
	  break;
	}
   }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereWithSpinSuperconductorOrderParameterOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetTargetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 1) | this->CreationSymmetryIndex2;
  double Sign = 1.0;
  if (this->Particle->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * Coefficient;      
	  }
      }
      break;
    case 3 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * Coefficient;      
	  }
	break;
      }
    case 1 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * (Sign * Coefficient);      
	  }
 	break;
     }
    case 2 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
	    if (Index != Dim)
	      vDestination[Index] += vSource[i] * Coefficient;      
	  }
	break;
      }
    }
  if (this->CombinationFlag == true)
    {
      SymmetryIndex = (this->CreationSymmetryIndex1SecondTerm << 1) | this->CreationSymmetryIndex2SecondTerm;
      switch (SymmetryIndex)
	{
	case 0 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AddAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * Coefficient;      
	      }
	  }
	  break;
	case 3 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AduAdu(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * Coefficient;      
	      }
	    break;
	  }
	case 1 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex2, this->CreationMomentumIndex1, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * (Sign * Coefficient);      
	      }
	    break;
	  }
	case 2 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		int Index = this->Particle->AduAdd(i, this->CreationMomentumIndex1, this->CreationMomentumIndex2, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += this->CombinationSign * vSource[i] * Coefficient;      
	      }
	break;
	  }
	}
    }
  return vDestination;
}

