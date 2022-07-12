////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of particle on sphere annihilation operator              //
//                                                                            //
//                        last modification : 22/07/2016                      //
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
#include "Operator/ParticleOnSphereAnnihilationOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default data
//
// particle = hilbert space associated to the particles
// index = index of the annihilation operator

ParticleOnSphereAnnihilationOperator::ParticleOnSphereAnnihilationOperator(ParticleOnSphere* particle, int index)
{
  this->Particle= (ParticleOnSphere*) (particle->Clone());
  this->OperatorIndex = index;
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnSphereAnnihilationOperator::ParticleOnSphereAnnihilationOperator(ParticleOnSphereAnnihilationOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->OperatorIndex = oper.OperatorIndex;
}

// destructor
//

ParticleOnSphereAnnihilationOperator::~ParticleOnSphereAnnihilationOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereAnnihilationOperator::Clone ()
{
  return new ParticleOnSphereAnnihilationOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereAnnihilationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereAnnihilationOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereAnnihilationOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnSphereAnnihilationOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  double Element = 0.0;
  double TmpCoefficient = 0.0;
  int TmpIndex ;
  int Dim = firstComponent + nbrComponent;
  for (int i = firstComponent; i < Dim; ++i)
    {
      TmpCoefficient = 0.0;
      TmpIndex = this->Particle->A(i, this->OperatorIndex, TmpCoefficient);
      if (TmpIndex < V1.GetVectorDimension())
	{
	  Element += V1[TmpIndex] * V2[i] * TmpCoefficient;
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

RealVector& ParticleOnSphereAnnihilationOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  double TmpCoefficient = 0.0;
  int TmpIndex ;
  for (int i = firstComponent; i < Last; ++i)
    {
      TmpCoefficient = 0.0;
      TmpIndex = this->Particle->A(i, this->OperatorIndex, TmpCoefficient);
      if (TmpIndex < vDestination.GetVectorDimension())
	{
	  vDestination[TmpIndex] += vSource[i] * TmpCoefficient;
	}
    }
  return vDestination;
}
  
// multiply a set of vectors by the current operator for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* ParticleOnSphereAnnihilationOperator::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  double TmpCoefficient = 0.0;
  int TmpIndex ;
  for (int i = firstComponent; i < Last; ++i)
    {
      TmpCoefficient = 0.0;
      TmpIndex = this->Particle->A(i, this->OperatorIndex, TmpCoefficient);
      if (TmpIndex < vDestinations[0].GetVectorDimension())
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][TmpIndex] += vSources[k][i] * TmpCoefficient;
	    }
	}
    }
  return vDestinations;
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereAnnihilationOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  Complex Element = 0.0;
  int Dim = firstComponent + nbrComponent;
  double TmpCoefficient = 0.0;
  int TmpIndex ;
  for (int i = firstComponent; i < Dim; ++i)
    {
      TmpCoefficient = 0.0;
      TmpIndex = this->Particle->A(i, this->OperatorIndex, TmpCoefficient);
      if (TmpIndex < V1.GetVectorDimension())
	{
	  Element += Conj(V1[TmpIndex]) * V2[i] * TmpCoefficient;
	}
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

ComplexVector& ParticleOnSphereAnnihilationOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  double TmpCoefficient = 0.0;
  int TmpIndex ;
  for (int i = firstComponent; i < Last; ++i)
    {
      TmpCoefficient = 0.0;
      TmpIndex = this->Particle->A(i, this->OperatorIndex, TmpCoefficient);
      if (TmpIndex < vDestination.GetVectorDimension())
	{
	  vDestination[TmpIndex] += vSource[i] * TmpCoefficient;
	}
    }
  return vDestination;
}
  
// multiply a set of vectors by the current operator for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnSphereAnnihilationOperator::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  double TmpCoefficient = 0.0;
  int TmpIndex ;
  for (int i = firstComponent; i < Last; ++i)
    {
      TmpCoefficient = 0.0;
      TmpIndex = this->Particle->A(i, this->OperatorIndex, TmpCoefficient);
      if (TmpIndex < vDestinations[0].GetVectorDimension())
	{
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][TmpIndex] += vSources[k][i] * TmpCoefficient;
	    }
	}
    }
  return vDestinations;
}



