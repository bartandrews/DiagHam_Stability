////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of S^- spin operator on the sphere                 //
//                                                                            //
//                        last modification : 20/02/2015                      //
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
#include "Operator/ParticleOnSphereWithSpinSMinusOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default data
//
// particle = hilbert space associated to the particles

ParticleOnSphereWithSpinSMinusOperator::ParticleOnSphereWithSpinSMinusOperator(ParticleOnSphereWithSpin* particle)
{
  this->Particle= (ParticleOnSphereWithSpin*) (particle->Clone());
  this->NbrOrbitals = this->Particle->GetNbrOrbitals();
  cout << "this->NbrOrbitals = " << this->NbrOrbitals << endl;
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnSphereWithSpinSMinusOperator::ParticleOnSphereWithSpinSMinusOperator(ParticleOnSphereWithSpinSMinusOperator& oper)
{
  this->Particle = (ParticleOnSphereWithSpin*) (oper.Particle->Clone());
  this->NbrOrbitals = oper.NbrOrbitals;
}

// destructor
//

ParticleOnSphereWithSpinSMinusOperator::~ParticleOnSphereWithSpinSMinusOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinSMinusOperator::Clone ()
{
  return new ParticleOnSphereWithSpinSMinusOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinSMinusOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinSMinusOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinSMinusOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereWithSpinSMinusOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int Index = 0;
  double Coefficient = 0.0;
  int TargetDim = this->Particle->GetTargetHilbertSpaceDimension();
  int NbrTranslations;
  for (int i = firstComponent; i < Last; ++i)
    {
      Complex& Tmp = vSource[i];
      for (int j = 0; j < this->NbrOrbitals; ++j)
	{
	  Index = this->Particle->AddAu(i, j, j, Coefficient);
	  if (Index < TargetDim)
	    {
	      vDestination[Index] += Tmp * Coefficient;		  
	    }
	}
    }
  return vDestination;
}
  


// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinSMinusOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int Index = 0;
  double Coefficient = 0.0;
  int TargetDim = this->Particle->GetTargetHilbertSpaceDimension();
  int NbrTranslations;
  for (int i = firstComponent; i < Last; ++i)
    {
      double& Tmp = vSource[i];
      for (int j = 0; j < this->NbrOrbitals; ++j)
	{
	  Index = this->Particle->AddAu(i, j, j, Coefficient);
	  if (Index < TargetDim)
	    {
	      vDestination[Index] += Tmp * Coefficient;		  
	    }
	}
    }
  return vDestination;
}
  


