////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class density operator                          // 
//              for particle on a torus with magentic translations            //
//                                                                            //
//                        last modification : 16/07/2015                      //
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
#include "Operator/ParticleOnTorusWithMagneticTranslationsDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// constructor
//
// particle = hilbert space associated to the particles
// m = index of the creation operator
// n = index of the annihilation operator

ParticleOnTorusWithMagneticTranslationsDensityOperator::ParticleOnTorusWithMagneticTranslationsDensityOperator(ParticleOnTorusWithMagneticTranslations* particle, 
													       int m, int n)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) (particle->Clone());
  this->AnnihilationIndex = n;
  this->CreationIndex = m;
}

// copy constructor
//

ParticleOnTorusWithMagneticTranslationsDensityOperator::ParticleOnTorusWithMagneticTranslationsDensityOperator(ParticleOnTorusWithMagneticTranslationsDensityOperator& oper)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) (oper.Particle->Clone());  
  this->AnnihilationIndex = oper.AnnihilationIndex;
  this->CreationIndex = oper.CreationIndex;
}

// destructor
//

ParticleOnTorusWithMagneticTranslationsDensityOperator::~ParticleOnTorusWithMagneticTranslationsDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusWithMagneticTranslationsDensityOperator::Clone ()
{
  return new ParticleOnTorusWithMagneticTranslationsDensityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithMagneticTranslationsDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusWithMagneticTranslationsDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusWithMagneticTranslationsDensityOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnTorusWithMagneticTranslationsDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int NbrFluxQuanta = this->Particle->GetNbrOrbitals();
  Complex* TmpPhases = new Complex [NbrFluxQuanta];
  for (int i = 0; i < NbrFluxQuanta; ++i)
    {
      TmpPhases[i] =  Phase (2.0 * M_PI * ((double) (i * this->Particle->GetKxMomentum())) / ((double) NbrFluxQuanta));
    }


  Complex Element = 0.0;
  double Coefficient1 = 0.0;
  double Coefficient2 = 0.0;
  int NbrTranslations;
  int Index;
  
  if (this->AnnihilationIndex == this->CreationIndex)
    {
      for (int i = (int) firstComponent; i < Last; ++i)
	{
	  Element += Conj(V1[i]) * V2[i] * this->Particle->AdA(i, this->AnnihilationIndex);      
	}
    }
  delete[] TmpPhases;
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

ComplexVector& ParticleOnTorusWithMagneticTranslationsDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											   int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int NbrFluxQuanta = this->Particle->GetNbrOrbitals();
  Complex* TmpPhases = new Complex [NbrFluxQuanta];
  for (int i = 0; i < NbrFluxQuanta; ++i)
    {
      TmpPhases[i] =  Phase (2.0 * M_PI * ((double) (i * this->Particle->GetKxMomentum())) / ((double) NbrFluxQuanta));
    }

  double Coefficient1 = 0.0;
  double Coefficient2 = 0.0;
  int NbrTranslations;
  int Index;

  if (this->AnnihilationIndex == this->CreationIndex)
    {
      for (int i = (int) firstComponent; i < Last; ++i)
	{
	  vDestination[Index] += vSource[i] * this->Particle->AdA(i, this->AnnihilationIndex);      
	}
    }

  delete[] TmpPhases;
  return vDestination;
}

