////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class density-density operator                      // 
//              for particle on a torus with magentic translations            //
//                                                                            //
//                        last modification : 14/07/2015                      //
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
#include "Operator/ParticleOnTorusWithMagneticTranslationsDensityDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// constructor
//
// particle = hilbert space associated to the particles
// m1 = index of the first creation operator
// m2 = index of the second creation operator
// n1 = index of the first annihilation operator
// n2 = index of the second annihilation operator

ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::ParticleOnTorusWithMagneticTranslationsDensityDensityOperator(ParticleOnTorusWithMagneticTranslations* particle, 
															     int m1, int m2, int n1, int n2)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) (particle->Clone());
  this->AnnihilationIndex1 = n1;
  this->AnnihilationIndex2 = n2;
  this->CreationIndex1 = m1;
  this->CreationIndex2 = m2;
}

// copy constructor
//

ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::ParticleOnTorusWithMagneticTranslationsDensityDensityOperator(ParticleOnTorusWithMagneticTranslationsDensityDensityOperator& oper)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) (oper.Particle->Clone());  
  this->AnnihilationIndex1 = oper.AnnihilationIndex1;
  this->AnnihilationIndex2 = oper.AnnihilationIndex2;
  this->CreationIndex1 = oper.CreationIndex1;
  this->CreationIndex2 = oper.CreationIndex2;
}

// destructor
//

ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::~ParticleOnTorusWithMagneticTranslationsDensityDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::Clone ()
{
  return new ParticleOnTorusWithMagneticTranslationsDensityDensityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
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
  
  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Coefficient1 = this->Particle->AA(i, this->AnnihilationIndex1, this->AnnihilationIndex2);
      if (Coefficient1 != 0.0)
	{
	  Index = this->Particle->AdAd(this->CreationIndex1, this->CreationIndex2, Coefficient2, NbrTranslations);
	  if (Index != Dim)
	    { 
	      Element += Conj(V1[Index]) * V2[i] * TmpPhases[NbrTranslations] * (Coefficient1 * Coefficient2);      
	    }
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

ComplexVector& ParticleOnTorusWithMagneticTranslationsDensityDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Coefficient1 = this->Particle->AA(i, this->AnnihilationIndex1, this->AnnihilationIndex2);
      if (Coefficient1 != 0.0)
	{
	  Index = this->Particle->AdAd(this->CreationIndex1, this->CreationIndex2, Coefficient2, NbrTranslations);
	  if (Index != Dim)
	    {
	      vDestination[Index] += vSource[i] * TmpPhases[NbrTranslations] * Coefficient1 * Coefficient2;      
	    }
	}
    }

  delete[] TmpPhases;
  return vDestination;
}

