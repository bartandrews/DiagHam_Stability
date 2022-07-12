////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class for the Sz parity operator                     // 
//             for particles with spin on a lattice using translations        //
//                                                                            //
//                        last modification : 04/07/2016                      //
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
#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// constructor from default data
//
// particle = hilbert space associated to the particles

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* particle)
{
  this->Particle = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) (particle->Clone());
}

// copy constructor
//

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator(ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator& oper)
{
  this->Particle = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) (oper.Particle->Clone());  
}

// destructor
//

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::~ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::Clone ()
{
  return new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::GetHilbertSpaceDimension ()
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

Complex ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();

  Complex** TmpPhases = new Complex* [this->Particle->GetMaxXMomentum()];
  for (int i = 0; i < this->Particle->GetMaxXMomentum(); ++i)
    {
      TmpPhases[i] = new Complex [this->Particle->GetMaxYMomentum()];
      for (int j = 0; j < this->Particle->GetMaxYMomentum(); ++j)
	TmpPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * this->Particle->GetKxMomentum())) / ((double) this->Particle->GetMaxXMomentum())
						+ ((double) (j * this->Particle->GetKyMomentum())) / ((double) this->Particle->GetMaxYMomentum())));
    }

  Complex Element = 0.0;
  double Coefficient;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;  

  
  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Index = this->Particle->ApplySzSymmetry(i, Coefficient, NbrTranslationsX, NbrTranslationsY);
      if (Index < V1.GetVectorDimension())
	Element += Conj(V1[Index]) * V2[i] * TmpPhases[NbrTranslationsX][NbrTranslationsY] * Coefficient;
      else
	{
	  cout << "error, " << i << " " << Index << " " << V1.GetVectorDimension() << " " << Dim << " " << Coefficient << endl;
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

ComplexVector& ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSzParityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												       int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();

  Complex** TmpPhases = new Complex* [this->Particle->GetMaxXMomentum()];
  for (int i = 0; i < this->Particle->GetMaxXMomentum(); ++i)
    {
      TmpPhases[i] = new Complex [this->Particle->GetMaxYMomentum()];
      for (int j = 0; j < this->Particle->GetMaxYMomentum(); ++j)
	TmpPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * this->Particle->GetKxMomentum())) / ((double) this->Particle->GetMaxXMomentum())
						+ ((double) (j * this->Particle->GetKyMomentum())) / ((double) this->Particle->GetMaxYMomentum())));
    }

  double Coefficient;
  int Index;
  int NbrTranslationsX;
  int NbrTranslationsY;  

  for (int i = (int) firstComponent; i < Last; ++i)
    {
      Index = this->Particle->ApplySzSymmetry(i, Coefficient, NbrTranslationsX, NbrTranslationsY);
      if (Index < vDestination.GetVectorDimension())
	vDestination[Index] += vSource[i] * TmpPhases[NbrTranslationsX][NbrTranslationsY] * Coefficient;      
    }

   return vDestination;
}

