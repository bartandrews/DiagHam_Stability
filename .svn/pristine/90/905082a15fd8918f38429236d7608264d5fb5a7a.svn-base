////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              class S^2 operator                            // 
//             for particle with spin on a lattice using translations         //
//                                                                            //
//                        last modification : 11/06/2015                      //
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
#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// fixedSzFlag = true if the Sz value is fixed for the Hilbert space

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* particle, bool fixedSzFlag)
{
  this->Particle = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) (particle->Clone());
  this->FixedSzFlag = fixedSzFlag;
}

// copy constructor
//

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator(ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator& oper)
{
  this->Particle = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) (oper.Particle->Clone());  
  this->FixedSzFlag = oper.FixedSzFlag;
}

// destructor
//

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::~ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::Clone ()
{
  return new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::GetHilbertSpaceDimension ()
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

Complex ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int NbrSites = this->Particle->GetNbrOrbitals();
  double S2ConstantShift = 0.5 * ((double) this->Particle->GetNbrParticles());
  if (this->FixedSzFlag == true)
    {
      S2ConstantShift += 0.25 * ((double) (this->Particle->GetTotalSpin() * this->Particle->GetTotalSpin()));
    }
  Complex** TmpPhases = new Complex* [this->Particle->GetMaxXMomentum()];
  for (int i = 0; i < this->Particle->GetMaxXMomentum(); ++i)
    {
      TmpPhases[i] = new Complex [this->Particle->GetMaxYMomentum()];
      for (int j = 0; j < this->Particle->GetMaxYMomentum(); ++j)
	TmpPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * this->Particle->GetKxMomentum())) / ((double) this->Particle->GetMaxXMomentum())
						+ ((double) (j * this->Particle->GetKyMomentum())) / ((double) this->Particle->GetMaxYMomentum())));
    }

  Complex Element = 0.0;
  double Coefficient1 = 0.0;
  double Coefficient2 = 0.0;
  int NbrTranslationsX;
  int NbrTranslationsY;  
  int Index;
  
  for (int i = (int) firstComponent; i < Last; ++i)
    {
      for (int m = 0; m < NbrSites; ++m)
	{
	  for (int n = 0; n < NbrSites; ++n)
	    {
	      Coefficient1 = this->Particle->AuAd(i, n, m);
	      if (Coefficient1 != 0.0)
		{
		  Index = this->Particle->AduAdd(m, n, Coefficient2, NbrTranslationsX, NbrTranslationsY);
		  if (Index != Dim)
		    {
		      Element += Conj(V1[Index]) * V2[i] * TmpPhases[NbrTranslationsX][NbrTranslationsY] * (Coefficient1 * Coefficient2);      
		    }
		}
	    }
	}
      Element += Conj(V1[i]) * V2[i] * S2ConstantShift;
      if (this->FixedSzFlag == false)
	{
	  double Tmp = 0.0;
	  for (int m= 0; m < NbrSites; ++m)
	    {
	      Tmp += (this->Particle->AduAu(i, m) - this->Particle->AddAd(i, m));
	    }
	  Element += Conj(V1[Index]) * V2[i] * (0.25 * Tmp * Tmp);
	}
    }

  for (int i = 0; i < this->Particle->GetMaxXMomentum(); ++i)
    delete[] TmpPhases[i];
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

ComplexVector& ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationS2Operator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
												 int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int NbrSites = this->Particle->GetNbrOrbitals();
  double S2ConstantShift = 0.5 * ((double) this->Particle->GetNbrParticles());
  if (this->FixedSzFlag == true)
    {
      S2ConstantShift += 0.25 * ((double) (this->Particle->GetTotalSpin() * this->Particle->GetTotalSpin()));
    }
  Complex** TmpPhases = new Complex* [this->Particle->GetMaxXMomentum()];
  for (int i = 0; i < this->Particle->GetMaxXMomentum(); ++i)
    {
      TmpPhases[i] = new Complex [this->Particle->GetMaxYMomentum()];
      for (int j = 0; j < this->Particle->GetMaxYMomentum(); ++j)
	TmpPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * this->Particle->GetKxMomentum())) / ((double) this->Particle->GetMaxXMomentum())
						+ ((double) (j * this->Particle->GetKyMomentum())) / ((double) this->Particle->GetMaxYMomentum())));
    }

  double Coefficient1 = 0.0;
  double Coefficient2 = 0.0;
  int NbrTranslationsX;
  int NbrTranslationsY;  
  int Index;

  for (int i = (int) firstComponent; i < Last; ++i)
    {
      for (int m = 0; m < NbrSites; ++m)
	{
	  for (int n = 0; n < NbrSites; ++n)
	    {
	      Coefficient1 = this->Particle->AuAd(i, n, m);
	      if (Coefficient1 != 0.0)
		{
		  Index = this->Particle->AduAdd(m, n, Coefficient2, NbrTranslationsX, NbrTranslationsY);
		  if (Index != Dim)
		    {
		      vDestination[Index] += vSource[i] * TmpPhases[NbrTranslationsX][NbrTranslationsY] * Coefficient1 * Coefficient2;      
		    }
		}
	    }
	}
      vDestination[i] += vSource[i] * S2ConstantShift;
      if (this->FixedSzFlag == false)
	{
	  double Tmp = 0.0;
	  for (int m= 0; m < NbrSites; ++m)
	    {
	      Tmp += (this->Particle->AduAu(i, m) - this->Particle->AddAd(i, m));
	    }
	  vDestination[i] += vSource[i] *(0.25 * Tmp * Tmp);
	}
    }

  for (int i = 0; i < this->Particle->GetMaxXMomentum(); ++i)
    delete[] TmpPhases[i];
  delete[] TmpPhases;
  return vDestination;
}

