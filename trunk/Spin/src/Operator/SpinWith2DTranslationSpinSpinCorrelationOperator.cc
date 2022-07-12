////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      Class author Cecile Repellin                          //
//                                                                            //
//                                                                            //
//            class of bond energy operator for spin-spin model with 2D       //
//                             translations                                   //
//                                                                            //
//                        last modification : 07/02/2017                      //
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


#include "Operator/SpinWith2DTranslationSpinSpinCorrelationOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/ComplexVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// chain = pointer to the Hilbert space
// nbrSpin = number of spins

SpinWith2DTranslationSpinSpinCorrelationOperator::SpinWith2DTranslationSpinSpinCorrelationOperator(Spin1_2ChainNewAnd2DTranslation* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, int siteIndex1, int siteIndex2)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->SiteIndex1 = siteIndex1;
  this->SiteIndex2 = siteIndex2;
  
  
  this->EvaluateExponentialFactors();
  
}

// destructor
//

SpinWith2DTranslationSpinSpinCorrelationOperator::~SpinWith2DTranslationSpinSpinCorrelationOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* SpinWith2DTranslationSpinSpinCorrelationOperator::Clone ()
{
  return new SpinWith2DTranslationSpinSpinCorrelationOperator(this->Chain, this->XMomentum, this->NbrSpinX, this->YMomentum, this->NbrSpinY, this->SiteIndex1, this->SiteIndex2);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinWith2DTranslationSpinSpinCorrelationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (Spin1_2ChainNewAnd2DTranslation*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinWith2DTranslationSpinSpinCorrelationOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int SpinWith2DTranslationSpinSpinCorrelationOperator::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex SpinWith2DTranslationSpinSpinCorrelationOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
						     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Element = 0.0;
  int pos;
  double TmpCoefficient;
  
  int nbrTranslationsX;
  int nbrTranslationsY;
    
  for (int i = (int) firstComponent; i < dim; ++i)
    {
	TmpCoefficient = this->Chain->SziSzj(this->SiteIndex1, this->SiteIndex2, i);
	Element += Conj(V1[i]) * TmpCoefficient * V2[i];
	
	
	pos = this->Chain->SmiSpj(this->SiteIndex2, this->SiteIndex1, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) * 0.5 * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpj(this->SiteIndex1, this->SiteIndex2, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) * 0.5 * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
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

ComplexVector& SpinWith2DTranslationSpinSpinCorrelationOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						    int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex TmpCoefficient;
  int pos;
  int pos2;
  double coef2;
  double coef;
  
  int nbrTranslationsX;
  int nbrTranslationsY;
  
 
  for (int i = (int) firstComponent; i < dim; ++i)
    {
// 	TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
// 	TmpIndex2 = this->GetLinearizedIndex(j, k);
      
	TmpCoefficient = this->Chain->SziSzj(this->SiteIndex1, this->SiteIndex2, i);
	vDestination[i] += TmpCoefficient * vSource[i];
	
	pos = this->Chain->SmiSpj(this->SiteIndex2, this->SiteIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += 0.5 * vSource[i] * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	pos = this->Chain->SmiSpj(this->SiteIndex1, this->SiteIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += 0.5 * vSource[i] * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
    }
  return vDestination;
}


// evaluate all exponential factors
//   

void SpinWith2DTranslationSpinSpinCorrelationOperator::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->NbrSpinY];
      for (int j = 0; j < this->NbrSpinY; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->NbrSpinX))
							       + (this->YMomentum * ((double) j) / ((double) this->NbrSpinY))));
	}
    }
}

