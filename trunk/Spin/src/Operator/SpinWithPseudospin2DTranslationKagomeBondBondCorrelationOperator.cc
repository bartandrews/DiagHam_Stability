////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      Class author Cecile Repellin                          //
//                                                                            //
//                                                                            //
//            class of bond energy operator for spin-pseudospin model         //
//                                                                            //
//                        last modification : 14/04/2017                      //
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


#include "Operator/SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator.h"
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

SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator(Spin1_2ChainWithPseudospinAnd2DTranslation* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4, int bondIndex1, int bondIndex2)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->Index1 = siteIndex1;
  this->Index2 = siteIndex2;
  this->Index3 = siteIndex3;
  this->Index4 = siteIndex4;
  this->BondIndex1 = bondIndex1;
  this->BondIndex2 = bondIndex2;
  
  if (((this->BondIndex1 & 1) != 0) && ((this->BondIndex2 & 1) == 0))
  {
    this->Index1 = siteIndex3;
    this->Index2 = siteIndex4;
    this->Index3 = siteIndex1;
    this->Index4 = siteIndex2;
    this->BondIndex1 = bondIndex2;
    this->BondIndex2 = bondIndex1;
  }
  
  // projection of kagome onto the s = 1/2 states on each triangle
  this->PseudospinCouplingElements = new double[3];    
  this->PseudospinCouplingElements[0] = 0.0;
  this->PseudospinCouplingElements[1] = 1.0/sqrt(3.0);
  this->PseudospinCouplingElements[2] = -1.0/sqrt(3.0);
  
  this-> PseudospinDiagCouplingElements = new double*[3];
  for (int i = 0; i < 3; ++i)
    this->PseudospinDiagCouplingElements[i] = new double[2];
  this->PseudospinDiagCouplingElements[0][0] = 1.0;
  this->PseudospinDiagCouplingElements[0][1] = -1.0/3.0;
  
  this->PseudospinDiagCouplingElements[1][0] = 0.0;
  this->PseudospinDiagCouplingElements[1][1] = 2.0/3.0;
  
  this->PseudospinDiagCouplingElements[2][0] = 0.0;
  this->PseudospinDiagCouplingElements[2][1] = 2.0/3.0;
  
  this->EvaluateExponentialFactors();
  
}

// destructor
//

SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::~SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::Clone ()
{
  return new SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator(this->Chain, this->XMomentum, this->NbrSpinX, this->YMomentum, this->NbrSpinY, this->Index1, this->Index2, this->Index3, this->Index4, this->BondIndex1, this->BondIndex2);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (Spin1_2ChainWithPseudospinAnd2DTranslation*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::GetHilbertSpaceDimension ()
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

Complex SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
						     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Element = 0.0;
  int pos;
  int pos2;
  double coef2;
  double coef;
  double TmpCoefficient;
  double TmpCoefficient1;
  
  int nbrTranslationsX;
  int nbrTranslationsY;
      
  double** TmpDiagCoupling = new double*[3];
  for (int l = 0; l < 3; ++l)
    TmpDiagCoupling[l] = new double[2];
  double* TmpOffDiagCoupling = new double[3];
  
  
  TmpDiagCoupling[0][0] = 0.0;
  TmpDiagCoupling[0][1] = -0.5;  
  TmpDiagCoupling[1][0] = -0.75;
  TmpDiagCoupling[1][1] = 0.25;
  TmpDiagCoupling[2][0] = 0.0;
  TmpDiagCoupling[2][1] = -0.5;
  
  TmpOffDiagCoupling[0] = sqrt(3) / 4.0;
  TmpOffDiagCoupling[1] = 0.0;
  TmpOffDiagCoupling[2] = -sqrt(3.0) / 4.0;
  
  int TmpBond1 = this->BondIndex1 / 2;
  int TmpBond2 = this->BondIndex2 / 2;
  
  int Tmp1;
  int Tmp2;
  if (TmpBond2 == 0)      //AB
  {
    Tmp1 = 0;
    Tmp2 = 1;
  }
  if (TmpBond2 == 1)      //BC
  {
    Tmp1 = 1;
    Tmp2 = 2;
  }
  if (TmpBond2 == 2)      //AC
  {
    Tmp1 = 0;
    Tmp2 = 2;
  }
  
  int Tmpa1;
  int Tmpa2;
  if (TmpBond1 == 0)      //AB
  {
    Tmpa1 = 0;
    Tmpa2 = 1;
  }
  if (TmpBond1 == 1)      //BC
  {
    Tmpa1 = 1;
    Tmpa2 = 2;
  }
  if (TmpBond1 == 2)      //AC
  {
    Tmpa1 = 0;
    Tmpa2 = 2;
  }
  
  
  if (((this->BondIndex1 & 1) != 0) && ((this->BondIndex2 & 1) != 0))
  {
    cout << "Warning: odd value (" << (this->BondIndex1) << ") not debugged. To be continued" << endl;
  }
	  
  for (int i = (int) firstComponent; i < dim; ++i)
    {
	if (((this->BondIndex1 & 1) == 0) && ((this->BondIndex2 & 1) == 0))
	{
	  pos = this->Chain->JoffiJoffj (this->Index1, this->Index3, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += TmpOffDiagCoupling[TmpBond1] * TmpOffDiagCoupling[TmpBond2] * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  coef = this->Chain->JDiagonali(this->Index1, i, TmpDiagCoupling[TmpBond1]);
	  pos = this->Chain->JOffDiagonali(this->Index3, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += TmpOffDiagCoupling[TmpBond2] * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, TmpDiagCoupling[TmpBond2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  
	  
	  coef2 = this->Chain->JDiagonali(this->Index3, i, TmpDiagCoupling[TmpBond2]);
	  Element += Conj(V1[i]) * V2[i] * coef2 * coef;
	}
	
	if (((this->BondIndex1 & 1) == 0) && ((this->BondIndex2 & 1) != 0))
	{
	  //OffDiag
	  //SzSz
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  ////SmSp
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  ////SpSm
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  
	  //Diag
// 	  SzSz
	  coef = this->Chain->JDiagonali(this->Index1, i, TmpDiagCoupling[TmpBond1]);
	  coef2 = this->Chain->SziSzj(this->Index3, this->Index4, i);
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, i, this->PseudospinDiagCouplingElements[Tmp1]) * this->Chain->JDiagonali(this->Index4, i, this->PseudospinDiagCouplingElements[Tmp2]);
	  Element += Conj(V1[i]) * V2[i] * coef * coef2 * TmpCoefficient;
	     
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, i, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->JOffDiagonali(this->Index4, i, TmpCoefficient1, nbrTranslationsX, nbrTranslationsY);
	  TmpCoefficient1 *= this->PseudospinCouplingElements[Tmp2];
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * TmpCoefficient1 * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  pos = this->Chain->JoffiJoffj(this->Index3, this->Index4, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1];
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2];
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  //////SmSp
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  ////SpSm
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	}
	if (((this->BondIndex1 & 1) != 0) && ((this->BondIndex2 & 1) != 0))
	{
	  //SzSzSzSz
	  coef = this->Chain->SziSzj(this->Index1, this->Index2, i);
	  coef *= this->Chain->SziSzj(this->Index3, this->Index4, i);
	  
	  //SzSzSzSz JdiagJdiagJJ
	  coef2 = this->Chain->JDiagonali(this->Index1, i, this->PseudospinDiagCouplingElements[Tmpa1]);
	  coef2 *= this->Chain->JDiagonali(this->Index2, i, this->PseudospinDiagCouplingElements[Tmpa2]);
	  
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, i, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, i, this->PseudospinDiagCouplingElements[Tmp2]);
	  Element += Conj(V1[i]) * V2[i] * coef * coef2 * TmpCoefficient;
	  
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, i, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->JOffDiagonali(this->Index4, i, TmpCoefficient1, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * TmpCoefficient1 * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  TmpCoefficient = this->Chain->JDiagonali(this->Index4, i, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->JOffDiagonali(this->Index3, i, TmpCoefficient1, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * TmpCoefficient1 * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  pos = this->Chain->JoffiJoffj(this->Index3, this->Index4, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  
	  //SzSzSzSz JdiagJoffJJ
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);  
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *=  this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]); 
	  TmpCoefficient = this->Chain->JDiagonali(this->Index4, i, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient =  this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSzSz JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  coef2 *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  coef2 *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);  
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *=  this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  coef2 *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]); 
	  TmpCoefficient = this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  coef2 *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient =  this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSzSz JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);  
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *=  this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1); 
	  TmpCoefficient = this->Chain->JDiagonali(this->Index4, i, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  coef2 = this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  coef2 *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient =  this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += Conj(V1[pos]) * V2[i] * coef * coef2 * TmpCoefficient * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSpSm
	  //SzSzSpSm JDiagJdiagJJ
	  coef = this->Chain->SziSzj(this->Index1, this->Index2, i);
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSpSm JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSpSm JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSmSp
	  //SzSzSmSp JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSpSm JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSpSm JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SzSzSpSm JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz
	  //SmSpSzSz JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SpsmSzSz
	  //SmSpSzSz JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.5 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
// 	  
// 	  
	  //SmSpSmSp
	  //SmSpSmSp JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSpSm
	  //SmSpSmSp JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index2, this->Index1);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSzSz JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JDiagJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JDiagJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index1, this->PseudospinDiagCouplingElements[Tmpa1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSmSp JOffJDiagJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index2, this->PseudospinDiagCouplingElements[Tmpa2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  //SmSpSpSm JOffJOffJJ
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	 
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->Chain->SmiSpj(this->Index1, this->Index2);
	  TmpCoefficient *= this->Chain->SmiSpj(this->Index4, this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa2] * this->Chain->JOffDiagonali(this->Index2);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmpa1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    Element += 0.25 * Conj(V1[pos]) * V2[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	}
    }
//   cout << (this->Index1) << " "<< (this->Index2) << " " << (this->Index3) << " " << (this->Index4) << " "  << Tmpa1 << " " << Tmpa2 << " " << Tmp1 << " " << Tmp2 << " -> " << Element << endl;
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

ComplexVector& SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						    int firstComponent, int nbrComponent)
{
  cout << "Warning: untested method SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::LowLevelMultiply" << endl;
  int dim = (int) (firstComponent + nbrComponent);
  double TmpCoefficient;
  double TmpCoefficient1;
  int pos;
  int pos2;
  double coef2;
  double coef;
  
  int nbrTranslationsX;
  int nbrTranslationsY;
  
  double** TmpDiagCoupling = new double*[3];
  for (int l = 0; l < 3; ++l)
    TmpDiagCoupling[l] = new double[2];
  double* TmpOffDiagCoupling = new double[3];
  
  
  TmpDiagCoupling[0][0] = 0.0;
  TmpDiagCoupling[0][1] = -0.5;  
  TmpDiagCoupling[1][0] = -0.75;
  TmpDiagCoupling[1][1] = 0.25;
  TmpDiagCoupling[2][0] = 0.0;
  TmpDiagCoupling[2][1] = -0.5;
  
  TmpOffDiagCoupling[0] = sqrt(3) / 4.0;
  TmpOffDiagCoupling[1] = 0.0;
  TmpOffDiagCoupling[2] = -sqrt(3.0) / 4.0;
  
  int TmpBond1 = this->BondIndex1 / 2;
  int TmpBond2 = this->BondIndex2 / 2;
  
  int Tmp1;
  int Tmp2;
  if (TmpBond2 == 0)      //AB
  {
    Tmp1 = 0;
    Tmp2 = 1;
  }
  if (TmpBond2 == 1)      //BC
  {
    Tmp1 = 1;
    Tmp2 = 2;
  }
  if (TmpBond2 == 2)      //AC
  {
    Tmp1 = 0;
    Tmp2 = 2;
  }
  if (((this->BondIndex1 & 1) != 0) && ((this->BondIndex2 & 1) != 0))
  {
    cout << "Error: odd value (" << (this->BondIndex1) << ") not implemented" << endl;
  }
	
	
  for (int i = (int) firstComponent; i < dim; ++i)
    {
	if (((this->BondIndex1 & 1) == 0) && ((this->BondIndex2 & 1) == 0))
	{
	  pos = this->Chain->JoffiJoffj (this->Index1, this->Index3, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += TmpOffDiagCoupling[TmpBond1] * TmpOffDiagCoupling[TmpBond2] * vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  coef = this->Chain->JDiagonali(this->Index1, i, TmpDiagCoupling[TmpBond1]);
	  pos = this->Chain->JOffDiagonali(this->Index3, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += TmpOffDiagCoupling[TmpBond2] * vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, TmpDiagCoupling[TmpBond2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  
	  
	  coef2 = this->Chain->JDiagonali(this->Index3, i, TmpDiagCoupling[TmpBond2]);
	  vDestination[i] += vSource[i] * coef2 * coef;
	}
	
	if (((this->BondIndex1 & 1) == 0) && ((this->BondIndex2 & 1) != 0))
	{
	  //OffDiag
	  //SzSz
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= this->Chain->SziSzj(this->Index3, this->Index4);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  ////SmSp
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  ////SpSm
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = TmpOffDiagCoupling[TmpBond1] * this->Chain->JOffDiagonali(this->Index1);
	  TmpCoefficient *= (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  
	  //Diag
// 	  SzSz
	  coef = this->Chain->JDiagonali(this->Index1, i, TmpDiagCoupling[TmpBond1]);
	  coef2 = this->Chain->SziSzj(this->Index3, this->Index4, i);
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, i, this->PseudospinDiagCouplingElements[Tmp1]) * this->Chain->JDiagonali(this->Index4, i, this->PseudospinDiagCouplingElements[Tmp2]);
	  vDestination[i] += vSource[i] * coef * coef2 * TmpCoefficient;
	     
	  TmpCoefficient = this->Chain->JDiagonali(this->Index3, i, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->JOffDiagonali(this->Index4, i, TmpCoefficient1, nbrTranslationsX, nbrTranslationsY);
	  TmpCoefficient1 *= this->PseudospinCouplingElements[Tmp2];
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * coef2 * TmpCoefficient * TmpCoefficient1 * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  pos = this->Chain->JoffiJoffj(this->Index3, this->Index4, i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1];
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2];
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * coef2 * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  //////SmSp
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index3, this->Index4));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  ////SpSm
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index4, this->PseudospinDiagCouplingElements[Tmp2]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  TmpCoefficient *= this->Chain->JDiagonali(this->Index3, this->PseudospinDiagCouplingElements[Tmp1]);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	  this->Chain->InitializeTransientState(i);
	  TmpCoefficient = (0.5 * this->Chain->SmiSpj(this->Index4, this->Index3));
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp1] * this->Chain->JOffDiagonali(this->Index3);
	  TmpCoefficient *= this->PseudospinCouplingElements[Tmp2] * this->Chain->JOffDiagonali(this->Index4);
	  pos = this->Chain->SymmetrizeResult(i, TmpCoefficient, nbrTranslationsX, nbrTranslationsY);
	  if (pos != dim)
	    vDestination[pos] += vSource[i] * coef * TmpCoefficient * this->ExponentialFactors[nbrTranslationsY][nbrTranslationsY];
	  
	}
    }
  return vDestination;
}


// evaluate all exponential factors
//   

void SpinWithPseudospin2DTranslationKagomeBondBondCorrelationOperator::EvaluateExponentialFactors()
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

