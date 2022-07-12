////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of periodic anisotropic magnetization operator          //
//                                                                            //
//                        last modification : 23/03/2002                      //
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


#include "Operator/PeriodicAnisotropicMagnetizationOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// sxFactor = array containing factors in front of Sx term (g mu B)
// syFactor = array containing factors in front of Sy term (g mu B)
// szFactor = array containing factors in front of Sz term (g mu B)

PeriodicAnisotropicMagnetizationOperator::PeriodicAnisotropicMagnetizationOperator(AbstractSpinChain* chain, int nbrSpin, 
										   double* sxFactor, double* syFactor, double* szFactor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SpFactorRe = new double [this->NbrSpin];
  this->SpFactorIm = new double [this->NbrSpin];
  this->SzFactor = new double [this->NbrSpin];
  for (int i = 0; i < (this->NbrSpin); i++)
    {
      this->SpFactorRe[i] = 0.5 * sxFactor[i];
      this->SpFactorIm[i] = - 0.5 * syFactor[i];
      this->SzFactor[i] = szFactor[i];
    }
}

// destructor
//

PeriodicAnisotropicMagnetizationOperator::~PeriodicAnisotropicMagnetizationOperator()
{
  delete[] this->SpFactorRe;
  delete[] this->SpFactorIm;
  delete[] this->SzFactor;
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator*  PeriodicAnisotropicMagnetizationOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicAnisotropicMagnetizationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (AbstractSpinChain*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PeriodicAnisotropicMagnetizationOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int PeriodicAnisotropicMagnetizationOperator::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicAnisotropicMagnetizationOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  return Complex();
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicAnisotropicMagnetizationOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  double x = 0.0;
  double y = 0.0;
  double coef;
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double TmpCoef1;
  double TmpCoef2;
  double Tmpx;
  double Tmpy;
  for (int j = 0; j < this->NbrSpin; ++j)
    {
      if (this->SzFactor[j] != 0)
	{
	  TmpCoef1 = this->SzFactor[j];
	  for (int i = 0; i < dim; ++i)
	    {	 	    
	      pos = this->Chain->Szi(j, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
 		  x += coef * (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
	  	  y += coef * (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		}
	    }
	}
      if ((this->SpFactorRe[j] != 0) || (this->SpFactorIm[j] != 0))
	{
	  TmpCoef1 = this->SpFactorRe[j];
	  TmpCoef2 = this->SpFactorIm[j];
	  for (int i = 0; i < dim; ++i)
	    {	 	    
	      pos = this->Chain->Spi(j, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 - Tmpy * TmpCoef2);
		  y += coef * (Tmpx * TmpCoef2 + Tmpy * TmpCoef1);
		}
	      pos = this->Chain->Smi(j, i, coef);
	      if (pos != dim)
		{
		  Tmpx = (V1.Re(pos) * V2.Re(i) + V1.Im(pos) * V2.Im(i));
		  Tmpy = (V1.Re(pos) * V2.Im(i) - V1.Im(pos) * V2.Re(i));
		  x += coef * (Tmpx * TmpCoef1 + Tmpy * TmpCoef2);
		  y += coef * (Tmpy * TmpCoef1 - Tmpx * TmpCoef2);
		}
	    }
	}
    }
  return Complex(x, y);
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicAnisotropicMagnetizationOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									  int firstComponent, int nbrComponent)
{
  double coef;
  int pos;
  int dim = this->Chain->GetHilbertSpaceDimension();
  int Max = firstComponent + nbrComponent;
  double TmpCoef1;
  double TmpCoef2;
  TmpCoef1 = this->SzFactor[0];
  for (int i = firstComponent; i < Max; ++i)
    {	 	    
      pos = this->Chain->Szi(0, i, coef);
      if (pos != dim)
	{
	  coef *= TmpCoef1;
	  vDestination.Re(pos) = coef * vSource.Re(i);
	  vDestination.Im(pos) = coef * vSource.Im(i);
	}
      else
	{
	  vDestination.Re(pos) = 0;
	  vDestination.Im(pos) = 0;
	}
    }
  if ((this->SpFactorRe[0] != 0) || (this->SpFactorIm[0] != 0))
    {
      TmpCoef1 = this->SpFactorRe[0];
      TmpCoef2 = this->SpFactorIm[0];
      for (int i = firstComponent; i < Max; ++i)
	{	 	    
	  pos = this->Chain->Spi(0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Re(i) * TmpCoef2 + vSource.Im(i) * TmpCoef1);
	    }
	  pos = this->Chain->Smi(0, i, coef);
	  if (pos != dim)
	    {
	      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
	      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
	    }
	}
    }
  for (int j = 1; j < this->NbrSpin; ++j)
    {
      if (this->SzFactor[j] != 0)
	{
	  TmpCoef1 = this->SzFactor[j];
	  for (int i = firstComponent; i < Max; ++i)
	    {	 	    
	      pos = this->Chain->Szi(j, i, coef);
	      if (pos != dim)
		{
		  coef *= TmpCoef1;
		  vDestination.Re(pos) = coef * vSource.Re(i);
		  vDestination.Im(pos) = coef * vSource.Im(i);
		}
	    }
	  if ((this->SpFactorRe[j] != 0) || (this->SpFactorIm[j] != 0))
	    {
	      TmpCoef1 = this->SpFactorRe[j];
	      TmpCoef2 = this->SpFactorIm[j];
	      for (int i = firstComponent; i < Max; ++i)
		{	 	    
		  pos = this->Chain->Spi(j, i, coef);
		  if (pos != dim)
		    {
		      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 - vSource.Im(i) * TmpCoef2);
		      vDestination.Im(pos) += coef * (vSource.Re(i) * TmpCoef2 + vSource.Im(i) * TmpCoef1);
		    }
		  pos = this->Chain->Smi(j, i, coef);
		  if (pos != dim)
		    {
		      vDestination.Re(pos) += coef * (vSource.Re(i) * TmpCoef1 + vSource.Im(i) * TmpCoef2);
		      vDestination.Im(pos) += coef * (vSource.Im(i) * TmpCoef1 - vSource.Re(i) * TmpCoef2);
		    }
		}
	    }
	}
    }
  return vDestination;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, PeriodicAnisotropicMagnetizationOperator& O)
{
  ComplexVector TmpV2 (O.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [O.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(O.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < O.Chain->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < O.Chain->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j].Re(i);
	  if (TmpV[j].Im(i) < 0.0)
	    {
	      Str << ((TmpV[j].Im(i))) << "i";
	    }
	  else
	    if (TmpV[j].Im(i) > 0.0)
	      {
		Str << "+" << TmpV[j].Im(i) << "i";
	      }
	  Str << "   ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, PeriodicAnisotropicMagnetizationOperator& O)
{
  ComplexVector TmpV2 (O.Chain->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [O.Chain->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Chain->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(O.Chain->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (O.Chain->GetHilbertSpaceDimension() - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (O.Chain->GetHilbertSpaceDimension() - 1); ++j)
	{
	  Str << TmpV[j].Re(i);
	  if (TmpV[j].Im(i) < 0)
	    {
	      Str << ((TmpV[j].Im(i))) << "I";
	    }
	  else
	    if (TmpV[j].Im(i) > 0)
	      {
		Str << "+" << ((TmpV[j].Im(i))) << "I";
	      }	  
	  Str << ",";
	}
      Str << TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Re(i);
      if (TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(i) < 0)
	{
	  Str << ((TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(i))) << "I";
	}
      else
	if (TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(i) > 0)
	  {
	    Str << "+" << ((TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(i))) << "I";
	  }	  
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (O.Chain->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j].Re(O.Chain->GetHilbertSpaceDimension() - 1);
      if (TmpV[j].Im(O.Chain->GetHilbertSpaceDimension() - 1) < 0)
	{
	  Str << ((TmpV[j].Im(O.Chain->GetHilbertSpaceDimension() - 1))) << "I";
	}
      else
	if (TmpV[j].Im(O.Chain->GetHilbertSpaceDimension() - 1) > 0)
	  {
	    Str << "+" << ((TmpV[j].Im(O.Chain->GetHilbertSpaceDimension() - 1))) << "I";
	  }	  
      Str << ",";
    }
  Str << TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Re(O.Chain->GetHilbertSpaceDimension() - 1);
  if (TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(O.Chain->GetHilbertSpaceDimension() - 1) < 0)
    {
      Str << ((TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(O.Chain->GetHilbertSpaceDimension() - 1))) << "I";
    }
  else
    if (TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(O.Chain->GetHilbertSpaceDimension() - 1) > 0)
      {
	Str << "+" << ((TmpV[O.Chain->GetHilbertSpaceDimension() - 1].Im(O.Chain->GetHilbertSpaceDimension() - 1))) << "I";
      }	  
  Str << "}}";
  return Str;
}


