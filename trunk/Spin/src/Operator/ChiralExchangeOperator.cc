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


#include "Operator/ChiralExchangeOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// collection = Hilbert-space to act upon
ChiralExchangeOperator::ChiralExchangeOperator(GenericSUNSpinCollection *collection)
{
  this->Collection=collection;
  this->S1=1;
  this->S2=2;
  this->S3=3;
}


// constructor from default datas
//
// collection = Hilbert-space to act upon
// s_i = spins being exchanged
ChiralExchangeOperator::ChiralExchangeOperator(GenericSUNSpinCollection *collection, int s1, int s2, int s3)
{
  this->Collection=collection;
  this->S1=s1;
  this->S2=s2;
  this->S3=s3;
}


// destructor
//

ChiralExchangeOperator::~ChiralExchangeOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator*  ChiralExchangeOperator::Clone ()
{
  return (AbstractOperator*) (new ChiralExchangeOperator(this->Collection, this->S1, this->S2, this->S3));
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ChiralExchangeOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Collection = (GenericSUNSpinCollection*) hilbertSpace;
}

// set indices of spins to be acted upon
//
// s1, s2, s3 : indices
//
void ChiralExchangeOperator::SetSpinIndices(int s1, int s2, int s3)
{
  this->S1=s1;
  this->S2=s2;
  this->S3=s3;
}


// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ChiralExchangeOperator::GetHilbertSpace ()
{
  return this->Collection;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ChiralExchangeOperator::GetHilbertSpaceDimension ()
{
  return this->Collection->GetHilbertSpaceDimension();
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ChiralExchangeOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  double Element = 0.0;
  int TmpIndex;
  int Dim = (int) (firstComponent + nbrComponent);
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      TmpIndex =  this->Collection->CyclicSpinPermutation(i, S1, S2, S3);
      Element += V1[TmpIndex] * V2[i];
      TmpIndex =  this->Collection->CyclicSpinPermutation(i, S3, S2, S1);
      Element += V1[TmpIndex] * V2[i];
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

RealVector& ChiralExchangeOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int TmpIndex;
  vDestination.ClearVector();
  for (int i = firstComponent; i < Last; ++i)
    {
      TmpIndex =  this->Collection->CyclicSpinPermutation(i, S1, S2, S3);
      vDestination[TmpIndex] += vSource[i];
      TmpIndex =  this->Collection->CyclicSpinPermutation(i, S3, S2, S1);
      vDestination[TmpIndex] += vSource[i];
    }
  return vDestination;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ChiralExchangeOperator& O)
{
  RealVector TmpV2 (O.Collection->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [O.Collection->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Collection->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(O.Collection->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2[i - 1] = 0.0;
	}
      TmpV2[i] = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < O.Collection->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < O.Collection->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i];
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

MathematicaOutput& operator << (MathematicaOutput& Str, ChiralExchangeOperator& O)
{
  RealVector TmpV2 (O.Collection->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [O.Collection->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Collection->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(O.Collection->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2[i - 1] = 0.0;
	}
      TmpV2[i] = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (O.Collection->GetHilbertSpaceDimension() - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (O.Collection->GetHilbertSpaceDimension() - 1); ++j)
	{
	  Str << TmpV[j][i];
	  Str << ",";
	}
      Str << TmpV[O.Collection->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (O.Collection->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][O.Collection->GetHilbertSpaceDimension() - 1];
      Str << ",";
    }
  Str << TmpV[O.Collection->GetHilbertSpaceDimension() - 1][O.Collection->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}


