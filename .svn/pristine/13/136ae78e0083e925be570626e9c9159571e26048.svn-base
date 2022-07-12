////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of non-local occupation number operator for the Potts models      //
//                                                                            //
//                        last modification : 05/12/2014                      //
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


#include "Operator/PottsNonLocalOccupationNumberOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// chain = pointer to the Hilbert space
// nbrSpin = number of spins

PottsNonLocalOccupationNumberOperator::PottsNonLocalOccupationNumberOperator(Spin1_2Chain* chain, int nbrSpin)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
}

// destructor
//

PottsNonLocalOccupationNumberOperator::~PottsNonLocalOccupationNumberOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* PottsNonLocalOccupationNumberOperator::Clone ()
{
  return new PottsNonLocalOccupationNumberOperator(this->Chain, this->NbrSpin);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PottsNonLocalOccupationNumberOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (Spin1_2Chain*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PottsNonLocalOccupationNumberOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int PottsNonLocalOccupationNumberOperator::GetHilbertSpaceDimension ()
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

Complex PottsNonLocalOccupationNumberOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, 
								     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      Tmp += V1[i] * V2[i] * (1.0 - 2.0 * ((double) this->Chain->Parity(i)));
    }
  return Complex(Tmp);
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& PottsNonLocalOccupationNumberOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      vDestination[i] = vSource[i] * (1.0 - 2.0 * ((double) this->Chain->Parity(i)));
    }
  return vDestination;
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex PottsNonLocalOccupationNumberOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
								     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Tmp = 0.0;
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      Tmp += Conj(V1[i]) * V2[i] * (1.0 - 2.0 * ((double) this->Chain->Parity(i)));
    }
  return Complex(Tmp);
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PottsNonLocalOccupationNumberOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								       int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  double Tmp = 0.0;
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
      vDestination[i] = vSource[i] * (1.0 - 2.0 * ((double) this->Chain->Parity(i)));
    }
  return vDestination;
}

