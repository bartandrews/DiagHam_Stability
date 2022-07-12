////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of Explicit hamiltonian                       //
//                                                                            //
//                        last modification : 23/05/2001                      //
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


#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Complex.h" 
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


// contructor from default datas without boundary operators
//
// hilbertSpace = Hilnbert associated to the Hamiltonian
// hamiltonian = matrix representation of the Hamiltonian

ExplicitHamiltonian::ExplicitHamiltonian(AbstractHilbertSpace* hilbertSpace, Matrix* hamiltonian)
{
  this->HilbertSpace = hilbertSpace->Clone();
  this->Hamiltonian = hamiltonian->Clone();
}

// contructor from default datas
//
// hilbertSpace = Hilnbert associated to the Hamiltonian
// hamiltonian = matrix representation of the Hamiltonian
// listRightInteractionOperator = List of operators used to describe right interaction
// listLeftInteractionOperator = List of operators used to describe left interaction

ExplicitHamiltonian::ExplicitHamiltonian(AbstractHilbertSpace* hilbertSpace, Matrix* hamiltonian, 
					 const List<Matrix*>& listRightInteractionOperator, 
					 const List<Matrix*>& listLeftInteractionOperator)
{
  this->HilbertSpace = hilbertSpace->Clone();
  this->Hamiltonian = hamiltonian->Clone();
  this->ListRightInteractionOperator = listRightInteractionOperator;
  this->ListLeftInteractionOperator = listLeftInteractionOperator;
}

// destructor
//

ExplicitHamiltonian::~ExplicitHamiltonian() 
{
  delete this->HilbertSpace;
  delete this->Hamiltonian;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ExplicitHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) 
{
  this->HilbertSpace = hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ExplicitHamiltonian::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ExplicitHamiltonian::GetHilbertSpaceDimension () 
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ExplicitHamiltonian::ShiftHamiltonian (double shift) 
{
  for (int i = 0; i < this->Hamiltonian->GetNbrRow(); i++)
    (*(this->Hamiltonian))(i, i) += shift;
}

// return matrix representation of current Hamiltonian
//
// return value = reference to representation

Matrix* ExplicitHamiltonian::GetHamiltonian ()
{
  return this->Hamiltonian;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ExplicitHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ExplicitHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vector where result has been stored

RealVector& ExplicitHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return vDestination.Multiply(*(this->Hamiltonian), vSource);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ExplicitHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						  int firstComponent, int nbrComponent) 
{
  return vDestination.Multiply(*(this->Hamiltonian), vSource, firstComponent, 1, nbrComponent, firstComponent, 1);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& ExplicitHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return vDestination.AddMultiply(*(this->Hamiltonian), vSource);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ExplicitHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return vDestination.AddMultiply(*(this->Hamiltonian), vSource, firstComponent, 1, nbrComponent, firstComponent, 1);
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& ExplicitHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination.Multiply(*(this->Hamiltonian), vSource);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ExplicitHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  return vDestination.Multiply(*(this->Hamiltonian), vSource, firstComponent, 1, nbrComponent, firstComponent, 1);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& ExplicitHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return vDestination.AddMultiply(*(this->Hamiltonian), vSource);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
ComplexVector& ExplicitHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							int firstComponent, int nbrComponent)
{
  return vDestination.AddMultiply(*(this->Hamiltonian), vSource, firstComponent, 1, nbrComponent, firstComponent, 1);
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ExplicitHamiltonian::LeftInteractionOperators() 
{
  return this->ListLeftInteractionOperator;
}


// return a list of right interaction operators 
//
// return value = list of right interaction operators

List<Matrix*> ExplicitHamiltonian::RightInteractionOperators() 
{
  return this->ListRightInteractionOperator;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ExplicitHamiltonian& H) 
{
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, ExplicitHamiltonian& H) 
{
  return Str;
}
