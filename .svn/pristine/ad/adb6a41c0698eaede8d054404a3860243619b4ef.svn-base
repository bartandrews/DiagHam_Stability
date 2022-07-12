////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of Mn_12 hamiltonian                        //
//                                                                            //
//                        last modification : 28/02/2001                      //
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


#include "Hamiltonian/Mn12Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::endl;
using std::ostream;


// contructor from default datas
//
// chain = reference on Hilbert space of the associated system
// j1 = coupling constant J1
// j2 = coupling constant J2
// j3 = coupling constant J3
// j3p = coupling constant J3p
// j4 = coupling constant J4

Mn12Hamiltonian::Mn12Hamiltonian(const ThierryChain& chain, double j1, double j2, double j3, double j3p, double j4) 
{
  this->Chain = chain;
  this->J1 = j1;
  this->HalfJ1 = 0.5 * this->J1;
  this->J2 = j2;
  this->HalfJ2 = 0.5 * this->J2;
  this->J3 = j3;
  this->HalfJ3 = 0.5 * this->J3;
  this->J3p = j3p;
  this->HalfJ3p = 0.5 * this->J3p;
  this->J4 = j4;
  this->HalfJ4 = 0.5 * this->J4;
  this->SzSzContributions = new double [this->Chain.GetHilbertSpaceDimension()];
  this->EvaluateSzSzContributions();
}

// destructor
//

Mn12Hamiltonian::~Mn12Hamiltonian() 
{
  delete[] this->SzSzContributions;
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

Mn12Hamiltonian& Mn12Hamiltonian::SetChain(const ThierryChain& chain)
{
  this->Chain = chain;
  delete[] this->SzSzContributions;
  this->SzSzContributions = new double [this->Chain.GetHilbertSpaceDimension()];
  this->EvaluateSzSzContributions();
  return *this;  
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* Mn12Hamiltonian::GetHilbertSpace ()
{
  return 0;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension
int Mn12Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain.GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void Mn12Hamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain.GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex Mn12Hamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  double coef;
  int pos;
  int dim = V1.GetVectorDimension();
  for (int i = 0; i < dim; i++)
    {
      x += V1[i] * this->SzSzContributions[i] * V2[i];

// J3 part of Hamiltonian      
      pos = this->Chain.SmiSpj(8, 9, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(9, 8, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(9, 10, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(10, 9, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(10, 11, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(11, 10, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(11, 8, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(8, 11, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3 * coef * V2[i];
      pos = this->Chain.SmiSpj(8, 10, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3p * coef * V2[i];
      pos = this->Chain.SmiSpj(10, 8, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3p * coef * V2[i];
      pos = this->Chain.SmiSpj(9, 11, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3p * coef * V2[i];
      pos = this->Chain.SmiSpj(11, 9, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ3p * coef * V2[i];

// J1 part of Hamiltonian      
      pos = this->Chain.SmiSpj(0, 8, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(8, 0, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(2, 9, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(9, 2, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(4, 10, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(10, 4, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(6, 11, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];
      pos = this->Chain.SmiSpj(11, 6, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ1 * coef * V2[i];

// J4 part of Hamiltonian      
      pos = this->Chain.SmiSpj(0, 1, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(1, 0, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(1, 2, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(2, 1, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(2, 3, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(3, 2, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(3, 4, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(4, 3, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
     pos = this->Chain.SmiSpj(4, 5, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(5, 4, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(5, 6, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(6, 5, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(6, 7, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(7, 6, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(7, 0, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];
      pos = this->Chain.SmiSpj(0, 7, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ4 * coef * V2[i];

// J2 part of Hamiltonian      
      pos = this->Chain.SmiSpj(1, 8, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(8, 1, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(1, 9, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(9, 1, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(3, 9, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(9, 3, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(3, 10, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(10, 3, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
     pos = this->Chain.SmiSpj(5, 10, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(10, 5, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(5, 11, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(11, 5, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(7, 11, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(11, 7, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(7, 8, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
      pos = this->Chain.SmiSpj(8, 7, i, &coef);
      if (pos != dim)
	x += V1[pos] * this->HalfJ2 * coef * V2[i];
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex Mn12Hamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& Mn12Hamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  int dim = this->Chain.GetHilbertSpaceDimension();
  double coef;
  int pos;
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {

// J3 part of Hamiltonian      
      pos = this->Chain.SmiSpj(8, 9, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(9, 8, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(9, 10, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(10, 9, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(10, 11, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(11, 10, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(11, 8, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(8, 11, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3 * coef * vSource[i];
      pos = this->Chain.SmiSpj(8, 10, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3p * coef * vSource[i];
      pos = this->Chain.SmiSpj(10, 8, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3p * coef * vSource[i];
      pos = this->Chain.SmiSpj(9, 11, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3p * coef * vSource[i];
      pos = this->Chain.SmiSpj(11, 9, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ3p * coef * vSource[i];

// J1 part of Hamiltonian      
      pos = this->Chain.SmiSpj(0, 8, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(8, 0, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(2, 9, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(9, 2, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(4, 10, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(10, 4, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(6, 11, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];
      pos = this->Chain.SmiSpj(11, 6, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ1 * coef * vSource[i];

// J4 part of Hamiltonian      
      pos = this->Chain.SmiSpj(0, 1, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(1, 0, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(1, 2, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(2, 1, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(2, 3, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(3, 2, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(3, 4, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(4, 3, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
     pos = this->Chain.SmiSpj(4, 5, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(5, 4, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(5, 6, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(6, 5, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(6, 7, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(7, 6, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(7, 0, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];
      pos = this->Chain.SmiSpj(0, 7, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ4 * coef * vSource[i];

// J2 part of Hamiltonian      
      pos = this->Chain.SmiSpj(1, 8, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(8, 1, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(1, 9, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(9, 1, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(3, 9, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(9, 3, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(3, 10, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(10, 3, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
     pos = this->Chain.SmiSpj(5, 10, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(10, 5, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(5, 11, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(11, 5, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(7, 11, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(11, 7, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(7, 8, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
      pos = this->Chain.SmiSpj(8, 7, i, &coef);
      if (pos != dim)
	vDestination[pos] += this->HalfJ2 * coef * vSource[i];
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vector where result has been stored

RealVector& Mn12Hamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
					      int firstComponent, int nbrComponent) 
{
  return this->LowLevelMultiply(vSource, vDestination);
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& Mn12Hamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// evaluate contribution of SziSzj terms
//

void Mn12Hamiltonian::EvaluateSzSzContributions()
{
  for (int i = 0; i < this->Chain.GetHilbertSpaceDimension(); i++)
    {

      // contribution of J3 terms
      this->SzSzContributions[i] = this->J3 * this->Chain.SziSzj(8, 9, i);
      this->SzSzContributions[i] += this->J3 * this->Chain.SziSzj(9, 10, i);
      this->SzSzContributions[i] += this->J3 * this->Chain.SziSzj(10, 11, i);
      this->SzSzContributions[i] += this->J3 * this->Chain.SziSzj(11, 8, i);
      this->SzSzContributions[i] += this->J3p * this->Chain.SziSzj(8, 10, i);
      this->SzSzContributions[i] += this->J3p * this->Chain.SziSzj(9, 11, i);
  
      // contribution of J1 terms
      this->SzSzContributions[i] += this->J1 * this->Chain.SziSzj(0, 8, i);
      this->SzSzContributions[i] += this->J1 * this->Chain.SziSzj(2, 9, i);
      this->SzSzContributions[i] += this->J1 * this->Chain.SziSzj(4, 10, i);
      this->SzSzContributions[i] += this->J1 * this->Chain.SziSzj(6, 11, i);

      // contribution of J4 terms
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(0, 1, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(1, 2, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(2, 3, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(3, 4, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(4, 5, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(5, 6, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(6, 7, i);
      this->SzSzContributions[i] += this->J4 * this->Chain.SziSzj(7, 0, i);

      // contribution of J2 terms
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(1, 8, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(1, 9, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(3, 9, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(3, 10, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(5, 10, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(5, 11, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(7, 11, i);
      this->SzSzContributions[i] += this->J2 * this->Chain.SziSzj(7, 8, i);
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, Mn12Hamiltonian& H) 
{
  RealVector TmpV2 (H.Chain.GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chain.GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain.GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chain.GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Chain.GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Chain.GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i] << "    ";
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

MathematicaOutput& operator << (MathematicaOutput& Str, Mn12Hamiltonian& H) 
{
  RealVector TmpV2 (H.Chain.GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Chain.GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Chain.GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Chain.GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Chain.GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Chain.GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Chain.GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Chain.GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Chain.GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Chain.GetHilbertSpaceDimension() - 1][H.Chain.GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

