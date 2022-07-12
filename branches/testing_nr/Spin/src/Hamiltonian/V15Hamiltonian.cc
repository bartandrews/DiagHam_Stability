////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of V_15 hamiltonian                         //
//                                                                            //
//                        last modification : 05/03/2001                      //
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


#include "Hamiltonian/V15Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "MathTools/Complex.h"
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
// j4 = coupling constant J4
// j5 = coupling constant J5

V15Hamiltonian::V15Hamiltonian(const Spin1_2Chain& chain, double j1, double j2, double j3, double j4, double j5) 
{
  this->Chain = chain;
  this->NbrLink = 36;
  this->J1 = j1;
  this->HalfJ1 = 0.5 * this->J1;
  this->J2 = j2;
  this->HalfJ2 = 0.5 * this->J2;
  this->J3 = j3;
  this->HalfJ3 = 0.5 * this->J3;
  this->J4 = j4;
  this->HalfJ4 = 0.5 * this->J4;
  this->J5 = j5;
  this->HalfJ5 = 0.5 * this->J5;
  this->SzSzContributions = new double [this->Chain.GetHilbertSpaceDimension()];
  this->MatrixElementIndex = new int* [this->NbrLink];
  this->MatrixElements = new double* [this->NbrLink];
  for (int i = 0; i < this->NbrLink; i++)
    {
      this->MatrixElementIndex[i] = new int [this->Chain.GetHilbertSpaceDimension()];
      this->MatrixElements[i] = new double [this->Chain.GetHilbertSpaceDimension()];
    }
  this->EvaluateMatrixElements();
}

// destructor
//

V15Hamiltonian::~V15Hamiltonian() 
{
  delete[] this->SzSzContributions;
  for (int i = 0; i < this->NbrLink; i++)
    {
      delete[] this->MatrixElementIndex[i];
      delete[] this->MatrixElements[i];
    }
  delete[] this->MatrixElementIndex;
  delete[] this->MatrixElements;  
}

// set chain
// 
// chain = reference on Hilbert space of the associated system
// return value = reference on current Hamiltonian

V15Hamiltonian& V15Hamiltonian::SetChain(const Spin1_2Chain& chain)
{  
  delete[] this->SzSzContributions;
  for (int i = 0; i < this->NbrLink; i++)
    {
      delete[] this->MatrixElementIndex[i];
      delete[] this->MatrixElements[i];
    }
  delete[] this->MatrixElementIndex;
  delete[] this->MatrixElements;  
  this->Chain = chain;
  this->SzSzContributions = new double [this->Chain.GetHilbertSpaceDimension()];
  this->MatrixElementIndex = new int* [this->NbrLink];
  this->MatrixElements = new double* [this->NbrLink];
  for (int i = 0; i < this->NbrLink; i++)
    {
      this->MatrixElementIndex[i] = new int [this->Chain.GetHilbertSpaceDimension()];
      this->MatrixElements[i] = new double [this->Chain.GetHilbertSpaceDimension()];
    }
  this->EvaluateMatrixElements();
  return *this;  
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* V15Hamiltonian::GetHilbertSpace ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void V15Hamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int V15Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain.GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void V15Hamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain.GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex V15Hamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Chain.GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {

      // J3 part of Hamiltonian      
      x += V1[MatrixElementIndex[1][i]] * MatrixElements[1][i] * V2[i];
      x += V1[MatrixElementIndex[3][i]] * MatrixElements[3][i] * V2[i];
      x += V1[MatrixElementIndex[5][i]] * MatrixElements[5][i] * V2[i];
      x += V1[MatrixElementIndex[19][i]] * MatrixElements[19][i] * V2[i];
      x += V1[MatrixElementIndex[21][i]] * MatrixElements[21][i] * V2[i];
      x += V1[MatrixElementIndex[23][i]] * MatrixElements[23][i] * V2[i];
      
      // J5 part of Hamiltonian      
      x += V1[MatrixElementIndex[7][i]] * MatrixElements[7][i] * V2[i];
      x += V1[MatrixElementIndex[11][i]] * MatrixElements[11][i] * V2[i];
      x += V1[MatrixElementIndex[15][i]] * MatrixElements[15][i] * V2[i];
      x += V1[MatrixElementIndex[8][i]] * MatrixElements[8][i] * V2[i];
      x += V1[MatrixElementIndex[12][i]] * MatrixElements[12][i] * V2[i];
      x += V1[MatrixElementIndex[16][i]] * MatrixElements[16][i] * V2[i];

      // J1 part of Hamiltonian      
      x += V1[MatrixElementIndex[0][i]] * MatrixElements[0][i] * V2[i];
      x += V1[MatrixElementIndex[2][i]] * MatrixElements[2][i] * V2[i];
      x += V1[MatrixElementIndex[4][i]] * MatrixElements[4][i] * V2[i];
      x += V1[MatrixElementIndex[18][i]] * MatrixElements[18][i] * V2[i];
      x += V1[MatrixElementIndex[20][i]] * MatrixElements[20][i] * V2[i];
      x += V1[MatrixElementIndex[22][i]] * MatrixElements[22][i] * V2[i];

      // J4 part of Hamiltonian      
      x += V1[MatrixElementIndex[6][i]] * MatrixElements[6][i] * V2[i];
      x += V1[MatrixElementIndex[10][i]] * MatrixElements[10][i] * V2[i];
      x += V1[MatrixElementIndex[14][i]] * MatrixElements[14][i] * V2[i];
      x += V1[MatrixElementIndex[9][i]] * MatrixElements[9][i] * V2[i];
      x += V1[MatrixElementIndex[13][i]] * MatrixElements[13][i] * V2[i];
      x += V1[MatrixElementIndex[17][i]] * MatrixElements[17][i] * V2[i];

      // J2 part of Hamiltonian      
      x += V1[MatrixElementIndex[24][i]] * MatrixElements[24][i] * V2[i];
      x += V1[MatrixElementIndex[25][i]] * MatrixElements[25][i] * V2[i];
      x += V1[MatrixElementIndex[26][i]] * MatrixElements[26][i] * V2[i];
      x += V1[MatrixElementIndex[27][i]] * MatrixElements[27][i] * V2[i];
      x += V1[MatrixElementIndex[28][i]] * MatrixElements[28][i] * V2[i];
      x += V1[MatrixElementIndex[29][i]] * MatrixElements[29][i] * V2[i];
      x += V1[MatrixElementIndex[30][i]] * MatrixElements[30][i] * V2[i];
      x += V1[MatrixElementIndex[31][i]] * MatrixElements[31][i] * V2[i];
      x += V1[MatrixElementIndex[32][i]] * MatrixElements[32][i] * V2[i];
      x += V1[MatrixElementIndex[33][i]] * MatrixElements[33][i] * V2[i];
      x += V1[MatrixElementIndex[34][i]] * MatrixElements[34][i] * V2[i];
      x += V1[MatrixElementIndex[35][i]] * MatrixElements[35][i] * V2[i];

      // SzSz contributions
      x += V1[i] * this->SzSzContributions[i] * V2[i];
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex V15Hamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& V15Hamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  int dim = this->Chain.GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    vDestination[i] = this->SzSzContributions[i] * vSource[i];
  for (int i = 0; i < dim; i++)
    {

      // J3 part of Hamiltonian      
      vDestination[MatrixElementIndex[1][i]] += MatrixElements[1][i] * vSource[i];
      vDestination[MatrixElementIndex[3][i]] += MatrixElements[3][i] * vSource[i];
      vDestination[MatrixElementIndex[5][i]] += MatrixElements[5][i] * vSource[i];
      vDestination[MatrixElementIndex[19][i]] += MatrixElements[19][i] * vSource[i];
      vDestination[MatrixElementIndex[21][i]] += MatrixElements[21][i] * vSource[i];
      vDestination[MatrixElementIndex[23][i]] += MatrixElements[23][i] * vSource[i];

      // J5 part of Hamiltonian      
      vDestination[MatrixElementIndex[7][i]] += MatrixElements[7][i] * vSource[i];
      vDestination[MatrixElementIndex[11][i]] += MatrixElements[11][i] * vSource[i];
      vDestination[MatrixElementIndex[15][i]] += MatrixElements[15][i] * vSource[i];
      vDestination[MatrixElementIndex[8][i]] += MatrixElements[8][i] * vSource[i];
      vDestination[MatrixElementIndex[12][i]] += MatrixElements[12][i] * vSource[i];
      vDestination[MatrixElementIndex[16][i]] += MatrixElements[16][i] * vSource[i];

      // J1 part of Hamiltonian      
      vDestination[MatrixElementIndex[0][i]] += MatrixElements[0][i] * vSource[i];
      vDestination[MatrixElementIndex[2][i]] += MatrixElements[2][i] * vSource[i];
      vDestination[MatrixElementIndex[4][i]] += MatrixElements[4][i] * vSource[i];
      vDestination[MatrixElementIndex[18][i]] += MatrixElements[18][i] * vSource[i];
      vDestination[MatrixElementIndex[20][i]] += MatrixElements[20][i] * vSource[i];
      vDestination[MatrixElementIndex[22][i]] += MatrixElements[22][i] * vSource[i];

      // J4 part of Hamiltonian      
      vDestination[MatrixElementIndex[6][i]] += MatrixElements[6][i] * vSource[i];
      vDestination[MatrixElementIndex[10][i]] += MatrixElements[10][i] * vSource[i];
      vDestination[MatrixElementIndex[14][i]] += MatrixElements[14][i] * vSource[i];
      vDestination[MatrixElementIndex[9][i]] += MatrixElements[9][i] * vSource[i];
      vDestination[MatrixElementIndex[13][i]] += MatrixElements[13][i] * vSource[i];
      vDestination[MatrixElementIndex[17][i]] += MatrixElements[17][i] * vSource[i];

      // J2 part of Hamiltonian      
      vDestination[MatrixElementIndex[24][i]] += MatrixElements[24][i] * vSource[i];
      vDestination[MatrixElementIndex[25][i]] += MatrixElements[25][i] * vSource[i];
      vDestination[MatrixElementIndex[26][i]] += MatrixElements[26][i] * vSource[i];
      vDestination[MatrixElementIndex[27][i]] += MatrixElements[27][i] * vSource[i];
      vDestination[MatrixElementIndex[28][i]] += MatrixElements[28][i] * vSource[i];
      vDestination[MatrixElementIndex[29][i]] += MatrixElements[29][i] * vSource[i];
      vDestination[MatrixElementIndex[30][i]] += MatrixElements[30][i] * vSource[i];
      vDestination[MatrixElementIndex[31][i]] += MatrixElements[31][i] * vSource[i];
      vDestination[MatrixElementIndex[32][i]] += MatrixElements[32][i] * vSource[i];
      vDestination[MatrixElementIndex[33][i]] += MatrixElements[33][i] * vSource[i];
      vDestination[MatrixElementIndex[34][i]] += MatrixElements[34][i] * vSource[i];
      vDestination[MatrixElementIndex[35][i]] += MatrixElements[35][i] * vSource[i];
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& V15Hamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// evaluate all matrix elements
// 

void V15Hamiltonian::EvaluateMatrixElements()
{
  int dim = this->Chain.GetHilbertSpaceDimension();
  int pos;

  // SzSz part
  for (int i = 0; i < dim; i++)
    {

      // SzSz part
      SzSzContributions[i] = this->J1 * this->Chain.SziSzj(0, 1, i);
      SzSzContributions[i] += this->J3 * this->Chain.SziSzj(1, 2, i);
      SzSzContributions[i] += this->J1 * this->Chain.SziSzj(2, 3, i);
      SzSzContributions[i] += this->J3 * this->Chain.SziSzj(3, 4, i);
      SzSzContributions[i] += this->J1 * this->Chain.SziSzj(4, 5, i);
      SzSzContributions[i] += this->J3 * this->Chain.SziSzj(5, 0, i);
      SzSzContributions[i] += this->J1 * this->Chain.SziSzj(9, 10, i);
      SzSzContributions[i] += this->J3 * this->Chain.SziSzj(10, 11, i);
      SzSzContributions[i] += this->J1 * this->Chain.SziSzj(11, 12, i);
      SzSzContributions[i] += this->J3 * this->Chain.SziSzj(12, 13, i);
      SzSzContributions[i] += this->J1 * this->Chain.SziSzj(13, 14, i);
      SzSzContributions[i] += this->J3 * this->Chain.SziSzj(14, 9, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(0, 2, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(2, 4, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(4, 0, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(1, 3, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(3, 5, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(5, 1, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(9, 11, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(11, 13, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(13, 9, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(10, 12, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(12, 14, i);
      SzSzContributions[i] += this->J2 * this->Chain.SziSzj(14, 10, i);
      SzSzContributions[i] += this->J4 * this->Chain.SziSzj(0, 6, i);
      SzSzContributions[i] += this->J5 * this->Chain.SziSzj(1, 6, i);
      SzSzContributions[i] += this->J4 * this->Chain.SziSzj(2, 7, i);
      SzSzContributions[i] += this->J5 * this->Chain.SziSzj(3, 7, i);
      SzSzContributions[i] += this->J4 * this->Chain.SziSzj(4, 8, i);
      SzSzContributions[i] += this->J5 * this->Chain.SziSzj(5, 8, i);
      SzSzContributions[i] += this->J5 * this->Chain.SziSzj(9, 6, i);
      SzSzContributions[i] += this->J4 * this->Chain.SziSzj(10, 6, i);
      SzSzContributions[i] += this->J5 * this->Chain.SziSzj(11, 7, i);
      SzSzContributions[i] += this->J4 * this->Chain.SziSzj(12, 7, i);
      SzSzContributions[i] += this->J5 * this->Chain.SziSzj(13, 8, i);
      SzSzContributions[i] += this->J4 * this->Chain.SziSzj(14, 8, i);

      // J1 contribution
      pos = this->Chain.Pij(0, 1, i);
      if (pos == dim)
	{
	  MatrixElementIndex[0][i] = 0;
	  MatrixElements[0][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[0][i] = pos;
	  MatrixElements[0][i] = this->HalfJ1;	  
	}
      pos = this->Chain.Pij(2, 3, i);
      if (pos == dim)
	{
	  MatrixElementIndex[2][i] = 0;
	  MatrixElements[2][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[2][i] = pos;
	  MatrixElements[2][i] = this->HalfJ1;	  
	}
      pos = this->Chain.Pij(4, 5, i);
      if (pos == dim)
	{
	  MatrixElementIndex[4][i] = 0;
	  MatrixElements[4][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[4][i] = pos;
	  MatrixElements[4][i] = this->HalfJ1;	  
	}
      pos = this->Chain.Pij(9, 10, i);
      if (pos == dim)
	{
	  MatrixElementIndex[18][i] = 0;
	  MatrixElements[18][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[18][i] = pos;
	  MatrixElements[18][i] = this->HalfJ1;	  
	}
      pos = this->Chain.Pij(11, 12, i);
      if (pos == dim)
	{
	  MatrixElementIndex[20][i] = 0;
	  MatrixElements[20][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[20][i] = pos;
	  MatrixElements[20][i] = this->HalfJ1;	  
	}
      pos = this->Chain.Pij(13, 14, i);
      if (pos == dim)
	{
	  MatrixElementIndex[22][i] = 0;
	  MatrixElements[22][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[22][i] = pos;
	  MatrixElements[22][i] = this->HalfJ1;	  
	}


      // J2 contribution
      pos = this->Chain.Pij(0, 2, i);
      if (pos == dim)
	{
	  MatrixElementIndex[24][i] = 0;
	  MatrixElements[24][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[24][i] = pos;
	  MatrixElements[24][i] = this->HalfJ2;	  
	}
      pos = this->Chain.Pij(2, 4, i);
      if (pos == dim)
	{
	  MatrixElementIndex[25][i] = 0;
	  MatrixElements[25][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[25][i] = pos;
	  MatrixElements[25][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(4, 0, i);
      if (pos == dim)
	{
	  MatrixElementIndex[26][i] = 0;
	  MatrixElements[26][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[26][i] = pos;
	  MatrixElements[26][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(1, 3, i);
      if (pos == dim)
	{
	  MatrixElementIndex[27][i] = 0;
	  MatrixElements[27][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[27][i] = pos;
	  MatrixElements[27][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(3, 5, i);
      if (pos == dim)
	{
	  MatrixElementIndex[28][i] = 0;
	  MatrixElements[28][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[28][i] = pos;
	  MatrixElements[28][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(5, 1, i);
      if (pos == dim)
	{
	  MatrixElementIndex[29][i] = 0;
	  MatrixElements[29][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[29][i] = pos;
	  MatrixElements[29][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(9, 11, i);
      if (pos == dim)
	{
	  MatrixElementIndex[30][i] = 0;
	  MatrixElements[30][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[30][i] = pos;
	  MatrixElements[30][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(11, 13, i);
      if (pos == dim)
	{
	  MatrixElementIndex[31][i] = 0;
	  MatrixElements[31][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[31][i] = pos;
	  MatrixElements[31][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(13, 9, i);
      if (pos == dim)
	{
	  MatrixElementIndex[32][i] = 0;
	  MatrixElements[32][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[32][i] = pos;
	  MatrixElements[32][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(10, 12, i);
      if (pos == dim)
	{
	  MatrixElementIndex[33][i] = 0;
	  MatrixElements[33][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[33][i] = pos;
	  MatrixElements[33][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(12, 14, i);
      if (pos == dim)
	{
	  MatrixElementIndex[34][i] = 0;
	  MatrixElements[34][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[34][i] = pos;
	  MatrixElements[34][i] = this->HalfJ2;	  
	}

      pos = this->Chain.Pij(14, 10, i);
      if (pos == dim)
	{
	  MatrixElementIndex[35][i] = 0;
	  MatrixElements[35][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[35][i] = pos;
	  MatrixElements[35][i] = this->HalfJ2;	  
	}

      // J4 contribution
      pos = this->Chain.Pij(0, 6, i);
      if (pos == dim)
	{
	  MatrixElementIndex[6][i] = 0;
	  MatrixElements[6][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[6][i] = pos;
	  MatrixElements[6][i] = this->HalfJ4;	  
	}

      pos = this->Chain.Pij(2, 7, i);
      if (pos == dim)
	{
	  MatrixElementIndex[10][i] = 0;
	  MatrixElements[10][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[10][i] = pos;
	  MatrixElements[10][i] = this->HalfJ4;	  
	}

      pos = this->Chain.Pij(4, 8, i);
      if (pos == dim)
	{
	  MatrixElementIndex[14][i] = 0;
	  MatrixElements[14][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[14][i] = pos;
	  MatrixElements[14][i] = this->HalfJ4;	  
	}

      pos = this->Chain.Pij(10, 6, i);
      if (pos == dim)
	{
	  MatrixElementIndex[9][i] = 0;
	  MatrixElements[9][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[9][i] = pos;
	  MatrixElements[9][i] = this->HalfJ4;	  
	}

      pos = this->Chain.Pij(12, 7, i);
      if (pos == dim)
	{
	  MatrixElementIndex[13][i] = 0;
	  MatrixElements[13][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[13][i] = pos;
	  MatrixElements[13][i] = this->HalfJ4;	  
	}

      pos = this->Chain.Pij(14, 8, i);
      if (pos == dim)
	{
	  MatrixElementIndex[17][i] = 0;
	  MatrixElements[17][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[17][i] = pos;
	  MatrixElements[17][i] = this->HalfJ4;	  
	}

      // J3 contribution
      pos = this->Chain.Pij(1, 2, i);
      if (pos == dim)
	{
	  MatrixElementIndex[1][i] = 0;
	  MatrixElements[1][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[1][i] = pos;
	  MatrixElements[1][i] = this->HalfJ3;	  
	}
      pos = this->Chain.Pij(3, 4, i);
      if (pos == dim)
	{
	  MatrixElementIndex[3][i] = 0;
	  MatrixElements[3][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[3][i] = pos;
	  MatrixElements[3][i] = this->HalfJ3;	  
	}
      pos = this->Chain.Pij(5, 0, i);
      if (pos == dim)
	{
	  MatrixElementIndex[5][i] = 0;
	  MatrixElements[5][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[5][i] = pos;
	  MatrixElements[5][i] = this->HalfJ3;	  
	}
      pos = this->Chain.Pij(10, 11, i);
      if (pos == dim)
	{
	  MatrixElementIndex[19][i] = 0;
	  MatrixElements[19][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[19][i] = pos;
	  MatrixElements[19][i] = this->HalfJ3;	  
	}
      pos = this->Chain.Pij(12, 13, i);
      if (pos == dim)
	{
	  MatrixElementIndex[21][i] = 0;
	  MatrixElements[21][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[21][i] = pos;
	  MatrixElements[21][i] = this->HalfJ3;	  
	}
      pos = this->Chain.Pij(14, 9, i);
      if (pos == dim)
	{
	  MatrixElementIndex[23][i] = 0;
	  MatrixElements[23][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[23][i] = pos;
	  MatrixElements[23][i] = this->HalfJ3;	  
	}

      // J5 contribution
      pos = this->Chain.Pij(1, 6, i);
      if (pos == dim)
	{
	  MatrixElementIndex[7][i] = 0;
	  MatrixElements[7][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[7][i] = pos;
	  MatrixElements[7][i] = this->HalfJ5;	  
	}
      pos = this->Chain.Pij(3, 7, i);
      if (pos == dim)
	{
	  MatrixElementIndex[11][i] = 0;
	  MatrixElements[11][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[11][i] = pos;
	  MatrixElements[11][i] = this->HalfJ5;	  
	}
      pos = this->Chain.Pij(5, 8, i);
      if (pos == dim)
	{
	  MatrixElementIndex[15][i] = 0;
	  MatrixElements[15][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[15][i] = pos;
	  MatrixElements[15][i] = this->HalfJ5;	  
	}
      pos = this->Chain.Pij(9, 6, i);
      if (pos == dim)
	{
	  MatrixElementIndex[8][i] = 0;
	  MatrixElements[8][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[8][i] = pos;
	  MatrixElements[8][i] = this->HalfJ5;	  
	}
      pos = this->Chain.Pij(11, 7, i);
      if (pos == dim)
	{
	  MatrixElementIndex[12][i] = 0;
	  MatrixElements[12][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[12][i] = pos;
	  MatrixElements[12][i] = this->HalfJ5;	  
	}
      pos = this->Chain.Pij(13, 8, i);
      if (pos == dim)
	{
	  MatrixElementIndex[16][i] = 0;
	  MatrixElements[16][i] = 0.0;	  
	}
      else
	{
	  MatrixElementIndex[16][i] = pos;
	  MatrixElements[16][i] = this->HalfJ5;	  
	}
    }
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, V15Hamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, V15Hamiltonian& H) 
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

