////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of Wigner small D matrix (i.e. Wigner D matrix restricted       //
//                        to rotations along the y axis)                      //
//                                                                            //
//                        last modification : 03/04/2009                      //
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
#include "MathTools/WignerSmallDMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"

// constructor 
//
// jValue =  twice the angular momentum of the sector the matrix is acting on

WignerSmallDMatrix::WignerSmallDMatrix(int jValue)
{
  this->JValue = jValue;
  this->EvaluateWignerSmallDMatrix();
}

// destructor
//

WignerSmallDMatrix::~WignerSmallDMatrix()
{
}

// evaluate all coefficients to compute the Wigner matrix
//

void WignerSmallDMatrix::EvaluateWignerSmallDMatrix()
{
  FactorialCoefficient TmpFactorial;
  BinomialCoefficients TmpBinomial(this->JValue);
  this->Coefficients = new double** [this->JValue + 1];
  this->JMin = new int* [this->JValue + 1];
  for (int i = 0; i <= this->JValue; ++i)
    {
      this->JMin[i] = new int[i + 1];	  
      this->Coefficients[i] = new double*[i + 1];
      for (int j = 0; j <= i; ++j)
	{
	  if (j <= (this->JValue - i))
	    this->JMin[i][j] = j;
	  else
	    this->JMin[i][j] = this->JValue - i;
	  this->Coefficients[i][j] = new double[this->JMin[i][j] + 1];
	  TmpFactorial.SetToOne();
	  TmpFactorial.FactorialMultiply(i);
	  TmpFactorial.FactorialDivide(j);
	  TmpFactorial.FactorialMultiply(this->JValue - i);
	  TmpFactorial.FactorialDivide(this->JValue - j);
	  double GlobalCoefficient = sqrt(TmpFactorial.GetNumericalValue());
	  if (((i -j) & 1) != 0)
	    GlobalCoefficient *= -1.0;
	  for (int k = 0; k <= this->JMin[i][j]; ++k)
	    {
	      this->Coefficients[i][j][k] = GlobalCoefficient * TmpBinomial(j, k)  * TmpBinomial(this->JValue - j, (i - j) + k) ;
	      GlobalCoefficient *= -1.0;
	    }
	}
    }
}
