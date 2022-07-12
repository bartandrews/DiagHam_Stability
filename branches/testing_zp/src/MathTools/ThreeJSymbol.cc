////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of clebsch gordan coefficients                   //
//                                                                            //
//                        last modification : 19/06/2002                      //
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
#include "MathTools/ThreeJSymbol.h"

#include <stdlib.h>
#include <math.h>


// default constructor
//

ThreeJSymbol::ThreeJSymbol() : ClebschGordanCoefficients()
{
}


// constructor 
//
// j1 = first angular momentum (twice the value to avoid half integer value)
// j2 = second angular momentum (twice the value to avoid half integer value)

ThreeJSymbol::ThreeJSymbol(int j1, int j2) : ClebschGordanCoefficients(j1,j2)
{
  this->ApplyThreeJPrefactors();
}

// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy

ThreeJSymbol::ThreeJSymbol (const ThreeJSymbol& coefficients)
{
  this->J1 = coefficients.J1;
  this->J2 = coefficients.J2;
  this->Coefficients = coefficients.Coefficients;
  this->JMin = coefficients.JMin;
  this->CurrentPosition = -1;
  this->Flag = coefficients.Flag;
}

// destructor
//

ThreeJSymbol::~ThreeJSymbol ()
{
//   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
//     {
//       int TotalM1 = (this->J1 + 1);
//       int TotalM2 = (this->J2 + 1);      
//       for (int i = 0; i < TotalM1; ++i)
// 	{
// 	  for (int j = 0; j < TotalM2; ++j)
// 	    delete[] this->Coefficients[i][j];
// 	  delete[] this->Coefficients[i];
// 	  delete[] this->JMin[i];
// 	}
//       delete[] this->JMin;
//       delete[] this->Coefficients;
//     }
}

// assignment (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to assign
// return value = reference on current Clebsch Gordan coefficients

ThreeJSymbol& ThreeJSymbol::operator = (const ThreeJSymbol& coefficients)
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      int TotalM1 = (this->J1 + 1);
      int TotalM2 = (this->J2 + 1);      
      for (int i = 0; i < TotalM1; ++i)
	{
	  for (int j = 0; j < TotalM2; ++j)
	    delete[] this->Coefficients[i][j];
	  delete[] this->Coefficients[i];
	  delete[] this->JMin[i];
	}
      delete[] this->JMin;
    }
  this->Flag = coefficients.Flag;
  this->J1 = coefficients.J1;
  this->J2 = coefficients.J2;
  this->Coefficients = coefficients.Coefficients;
  this->JMin = coefficients.JMin;
  this->CurrentPosition = -1;
  return *this;
}

// get a particular coefficient (without testing if m1, m2 and j are valid)
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum 
// j = resulting angular momentum
// return value = corresponding Clebsch Gordan coefficient

double ThreeJSymbol::GetCoefficient (int m1, int m2, int j)
{
  int TmpPos1 = (this->J1 + m1) >> 1;
  int TmpPos2 = (this->J2 + m2) >> 1;
  return this->Coefficients[TmpPos1][TmpPos2][(j - this->JMin[TmpPos1][TmpPos2]) >> 1];
}

// initial an iterator on all Clebsch Gordan coefficients for fixed m1 and m2 values
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum   

void ThreeJSymbol::InitializeCoefficientIterator(int m1, int m2)
{
  this->M1 = (this->J1 + m1) >> 1;
  this->M2 = (this->J2 + m2) >> 1;
  this->CurrentPosition = 0;
  this->J = this->JMin[this->M1][this->M2];
}

// return next coefficient associated with current iterator (with increasing j value)
//
// j = reference on integer where resulting angular momentum has to be stored
// coefficient = reference on double where Clebsch Gordan coefficient has to be stored
// return value = false if no coefficient has been returned

bool ThreeJSymbol::Iterate(int& j, double& coefficient)
{
  if (this->CurrentPosition == -1)
    return false;
  coefficient = this->Coefficients[this->M1][this->M2][this->CurrentPosition];
  j = this->J;
  this->J += 2;
  if (this->J > (this->J1 + this->J2))
    this->CurrentPosition = - 1;
  else
    ++this->CurrentPosition;
  return true;
}

// print a particular coefficient (without testing if m1, m2 and j are valid)
//
// str = reference on output stream
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum 
// j = resulting angular momentum
// return value = reference on output stream

ostream& ThreeJSymbol::PrintCoefficient (ostream& str, int m1, int m2, int j)
{
  int TmpPos1 = (this-> J1 + m1) >> 1;
  int TmpPos2 = (this-> J2 + m2) >> 1;
  str << " < ";
  if ((this->J1 & 0x1) == 0)
    str << (this->J1 >> 1) << " " << (m1 / 2) << " ; ";
  else
    str << this->J1 << "/2 " << m1 << "/2 ; ";
  if ((this->J2 & 0x1) == 0)
    str << (this->J2 >> 1) << " " << (m2 / 2) << " | ";
  else
    str << this->J2 << "/2 " << m2 << "/2 | ";
  if ((j & 0x1) == 0)
    str << (j >> 1) << " " << ((m1 + m2) / 2);
  else
    str << j << "/2 " << (m1 + m2) << "/2";
  str << " > = " << this->Coefficients[TmpPos1][TmpPos2][(j - this->JMin[TmpPos1][TmpPos2]) >> 1];
  return str;
}

// Apply prefactors to convert Clebsch-Gordon Coefficient to a 3J symbol
void ThreeJSymbol::ApplyThreeJPrefactors()
{
  int SumJ = this->J1 + this->J2;
  int DiffJ = abs(this->J1 - this->J2);
  int M2Max;
  int M2Min;
  for (int m1 = -this->J1; m1 <= this->J1; m1 += 2)
    {
      for (int j = DiffJ; j <= SumJ; j += 2)
	{
	  double J3Factor = 1.0/sqrt(1.0+(double)j);
	  M2Max = j - m1;
	  if(M2Max > this->J2)
	    M2Max = this->J2;
	  M2Min = -j - m1;
	  if(M2Min < (-this->J2))
	    M2Min = -this->J2;
	  int M1Pos = (m1 + this->J1) >> 1;
	  for (int m2 = M2Min; m2<=M2Max; m2 += 2)
	    {	      
	      int M2Pos = (m2 + this->J2) >> 1;
	      int Parity = (this->J1-this->J2+m1+m2)>>1;
	      if (Parity & 0x1)
		this->Coefficients[M1Pos][M2Pos][(j - this->JMin[M1Pos][M2Pos]) >> 1] *= -J3Factor;
	      else
		this->Coefficients[M1Pos][M2Pos][(j - this->JMin[M1Pos][M2Pos]) >> 1] *= J3Factor;
	    }
	}
    }
}
