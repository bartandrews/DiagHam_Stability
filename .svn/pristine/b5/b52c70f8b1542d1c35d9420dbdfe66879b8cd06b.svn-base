////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Wigner small D matrix (i.e. Wigner D matrix restricted      //
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


#ifndef WIGNERSMALLDMATRIX_H
#define WIGNERSMALLDMATRIX_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <math.h>

#include <iostream>


using std::cout;
using std::endl;




class WignerSmallDMatrix
{

 private:

  // twice the angular momentum of the sector the matrix is acting on
  int JValue;

  // array of coefficients using to compute one matrix element
  double*** Coefficients;
  // number of coefficients per component to compute
  int** JMin;
  // garbage flag associated to coefficient array
  GarbageFlag Flag;


 public:

  // constructor 
  //
  // jValue =  twice the angular momentum of the sector the matrix is acting on
  WignerSmallDMatrix(int jValue);

  // destructor
  //
  ~WignerSmallDMatrix();

  // get a particular coefficient (without testing if m1 and  m2 are valid)
  //
  // m1 = first angular momentum projection
  // m2 = second angular momentum projection
  // angle = angle of the rotaion along the y axis
  // return value = corresponding Clebsch Gordan coefficient
  double operator () (int m1, int m2, double angle);

 private:
  
  // evaluate all coefficients to compute the Wigner matrix
  //
  void EvaluateWignerSmallDMatrix();

};

// get a particular coefficient (without testing if m1 and  m2 are valid)
//
// m1 = first angular momentum projection
// m2 = second angular momentum projection
// angle = angle of the rotaion along the y axis
// return value = corresponding Clebsch Gordan coefficient

inline double WignerSmallDMatrix::operator () (int m1, int m2, double angle)
{
  double Sign = 1.0;
  m1 += this->JValue;
  m1 >>= 1;
  m2 += this->JValue;
  m2 >>= 1;
  if (m2 > m1)
    {
      int Tmp = m1;
      m1 = m2;
      m2 = Tmp;
      if (((m1 - m2) & 1) != 0)
	Sign = -1.0;
    }
  angle *= 0.5;
  double CosValue = cos(angle);
  double SinValue = sin(angle);
  int Diff = m1 - m2;;
  double Coefficient = 0.0;
  double* TmpCoefficients = this->Coefficients[m1][m2];
  double CurrentSin = 1.0;
  double CurrentCos = 1.0;  
  for (int i = 0; i < Diff; ++i)
    CurrentSin *= SinValue;
  int TmpMinCos = -Diff;
  TmpMinCos += this->JValue;
  TmpMinCos -= 2 * JMin[m1][m2];
  for (int i = 0; i < TmpMinCos; ++i)
    CurrentCos *= CosValue;
  CosValue *= CosValue;
  SinValue *= SinValue;
  int Max = JMin[m1][m2];
//  cout << CurrentSin << " " << CurrentCos << " " << CosValue << " " << SinValue << " " << Max << endl;
  for (int i = 0; i <= Max; ++i)
    {
      double TmpCos = CurrentCos;
      for (int j = Max; j > i; --j)
	TmpCos *= CosValue;
//      Coefficient += CurrentSin * TmpCos;
      Coefficient += TmpCoefficients[i] * CurrentSin * TmpCos;
      CurrentSin *= SinValue;
    }
  return (Sign * Coefficient);
}

#endif
