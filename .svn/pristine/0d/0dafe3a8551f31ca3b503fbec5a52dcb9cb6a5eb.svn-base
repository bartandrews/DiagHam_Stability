////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author : Cecile Repellin                      //
//                                                                            //
//                 class of SU(3) irreducible representations                 //
//                                                                            //
//                        last modification : 17/02/2013                      //
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
#include "SU3IrreducibleRepresentations.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// constructor 
//
// (p,q) = indices that define an irreducible SU(3) representation
SU3IrreducibleRepresentations::SU3IrreducibleRepresentations(int p, int q)
{
  this->P = p;
  this->Q = q;
  this->RepresentationDimension = (this->P + 1)*(this->Q + 1)*(this->P + this->Q + 2)/2;
  this->States = new int* [this->RepresentationDimension];
  this->QuantumNumbers = new int*[this->RepresentationDimension];
  for (int i = 0; i < this->RepresentationDimension; ++i)
  {
    this->States[i] = new int [3];
    this->QuantumNumbers[i] = new int[2];
  }
  this->IrreducibleRepresentationGenerateStates();
  int shiftedTz;
  int shiftedY;
  this->NbrTzValuesPerShiftedY = new int[this->P + this->Q + 1];
  for (int i = 0; i < this->P ; ++i)
     this->NbrTzValuesPerShiftedY[i] = i + this->Q + 1;
  for (int i = this->P; i <= this->P + this->Q; ++i)
    this->NbrTzValuesPerShiftedY[i] = 2*this->P + this->Q - i + 1;
  for (int i = 0; i < this->RepresentationDimension; ++i)
    this->GetQuantumNumbers(i, this->QuantumNumbers[i][0], this->QuantumNumbers[i][1]); 
  this->ComputeDegeneracy();  
}


//destructor
SU3IrreducibleRepresentations::~SU3IrreducibleRepresentations()
{
}

// generate all states in an irreducible representation (p, q) of SU(3)
//
// p = integer that defines the irrep together with q
// q = integer 
//return value = pointer to an array that gives each point in the representation. 
void SU3IrreducibleRepresentations::IrreducibleRepresentationGenerateStates()
{
  int TmpM12 = this->Q;
  int TmpM22;
  int TmpM11;
  int index = 0;
  while (TmpM12 <= (this->P + this->Q))
  {
   TmpM22 = 0;
   while (TmpM22 <= this->Q)
   {
     TmpM11 = TmpM22;
     while (TmpM11 <= TmpM12)
     {
	this->States[index][0] = TmpM12;
	this->States[index][1] = TmpM22;
	this->States[index][2] = TmpM11;
	TmpM11 += 1;
	index += 1;
      }
      TmpM22 += 1;
     }
     TmpM12 += 1;
   }
}


//Compute the degeneracy for all values of (tz,y) accessible in representation (p,q)
//
void SU3IrreducibleRepresentations::ComputeDegeneracy()
{
  int tz;
  int y;
  int TmpTz;
  int TmpY;
  this->Degeneracy = new int*[this->P + this->Q + 1];
  this->QIndices = new int**[this->P + this->Q + 1];
  for (int shiftedY = 0; shiftedY <= this->P + this->Q; ++shiftedY)
  {
   this->Degeneracy[shiftedY] = new int[this->NbrTzValuesPerShiftedY[shiftedY]];
   this->QIndices[shiftedY] = new int*[this->NbrTzValuesPerShiftedY[shiftedY]];
   for (int shiftedTz = 0; shiftedTz < this->NbrTzValuesPerShiftedY[shiftedY] ; ++shiftedTz)
   {
     this->GetQuantumNumbersFromShiftedIndices(tz, y, shiftedTz, shiftedY);
     this->Degeneracy[shiftedY][shiftedTz] = 0;
     for (int index = 0; index < this->RepresentationDimension; ++index)
     {
      if ((QuantumNumbers[index][0] == tz ) && (QuantumNumbers[index][1]== y))
      {
	this->Degeneracy[shiftedY][shiftedTz] += 1;
      }
     }
     this->QIndices[shiftedY][shiftedTz] = new int[this->Degeneracy[shiftedY][shiftedTz]];
     this->Degeneracy[shiftedY][shiftedTz] = 0;
     for (int index = 0; index < this->RepresentationDimension; ++index)
     {
      if ((QuantumNumbers[index][0] == tz ) && (QuantumNumbers[index][1]== y))
      {
	this->QIndices[shiftedY][shiftedTz][this->Degeneracy[shiftedY][shiftedTz]] = index;
	this->Degeneracy[shiftedY][shiftedTz] += 1;
      }
     }
   }
  } 
}