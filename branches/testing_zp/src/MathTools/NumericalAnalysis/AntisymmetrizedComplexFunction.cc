////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//       class implementing antisymmetrizer for complex function //
//                                                                            //
//                        last modification : 06/12/2007                      //
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

#include "MathTools/NumericalAnalysis/AntisymmetrizedComplexFunction.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
using std::cout;
using std::endl;

// default constructor
AntisymmetrizedComplexFunction::AntisymmetrizedComplexFunction()
{
  this->OriginalFunction=0;
  this->NbrCoordinates=0;
}

// construct antisymmetric form of given function
// originalFunction = function to be AS
// nbrCoordinates = number of coordinates to be permuted
// dimCoordinate = number of elements per coordinate
AntisymmetrizedComplexFunction::AntisymmetrizedComplexFunction(Abstract1DComplexFunction *originalFunction, int nbrCoordinates, int dimCoordinate)
{
  this->OriginalFunction=originalFunction;
  this->NbrCoordinates=nbrCoordinates;
  this->CoordinateDimension=dimCoordinate;
  FactorialCoefficient Tmp;
  Tmp.SetToOne();
  Tmp.FactorialDivide(nbrCoordinates);
  this->FactorialPrefactor=Tmp.GetNumericalValue();
  this->InternalCoordinates.Resize(this->CoordinateDimension*this->NbrCoordinates);
  this->Permutations = new int[this->NbrCoordinates];
}


// copy constructor
AntisymmetrizedComplexFunction::AntisymmetrizedComplexFunction(AntisymmetrizedComplexFunction &F)
{
  this->OriginalFunction=F.OriginalFunction;
  this->NbrCoordinates=F.NbrCoordinates;
  this->CoordinateDimension=F.CoordinateDimension;
  this->FactorialPrefactor=F.FactorialPrefactor;
  this->InternalCoordinates=RealVector(F.InternalCoordinates,true);
  this->Permutations = new int[this->NbrCoordinates];
}

// destructor
AntisymmetrizedComplexFunction::~AntisymmetrizedComplexFunction()
{
  if (this->NbrCoordinates!=0)
    delete [] Permutations;
}

  // clone function 
  //
  // return value = clone of the function 
Abstract1DComplexFunction* AntisymmetrizedComplexFunction::Clone ()
{
  return new AntisymmetrizedComplexFunction(*this);
}

 
// routines used to generate permutations

void AntisymmetrizedComplexFunction::AddASTerm(int sign)
{
  //cout << "AddTerm: " << NbrCoordinates << " sign: " << sign << endl;
  for (int i = 0; i < NbrCoordinates; i++)
    for (int j=0; j< this->CoordinateDimension; ++j)
      this->InternalCoordinates[this->CoordinateDimension*i+j]
	=(*(this->GivenCoordinates))[this->CoordinateDimension*Permutations[i]+j];
  ASSum+=(-1.0+2.0*sign)*(*OriginalFunction)(InternalCoordinates);
} // AddASTerm


void AntisymmetrizedComplexFunction::Swap(const int i, const int j, int &sign)
{
  //  cout << "Swap ("<<i<<", "<<j<<" sign: " << sign << endl;
  int t;
  t = Permutations[i];
  Permutations[i] = Permutations[j];
  Permutations[j] = t;
  sign ^= 1;
} // Swap


void AntisymmetrizedComplexFunction::RotateLeft(const int start, const int n, int &sign)
{
  //cout << "RotateLeft ("<<start<<", "<<n<<" sign: " << sign << endl;
  int tmp = Permutations[start];
  for (int i = start; i < n-1; i++) {
    Permutations[i] = Permutations[i+1];
  }
  Permutations[n-1] = tmp;
  sign ^= (((n-start)^1)&1);
} // RotateLeft


void AntisymmetrizedComplexFunction::Permute(const int start, const int n, int sign)
{
  AddASTerm(sign);
  //cout << "Permute ("<<start<<", "<<n<<" sign: " << sign << endl;
  if (start < n) {
    int i, j;
    for (i = n-2; i >= start; i--) {
      for (j = i + 1; j < n; j++) {
	this->Swap(i, j,sign);
	Permute(i+1, n,sign);
      } // for j
      RotateLeft(i, n,sign);
    } // for i
  }
} // Permute



// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  
Complex AntisymmetrizedComplexFunction::operator ()(RealVector& x)
{
  int sign=1;
  this->ASSum=0.0;
  this->GivenCoordinates=&x;
  for (int i = 0; i < this->NbrCoordinates; i++) this->Permutations[i] = i;
  this->Permute(0, this->NbrCoordinates, sign);
  return (this->ASSum*this->FactorialPrefactor);
}


