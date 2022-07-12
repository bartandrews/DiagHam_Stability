////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//       class implementing Symmetrizer for complex function //
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

#include "MathTools/NumericalAnalysis/SymmetrizedComplexFunction.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
using std::cout;
using std::endl;

// default constructor
SymmetrizedComplexFunction::SymmetrizedComplexFunction()
{
  this->OriginalFunction=0;
  this->NbrCoordinates=0;
}

// construct Symmetric form of given function
// originalFunction = function to be AS
// nbrCoordinates = number of coordinates to be permuted
// dimCoordinate = number of elements per coordinate
SymmetrizedComplexFunction::SymmetrizedComplexFunction(Abstract1DComplexFunction *originalFunction, int nbrCoordinates, int dimCoordinate)
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
SymmetrizedComplexFunction::SymmetrizedComplexFunction(SymmetrizedComplexFunction &F)
{
  this->OriginalFunction=F.OriginalFunction;
  this->NbrCoordinates=F.NbrCoordinates;
  this->CoordinateDimension=F.CoordinateDimension;
  this->FactorialPrefactor=F.FactorialPrefactor;
  this->InternalCoordinates=RealVector(F.InternalCoordinates,true);
  this->Permutations = new int[this->NbrCoordinates];
}

// destructor
SymmetrizedComplexFunction::~SymmetrizedComplexFunction()
{
  if (this->NbrCoordinates!=0)
    delete [] Permutations;
}

  // clone function 
  //
  // return value = clone of the function 
Abstract1DComplexFunction* SymmetrizedComplexFunction::Clone ()
{
  return new SymmetrizedComplexFunction(*this);
}

 
// routines used to generate permutations

void SymmetrizedComplexFunction::AddSTerm(int sign)
{
//   cout << "Add term " << NbrTerms++ << " ["<<Permutations[0];
//   for (int i = 1; i < NbrCoordinates; i++)
//     cout << ", " << Permutations[i];
//   cout <<"]"<<endl;
  for (int i = 0; i < NbrCoordinates; i++)
    for (int j=0; j< this->CoordinateDimension; ++j)
      this->InternalCoordinates[this->CoordinateDimension*i+j]
	=(*(this->GivenCoordinates))[this->CoordinateDimension*Permutations[i]+j];
  cout << "ExContribution:" <<(*OriginalFunction)(InternalCoordinates)<<" P=["<<Permutations[0];
  for (int i = 1; i < NbrCoordinates; i++)
     cout << ", " << Permutations[i];
  cout <<"]"<<endl;
  SSum+=(*OriginalFunction)(InternalCoordinates);
} // AddASTerm


void SymmetrizedComplexFunction::Swap(const int i, const int j, int &sign)
{
  //  cout << "Swap ("<<i<<", "<<j<<" sign: " << sign << endl;
  int t;
  t = Permutations[i];
  Permutations[i] = Permutations[j];
  Permutations[j] = t;
  sign ^= 1;
} // Swap


void SymmetrizedComplexFunction::RotateLeft(const int start, const int n, int &sign)
{
  //cout << "RotateLeft ("<<start<<", "<<n<<" sign: " << sign << endl;
  int tmp = Permutations[start];
  for (int i = start; i < n-1; i++) {
    Permutations[i] = Permutations[i+1];
  }
  Permutations[n-1] = tmp;
  sign ^= (((n-start)^1)&1);
} // RotateLeft


void SymmetrizedComplexFunction::Permute(const int start, const int n, int sign)
{
  AddSTerm(sign);
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
Complex SymmetrizedComplexFunction::operator ()(RealVector& x)
{
  int sign=1;
  this->SSum=0.0;
//   cout << "Evaluating Function"<<endl;
//   this->NbrTerms=1;
  this->GivenCoordinates=&x;
  for (int i = 0; i < this->NbrCoordinates; i++) this->Permutations[i] = i;
  this->Permute(0, this->NbrCoordinates, sign);
//   cout << "Function Value="<<this->SSum*this->FactorialPrefactor<<endl;
  return (this->SSum*this->FactorialPrefactor);
}


