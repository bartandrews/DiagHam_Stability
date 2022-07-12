////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                a set of functions usefull for integer algebra              //
//                                                                            //
//                        last modification : 11/09/2003                      //
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


#include "MathTools/IntegerAlgebraTools.h"


// find greatest common divider
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

int FindGCD(int m, int n)
{
  if (m < n)
    return RecursiveFindGCD (m, n);
  else
    return RecursiveFindGCD (n, m);
  return n;
}

// find greatest common divider (recurisive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

int RecursiveFindGCD(int m, int n)
{
  if (m == 0)
    return n;
  else
    return RecursiveFindGCD ((n % m), m);
}

// get all binomial coefficients up to a given number of element
//
// n = maximum number of elements
// return value = array containing the binomial coefficients (first index (i) corresponding to the number of elements, the second index going from 0 to i)

long** GetBinomialCoefficients (int n)
{
  if (n == 0)
    n = 1;
  long** TmpBinomialCoefficients =  new long* [n + 1]; 
  TmpBinomialCoefficients[0] = new long [1];
  TmpBinomialCoefficients[1] = new long [2];
  TmpBinomialCoefficients[0][0] = 1l;
  TmpBinomialCoefficients[1][0] = 1l;
  TmpBinomialCoefficients[1][1] = 1l;
  for (int i = 2; i <= n; ++i)
    {
      long*& Tmp1 = TmpBinomialCoefficients[i - 1]; 
      TmpBinomialCoefficients[i] = new long [i + 1];
      long*& Tmp2 = TmpBinomialCoefficients[i]; 
      Tmp2[0] = 1l;
      Tmp2[i] = 1l;
      for (int j = 1; j < i; ++j)
	{
	  Tmp2[j] = Tmp1[j - 1] + Tmp1[j];
	}
    }
  return TmpBinomialCoefficients;
}

// get all dimensions of irreducible representations of the symmetric group
//
// n = maximum number of elements
// return value = array containing binomial coeffcient (first index corresponds to the number of elements, the second is the number of indices,
//                if the number of indices is zero, the dimension is then equal to one)

long** GetIrreducibleRepresentationDimensionSymmetricGroup (int n)
{
  if (n == 0)
    n = 1;
  long** TmpBinomialCoefficients =  new long* [n + 1]; 
  for (int i = 0; i <= n; ++i)
    {
      TmpBinomialCoefficients[i] = new long [n + 1];
    }
  TmpBinomialCoefficients[0][0] = 1l;
  TmpBinomialCoefficients[1][0] = 1l;
  for (int i = 1; i <= n; ++i)
    {
      TmpBinomialCoefficients[0][i] = 1l;
      TmpBinomialCoefficients[1][i] = 1l;
    }

  for (int i = 2; i <= n; ++i)
    {
      TmpBinomialCoefficients[i][0] = 1l;
      TmpBinomialCoefficients[i][1] = (long) i;
      for (int j = 2; j <= n; ++j)
	TmpBinomialCoefficients[i][j] = TmpBinomialCoefficients[i - 1][j] + TmpBinomialCoefficients[i][j - 1];	  
    }
  return TmpBinomialCoefficients;
}
