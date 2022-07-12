////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract quantum Hall hamiltonian associated          //
//             to particles on a torus with n-body interaction terms          //
//                                                                            //
//                        last modification : 04/06/2010                      //
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
#include "Hamiltonian/AbstractQHEOnTorusNBodyInteractionHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"


#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

AbstractQHEOnTorusNBodyInteractionHamiltonian::~AbstractQHEOnTorusNBodyInteractionHamiltonian()
{
}


// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnTorusNBodyInteractionHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
										 int**& sortedIndicesPerSum)
{
  long** BinomialCoefficients = GetBinomialCoefficients(nbrValues);
  long NbrElements = BinomialCoefficients[nbrValues][nbrIndices];
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Min = nbrIndices - 1;
  int Max;
  int Step;
  int Pos = 0;
  for (int i = nbrValues - 1; i >= Min; --i)
    {
      Step = BinomialCoefficients[i][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      Min = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1] - 1;
	  for (; Max >= Min; --Max)
	    {
	      Step = BinomialCoefficients[Max][Min];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }
  
  for (int i = 0; i < NbrElements; ++i)
    {
      Sum[i] = (Sum[i] % nbrValues);
    }
  
  int MaxSum = (nbrValues - 1);
  int MinSum = 0;
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    ++nbrSortedIndicesPerSum[Sum[i]];
  long* TmpPos = new long [MaxSum + 1];
  for (int i = MinSum; i <= MaxSum; ++i)
    {
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      nbrSortedIndicesPerSum[i] = 0;
      TmpPos[i] = 0l;      
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = Indices[i][j];
	  ++TmpPos[Pos];
	}
      ++nbrSortedIndicesPerSum[Pos];
      delete[] Indices[i];
    }
  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] BinomialCoefficients[i];
  delete[] BinomialCoefficients;
  return NbrElements;
}

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnTorusNBodyInteractionHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
									    int**& sortedIndicesPerSum, double**& sortedIndicesPerSumSymmetryFactor)
{
  long** DimensionSymmetricGroup;
  if (nbrValues >= nbrIndices)
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrValues);
  else
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrIndices);
  long NbrElements = DimensionSymmetricGroup[nbrValues][nbrIndices];
  
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (int i = nbrValues - 1; i >= 0; --i)
    {
      Step = DimensionSymmetricGroup[i + 1][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      TmpNbrIndices = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1];
	  Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	  for (int j = 0; j < Step; ++j)
	    {
	      Indices[Pos][i] = Max;
	      Sum[Pos] += Max;
	      ++Pos;
	    }
	  --Max;
	  for (; Max >= 0; --Max)
	    {
	      Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }
  
  for (int i = 0; i < NbrElements; ++i)
    {
      Sum[i] = (Sum[i] % nbrValues);
    }
		
  int MaxSum = nbrValues - 1;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    {
      ++nbrSortedIndicesPerSum[Sum[i]];
    }
  for (int i = 0; i <= MaxSum; ++i)
    {
      TmpPos[i] = 0l;
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
      nbrSortedIndicesPerSum[i] = 0;
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      int* TmpIndices = Indices[i];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = TmpIndices[j];
	  ++TmpPos[Pos];
	}
      
      double& SymmetryFactor = sortedIndicesPerSumSymmetryFactor[Pos][Max];
      SymmetryFactor = 1.0;
      for (int j = 1; j < nbrIndices; ++j)
	{
	  int TmpSymmetryFactor = 1;
	  while ((j < nbrIndices) && (TmpIndices[j - 1] == TmpIndices[j]))
	    {
	      ++TmpSymmetryFactor;
	      ++j;
	    }
	  if (TmpSymmetryFactor != 1)
	    for (int k = 2; k <= TmpSymmetryFactor; ++k)
	      SymmetryFactor *= (double) k;
	}
      delete[] TmpIndices;
      SymmetryFactor = 1.0 / SymmetryFactor;
      ++nbrSortedIndicesPerSum[Pos];
    }
			       
  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] DimensionSymmetricGroup[i];
  delete[] DimensionSymmetricGroup;
  return NbrElements;
}
