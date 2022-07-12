////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                       class author: Nicolas Regnault                       //
//                                                                            //
//              class of quantum Hall many-body hamiltonian associated        //
//           to particles on a torus with spin and magnetic translations      //
//                                                                            //
//                        last modification : 23/10/2015                      //
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
#include "Hamiltonian/AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "MathTools/Complex.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::hex;
using std::dec;
using std::endl;
using std::ostream;


// destructor
//

AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::~AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian()
{  
  delete[] this->ExponentialFactors;
}

// evaluate all exponential factors
//   

void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex[this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
}

// get all the indices that should appear in the annihilation/creation operators
//
// onlyIntraFlag = true if only the intra component indices have to be computed

void AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::GetIndices(bool onlyIntraFlag)
{
  this->NBodySign = new double* [this->MaxNBody + 1];
  this->SpinIndices = new int**[this->MaxNBody + 1];
  this->SpinIndicesShort = new int*[this->MaxNBody + 1];
  this->NbrNBodySpinMomentumSectorSum = new int* [this->MaxNBody + 1];
  this->NbrNBodySpinMomentumSectorIndicesPerSum = new int** [this->MaxNBody + 1];
  this->NBodySpinMomentumSectorIndicesPerSum = new int*** [this->MaxNBody + 1];
  for (int TmpNBody = 0; TmpNBody <= this->MaxNBody; ++TmpNBody)
    {
      if (this->NBodyFlags[TmpNBody] == true)
	{
	  if (onlyIntraFlag == false)
	    this->NbrSpinSectors[TmpNBody] = TmpNBody + 1;
	  else
	    this->NbrSpinSectors[TmpNBody] = 2;
	  this->NBodySign[TmpNBody] = new double [this->NbrSpinSectors[TmpNBody]];
	  this->SpinIndices[TmpNBody] = new int*[this->NbrSpinSectors[TmpNBody]];
	  this->SpinIndicesShort[TmpNBody] = new int[this->NbrSpinSectors[TmpNBody]];
	  this->NbrNBodySpinMomentumSectorSum[TmpNBody] = new int[this->NbrSpinSectors[TmpNBody]];
	  this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody] = new int*[this->NbrSpinSectors[TmpNBody]];
	  this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody] = new int**[this->NbrSpinSectors[TmpNBody]];
	  int Index = 0;
	  for (int i = 0; i <= TmpNBody; ++i)
	    {
	      if ((onlyIntraFlag == false) || (i == 0) || (i == TmpNBody))
		{
		  this->NBodySign[TmpNBody][Index] = 1.0;
		  this->SpinIndicesShort[TmpNBody][Index] = ((1 << (TmpNBody - i)) - 1);
		  this->SpinIndicesShort[TmpNBody][Index] |= this->SpinIndicesShort[TmpNBody][Index] << TmpNBody;
		  this->SpinIndices[TmpNBody][Index] = new int [2 * TmpNBody];
		  for (int j = 0; j < (TmpNBody - i); ++j)
		    {
		      this->SpinIndices[TmpNBody][Index][j] = 1;
		      this->SpinIndices[TmpNBody][Index][j + TmpNBody] = 1;
		    }
		  for (int j = TmpNBody - i; j < TmpNBody; ++j)
		    {
		      this->SpinIndices[TmpNBody][Index][j] = 0;     
		      this->SpinIndices[TmpNBody][Index][j + TmpNBody] = 0;     
		    }
		  this->NbrNBodySpinMomentumSectorSum[TmpNBody][Index] = this->NbrLzValue;
		  if ((i == 0) || (i == TmpNBody))
		    {
		      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
			{
			  this->GetAllSkewSymmetricIndices(this->NbrLzValue, TmpNBody, this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index],  
							   this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index]);
			}
		      else
			{
			  double** TmpSymmetryFactors;
			  this->GetAllSymmetricIndices(this->NbrLzValue, TmpNBody, this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index],  
						       this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index], TmpSymmetryFactors);
			}
		    }
		  else
		    {
		      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
			{
			  this->GetAllTwoSetSkewSymmetricIndices(this->NbrLzValue, TmpNBody - i, i, this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index],  
								 this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index]);
			}
		      else
			{
//			  double** TmpSymmetryFactors;
// 			  this->GetAllTwoSetSymmetricIndices(this->NbrLzValue, TmpNBody - i, i, this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index],  
// 							     this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index], TmpSymmetryFactors);
			  this->GetAllIndices(this->NbrLzValue, TmpNBody, this->NbrNBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index],  
					      this->NBodySpinMomentumSectorIndicesPerSum[TmpNBody][Index]);
			}
		    }
		  ++Index;
		}
	    }
	}
    }


  if (this->FullTwoBodyFlag == true)
    {
      
      this->NbrInterSectorSums = this->NbrLzValue;
      this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	this->NbrInterSectorIndicesPerSum[i] = 0;      
      this->NbrIntraSectorSums = this->NbrLzValue;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      

      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  ++this->NbrInterSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  if (this->NbrInterSectorIndicesPerSum[i]  > 0)
	    {
	      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	      this->NbrInterSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 + m2) % this->NbrLzValue;
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
      
      
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
	  this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
		{
		  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
		  this->NbrIntraSectorIndicesPerSum[i] = 0;
		}
	    }
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      {
		int TmpSum = (m1 + m2) % this->NbrLzValue;
		this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = m1;
		this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = m2;
		++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	      }
	}
      else
	{
	  
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
	  this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
		{
		  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
		  this->NbrIntraSectorIndicesPerSum[i] = 0;
		}
	    }
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      {
		int TmpSum = (m1 + m2) % this->NbrLzValue;
		this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = m1;
		this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = m2;
		++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	      }
	}
    }
  else
    {
      this->NbrInterSectorSums = 0;
      this->NbrInterSectorIndicesPerSum = 0;
      this->NbrIntraSectorSums = 0;
      this->NbrIntraSectorIndicesPerSum = 0;
    }
}


// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
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

long AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
												int**& sortedIndicesPerSum,
												double**& sortedIndicesPerSumSymmetryFactor)
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

  int MaxSum = (nbrValues - 1);
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

// get all indices needed to characterize a  tensor made of two completly skew symmetric  sets of indices, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
// nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::GetAllTwoSetSkewSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, 
													   int**& sortedIndicesPerSum)
{
  int* TmpNbrSortedIndicesPerSumUp;
  int** TmpSortedIndicesPerSumUp;
  this->GetAllSkewSymmetricIndices(nbrValues, nbrIndicesUp, TmpNbrSortedIndicesPerSumUp, TmpSortedIndicesPerSumUp);
  int MaxSumUp = (nbrValues - 1);
  int MinSumUp = 0;
  long NbrElements = 0l;
  int* TmpNbrSortedIndicesPerSumDown;
  int** TmpSortedIndicesPerSumDown;
  this->GetAllSkewSymmetricIndices(nbrValues, nbrIndicesDown, TmpNbrSortedIndicesPerSumDown, TmpSortedIndicesPerSumDown);
  int MaxSumDown = (nbrValues - 1);
  int MinSumDown = 0;
  int MaxSum = (nbrValues - 1);
  int MinSum = 0;
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = MinSumUp; i <= MaxSumUp; ++i)
    {
      for (int j = MinSumDown; j <= MaxSumDown; ++j)
	nbrSortedIndicesPerSum[(i + j) % nbrValues] += TmpNbrSortedIndicesPerSumUp[i] * TmpNbrSortedIndicesPerSumDown[j];
    }
  
  for (int i = MinSum; i <= MaxSum; ++i)
    {
      NbrElements += (long) nbrSortedIndicesPerSum[i];
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + nbrIndicesDown)];
      nbrSortedIndicesPerSum[i] = 0;
    }

  for (int i = MinSumDown; i <= MaxSumDown; ++i)
    {
      int* TmpSortedIndicesPerSumDown2 = TmpSortedIndicesPerSumDown[i];
      int Lim2 = TmpNbrSortedIndicesPerSumDown[i];
      int Pos3 = 0;
      for (int m = 0; m < Lim2; ++m)
	{
	  for (int j = MinSumUp; j <= MaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int TmpSum = (i + j) % nbrValues;
	      int* TmpSortedIndicesPerSum = sortedIndicesPerSum[TmpSum];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  int Pos2 = nbrSortedIndicesPerSum[TmpSum] * (nbrIndicesUp + nbrIndicesDown);
		  for (int l = 0; l < nbrIndicesUp; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
		  for (int l = 0; l < nbrIndicesDown; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumDown2[Pos3 + l];		    
		  ++nbrSortedIndicesPerSum[TmpSum];
		}
	    }
	  Pos3 += nbrIndicesDown;
	}
    }
  for (long i = MinSumDown; i <= MaxSumDown; ++i)
    delete[] TmpSortedIndicesPerSumDown[i];
  delete[] TmpNbrSortedIndicesPerSumDown;
  delete[] TmpSortedIndicesPerSumDown;
  for (long i = MinSumUp; i <= MaxSumUp; ++i)
    delete[] TmpSortedIndicesPerSumUp[i];
  delete[] TmpNbrSortedIndicesPerSumUp;
  delete[] TmpSortedIndicesPerSumUp;
  return NbrElements;
}

// get all indices needed to characterize a  tensor made of two completly symmetric  sets of indices, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
// nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::GetAllTwoSetSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum,
												      double**& sortedIndicesPerSumSymmetryFactor)
{
  int* TmpNbrSortedIndicesPerSumUp;
  int** TmpSortedIndicesPerSumUp;
  double** TmpSortedIndicesPerSumSymmetryFactorUp;
  this->GetAllSymmetricIndices(nbrValues, nbrIndicesUp, TmpNbrSortedIndicesPerSumUp, TmpSortedIndicesPerSumUp,
			       TmpSortedIndicesPerSumSymmetryFactorUp);
  int MaxSumUp = (nbrValues - 1);
  int MinSumUp = 0;
  long NbrElements = 0l;
  int* TmpNbrSortedIndicesPerSumDown;
  int** TmpSortedIndicesPerSumDown;
  double** TmpSortedIndicesPerSumSymmetryFactorDown;
  this->GetAllSymmetricIndices(nbrValues, nbrIndicesDown, TmpNbrSortedIndicesPerSumDown, TmpSortedIndicesPerSumDown,
			       TmpSortedIndicesPerSumSymmetryFactorDown);
  int MaxSumDown = (nbrValues - 1);
  int MinSumDown = 0;
  int MaxSum = (nbrValues - 1);
  int MinSum = 0;
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = MinSum; i <= MaxSum; ++i)
    {
      int TmpMinSumUp = i - MaxSumDown;
      if (TmpMinSumUp < MinSumUp)
	TmpMinSumUp = MinSumUp;
      int TmpMaxSumUp = i - MinSumDown;
      if (TmpMaxSumUp > MaxSumUp)
	TmpMaxSumUp = MaxSumUp;
      for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j] * TmpNbrSortedIndicesPerSumDown[i - j];
      NbrElements += (long) nbrSortedIndicesPerSum[i];
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + nbrIndicesDown)];
      sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
      int Pos2 = 0;
      int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
      double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[i];
      for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)	    
	{
	  int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	  double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[j];
	  int Lim = TmpNbrSortedIndicesPerSumUp[j];
	  int Pos = 0;
	  int Pos4 = 0;
	  for (int k = 0; k < Lim; ++k)
	    {
	      int Pos3 = 0;
	      int* TmpSortedIndicesPerSumDown2 = TmpSortedIndicesPerSumDown[i - j];
	      double* TmpSortedIndicesPerSumSymmetryFactorDown2 = TmpSortedIndicesPerSumSymmetryFactorDown[i - j];
	      int Lim2 = TmpNbrSortedIndicesPerSumDown[i - j];
	      for (int m = 0; m < Lim2; ++m)
		{
		  for (int l = 0; l < nbrIndicesUp; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos + l];
		  for (int l = 0; l < nbrIndicesDown; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumDown2[Pos3++];		      
		  TmpSortedIndicesPerSumSymmetryFactor[Pos4] = (TmpSortedIndicesPerSumSymmetryFactorUp2[k] * 
								TmpSortedIndicesPerSumSymmetryFactorDown2[m]);
		  ++Pos4;
		}
	      Pos += nbrIndicesUp;
	    }
	}
    }
  for (long i = MinSumDown; i <= MaxSumDown; ++i)
    delete[] TmpSortedIndicesPerSumDown[i];
  delete[] TmpNbrSortedIndicesPerSumDown;
  delete[] TmpSortedIndicesPerSumDown;
  for (long i = MinSumUp; i <= MaxSumUp; ++i)
    delete[] TmpSortedIndicesPerSumUp[i];
  delete[] TmpNbrSortedIndicesPerSumUp;
  delete[] TmpSortedIndicesPerSumUp;
  return NbrElements;
}

// get all indices (without assuming any symmetry), sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnTorusWithSpinAndMagneticTranslationsNBodyHamiltonian::GetAllIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum)
{
  long NbrElements = nbrValues;
  for (int i = 1; i < nbrIndices; ++i)
    NbrElements *= (long) nbrValues;
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (long i = 0l; i < NbrElements; ++i)
    {
      Indices[Pos] = new int [nbrIndices];
      Sum[Pos] = 0;
      long TmpIndex = i;
      for (int j = 0; j < nbrIndices; ++j)
	{	  
	  Indices[Pos][j] = (int) (TmpIndex % nbrValues);
	  TmpIndex /= (long) nbrValues;
	  Sum[Pos] += Indices[Pos][j];
	}      
      ++Pos;
    }
  for (int i = 0; i < NbrElements; ++i)
    {
      Sum[i] = (Sum[i] % nbrValues);
    }
  int MaxSum = nbrValues - 1;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
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
      nbrSortedIndicesPerSum[i] = 0;
    }
  for (long i = 0l; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      int* TmpIndices = Indices[i];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = TmpIndices[j];
	  ++TmpPos[Pos];
	}
      delete[] TmpIndices;
      ++nbrSortedIndicesPerSum[Pos];
    }

  delete[] TmpPos;
  delete[] Sum;
  delete[] Indices;
  return NbrElements;
}

