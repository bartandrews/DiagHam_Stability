////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of hamiltonian with particles on                   //
//               Chern insulator in the single band approximation             //
//                           and n-body interaction                           //
//                                                                            //
//                        last modification : 03/08/2011                      //
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
#include "MathTools/FactorialCoefficient.h"
#include "Hamiltonian/ParticleOnCylinderNBodyHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>

#define Pi_L 3.14159265358979323846264338328L

using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnCylinderNBodyHamiltonian::ParticleOnCylinderNBodyHamiltonian()
{
  this->HermitianSymmetryFlag = true;
  this->TwoBodyFlag = false;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum momentum
// nBody = order of n-body interaction
// ratio = aspect ratio of the cylinder// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCylinderNBodyHamiltonian::ParticleOnCylinderNBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum, int nBody, double ratio, AbstractArchitecture* architecture, long memory)
{
  cout<<"Nbody constructor "<<endl;
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->LzMax = maxMomentum;
  this->Ratio = ratio;
  this->HamiltonianShift = 0.0;
  this->NBodyValue = nBody;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->NbrSectorSums = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  this->HermitianSymmetryFlag = true;
  this->TwoBodyFlag = false;
  cout<<"Nbr sums "<<this->NbrNBodySectorSums<<endl;

  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

//ParticleOnCylinderNBodyHamiltonian::~ParticleOnCylinderNBodyHamiltonian()
//{
/*
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      if (this->NbrNBodySectorIndicesPerSum[i] > 0)
	{
	  delete[] this->NBodySectorIndicesPerSum[i];
	  delete[] this->NBodyInteractionFactors[i];
	}
    }
  if (this->NbrNBodySectorSums>0)
    {
      delete[] this->NBodyInteractionFactors;
      delete[] this->NbrNBodySectorIndicesPerSum;
      delete[] this->NBodySectorIndicesPerSum;
    }
*/
//}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderNBodyHamiltonian::EvaluateInteractionFactors()
{



  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
            cout<<"Evaluating interaction elements"<<endl;

/*
	    this->MinSumIndices = 0;
	    this->MaxSumIndices = this->NBodyValue * this->LzMax;
            this->NbrNBodySectorSums = this->NBodyValue * this->LzMax + 1; 


	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->LzMax + 1, this->NBodyValue, this->NbrNBodySectorIndicesPerSum, this->NBodySectorIndicesPerSum,
				   SortedIndicesPerSumSymmetryFactor);
	    this->InteractionFactors = new Complex* [this->NbrNBodySectorSums];
 	    int Lim;
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices; ++MinSum)
	      {
		Lim = this->NbrNBodySectorIndicesPerSum[MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		this->InteractionFactors[MinSum] = new Complex [Lim];
		Complex* TmpNBodyInteractionFactors = this->InteractionFactors[MinSum];		  
		int* TmpMIndices = this->NBodySectorIndicesPerSum[MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    double Coefficient = TmpSymmetryFactors[i];
		    TmpNBodyInteractionFactors[i] = Coefficient * this->EvaluateInteractionCoefficient(TmpMIndices[i]);
		    TmpMIndices += this->NBodyValue;
		  }
	      }
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;

*/

      this->NbrNBodySectorSums = this->NBodyValue * this->LzMax + 1;
      //this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      
      //for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	//this->NbrNBodySectorIndicesPerSum[i] = 0;
      
      int** TmpSortedIndicesPerSumSymmetryFactor;
      
      GetAllSkewSymmetricIndices (this->LzMax + 1, this->NBodyValue, this->NbrNBodySectorIndicesPerSum, this->NBodySectorIndicesPerSum);
      //for (int j = 0; j < this->NbrNBodySectorSums; ++j)
      //  cout<<"Sum= "<<j<<" "<<this->NbrNBodySectorIndicesPerSum[j]<<endl;

/*
  for (int j = 0; j < this->NbrNBodySectorSums; ++j)
   {
      cout<<"Sum= "<<j<<endl;
      int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
      int* TmpIndices = this->NBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	{
	  int* Tmp1 = (TmpIndices + i1);
          for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
	       {
		 int* Tmp2 = (TmpIndices + i2);

		  for (int n1=0; n1<this->NBodyValue; n1++)
                     cout<<Tmp1[n1]<<" ";
		  for (int n1=0; n1<this->NBodyValue; n1++)
                     cout<<Tmp2[n1]<<" ";
                  cout<<endl;

		}
	    
	}
    }

 exit(1);
*/
 
     
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
       
      int TotalNbrInteractionFactors=0;
      int* TmpIndices;
      double MaxMatEl = 0.0;
      double MinMatEl = (double)1e8;

      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
          int Index = 0;
          int* TmpIndices = this->NBodySectorIndicesPerSum[i];
          int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i];

          for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	    { 
              
	      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		{
                  this->NBodyInteractionFactors[i][Index] =  this->EvaluateInteractionCoefficient(TmpIndices + i1, TmpIndices + i2);
                  
                  if (Norm(this->NBodyInteractionFactors[i][Index])>MaxMatEl)
                     MaxMatEl = Norm(this->NBodyInteractionFactors[i][Index]);
                  else if (Norm(this->NBodyInteractionFactors[i][Index])<MinMatEl)
                     MinMatEl = Norm(this->NBodyInteractionFactors[i][Index]);

                  if (Norm(this->NBodyInteractionFactors[i][Index])>1e-14)
                  {
                  int* Tmp1=TmpIndices+i1;
                  for (int n=0; n<this->NBodyValue; ++n)
                     cout<<Tmp1[n]<<" ";

                  int* Tmp2=TmpIndices+i2;
                  for (int n=0; n<this->NBodyValue; ++n)
                     cout<<Tmp2[n]<<" ";
                  cout<<this->EvaluateInteractionCoefficient(TmpIndices + i1, TmpIndices + i2) <<endl;
                  }

 
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }

          //cout<<"Totalnbr = "<<TotalNbrInteractionFactors<<endl;
	}
      
      cout<<"Max matrix el = "<<MaxMatEl<<" min= "<<MinMatEl<<endl;

      cout<<"Done evaluating interaction "<<TotalNbrInteractionFactors<<endl;

    }
  else
    {
            cout<<"Evaluating interaction elements"<<endl;

/*
	    this->MinSumIndices = 0;
	    this->MaxSumIndices = this->NBodyValue * this->LzMax;
            this->NbrNBodySectorSums = this->NBodyValue * this->LzMax + 1; 


	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->LzMax + 1, this->NBodyValue, this->NbrNBodySectorIndicesPerSum, this->NBodySectorIndicesPerSum,
				   SortedIndicesPerSumSymmetryFactor);
	    this->InteractionFactors = new Complex* [this->NbrNBodySectorSums];
 	    int Lim;
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices; ++MinSum)
	      {
		Lim = this->NbrNBodySectorIndicesPerSum[MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		this->InteractionFactors[MinSum] = new Complex [Lim];
		Complex* TmpNBodyInteractionFactors = this->InteractionFactors[MinSum];		  
		int* TmpMIndices = this->NBodySectorIndicesPerSum[MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    double Coefficient = TmpSymmetryFactors[i];
		    TmpNBodyInteractionFactors[i] = Coefficient * this->EvaluateInteractionCoefficient(TmpMIndices[i]);
		    TmpMIndices += this->NBodyValue;
		  }
	      }
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;

*/

      this->NbrNBodySectorSums = this->NBodyValue * this->LzMax + 1;
      //this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      
      //for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	//this->NbrNBodySectorIndicesPerSum[i] = 0;
      
      int** TmpSortedIndicesPerSumSymmetryFactor;
      
      GetAllSymmetricIndices (this->LzMax + 1, this->NBodyValue, this->NbrNBodySectorIndicesPerSum, this->NBodySectorIndicesPerSum, TmpSortedIndicesPerSumSymmetryFactor);
      //for (int j = 0; j < this->NbrNBodySectorSums; ++j)
      //  cout<<"Sum= "<<j<<" "<<this->NbrNBodySectorIndicesPerSum[j]<<endl;
/*
  for (int j = 0; j < this->NbrNBodySectorSums; ++j)
   {
      cout<<"Sum= "<<j<<endl;
      int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[j];
      int* TmpIndices = this->NBodySectorIndicesPerSum[j];
      for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	{
	  int* Tmp1 = (TmpIndices + i1);
          for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
	       {
		 int* Tmp2 = (TmpIndices + i2);

		  for (int n1=0; n1<this->NBodyValue; n1++)
                     cout<<Tmp1[n1]<<" ";
		  for (int n1=0; n1<this->NBodyValue; n1++)
                     cout<<Tmp2[n1]<<" ";
                  cout<<endl;

		}
	    
	}
    }
*/
      
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      

     
      int TotalNbrInteractionFactors=0;
      int* TmpIndices;
      double MaxMatEl = 0.0;
      double MinMatEl = (double)1e8;

      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
          int Index = 0;
          int* TmpIndices = this->NBodySectorIndicesPerSum[i];
          int Lim = this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i];

          for (int i1 = 0; i1 < Lim; i1 += this->NBodyValue)
	    { 
              
	      for (int i2 = 0; i2 < Lim; i2 += this->NBodyValue)
		{
                  FactorialCoefficient TmpFact;
                  TmpFact.FactorialMultiply(this->NBodyValue);
                  double NBodyFact = TmpFact.GetNumericalValue();

                  this->NBodyInteractionFactors[i][Index] =  this->EvaluateInteractionCoefficient(TmpIndices + i1, TmpIndices + i2) * (NBodyFact/TmpSortedIndicesPerSumSymmetryFactor[i][i1/this->NBodyValue]) * (NBodyFact/TmpSortedIndicesPerSumSymmetryFactor[i][i2/this->NBodyValue]);
                  
                  if (Norm(this->NBodyInteractionFactors[i][Index])>MaxMatEl)
                     MaxMatEl = Norm(this->NBodyInteractionFactors[i][Index]);
                  else if (Norm(this->NBodyInteractionFactors[i][Index])<MinMatEl)
                     MinMatEl = Norm(this->NBodyInteractionFactors[i][Index]);

                  //int* Tmp1=TmpIndices+i1;
                  //for (int n=0; n<this->NBodyValue; ++n)
                  //   cout<<Tmp1[n]<<" ";

                  //int* Tmp2=TmpIndices+i2;
                  //for (int n=0; n<this->NBodyValue; ++n)
                  //   cout<<Tmp2[n]<<" ";
                  //cout<<this->EvaluateInteractionCoefficient(TmpIndices + i1, TmpIndices + i2) <<" sym " << TmpSortedIndicesPerSumSymmetryFactor[i][i1/this->NBodyValue] <<" "<< TmpSortedIndicesPerSumSymmetryFactor[i][i2/this->NBodyValue]<<endl;

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }

          //cout<<"Totalnbr = "<<TotalNbrInteractionFactors<<endl;
	}
      
      cout<<"Max matrix el = "<<MaxMatEl<<" min= "<<MinMatEl<<endl;

      for (int Sum = 0; Sum < this->NbrNBodySectorSums; Sum++)
	{
	  delete [] TmpSortedIndicesPerSumSymmetryFactor[Sum];
	}
      
      delete [] TmpSortedIndicesPerSumSymmetryFactor;

           cout<<"Done evaluating interaction "<<TotalNbrInteractionFactors<<endl;

    }

}

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long ParticleOnCylinderNBodyHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
									     int**& sortedIndicesPerSum,
									     int**& sortedIndicesPerSumSymmetryFactor)
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

  int MaxSum = (nbrValues - 1) * nbrIndices;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new int* [MaxSum + 1];
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
      sortedIndicesPerSumSymmetryFactor[i] = new int [nbrSortedIndicesPerSum[i]];
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
      int& SymmetryFactor = sortedIndicesPerSumSymmetryFactor[Pos][Max];
      SymmetryFactor = 1;
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
	      SymmetryFactor *= k;
	}
      //SymmetryFactor = 1.0 / SymmetryFactor;
      ++nbrSortedIndicesPerSum[Pos];
    }


  for (int i = 0; i < NbrElements; ++i)
    {
      if (Indices[i] != 0)
        delete[] Indices[i];
    }
  delete[]Indices;

  delete[] TmpPos;
  delete[] Sum;

  for (int i = 0; i <= nbrValues; ++i)
    delete[] DimensionSymmetricGroup[i];
  delete[] DimensionSymmetricGroup;
  return NbrElements;
}

// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  ParticleOnCylinderNBodyHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
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

  int MaxSum = (((nbrValues - 1) * nbrValues) - ((nbrIndices - 1) * (nbrIndices - 2)))/ 2;
  int MinSum = (nbrIndices * (nbrIndices - 1)) / 2;
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

// evaluate the numerical coefficient  in front of the Prod a_m^+ Prod a_m coupling term for bosons
// Indices1,Indices2 = array containing indices
// return value = numerical coefficient

Complex ParticleOnCylinderNBodyHamiltonian::EvaluateInteractionCoefficient(int* Indices1, int* Indices2)
{
  long double Length = sqrt(2.0L * Pi_L * (long double)this->Ratio * (long double)(this->LzMax + 1));
  long double kappa = 2.0L * Pi_L/Length;
  long double Normalization = 1.0L;
  Complex Coefficient;
  long double MatrixElement;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      long double Sum = 0.0L;
      long double SumSq = 0.0L;
      for (int i = 0; i < this->NBodyValue; ++i)
       {
         Sum += (long double)(Indices1[i] + Indices2[i]);
         SumSq += (long double)(Indices1[i] * Indices1[i] + Indices2[i] * Indices2[i]);
       }


      long double Pairwise = 1.0L;
      for (int i1 = 0; i1 < this->NBodyValue; ++i1)
        for (int i2 = i1 + 1; i2 < this->NBodyValue; ++i2)
          {
            Pairwise *= (long double)(Indices1[i1] - Indices1[i2]);
            Pairwise *= (long double)(Indices2[i1] - Indices2[i2]);
          }

     // for (int i1 = this->NBodyValue - 1; i1 >= 0; --i1)
     //   for (int i2 = i1 - 1; i2 >= 0; --i2)
     //     {
     //       Pairwise *= (double)(Indices2[i1] - Indices2[i2]);
     //     }

      Pairwise *= powl(-1.0L, this->NBodyValue * (this->NBodyValue + 1)/2) * powl(kappa, this->NBodyValue * (this->NBodyValue - 1)/2) * powl(kappa, this->NBodyValue * (this->NBodyValue - 1)/2);

     //if ((0.5 * (SumSq - Sum * Sum/(2.0 * this->NBodyValue) )) <= 100000000.0)
     //  {
          MatrixElement = Pairwise * expl(- 0.5L * kappa * kappa * (SumSq - Sum * Sum/(2.0L * this->NBodyValue) ) );   
     //  }

     //cout<<"Mat el ";
     //for (int i=0; i<this->NBodyValue; ++i)
     //  cout<<Indices1[i]<<" ";

     //for (int i=0; i<this->NBodyValue; ++i)
     //  cout<<Indices2[i]<<" ";

     //cout<<"    "<<Coefficient/sqrt(this->Ratio * (this->LzMax+1))<<endl;

     if (this->NBodyValue == 2)
       //fermions 1/3
       Normalization = 4.0L * sqrtl(2.0L * Pi_L)/Length;
     else if (this->NBodyValue == 3)
     //fermions 1/2
      Normalization = 24.0L * (1.0L/sqrtl(3.0L)) * powl(2.0L * Pi_L, 1.0L)/powl(Length, 2.0L);
     else if (this->NBodyValue == 4) 
     //fermions 3/5
       Normalization = 128.0L * 1.0L/(3.0L * sqrtl(2.0l)) * powl(2.0l * Pi_L, 1.5l)/powl(Length, 3.0l);
     else if (this->NBodyValue == 5)
     //fermions 2/3  
       Normalization = (1280.0l/3.0l)* 1.0l/(6.0l * sqrtl(5.0l)) * powl(2.0l * Pi_L, 2.0l)/pow(Length, 4.0l);
     else if (this->NBodyValue == 6)
     //fermions 5/7
       Normalization = (2048.0l/3.0l) * 1.0l/(30.0l * sqrtl(3.0l)) * powl(2.0l * Pi_L, 2.5l)/powl(Length, 5.0l);
     else if (this->NBodyValue == 7)
     //fermions 3/4
       Normalization =  (425.0l) * 1.0l/(90.0l * sqrtl(7.0l)) * powl(2.0l * Pi_L, 3.0l)/powl(Length, 6.0l);
     else if (this->NBodyValue == 8)
     //fermions 7/9 -- not sure if correct
       Normalization =  (86.5l) * 1.0l/(180.0l * 7.0l) * powl(2.0l * Pi_L, 3.5l)/powl(Length, 7.0l);

     MatrixElement *= Normalization;

    }
  else //Bosons
    {
      long double Sum = 0.0l;
      long double SumSq = 0.0l;
      for (int i = 0; i < this->NBodyValue; ++i)
       {
         Sum += (long double)(Indices1[i] + Indices2[i]);
         SumSq += (long double)(Indices1[i] * Indices1[i] + Indices2[i] * Indices2[i]);
       }

     //if ((0.5 * (SumSq - Sum * Sum/(2.0 * this->NBodyValue) )) <= 100000000.0)
     //  {
          MatrixElement = expl(- 0.5l * kappa * kappa * (SumSq - Sum * Sum/(2.0l * this->NBodyValue) ) );   
     //  }

     //cout<<"Mat el ";
     //for (int i=0; i<this->NBodyValue; ++i)
     //  cout<<Indices1[i]<<" ";

     //for (int i=0; i<this->NBodyValue; ++i)
     //  cout<<Indices2[i]<<" ";

     //cout<<"    "<<Coefficient/sqrt(this->Ratio * (this->LzMax+1))<<endl;


     if (this->NBodyValue == 2)
       Normalization = sqrtl(2.0L * Pi_L)/Length;
     else if (this->NBodyValue == 3)
      Normalization = (1.0L/sqrtl(3.0L)) * powl(2.0L * Pi_L, 1.0L)/powl(Length, 2.0L);
     else if (this->NBodyValue == 4) 
       Normalization = 1.0L/(3.0L * sqrtl(2.0l)) * powl(2.0l * Pi_L, 1.5l)/powl(Length, 3.0l);
     else if (this->NBodyValue == 5)
       Normalization = 1.0l/(6.0l * sqrtl(5.0l)) * powl(2.0l * Pi_L, 2.0l)/pow(Length, 4.0l);
     else if (this->NBodyValue == 6)
       Normalization = 1.0l/(30.0l * sqrtl(3.0l)) * powl(2.0l * Pi_L, 2.5l)/powl(Length, 5.0l);
     else if (this->NBodyValue == 7)
       Normalization =  1.0l/(90.0l * sqrtl(7.0l)) * powl(2.0l * Pi_L, 3.0l)/powl(Length, 6.0l);
     else if (this->NBodyValue == 8)
       Normalization =  1.0l/(180.0l * 7.0l) * powl(2.0l * Pi_L, 3.5l)/powl(Length, 7.0l);
     else if (this->NBodyValue == 9)
       Normalization = (1.0l/(180.0l*42.0l)) * powl(2.0l * Pi_L, 4.0l)/powl(Length, 8.0l);

     MatrixElement *= Normalization;

   }

  Coefficient.Re = MatrixElement;
  Coefficient.Im = 0.0;

  return Coefficient;
}
