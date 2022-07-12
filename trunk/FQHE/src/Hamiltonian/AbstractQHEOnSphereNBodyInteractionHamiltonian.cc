////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract quantum Hall hamiltonian associated          //
//            to particles on a sphere with n-body interaction terms          //
//                                                                            //
//                        last modification : 22/09/2004                      //
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
#include "Hamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

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

// default constructor
//

AbstractQHEOnSphereNBodyInteractionHamiltonian::AbstractQHEOnSphereNBodyInteractionHamiltonian()
{
  this->M4Value = 0;
}

// destructor
//

AbstractQHEOnSphereNBodyInteractionHamiltonian::~AbstractQHEOnSphereNBodyInteractionHamiltonian()
{
  if (this->M4Value != 0)
    {
      delete[] this->M4Value;
    }
}

// ask if Hamiltonian implements hermitian symmetry operations
//
bool AbstractQHEOnSphereNBodyInteractionHamiltonian::IsHermitian()
{
  return false;
}

// ask if Hamiltonian implements conjugate methods
//
bool AbstractQHEOnSphereNBodyInteractionHamiltonian::IsConjugate()
{
  return false;
}

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
//
//bool AbstractQHEOnSphereNBodyInteractionHamiltonian::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
//{
//  return false;
//}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereNBodyInteractionHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
										int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;

  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int* TmpMIndices;
      int* TmpNIndices;
      double* TmpInteraction;
      double Coefficient2;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
 	      int m1 = this->OneBodyMValues[j];
	      int m2 = this->OneBodyNValues[j];
	      double TmpInteraction2 = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = this->Particles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * TmpInteraction2 * vSource[i];		  
		}
	    }
	}
     if (this->NBodyFlags[1] == true)
	{
	  int TmpMinSumIndices = this->MinSumIndices[1];
	  int TmpMaxSumIndices = this->MaxSumIndices[1];	      
	  double Sign = this->NBodySign[1];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
		{
		  int Lim = NbrSortedIndicesPerSum[1][j];
		  TmpNIndices = this->SortedIndicesPerSum[1][j];
		  TmpInteraction = this->NBodyInteractionFactors[1][j];
		  for (int i1 = 0; i1 < Lim; ++i1)
		    {
		      Coefficient2 = Sign * vSource[i] * TmpInteraction[i1];
		      TmpMIndices = this->SortedIndicesPerSum[1][j];
		      for (int i2 = 0; i2 < Lim; ++i2)
			{
			  Index = this->Particles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
			  if (Index < Dim)
			      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction[i2];
			  ++TmpMIndices;
			}
		    }
		  ++TmpNIndices;
		}
	    }
	}
     if (this->FullTwoBodyFlag == true)
       {
	 int m1;
	 int m3;
	 if (this->NbrM12Indices == 0)
	   {
	     double TmpTwoBodyInteraction;
	     int m2;
	     int m4;	     
	     for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	       {
		 m1 = this->M1Value[j];
		 m2 = this->M2Value[j];
		 m3 = this->M3Value[j];
		 if (this->M4Value != 0)
		   m4 = this->M4Value[j];
		 else
		   m4 = m1 + m2 - m3;
		 TmpTwoBodyInteraction = this->InteractionFactors[j];
		 for (int i = firstComponent; i < LastComponent; ++i)
		   {
		     Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		     if (Index < Dim)
		       vDestination[Index] += Coefficient * TmpTwoBodyInteraction * vSource[i];
		   }
	       }
	   }
	 else
	   {
	     double Coefficient2;
	     int SumIndices;
	     int TmpNbrM3Values;
	     int* TmpM3Values;
	     int ReducedNbrInteractionFactors;
	     for (int i = firstComponent; i < LastComponent; ++i)
	       {
		 ReducedNbrInteractionFactors = 0;
		 for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		   {
		     Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		     if (Coefficient != 0.0)
		       {
			 SumIndices = this->M1Value[m1] + this->M2Value[m1];
			 Coefficient *= vSource[i];
			 TmpNbrM3Values = this->NbrM3Values[m1];
			 TmpM3Values = this->M3Values[m1];
			 for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			   {
			     Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			     if (Index < Dim)			
			       vDestination[Index] += Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			     ++ReducedNbrInteractionFactors;
			   }
		       }
		     else
		       ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		   }
	       }
	    }	   
       }
     for (int k = 2; k <= this->MaxNBody; ++k)
       if (this->NBodyFlags[k] == true)
	 {
	   for (int i = firstComponent; i < LastComponent; ++i)
	     this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
	 }	
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
                {
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
                }
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	      int* TmpIndexArray;
	      double* TmpCoefficientArray; 
	      int j;
	      int TmpNbrInteraction;
	      firstComponent -= this->PrecalculationShift;
	      LastComponent -= this->PrecalculationShift;
	      int Pos = firstComponent / this->FastMultiplicationStep; 
	      int PosMod = firstComponent % this->FastMultiplicationStep;
	      if (PosMod != 0)
		{
		  ++Pos;
		  PosMod = this->FastMultiplicationStep - PosMod;
		}
	      int l =  PosMod + firstComponent + this->PrecalculationShift;
	      for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
		  TmpIndexArray = this->InteractionPerComponentIndex[Pos];
		  TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
		  Coefficient = vSource[l];
		  for (j = 0; j < TmpNbrInteraction; ++j)
		    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  vDestination[l] += this->HamiltonianShift * Coefficient;
		  l += this->FastMultiplicationStep;
		  ++Pos;
		}
	      int Index;
	      int* TmpMIndices;
	      int* TmpNIndices;
	      double* TmpInteraction;
	      double Coefficient2;
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    if (this->OneBodyTermFlag == true)
		      {
			for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
			  {
			    int m1 = this->OneBodyMValues[j];
			    int m2 = this->OneBodyNValues[j];
			    double TmpInteraction2 = this->OneBodyInteractionFactors[j];
			    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			      {
				Index = this->Particles->AdA(i, m1, m2, Coefficient);
				if (Index < Dim)
				  vDestination[Index] += Coefficient * TmpInteraction2 * vSource[i];		  
			      }
			  }
		      }
		    if (this->NBodyFlags[1] == true)
		      {
			double Sign = this->NBodySign[1];
			int TmpMinSumIndices = this->MinSumIndices[1];
			int TmpMaxSumIndices = this->MaxSumIndices[1];	      
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  {
			    for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
			      {
				int Lim = NbrSortedIndicesPerSum[1][j];
				TmpNIndices = this->SortedIndicesPerSum[1][j];
                                TmpInteraction = this->NBodyInteractionFactors[1][j];
				for (int i1 = 0; i1 < Lim; ++i1)
				  {
				    Coefficient2 = Sign * vSource[i] * TmpInteraction[i1];
				    TmpMIndices = this->SortedIndicesPerSum[1][j];
				    for (int i2 = 0; i2 < Lim; ++i2)
				      {
					Index = this->Particles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
					if (Index < Dim)
					  {
					    vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction[i2];
					  }
					++TmpMIndices;
				      }
				  }
				++TmpNIndices;
			      }
			  }
		      }
		    if (this->FullTwoBodyFlag == true)
		      {
			int m1;
			int m3;
			if (this->NbrM12Indices == 0)
			  {
			    double TmpTwoBodyInteraction;
			    int m2;
			    int m4;
			    for (int j = 0; j < this->NbrInteractionFactors; ++j) 
			      {
				m1 = this->M1Value[j];
				m2 = this->M2Value[j];
				m3 = this->M3Value[j];
				TmpTwoBodyInteraction = this->InteractionFactors[j];
				if (this->M4Value != 0)
				  m4 = this->M4Value[j];
				else
				  m4 = m1 + m2 - m3;
				for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
				  {
				    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
				    if (Index < Dim)
				      vDestination[Index] += Coefficient * TmpTwoBodyInteraction * vSource[i];
				  }
			      }
			  }
			else
			  {
			    double Coefficient2;
			    int SumIndices;
			    int TmpNbrM3Values;
			    int* TmpM3Values;
			    int ReducedNbrInteractionFactors;
			    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			      {
				ReducedNbrInteractionFactors = 0;
				for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
				  {
				    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
				    if (Coefficient != 0.0)
				      {
					SumIndices = this->M1Value[m1] + this->M2Value[m1];
					Coefficient *= vSource[i];
					TmpNbrM3Values = this->NbrM3Values[m1];
					TmpM3Values = this->M3Values[m1];
					for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
					  {
					    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					    if (Index < Dim)			
					      vDestination[Index] += Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
					    ++ReducedNbrInteractionFactors;
					  }
				      }
				    else
				      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
				  }
			      }			    
			  }
		      }
		    for (int k = 2; k <= this->MaxNBody; ++k)
		      if (this->NBodyFlags[k] == true)
			{
			  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			    this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
			}
		    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		      vDestination[i] += this->HamiltonianShift * vSource[i];					
		  }
	      delete TmpParticles;
	    }
	  else
	    {
	      int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      int TmpNbrIteration = nbrComponent / this->BufferSize;
	      int* TmpIndexArray;
	      double* TmpCoefficientArray;
	      int TmpNbrInteraction;
	      int k = firstComponent;
	      int EffectiveHilbertSpaceDimension;
	      firstComponent -= this->PrecalculationShift;

	      ifstream File;
	      File.open(this->DiskStorageFileName, ios::binary | ios::in);
	      File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
	      long FileJump = 0;
	      for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
		FileJump += (long) this->NbrInteractionPerComponent[i];
	      FileJump *= sizeof(int);
	      long FileOffset = 0;
	      for (int i = this->DiskStorageStart; i < firstComponent; ++i)
		FileOffset += this->NbrInteractionPerComponent[i];
	      File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
	      FileJump += (sizeof(double) - sizeof(int)) * FileOffset;

	      for (int i = 0; i < TmpNbrIteration; ++i)
		{
		  int TmpPos = firstComponent;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < this->BufferSize; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);
		  
		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int l = 0; l < this->BufferSize; ++l)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      Coefficient = vSource[k];
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      vDestination[k] += this->HamiltonianShift * Coefficient;
		      ++k;
		      ++firstComponent;
		    }
		}

	      if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
		{
		  int TmpPos = firstComponent;
		  int Lim =  nbrComponent % this->BufferSize;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < Lim; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);

		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int i = 0; i < Lim; ++i)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      Coefficient = vSource[k];
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      vDestination[k] += this->HamiltonianShift * Coefficient;
		      ++k;
		      ++firstComponent;
		    }
		}

	      File.close();
	      delete[] BufferIndexArray;
	      delete[] BufferCoefficientArray;
	    }
	  
	}
    }
  if (this->L2Operator != 0)
    this->L2Operator->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereNBodyInteractionHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int* TmpMIndices;
      int* TmpNIndices;
      double* TmpInteraction;
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      int m1 = this->OneBodyMValues[j];
	      int m2 = this->OneBodyNValues[j];
	      double TmpInteraction2 = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = this->Particles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction2;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += Coefficient * vSources[l][i];
		    }
		}
	    }
	}
      if (this->NBodyFlags[1] == true)
	{
	  double Sign = this->NBodySign[1];
	  int TmpMinSumIndices = this->MinSumIndices[1];
	  int TmpMaxSumIndices = this->MaxSumIndices[1];	      
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
		{
		  TmpInteraction = this->NBodyInteractionFactors[1][j];
		  int Lim = NbrSortedIndicesPerSum[1][j];
		  TmpNIndices = this->SortedIndicesPerSum[1][j];
		  for (int i1 = 0; i1 < Lim; ++i1)
		    {
		      for (int l = 0; l < nbrVectors; ++l)
			Coefficient2[l] = Sign * vSources[l][i] * TmpInteraction[i1];
		      TmpMIndices = this->SortedIndicesPerSum[1][j];
		      for (int i2 = 0; i2 < Lim; ++i2)
			{
			  Index = this->Particles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
			  if (Index < Dim)
			    {
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += Coefficient * TmpInteraction[i2] * Coefficient2[l];
			    }
			  ++TmpMIndices;
			}
		    }
		  ++TmpNIndices;
		}
	    }
	}
     if (this->FullTwoBodyFlag == true)
       {
	 int m1;
	 int m3;
	 if (this->NbrM12Indices == 0)
	   {	     
	     double TmpTwoBodyInteraction;
	     int m2;
	     int m4;
	     for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	       {
		 m1 = this->M1Value[j];
		 m2 = this->M2Value[j];
		 m3 = this->M3Value[j];
		 if (this->M4Value != 0)
		   m4 = this->M4Value[j];
		 else
		   m4 = m1 + m2 - m3;
		 TmpTwoBodyInteraction = this->InteractionFactors[j];
		 for (int i = firstComponent; i < LastComponent; ++i)
		   {
		     Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		     if (Index < Dim)
		       {
			 Coefficient *= TmpTwoBodyInteraction;
			 for (int l = 0; l < nbrVectors; ++l)
			   vDestinations[l][Index] += Coefficient * vSources[l][i];
		       }
		   }
	       }
	   }
	 else
	   {
	     double Coefficient2;
	     int SumIndices;
	     int TmpNbrM3Values;
	     int* TmpM3Values;
	     double* TmpCoefficients = new double[nbrVectors];
	     int ReducedNbrInteractionFactors;
	     for (int i = firstComponent; i < LastComponent; ++i)
	       {
		 ReducedNbrInteractionFactors = 0;
		 for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		   {
		     Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		     if (Coefficient != 0.0)
		       {
			 SumIndices = this->M1Value[m1] + this->M2Value[m1];
			 TmpNbrM3Values = this->NbrM3Values[m1];
			 TmpM3Values = this->M3Values[m1];
			 for (int l = 0; l < nbrVectors; ++l)
			   TmpCoefficients[l] = Coefficient * vSources[l][i];
			 for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			   {
			     Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			     if (Index < Dim)
			       for (int l = 0; l < nbrVectors; ++l)
				 vDestinations[l][Index] += TmpCoefficients[l] * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			     ++ReducedNbrInteractionFactors;
			   }
		       }
		     else
		       ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		   }	     
	       }
	     delete[] TmpCoefficients;
	   }
       }
     for (int k = 2; k <= this->MaxNBody; ++k)
       if (this->NBodyFlags[k] == true)
	 {
	   for (int i = firstComponent; i < LastComponent; ++i)
	     this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2) ;
	  }	    
      for (int l = 0; l < nbrVectors; ++l)
	{
	  RealVector& TmpSourceVector = vSources[l];
	  RealVector& TmpDestinationVector = vDestinations[l];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      delete TmpParticles;
      delete[] Coefficient2;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  double* Coefficient2 = new double [nbrVectors];
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	      int* TmpIndexArray;
	      double* TmpCoefficientArray; 
	      double* Coefficient2 = new double [nbrVectors];
	      int j;
	      int Pos2;
	      int TmpNbrInteraction;
	      firstComponent -= this->PrecalculationShift;
	      LastComponent -= this->PrecalculationShift;
	      int Pos = firstComponent / this->FastMultiplicationStep; 
	      int PosMod = firstComponent % this->FastMultiplicationStep;
	      if (PosMod != 0)
		{
		  ++Pos;
		  PosMod = this->FastMultiplicationStep - PosMod;
		}
	      int l =  PosMod + firstComponent + this->PrecalculationShift;
	      for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
		  TmpIndexArray = this->InteractionPerComponentIndex[Pos];
		  TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      Coefficient2[k] = vSources[k][l];
		    }
		  for (j = 0; j < TmpNbrInteraction; ++j)
		    {
		      Pos2 = TmpIndexArray[j];
		      Coefficient = TmpCoefficientArray[j];
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
			}
		    }
		  l += this->FastMultiplicationStep;
		  ++Pos;
		}
	      int Index;
	      int* TmpMIndices;
	      int* TmpNIndices;
	      double* TmpInteraction;
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    if (this->OneBodyTermFlag == true)
		      {
			for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
			  {
			    int m1 = this->OneBodyMValues[j];
			    int m2 = this->OneBodyNValues[j];
			    double TmpInteraction2 = this->OneBodyInteractionFactors[j];
			    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			      {
				Index = this->Particles->AdA(i, m1, m2, Coefficient);
				if (Index < Dim)
				  {
				    Coefficient *= TmpInteraction2;
				    for (int p = 0; p < nbrVectors; ++p)
				      vDestinations[p][Index] += Coefficient * vSources[p][i];
				  }
			      }
			  }
		      }
		    if (this->NBodyFlags[1] == true)
		      {
			double Sign = this->NBodySign[1];
			int TmpMinSumIndices = this->MinSumIndices[1];
			int TmpMaxSumIndices = this->MaxSumIndices[1];	      
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  {
			    for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
			      {
				int Lim = NbrSortedIndicesPerSum[1][j];
				TmpInteraction = this->NBodyInteractionFactors[1][j];
				TmpNIndices = this->SortedIndicesPerSum[1][j];
				for (int i1 = 0; i1 < Lim; ++i1)
				  {
				    for (int l = 0; l < nbrVectors; ++l)
				      Coefficient2[l] = Sign * vSources[l][i] * TmpInteraction[i1];
				    TmpMIndices = this->SortedIndicesPerSum[1][j];
				    for (int i2 = 0; i2 < Lim; ++i2)
				      {
					Index = this->Particles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
					if (Index < Dim)
					  {
					    for (int p = 0; p < nbrVectors; ++p)
					      vDestinations[p][Index] += Coefficient * TmpInteraction[i2] * Coefficient2[p];
					  }
					++TmpMIndices;
				      }
				  }
				++TmpNIndices;
			      }
			  }
		      }
		    if (this->FullTwoBodyFlag == true)
		      {
			int m1;
			int m3;
			if (this->NbrM12Indices == 0)
			  {
			    double TmpTwoBodyInteraction;
			    int m2;
			    int m4;
			    for (int j = 0; j < this->NbrInteractionFactors; ++j) 
			      {
				m1 = this->M1Value[j];
				m2 = this->M2Value[j];
				m3 = this->M3Value[j];
				TmpTwoBodyInteraction = this->InteractionFactors[j];
				if (this->M4Value != 0)
				  m4 = this->M4Value[j];
				else
				  m4 = m1 + m2 - m3;
				for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
				  {
				    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
				    if (Index < Dim)
				      {
					Coefficient *= TmpTwoBodyInteraction;
					for (int p = 0; p < nbrVectors; ++p)
					  vDestinations[p][Index] += Coefficient * vSources[p][i];
				      }
				  }
			      }
			  }
			else
			  {
			    double Coefficient2;
			    int SumIndices;
			    int TmpNbrM3Values;
			    int* TmpM3Values;
			    double* TmpCoefficients = new double[nbrVectors];
			    int ReducedNbrInteractionFactors;
			    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			      {
				ReducedNbrInteractionFactors = 0;
				for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
				  {
				    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
				    if (Coefficient != 0.0)
				      {
					SumIndices = this->M1Value[m1] + this->M2Value[m1];
					TmpNbrM3Values = this->NbrM3Values[m1];
					TmpM3Values = this->M3Values[m1];
					for (int l = 0; l < nbrVectors; ++l)
					  TmpCoefficients[l] = Coefficient * vSources[l][i];
					for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
					  {
					    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					    if (Index < Dim)
					      for (int p = 0; p < nbrVectors; ++p)
						vDestinations[p][Index] += TmpCoefficients[p] * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
					    ++ReducedNbrInteractionFactors;
					  }
				      }
				    else
				      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
				  }
			      }
			    delete[] TmpCoefficients;
			  }
		      }
		    for (int k = 2; k <= this->MaxNBody; ++k)
		      if (this->NBodyFlags[k] == true)
			{
			  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			    this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2) ;
			}
		  }
	      for (int p = 0; p < nbrVectors; ++p)
		{
		  RealVector& TmpSourceVector = vSources[p];
		  RealVector& TmpDestinationVector = vDestinations[p];
		  for (int i = firstComponent; i < LastComponent; ++i)
		    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		}
	      delete[] Coefficient2;
	      delete TmpParticles;
	    }
	  else
	    {
	      int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      double* Coefficient2 = new double [nbrVectors];
	      int TmpNbrIteration = nbrComponent / this->BufferSize;
	      int* TmpIndexArray;
	      double* TmpCoefficientArray;
	      int TmpNbrInteraction;
	      int k = firstComponent;
	      int EffectiveHilbertSpaceDimension;
	      int Pos;
	      firstComponent -= this->PrecalculationShift;

	      ifstream File;
	      File.open(this->DiskStorageFileName, ios::binary | ios::in);
	      File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
	      long FileJump = 0;
	      for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
		FileJump += (long) this->NbrInteractionPerComponent[i];
	      FileJump *= sizeof(int);
	      long FileOffset = 0;
	      for (int i = this->DiskStorageStart; i < firstComponent; ++i)
		FileOffset += this->NbrInteractionPerComponent[i];
	      File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
	      FileJump += (sizeof(double) - sizeof(int)) * FileOffset;

	      for (int i = 0; i < TmpNbrIteration; ++i)
		{
		  int TmpPos = firstComponent;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < this->BufferSize; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);
		  
		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int m = 0; m < this->BufferSize; ++m)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      for (int l = 0; l < nbrVectors; ++l)
			{
			  Coefficient2[l] = vSources[l][k];
			  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
			}
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    {
			      Pos = TmpIndexArray[j];
			      Coefficient = TmpCoefficientArray[j];
			      for (int l = 0; l < nbrVectors; ++l)
				{
				  vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
				}
			    }
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      ++k;
		      ++firstComponent;
		    }
		}

	      if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
		{
		  int TmpPos = firstComponent;
		  int Lim =  nbrComponent % this->BufferSize;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < Lim; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);

		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int m = 0; m < Lim; ++m)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      for (int l = 0; l < nbrVectors; ++l)
			{
			  Coefficient2[l] = vSources[l][k];
			  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
			}
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    {
			      Pos = TmpIndexArray[j];
			      Coefficient = TmpCoefficientArray[j];
			      for (int l = 0; l < nbrVectors; ++l)
				{
				  vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
				}
			    }
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      ++k;
		      ++firstComponent;
		    }
		}

	      File.close();
	      delete[] BufferIndexArray;
	      delete[] BufferCoefficientArray;
	      delete[] Coefficient2;
	    }
	}
    }
  if (this->L2Operator != 0)
    for (int l = 0; l < nbrVectors; ++l)
      this->L2Operator->LowLevelAddMultiply(vSources[l], vDestinations[l], firstComponent, nbrComponent);
  return vDestinations;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereNBodyInteractionHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  double Coefficient2;
  long Memory = 0;
  int* TmpMIndices;
  int* TmpNIndices;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;

  for (int i = firstComponent; i < LastComponent; ++i)
    {	      
      if (this->OneBodyTermFlag == true)
      {
        for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	  {
	    int m1 = this->OneBodyMValues[j];
	    int m2 = this->OneBodyNValues[j];
	    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
	    if (Index < this->Particles->GetHilbertSpaceDimension())
	     {
	       ++Memory;
	       ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
	     }
	  }
      }


     if (this->NBodyFlags[1] == true)
      {
        int TmpMinSumIndices = this->MinSumIndices[1];
        int TmpMaxSumIndices = this->MaxSumIndices[1];	      
        for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
	 {
	  int Lim = this->NbrSortedIndicesPerSum[1][j];
	  TmpNIndices = this->SortedIndicesPerSum[1][j];
	  for (int i1 = 0; i1 < Lim; ++i1)
	   {
	     TmpMIndices = this->SortedIndicesPerSum[1][j];
	     for (int i2 = 0; i2 < Lim; ++i2)
	      {
	        Index = TmpParticles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
	        if (Index < this->Particles->GetHilbertSpaceDimension())
	  	  {
		    ++Memory;
		    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		  }
		++TmpMIndices;
	       }
	   }
	  ++TmpNIndices;
         }
      }



    if (this->FullTwoBodyFlag == true)
      {
        int m1;
        int m3;
        if (this->NbrM12Indices == 0)
	  {
	    int m2;
	    int m4;      
	    for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	      {
	        m1 = this->M1Value[j];
	        m2 = this->M2Value[j];
	        m3 = this->M3Value[j];
		if (this->M4Value != 0)
		  m4 = this->M4Value[j];
		else
		  m4 = m1 + m2 - m3;
	        Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		if (Index < this->Particles->GetHilbertSpaceDimension())
		 {
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		 }		
	       }
	   }
        else
	  {
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
            for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
	     {
	       Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
	       if (Coefficient != 0.0)
		{
		  SumIndices = this->M1Value[m1] + this->M2Value[m1];
		  TmpM3Values = this->M3Values[m1];
		  TmpNbrM3Values = this->NbrM3Values[m1];
		  for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
		   {
		     if (TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
		       {
		         ++Memory;
		         ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		       }
		    }    
		 }
	      }
	   }
        }


   for (int k = 2; k <= this->MaxNBody; ++k)
      if (this->NBodyFlags[k] == true)
        {
	  if (this->MNNBodyInteractionFactors == 0)
	    {
	      int TmpMinSumIndices = this->MinSumIndices[k];
	      int TmpMaxSumIndices = this->MaxSumIndices[k];	      
 	      for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
		{
		  int Lim = NbrSortedIndicesPerSum[k][j];
		  TmpNIndices = this->SortedIndicesPerSum[k][j];
		  for (int i1 = 0; i1 < Lim; ++i1)
		   {
		     Coefficient2 = TmpParticles->ProdA(i, TmpNIndices, k);
		     if (Coefficient2 != 0.0)
		       {
			 TmpMIndices = this->SortedIndicesPerSum[k][j];
			 for (int i2 = 0; i2 < Lim; ++i2)
			  {
			    Index = TmpParticles->ProdAd(TmpMIndices, k, Coefficient);
			    if (Index < this->Particles->GetHilbertSpaceDimension())
			      {
				++Memory;
				++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			      }
			    TmpMIndices += k;
			  }
		        }
		      TmpNIndices += k;
		  }
	        }
	      }
	 else
	  {
	    long TmpNbrNIndices = this->NbrNIndices[k];
 	    TmpNIndices = this->NIndices[k];
	    for (long j = 0; j < TmpNbrNIndices; ++j)
	     {
	       Coefficient2 = TmpParticles->ProdA(i, TmpNIndices, k);
		if (Coefficient2 != 0.0)
		 {
		   long TmpNbrMIndices = this->NbrMIndices[k][j];
		   TmpMIndices = this->MIndices[k][j];
		   for (long i2 = 0; i2 < TmpNbrMIndices; ++i2)
		    {
		      int Index = TmpParticles->ProdAd(TmpMIndices, k, Coefficient);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
		        {
		   	  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		       TmpMIndices += k;	      
	             }
		  }
	        TmpNIndices += k;
	     }
	  }
      }


   } //i ranges over Hilbert space

  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnSphereNBodyInteractionHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];
  for (int i = 0; i < ReducedSpaceDimension; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
    }

  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);

  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnSphereNBodyInteractionHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long TotalPos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++TotalPos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }

  int Index;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  double* TmpInteraction;
  int* TmpMIndices;
  int* TmpNIndices;
  int m1;
  int m2;
  int m3;
  int m4;
  
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];

      int Pos = 0;
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      int m1 = this->OneBodyMValues[j];
	      int m2 = this->OneBodyNValues[j];
	      Index = TmpParticles->AdA(i, m1, m2, Coefficient);
	      if (Index < TmpParticles->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Pos] = Index;
		  TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		  ++Pos;
		}
	    }
	}
      if (this->NBodyFlags[1] == true)
	{
	  double Sign = this->NBodySign[1];
	  int TmpMinSumIndices = this->MinSumIndices[1];
	  int TmpMaxSumIndices = this->MaxSumIndices[1];	      
	  for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
	    {
	      int Lim = this->NbrSortedIndicesPerSum[1][j];
	      TmpNIndices = this->SortedIndicesPerSum[1][j];
	      TmpInteraction = this->NBodyInteractionFactors[1][j];
	      for (int i1 = 0; i1 < Lim; ++i1)
		{
		  TmpMIndices = this->SortedIndicesPerSum[1][j];
		  for (int i2 = 0; i2 < Lim; ++i2)
		    {
		      Index = TmpParticles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
		      if (Index < TmpParticles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Sign * Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
			  ++Pos;
			}
		      ++TmpMIndices;
		    }
		}
	      ++TmpNIndices;
	    }
	}
      if (this->FullTwoBodyFlag == true)
	{
	  if (this->NbrM12Indices == 0)
	    {
	      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		{
		  m1 = this->M1Value[j];
		  m2 = this->M2Value[j];
		  m3 = this->M3Value[j];
		  if (this->M4Value != 0)
		    m4 = this->M4Value[j];
		  else
		    m4 = m1 + m2 - m3;
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < TmpParticles->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
	  else
	    {
	      double Coefficient2;
	      int SumIndices;
	      int TmpNbrM3Values;
	      int* TmpM3Values;
	      int ReducedNbrInteractionFactors = 0;
	      //Pos = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      TmpM3Values = this->M3Values[m1];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < TmpParticles->GetHilbertSpaceDimension())
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->InteractionFactors[ReducedNbrInteractionFactors];
			      ++Pos;
			    }		      
			  ++ReducedNbrInteractionFactors;
			}    
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	    }
	}
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  this->EvaluateMNNBodyFastMultiplicationComponent(TmpParticles, i + this->PrecalculationShift, k, TmpIndexArray, TmpCoefficientArray, Pos);
      ++TotalPos;
    }
  delete TmpParticles;
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnSphereNBodyInteractionHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  if (this->FastMultiplicationStep == 1)
    {
      this->DiskStorageFlag = false;
      this->DiskStorageFileName = 0;
      this->EnableFastMultiplication();
      return;
    }

  this->DiskStorageFlag = true;
  this->DiskStorageFileName = new char [strlen(fileName) + 8];
  sprintf (this->DiskStorageFileName, "%s.ham", fileName);

  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->DiskStorageStart = (int) MinIndex;
  int DiskStorageEnd = 1 + (int) MaxIndex;

  int Index;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  double* TmpInteraction;
  int* TmpMIndices;
  int* TmpNIndices;
  int m1;
  int m2;
  int m3;
  int m4;
  this->MaxNbrInteractionPerComponent = 0;

  long TotalPos = 0l;
  ofstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::out);
 
  File.write((char*) &(EffectiveHilbertSpaceDimension), sizeof(int));
  File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
  File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * EffectiveHilbertSpaceDimension);

  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      FileJump += (long) this->NbrInteractionPerComponent[i];
      if (this->MaxNbrInteractionPerComponent < this->NbrInteractionPerComponent[i])
	this->MaxNbrInteractionPerComponent = this->NbrInteractionPerComponent[i];
    }
  FileJump *= sizeof(int);

  TmpIndexArray = new int [this->MaxNbrInteractionPerComponent];
  TmpCoefficientArray = new double [this->MaxNbrInteractionPerComponent];      

  for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
    {
      if (this->NbrInteractionPerComponent[TotalPos] > 0)
	{
	  Pos = 0;
	  if (this->OneBodyTermFlag == true)
	    {
	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
		{
		  int m1 = this->OneBodyMValues[j];
		  int m2 = this->OneBodyNValues[j];
		  Index = this->Particles->AdA(i, m1, m2, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
	  if (this->NBodyFlags[1] == true)
	    {
	      double Sign = this->NBodySign[1];
	      int TmpMinSumIndices = this->MinSumIndices[1];
	      int TmpMaxSumIndices = this->MaxSumIndices[1];	      
	      for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
		{
		  int Lim = NbrSortedIndicesPerSum[1][j];
		  TmpNIndices = this->SortedIndicesPerSum[1][j];
		  TmpInteraction = this->NBodyInteractionFactors[1][j];
		  for (int i1 = 0; i1 < Lim; ++i1)
		    {
		      TmpMIndices = this->SortedIndicesPerSum[1][j];
		      for (int i2 = 0; i2 < Lim; ++i2)
			{
			  Index = this->Particles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
			  if (Index < this->Particles->GetHilbertSpaceDimension())
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Sign * Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
			      ++Pos;
			    }
			  ++TmpMIndices;
			}
		    }
		  ++TmpNIndices;
		}
	    }
	  if (this->FullTwoBodyFlag == true)
	    {
	      if (this->NbrM12Indices == 0)
		{
		  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		    {
		      m1 = this->M1Value[j];
		      m2 = this->M2Value[j];
		      m3 = this->M3Value[j];
		      if (this->M4Value != 0)
			m4 = this->M4Value[j];
		      else
			m4 = m1 + m2 - m3;
		      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
			  ++Pos;
			}
		    }
		}
	      else
		{
		  double Coefficient2;
		  int SumIndices;
		  int TmpNbrM3Values;
		  int* TmpM3Values;
		  int ReducedNbrInteractionFactors = 0;
		  for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		    {
		      Coefficient = this->Particles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		      if (Coefficient != 0.0)
			{
			  SumIndices = this->M1Value[m1] + this->M2Value[m1];
			  TmpM3Values = this->M3Values[m1];
			  TmpNbrM3Values = this->NbrM3Values[m1];
			  for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			    {
			      Index = this->Particles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			      if (Index < this->Particles->GetHilbertSpaceDimension())
				{
				  TmpIndexArray[Pos] = Index;
				  TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
				  ++Pos;
				}
			      ++ReducedNbrInteractionFactors;
			    }    
			}
		      else
			ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		    }
		}
	    }
	  for (int k = 2; k <= this->MaxNBody; ++k)
	    if (this->NBodyFlags[k] == true)
	      this->EvaluateMNNBodyFastMultiplicationComponent(this->Particles, i, k, TmpIndexArray, TmpCoefficientArray, Pos);
	  File.write((char*) TmpIndexArray, sizeof(int) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump -= sizeof(int) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(FileJump, ios::cur);
	  File.write((char*) TmpCoefficientArray, sizeof(double) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump += sizeof(double) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(-FileJump, ios::cur);	  
	}
      ++TotalPos;
    }
  delete[] TmpIndexArray;
  delete[] TmpCoefficientArray;
  File.close();

  this->FastMultiplicationFlag = true;
  this->BufferSize = this->Memory / ((this->MaxNbrInteractionPerComponent * (sizeof(int) + sizeof(double))) + sizeof(int*) + sizeof(double*));

  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnSphereNBodyInteractionHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
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

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnSphereNBodyInteractionHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
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

  int MaxSum = (nbrValues - 1) * nbrIndices;
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

