////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(3) spin on a sphere           //
//                                                                            //
//                        last modification : 20/01/2008                      //
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
#include "Hamiltonian/AbstractQHEOnSphereWithSU3SpinHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <fstream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

AbstractQHEOnSphereWithSU3SpinHamiltonian::~AbstractQHEOnSphereWithSU3SpinHamiltonian()
{
  if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
    {
      delete[] this->NbrIntraSectorIndicesPerSum;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->IntraSectorIndicesPerSum[i];
	  delete[] this->InteractionFactors11[i];
	  delete[] this->InteractionFactors22[i];
	  delete[] this->InteractionFactors33[i];
	}
      delete[] this->IntraSectorIndicesPerSum;
      delete[] this->NbrInterSectorIndicesPerSum;
      delete[] this->InteractionFactors11;
      delete[] this->InteractionFactors22;
      delete[] this->InteractionFactors33;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  delete[] this->InterSectorIndicesPerSum[i];      
	  delete[] this->InteractionFactors12[i];
	  delete[] this->InteractionFactors13[i];
	  delete[] this->InteractionFactors23[i];
	}
      delete[] this->InteractionFactors12;
      delete[] this->InteractionFactors13;
      delete[] this->InteractionFactors23;
      delete[] this->InterSectorIndicesPerSum;
    }
  else
    {
      delete[] this->M1IntraValue;
      delete[] this->M2IntraValue;
      for (int i = 0; i < this->NbrM12IntraIndices; ++i)
	delete[] this->M3IntraValues[i];
      delete[] this->M3IntraValues;
      delete[] this->NbrM3IntraValues;

      delete[] this->M1InterValue;
      delete[] this->M2InterValue;
      for (int i = 0; i < this->NbrM12InterIndices; ++i)
	delete[] this->M3InterValues[i];
      delete[] this->M3InterValues;
      delete[] this->NbrM3InterValues;

      delete[] this->M12InteractionFactors11;
      delete[] this->M12InteractionFactors22;
      delete[] this->M12InteractionFactors33;
      delete[] this->M12InteractionFactors12;
      delete[] this->M12InteractionFactors13;
      delete[] this->M12InteractionFactors23;
    }
  if (this->OneBodyInteractionFactors11 != 0)
    delete[] this->OneBodyInteractionFactors11;
  if (this->OneBodyInteractionFactors22 != 0)
    delete[] this->OneBodyInteractionFactors22;
  if (this->OneBodyInteractionFactors33 != 0)
    delete[] this->OneBodyInteractionFactors33;  
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereWithSU3SpinHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int* TmpIndices;
      double* TmpInteractionFactor;
      double Coefficient3;
      ParticleOnSphereWithSU3Spin* TmpParticles = (ParticleOnSphereWithSU3Spin*) this->Particles->Clone();
      if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		{
		  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		  TmpIndices = this->IntraSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      Coefficient3 = TmpParticles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
			  Coefficient3 *= vSource[i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
			  Coefficient3 *= vSource[i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
			  Coefficient3 *= vSource[i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int j = 0; j < this->NbrInterSectorSums; ++j)
		{
		  int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
		  TmpIndices = this->InterSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      Coefficient3 = TmpParticles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
			  Coefficient3 *= vSource[i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
			  Coefficient3 *= vSource[i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
			  Coefficient3 *= vSource[i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
			      ++TmpInteractionFactor;
			    }
			}
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
	  int ReducedNbrInteractionFactors;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		{
		  Coefficient = TmpParticles->A1A1(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3IntraValues[m1];
		      TmpM3Values = this->M3IntraValues[m1];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad1Ad1(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->M12InteractionFactors11[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		{
		  Coefficient = TmpParticles->A2A2(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3IntraValues[m1];
		      TmpM3Values = this->M3IntraValues[m1];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad2Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->M12InteractionFactors22[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		{
		  Coefficient = TmpParticles->A3A3(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3IntraValues[m1];
		      TmpM3Values = this->M3IntraValues[m1];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad3Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->M12InteractionFactors33[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		{
		  Coefficient = TmpParticles->A1A2(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3InterValues[m1];
		      TmpM3Values = this->M3InterValues[m1];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad1Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->M12InteractionFactors12[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		{
		  Coefficient = TmpParticles->A1A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3InterValues[m1];
		      TmpM3Values = this->M3InterValues[m1];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad1Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->M12InteractionFactors13[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		{
		  Coefficient = TmpParticles->A2A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3InterValues[m1];
		      TmpM3Values = this->M3InterValues[m1];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad2Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->M12InteractionFactors23[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		}
	    }	  
	}
      if (this->OneBodyInteractionFactors11 != 0) 
	{
	  double TmpDiagonal = 0.0;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    { 
	      TmpDiagonal = 0.0;
	      for (int j = 0; j <= this->LzMax; ++j) 
		{
		  TmpDiagonal += this->OneBodyInteractionFactors11[j] * TmpParticles->Ad1A1(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactors22[j] * TmpParticles->Ad2A2(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactors33[j] * TmpParticles->Ad3A3(i, j);
		}
	      vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
	    }
	}
      else
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
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphereWithSU3Spin* TmpParticles = (ParticleOnSphereWithSU3Spin*) this->Particles->Clone();
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
	      int* TmpIndices;
	      double* TmpInteractionFactor;
	      double Coefficient3;
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
		      {
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  {
			    for (int j = 0; j < this->NbrIntraSectorSums; ++j)
			      {
				int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
				TmpIndices = this->IntraSectorIndicesPerSum[j];
				for (int i1 = 0; i1 < Lim; i1 += 2)
				  {
				    Coefficient3 = TmpParticles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
				    if (Coefficient3 != 0.0)
				      {
					TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
					Coefficient3 *= vSource[i];
					for (int i2 = 0; i2 < Lim; i2 += 2)
					  {
					    Index = TmpParticles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
					    if (Index < Dim)
					      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
					    ++TmpInteractionFactor;
					  }
				      }
				    Coefficient3 = TmpParticles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
				    if (Coefficient3 != 0.0)
				      {
					TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
					Coefficient3 *= vSource[i];
					for (int i2 = 0; i2 < Lim; i2 += 2)
					  {
					    Index = TmpParticles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
					    if (Index < Dim)
					      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
					    ++TmpInteractionFactor;
					  }
				      }
				    Coefficient3 = TmpParticles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
				    if (Coefficient3 != 0.0)
				      {
					TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
					Coefficient3 *= vSource[i];
					for (int i2 = 0; i2 < Lim; i2 += 2)
					  {
					    Index = TmpParticles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
					    if (Index < Dim)
					      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
					    ++TmpInteractionFactor;
					  }
				      }
				  }
			      }
			    for (int j = 0; j < this->NbrInterSectorSums; ++j)
			      {
				int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
				TmpIndices = this->InterSectorIndicesPerSum[j];
				for (int i1 = 0; i1 < Lim; i1 += 2)
				  {
				    Coefficient3 = TmpParticles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
				    if (Coefficient3 != 0.0)
				      {
					TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
					Coefficient3 *= vSource[i];
					for (int i2 = 0; i2 < Lim; i2 += 2)
					  {
					    Index = TmpParticles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
					    if (Index < Dim)
					      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
					    ++TmpInteractionFactor;
					  }
				      }
				    Coefficient3 = TmpParticles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
				    if (Coefficient3 != 0.0)
				      {
					TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
					Coefficient3 *= vSource[i];
					for (int i2 = 0; i2 < Lim; i2 += 2)
					  {
					    Index = TmpParticles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
					    if (Index < Dim)
					      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
					    ++TmpInteractionFactor;
					  }
				  }
				    Coefficient3 = TmpParticles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
				    if (Coefficient3 != 0.0)
				      {
					TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
					Coefficient3 *= vSource[i];
					for (int i2 = 0; i2 < Lim; i2 += 2)
					  {
					    Index = TmpParticles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
					    if (Index < Dim)
					      vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
					    ++TmpInteractionFactor;
					  }
				      }
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
			int ReducedNbrInteractionFactors;
			int m1;
			int m3;
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  {
			    ReducedNbrInteractionFactors = 0;
			    for (m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
			      {
				Coefficient = TmpParticles->A1A1(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
				if (Coefficient != 0.0)
				  {
				    SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
				    Coefficient *= vSource[i];
				    TmpNbrM3Values = this->NbrM3IntraValues[m1];
				    TmpM3Values = this->M3IntraValues[m1];
				    for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
				      {
					Index = TmpParticles->Ad1Ad1(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					if (Index < Dim)			
					  vDestination[Index] += Coefficient * this->M12InteractionFactors11[ReducedNbrInteractionFactors] * Coefficient2;
					++ReducedNbrInteractionFactors;
				      }
				  }
				else
				  ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
			      }
			    ReducedNbrInteractionFactors = 0;
			    for (m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
			      {
				Coefficient = TmpParticles->A2A2(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
				if (Coefficient != 0.0)
				  {
				    SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
				    Coefficient *= vSource[i];
				    TmpNbrM3Values = this->NbrM3IntraValues[m1];
				    TmpM3Values = this->M3IntraValues[m1];
				    for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
				      {
					Index = TmpParticles->Ad2Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					if (Index < Dim)			
					  vDestination[Index] += Coefficient * this->M12InteractionFactors22[ReducedNbrInteractionFactors] * Coefficient2;
					++ReducedNbrInteractionFactors;
				      }
				  }
				else
				  ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
			      }
			    ReducedNbrInteractionFactors = 0;
			    for (m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
			      {
				Coefficient = TmpParticles->A3A3(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
				if (Coefficient != 0.0)
				  {
				    SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
				    Coefficient *= vSource[i];
				    TmpNbrM3Values = this->NbrM3IntraValues[m1];
				    TmpM3Values = this->M3IntraValues[m1];
				    for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
				      {
					Index = TmpParticles->Ad3Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					if (Index < Dim)			
					  vDestination[Index] += Coefficient * this->M12InteractionFactors33[ReducedNbrInteractionFactors] * Coefficient2;
					++ReducedNbrInteractionFactors;
				      }
				  }
				else
				  ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
			      }
			    ReducedNbrInteractionFactors = 0;
			    for (m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
			      {
				Coefficient = TmpParticles->A1A2(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
				if (Coefficient != 0.0)
				  {
				    SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
				    Coefficient *= vSource[i];
				    TmpNbrM3Values = this->NbrM3InterValues[m1];
				    TmpM3Values = this->M3InterValues[m1];
				    for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
				      {
					Index = TmpParticles->Ad1Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					if (Index < Dim)			
					  vDestination[Index] += Coefficient * this->M12InteractionFactors12[ReducedNbrInteractionFactors] * Coefficient2;
					++ReducedNbrInteractionFactors;
				      }
				  }
				else
				  ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
			      }
			    ReducedNbrInteractionFactors = 0;
			    for (m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
			      {
				Coefficient = TmpParticles->A1A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
				if (Coefficient != 0.0)
				  {
				    SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
				    Coefficient *= vSource[i];
				    TmpNbrM3Values = this->NbrM3InterValues[m1];
				    TmpM3Values = this->M3InterValues[m1];
				    for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
				      {
					Index = TmpParticles->Ad1Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					if (Index < Dim)			
					  vDestination[Index] += Coefficient * this->M12InteractionFactors13[ReducedNbrInteractionFactors] * Coefficient2;
					++ReducedNbrInteractionFactors;
				      }
				  }
				else
				  ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
			      }
			    ReducedNbrInteractionFactors = 0;
			    for (m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
			      {
				Coefficient = TmpParticles->A2A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
				if (Coefficient != 0.0)
				  {
				    SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
				    Coefficient *= vSource[i];
				    TmpNbrM3Values = this->NbrM3InterValues[m1];
				    TmpM3Values = this->M3InterValues[m1];
				    for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
				      {
					Index = TmpParticles->Ad2Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
					if (Index < Dim)			
					  vDestination[Index] += Coefficient * this->M12InteractionFactors23[ReducedNbrInteractionFactors] * Coefficient2;
					++ReducedNbrInteractionFactors;
				      }
				  }
				else
				  ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
			      }
			  }
		      }
		    if (this->OneBodyInteractionFactors11 != 0) 
		      {
			double TmpDiagonal = 0.0;
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  { 
			    TmpDiagonal = 0.0;
			    for (int j = 0; j <= this->LzMax; ++j) 
			      {
				TmpDiagonal += this->OneBodyInteractionFactors11[j] * TmpParticles->Ad1A1(i, j);
				TmpDiagonal += this->OneBodyInteractionFactors22[j] * TmpParticles->Ad2A2(i, j);
				TmpDiagonal += this->OneBodyInteractionFactors33[j] * TmpParticles->Ad3A3(i, j);
			      }
			    vDestination[i] += (this->HamiltonianShift + TmpDiagonal)* vSource[i];
			  }
		      }
		    else
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

  return vDestination;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereWithSU3SpinHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int* TmpIndices;
      double* TmpInteractionFactor;
      double Coefficient3;
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphereWithSU3Spin* TmpParticles = (ParticleOnSphereWithSU3Spin*) this->Particles->Clone();
      if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		{
		  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		  TmpIndices = this->IntraSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      Coefficient3 = TmpParticles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
			  for (int p = 0; p < nbrVectors; ++p)
			    Coefficient2[p] = Coefficient3 * vSources[p][i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
			  for (int p = 0; p < nbrVectors; ++p)
			    Coefficient2[p] = Coefficient3 * vSources[p][i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
			  ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
			  for (int p = 0; p < nbrVectors; ++p)
			    Coefficient2[p] = Coefficient3 * vSources[p][i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int j = 0; j < this->NbrInterSectorSums; ++j)
		{
		  int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
		  TmpIndices = this->InterSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      Coefficient3 = TmpParticles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
			  for (int p = 0; p < nbrVectors; ++p)
			    Coefficient2[p] = Coefficient3 * vSources[p][i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
			  for (int p = 0; p < nbrVectors; ++p)
			    Coefficient2[p] = Coefficient3 * vSources[p][i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
			      ++TmpInteractionFactor;
			    }
			}
		      Coefficient3 = TmpParticles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		      if (Coefficient3 != 0.0)
			{
			  TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
			  for (int p = 0; p < nbrVectors; ++p)
			    Coefficient2[p] = Coefficient3 * vSources[p][i];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Index = TmpParticles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			      if (Index < Dim)
				for (int p = 0; p < nbrVectors; ++p)
				  vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
			      ++TmpInteractionFactor;
			    }
			}
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
	      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		{
		  Coefficient = TmpParticles->A1A1(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		      TmpNbrM3Values = this->NbrM3IntraValues[m1];
		      TmpM3Values = this->M3IntraValues[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad1Ad1(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->M12InteractionFactors11[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		{
		  Coefficient = TmpParticles->A2A2(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		      TmpNbrM3Values = this->NbrM3IntraValues[m1];
		      TmpM3Values = this->M3IntraValues[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad2Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->M12InteractionFactors22[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		{
		  Coefficient = TmpParticles->A3A3(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		      TmpNbrM3Values = this->NbrM3IntraValues[m1];
		      TmpM3Values = this->M3IntraValues[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad3Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->M12InteractionFactors33[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		{
		  Coefficient = TmpParticles->A1A2(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		      TmpNbrM3Values = this->NbrM3InterValues[m1];
		      TmpM3Values = this->M3InterValues[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad1Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->M12InteractionFactors12[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		{
		  Coefficient = TmpParticles->A1A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		      TmpNbrM3Values = this->NbrM3InterValues[m1];
		      TmpM3Values = this->M3InterValues[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad1Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->M12InteractionFactors13[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		}
	      ReducedNbrInteractionFactors = 0;
	      for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		{
		  Coefficient = TmpParticles->A2A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		      TmpNbrM3Values = this->NbrM3InterValues[m1];
		      TmpM3Values = this->M3InterValues[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->Ad2Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->M12InteractionFactors23[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		}
	    }	  
	}
      if (this->OneBodyInteractionFactors11 != 0) 
	{
	  double TmpDiagonal = 0.0;
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      RealVector& TmpSourceVector = vSources[l];
	      RealVector& TmpDestinationVector = vDestinations[l];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{ 
		  TmpDiagonal = 0.0;
		  for (int j = 0; j <= this->LzMax; ++j) 
		    {
		      TmpDiagonal += this->OneBodyInteractionFactors11[j] * TmpParticles->Ad1A1(i, j);
		      TmpDiagonal += this->OneBodyInteractionFactors22[j] * TmpParticles->Ad2A2(i, j);
		      TmpDiagonal += this->OneBodyInteractionFactors33[j] * TmpParticles->Ad3A3(i, j);
		    }
		  TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
		}
	    }
	}
      else
	for (int l = 0; l < nbrVectors; ++l)
	  {
	    RealVector& TmpSourceVector = vSources[l];
	    RealVector& TmpDestinationVector = vDestinations[l];
	    for (int i = firstComponent; i < LastComponent; ++i)
	      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	  }
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* Coefficient2 = new double [nbrVectors];
	  double* TmpCoefficientArray; 
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
		    vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	  else
	    this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereWithSU3SpinHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
												   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphereWithSU3Spin* TmpParticles = (ParticleOnSphereWithSU3Spin*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  double Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  double* Coefficient2 = new double [nbrVectors];
  int Dim = this->Particles->GetHilbertSpaceDimension();
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
	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int* TmpIndices;
  double* TmpInteractionFactor;
  double Coefficient3;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (l = 0; l < this->FastMultiplicationStep; ++l)
    if (PosMod != l)
      {	
	if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
	  {
	    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		  {
		    int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		    TmpIndices = this->IntraSectorIndicesPerSum[j];
		    for (int i1 = 0; i1 < Lim; i1 += 2)
		      {
			Coefficient3 = TmpParticles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
			if (Coefficient3 != 0.0)
			  {
			    TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
			    for (int p = 0; p < nbrVectors; ++p)
			      Coefficient2[p] = Coefficient3 * vSources[p][i];
			    for (int i2 = 0; i2 < Lim; i2 += 2)
			      {
				Index = TmpParticles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
				if (Index < Dim)
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
				++TmpInteractionFactor;
			      }
			  }
			Coefficient3 = TmpParticles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
			if (Coefficient3 != 0.0)
			  {
			    TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
			    for (int p = 0; p < nbrVectors; ++p)
			      Coefficient2[p] = Coefficient3 * vSources[p][i];
			    for (int i2 = 0; i2 < Lim; i2 += 2)
			      {
				Index = TmpParticles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
				if (Index < Dim)
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
				++TmpInteractionFactor;
			      }
			  }
			Coefficient3 = TmpParticles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
			if (Coefficient3 != 0.0)
			  {
			    TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
			    for (int p = 0; p < nbrVectors; ++p)
			      Coefficient2[p] = Coefficient3 * vSources[p][i];
			    for (int i2 = 0; i2 < Lim; i2 += 2)
			      {
				Index = TmpParticles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
				if (Index < Dim)
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
				++TmpInteractionFactor;
			      }
			  }
		      }
		  }
		for (int j = 0; j < this->NbrInterSectorSums; ++j)
		  {
		    int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
		    TmpIndices = this->InterSectorIndicesPerSum[j];
		    for (int i1 = 0; i1 < Lim; i1 += 2)
		      {
			Coefficient3 = TmpParticles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
			if (Coefficient3 != 0.0)
			  {
			    TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
			    for (int p = 0; p < nbrVectors; ++p)
			      Coefficient2[p] = Coefficient3 * vSources[p][i];
			    for (int i2 = 0; i2 < Lim; i2 += 2)
			      {
				Index = TmpParticles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
				if (Index < Dim)
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
				++TmpInteractionFactor;
			      }
			  }
			Coefficient3 = TmpParticles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
			if (Coefficient3 != 0.0)
			  {
			    TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
			    for (int p = 0; p < nbrVectors; ++p)
			      Coefficient2[p] = Coefficient3 * vSources[p][i];
			    for (int i2 = 0; i2 < Lim; i2 += 2)
			      {
				Index = TmpParticles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
				if (Index < Dim)
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
				++TmpInteractionFactor;
			      }
			  }
			Coefficient3 = TmpParticles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
			if (Coefficient3 != 0.0)
			  {
			    TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
			    for (int p = 0; p < nbrVectors; ++p)
			      Coefficient2[p] = Coefficient3 * vSources[p][i];
			    for (int i2 = 0; i2 < Lim; i2 += 2)
			      {
				Index = TmpParticles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
				if (Index < Dim)
				  for (int p = 0; p < nbrVectors; ++p)
				    vDestinations[p][Index] += Coefficient * (*TmpInteractionFactor) * Coefficient2[p];
				++TmpInteractionFactor;
			      }
			  }
		      }
		  }
	      }
	  }
	else
	  {
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    int ReducedNbrInteractionFactors;
	    int m1;
	    int m3;
	    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		  {
		    Coefficient3 = TmpParticles->A1A1(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		    if (Coefficient3 != 0.0)
		      {
			SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
			TmpNbrM3Values = this->NbrM3IntraValues[m1];
			TmpM3Values = this->M3IntraValues[m1];
			for (int p = 0; p < nbrVectors; ++p)
			  Coefficient2[p] = Coefficient3 * vSources[p][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->Ad1Ad1(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient);
			    if (Index < Dim)
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->M12InteractionFactors11[ReducedNbrInteractionFactors] * Coefficient2[p];
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		  }
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		  {
		    Coefficient3 = TmpParticles->A2A2(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		    if (Coefficient3 != 0.0)
		      {
			SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
			TmpNbrM3Values = this->NbrM3IntraValues[m1];
			TmpM3Values = this->M3IntraValues[m1];
			for (int p = 0; p < nbrVectors; ++p)
			  Coefficient2[p] = Coefficient3 * vSources[p][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->Ad2Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient);
			    if (Index < Dim)
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->M12InteractionFactors22[ReducedNbrInteractionFactors] * Coefficient2[p];
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		  }
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
		  {
		    Coefficient3 = TmpParticles->A3A3(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
		    if (Coefficient3 != 0.0)
		      {
			SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
			TmpNbrM3Values = this->NbrM3IntraValues[m1];
			TmpM3Values = this->M3IntraValues[m1];
			for (int p = 0; p < nbrVectors; ++p)
			  Coefficient2[p] = Coefficient3 * vSources[p][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->Ad3Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient);
			    if (Index < Dim)
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->M12InteractionFactors33[ReducedNbrInteractionFactors] * Coefficient2[p];
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3IntraValues[m1];
		  }
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		  {
		    Coefficient3 = TmpParticles->A1A2(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		    if (Coefficient3 != 0.0)
		      {
			SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
			TmpNbrM3Values = this->NbrM3InterValues[m1];
			TmpM3Values = this->M3InterValues[m1];
			for (int p = 0; p < nbrVectors; ++p)
			  Coefficient2[p] = Coefficient3 * vSources[p][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->Ad1Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient);
			    if (Index < Dim)
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->M12InteractionFactors12[ReducedNbrInteractionFactors] * Coefficient2[p];
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		  }
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		  {
		    Coefficient3 = TmpParticles->A1A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		    if (Coefficient3 != 0.0)
		      {
			SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
			TmpNbrM3Values = this->NbrM3InterValues[m1];
			TmpM3Values = this->M3InterValues[m1];
			for (int p = 0; p < nbrVectors; ++p)
			  Coefficient2[p] = Coefficient3 * vSources[p][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->Ad1Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient);
			    if (Index < Dim)
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->M12InteractionFactors13[ReducedNbrInteractionFactors] * Coefficient2[p];
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		  }
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
		  {
		    Coefficient3 = TmpParticles->A2A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
		    if (Coefficient3 != 0.0)
		      {
			SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
			TmpNbrM3Values = this->NbrM3InterValues[m1];
			TmpM3Values = this->M3InterValues[m1];
			for (int p = 0; p < nbrVectors; ++p)
			  Coefficient2[p] = Coefficient3 * vSources[p][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->Ad2Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient);
			    if (Index < Dim)
			      for (int p = 0; p < nbrVectors; ++p)
				vDestinations[p][Index] += Coefficient * this->M12InteractionFactors23[ReducedNbrInteractionFactors] * Coefficient2[p];
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3InterValues[m1];
		  }
	      }
	  }
	if (this->OneBodyInteractionFactors11 != 0) 
	  {
	    double TmpDiagonal = 0.0;
	    for (int p = 0; p < nbrVectors; ++p)
	      {
		RealVector& TmpSourceVector = vSources[p];
		RealVector& TmpDestinationVector = vDestinations[p];
		for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		  { 
	      TmpDiagonal = 0.0;
	      for (int j = 0; j <= this->LzMax; ++j) 
		{
		  TmpDiagonal += this->OneBodyInteractionFactors11[j] * TmpParticles->Ad1A1(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactors22[j] * TmpParticles->Ad2A2(i, j);
		  TmpDiagonal += this->OneBodyInteractionFactors33[j] * TmpParticles->Ad2A2(i, j);
		}
	      TmpDestinationVector[i] += (this->HamiltonianShift + TmpDiagonal)* TmpSourceVector[i];
		  }
	      }
	  }
	else
	  for (int p = 0; p < nbrVectors; ++p)
	    {
	      RealVector& TmpSourceVector = vSources[p];
	      RealVector& TmpDestinationVector = vDestinations[p];
	      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	    }
      }
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereWithSU3SpinHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  long Memory = 0;
  int* TmpIndices;
  ParticleOnSphereWithSU3Spin* TmpParticles = (ParticleOnSphereWithSU3Spin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();

  if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = TmpParticles->A1A1(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		  Coefficient2 = TmpParticles->A2A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++Memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		  Coefficient2 = TmpParticles->A3A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++Memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		}
	    }
	  
	  for (int j = 0; j < this->NbrInterSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = TmpParticles->A1A2(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++Memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		  Coefficient2 = TmpParticles->A1A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++Memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		  Coefficient2 = TmpParticles->A2A3(i, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++Memory;
			      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			    }
			}
		    }
		}
	      if (this->OneBodyInteractionFactors11 != 0)
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
		}
	    }	
	}
    }
  else
    {
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	    {
	      Coefficient = TmpParticles->A1A1(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->Ad1Ad1(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	      Coefficient = TmpParticles->A2A2(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->Ad2Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	      Coefficient = TmpParticles->A3A3(i, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->Ad3Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	    }
	  for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	    {
	      Coefficient = TmpParticles->A1A2(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->Ad1Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	      Coefficient = TmpParticles->A1A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->Ad1Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	      Coefficient = TmpParticles->A2A3(i, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->Ad2Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	    }
	  if ((this->OneBodyInteractionFactors11 != 0) || (this->OneBodyInteractionFactors22 != 0) || (this->OneBodyInteractionFactors33 != 0))
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	    }
	}
    }
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnSphereWithSU3SpinHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int Index;
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
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
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double* TmpInteractionFactor;
  int* TmpIndices;
  ParticleOnSphereWithSU3Spin* TmpParticles = (ParticleOnSphereWithSU3Spin*) this->Particles;

  int TotalPos = 0;

  if ((this->NbrIntraSectorSums != 0) || (this->NbrInterSectorSums != 0))
    {
      for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	{
	  this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
	  this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
	  TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
	  TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
	  Pos = 0;
	  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
	      TmpIndices = this->IntraSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = TmpParticles->A1A1(i + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors11[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad1Ad1(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			      ++Pos;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient2 = TmpParticles->A2A2(i + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors22[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad2Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			      ++Pos;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient2 = TmpParticles->A3A3(i + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors33[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad3Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			      ++Pos;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	  for (int j = 0; j < this->NbrInterSectorSums; ++j)
	    {
	      int Lim = 2 * this->NbrInterSectorIndicesPerSum[j];
	      TmpIndices = this->InterSectorIndicesPerSum[j];
	      for (int i1 = 0; i1 < Lim; i1 += 2)
		{
		  Coefficient2 = TmpParticles->A1A2(i + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors12[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad1Ad2(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			      ++Pos;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient2 = TmpParticles->A1A3(i + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors13[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad1Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			      ++Pos;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		  Coefficient2 = TmpParticles->A2A3(i + this->PrecalculationShift, TmpIndices[i1], TmpIndices[i1 + 1]);
		  if (Coefficient2 != 0.0)
		    {
		      TmpInteractionFactor = &(this->InteractionFactors23[j][(i1 * Lim) >> 2]);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Index = TmpParticles->Ad2Ad3(TmpIndices[i2], TmpIndices[i2 + 1], Coefficient);
			  if (Index < Dim)
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * (*TmpInteractionFactor);
			      ++Pos;
			    }
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	  if (this->OneBodyInteractionFactors11 != 0)
	    {
	      double TmpDiagonal = 0.0;
	      if (this->OneBodyInteractionFactors11 != 0)
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors11[j] * TmpParticles->Ad1A1(i + this->PrecalculationShift, j);
	      if (this->OneBodyInteractionFactors22 != 0)
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors22[j] * TmpParticles->Ad2A2(i + this->PrecalculationShift, j);
	      if (this->OneBodyInteractionFactors22 != 0)
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors33[j] * TmpParticles->Ad3A3(i + this->PrecalculationShift, j);
	      TmpIndexArray[Pos] = i + this->PrecalculationShift;
	      TmpCoefficientArray[Pos] = TmpDiagonal;
	      ++Pos;	  
	    }
	  ++TotalPos;
	}
    }
  else
    {
      double Coefficient2;
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      int ReducedNbrInteractionFactors11;
      int ReducedNbrInteractionFactors22;
      int ReducedNbrInteractionFactors33;
      int ReducedNbrInteractionFactors12;
      int ReducedNbrInteractionFactors13;
      int ReducedNbrInteractionFactors23;
      for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	{
	  this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
	  this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
	  TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
	  TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
	  Pos = 0;
	  ReducedNbrInteractionFactors11 = 0;
	  ReducedNbrInteractionFactors22 = 0;
	  ReducedNbrInteractionFactors33 = 0;
	  for (int m1 = 0; m1 < this->NbrM12IntraIndices; ++m1)
	    {
	      Coefficient = TmpParticles->A1A1(i + this->PrecalculationShift, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = TmpParticles->Ad1Ad1(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->M12InteractionFactors11[ReducedNbrInteractionFactors11];
			  ++Pos;
			}		      
		      ++ReducedNbrInteractionFactors11;
		    }    
		}
	      else
		ReducedNbrInteractionFactors11 += this->NbrM3IntraValues[m1];
	      Coefficient = TmpParticles->A2A2(i + this->PrecalculationShift, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = TmpParticles->Ad2Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->M12InteractionFactors22[ReducedNbrInteractionFactors22];
			  ++Pos;
			}		      
		      ++ReducedNbrInteractionFactors22;
		    }    
		}
	      else
		ReducedNbrInteractionFactors22 += this->NbrM3IntraValues[m1];
	      Coefficient = TmpParticles->A3A3(i + this->PrecalculationShift, this->M1IntraValue[m1], this->M2IntraValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1IntraValue[m1] + this->M2IntraValue[m1];
		  TmpM3Values = this->M3IntraValues[m1];
		  TmpNbrM3Values = this->NbrM3IntraValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = TmpParticles->Ad3Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->M12InteractionFactors33[ReducedNbrInteractionFactors33];
			  ++Pos;
			}		      
		      ++ReducedNbrInteractionFactors33;
		    }    
		}
	      else
		ReducedNbrInteractionFactors33 += this->NbrM3IntraValues[m1];
	    }
	  ReducedNbrInteractionFactors12 = 0;
	  for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	    {
	      Coefficient = TmpParticles->A1A2(i + this->PrecalculationShift, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = TmpParticles->Ad1Ad2(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->M12InteractionFactors12[ReducedNbrInteractionFactors12];
			  ++Pos;
			}		      
		      ++ReducedNbrInteractionFactors12;
		    }    
		}
	      else
		ReducedNbrInteractionFactors12 += this->NbrM3InterValues[m1];
	    }
	  ReducedNbrInteractionFactors13 = 0;
	  for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	    {
	      Coefficient = TmpParticles->A1A3(i + this->PrecalculationShift, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = TmpParticles->Ad1Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->M12InteractionFactors13[ReducedNbrInteractionFactors13];
			  ++Pos;
			}		      
		      ++ReducedNbrInteractionFactors13;
		    }    
		}
	      else
		ReducedNbrInteractionFactors13 += this->NbrM3InterValues[m1];
	    }
	  ReducedNbrInteractionFactors23 = 0;
	  for (int m1 = 0; m1 < this->NbrM12InterIndices; ++m1)
	    {
	      Coefficient = TmpParticles->A2A3(i + this->PrecalculationShift, this->M1InterValue[m1], this->M2InterValue[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1InterValue[m1] + this->M2InterValue[m1];
		  TmpM3Values = this->M3InterValues[m1];
		  TmpNbrM3Values = this->NbrM3InterValues[m1];
		  for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      Index = TmpParticles->Ad2Ad3(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * Coefficient2 * this->M12InteractionFactors23[ReducedNbrInteractionFactors23];
			  ++Pos;
			}		      
		      ++ReducedNbrInteractionFactors23;
		    }    
		}
	      else
		ReducedNbrInteractionFactors23 += this->NbrM3InterValues[m1];
	    }
	  if ((this->OneBodyInteractionFactors11 != 0) || (this->OneBodyInteractionFactors22 != 0) || (this->OneBodyInteractionFactors33 != 0))
	    {
	      double TmpDiagonal = 0.0;
	      if (this->OneBodyInteractionFactors11 != 0)
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors11[j] * TmpParticles->Ad1A1(i + this->PrecalculationShift, j);
	      if (this->OneBodyInteractionFactors22 != 0)
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors22[j] * TmpParticles->Ad2A2(i + this->PrecalculationShift, j);
	      if (this->OneBodyInteractionFactors33 != 0)
		for (int j = 0; j <= this->LzMax; ++j) 
		  TmpDiagonal += this->OneBodyInteractionFactors33[j] * TmpParticles->Ad3A3(i + this->PrecalculationShift, j);
	      TmpIndexArray[Pos] = i + this->PrecalculationShift;
	      TmpCoefficientArray[Pos] = TmpDiagonal;
	      ++Pos;	  
	    }
	  ++TotalPos;
	}
    }
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
// lastComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnSphereWithSU3SpinHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnSphereWithSU3SpinHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
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
  double Coefficient2;
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
  int* MIndices;
  int* NIndices;
  this->MaxNbrInteractionPerComponent = 0;

  int TotalPos = 0;
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
// 	  for (int k = 2; k <= this->MaxNBody; ++k)
// 	    if (this->NBodyFlags[k] == true)
// 	      {
// 		double Sign = this->NBodySign[k];
// 		int TmpMinSumIndices = this->MinSumIndices[k];
// 		int TmpMaxSumIndices = this->MaxSumIndices[k];	      
// 		for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
// 		  {
// 		    TmpInteraction = this->NBodyInteractionFactors[k][j];
// 		    int Lim = this->NbrSortedIndicesPerSum[k][j];
// 		    NIndices = this->SortedIndicesPerSum[k][j];
// 		    for (int i1 = 0; i1 < Lim; ++i1)
// 		      {
// 			Coefficient2 = Sign * this->Particles->ProdA(i, NIndices, k);
// 			if (Coefficient2 != 0.0)
// 			  {
// 			    MIndices = this->SortedIndicesPerSum[k][j];
// 			    for (int i2 = 0; i2 < Lim; ++i2)
// 			      {
// 				Index = this->Particles->ProdAd(MIndices, k, Coefficient);
// 				if (Index < this->Particles->GetHilbertSpaceDimension())
// 				  {
// 				    TmpIndexArray[Pos] = Index;
// 				    TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
// 				    ++Pos;
// 				  }
// 				MIndices += k;
// 			      }
// 			  }
// 			NIndices += k;
// 		      }
// 		  }
// 	      }
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

