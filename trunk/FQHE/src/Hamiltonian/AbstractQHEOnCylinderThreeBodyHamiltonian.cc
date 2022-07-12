////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                         to particles on a torus with                       //
//                                                                            //
//                        last modification : 27/06/2003                      //
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
#include "Hamiltonian/AbstractQHEOnCylinderThreeBodyHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
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
using std::endl;
using std::ostream;

// default constructor
//
AbstractQHEOnCylinderThreeBodyHamiltonian::AbstractQHEOnCylinderThreeBodyHamiltonian()
{
  this->HermitianSymmetryFlag=true;
}


// destructor
//

AbstractQHEOnCylinderThreeBodyHamiltonian::~AbstractQHEOnCylinderThreeBodyHamiltonian()
{
}

// ask if Hamiltonian implements hermitian symmetry operations
//
bool AbstractQHEOnCylinderThreeBodyHamiltonian::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnCylinderThreeBodyHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnCylinderThreeBodyHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnCylinderThreeBodyHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnCylinderThreeBodyHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnCylinderThreeBodyHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnCylinderThreeBodyHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnCylinderThreeBodyHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
RealVector& AbstractQHEOnCylinderThreeBodyHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  return vDestination;
}
 

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnCylinderThreeBodyHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->HermitianLowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// WORKING VERSION -- FULL
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
/*
ComplexVector& AbstractQHEOnCylinderThreeBodyHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Shift = this->EnergyShift;
  double Coefficient, Coefficient2;
  if (this->FastMultiplicationFlag == false)
    {
      //cout<<"Flagfalse"<<endl;
      double Coefficient;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      int m5;
      int m6;
      int TmpA[3];
      int TmpC[3];
      Complex TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	{
	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];
          
	  TmpInteraction = this->InteractionFactors[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	           if (Index < Dim)
		     vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                } 
	    }
	}
      TmpC[0] = this->M1Value[ReducedNbrInteractionFactors];
      TmpC[1] = this->M2Value[ReducedNbrInteractionFactors];
      TmpC[2] = this->M3Value[ReducedNbrInteractionFactors];
      TmpA[0] = this->M4Value[ReducedNbrInteractionFactors];
      TmpA[1] = this->M5Value[ReducedNbrInteractionFactors];
      TmpA[2] = this->M6Value[ReducedNbrInteractionFactors];

      TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
          Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	    if (Index < Dim)
	      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
           } 
	 vDestination[i] += this->EnergyShift * vSource[i];
	}

      if (this->OneBodyInteractionFactors != 0)
        for (int i = firstComponent; i < LastComponent; ++i)
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
		vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];

      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
          ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * vSource[i];
	      vDestination[i] += Shift * vSource[i];
	    }
        
          if (this->OneBodyInteractionFactors != 0)
           for (int i = firstComponent; i < LastComponent; ++i)
	     for (int j = 0; j <= this->MaxMomentum; ++j) 
		vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];
          delete TmpParticles;
	}
      else
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * vSource[i];
	      vDestination[i] += Shift * vSource[i];
	      ++Pos;
	    }
           if (this->OneBodyInteractionFactors != 0)
            for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	     for (int j = 0; j <= this->MaxMomentum; ++j) 
	       vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];

	  int Index;
	  int m1;
	  int m2;
	  int m3;
	  int m4;
          int m5;
          int m6;
          int TmpA[3];
          int TmpC[3];
	  Complex TmpInteraction;
	  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;

	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
		  {
	            TmpC[0] = this->M1Value[j];
	            TmpC[1] = this->M2Value[j];
	            TmpC[2] = this->M3Value[j];
	            TmpA[0] = this->M4Value[j];
	            TmpA[1] = this->M5Value[j];
	            TmpA[2] = this->M6Value[j];

		    TmpInteraction = this->InteractionFactors[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	                   if (Index < Dim)
		             vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                         } 
		      }
		  }
                TmpC[0] = this->M1Value[ReducedNbrInteractionFactors];
                TmpC[1] = this->M2Value[ReducedNbrInteractionFactors];
                TmpC[2] = this->M3Value[ReducedNbrInteractionFactors];
                TmpA[0] = this->M4Value[ReducedNbrInteractionFactors];
                TmpA[1] = this->M5Value[ReducedNbrInteractionFactors];
                TmpA[2] = this->M6Value[ReducedNbrInteractionFactors];

		TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
	            Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
                    if (Coefficient2 != 0.0)
                     {
                      Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	              if (Index < Dim)
		         vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      } 
		    vDestination[i] += this->EnergyShift * vSource[i];
		  }
                if (this->OneBodyInteractionFactors != 0)
                 for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	          for (int j = 0; j <= this->MaxMomentum; ++j) 
		     vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];
	      }
	  delete TmpParticles;
	}
   }
  return vDestination;
}
*/

ComplexVector& AbstractQHEOnCylinderThreeBodyHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Shift = this->EnergyShift;
  double Coefficient, Coefficient2;
  if (this->FastMultiplicationFlag == false)
    {
      //cout<<"Flagfalse"<<endl;
      double Coefficient;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      int m5;
      int m6;
      int TmpA[3];
      int TmpC[3];
      Complex TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	{
	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];
          
	  TmpInteraction = this->InteractionFactors[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }
	}
      TmpC[0] = this->M1Value[ReducedNbrInteractionFactors];
      TmpC[1] = this->M2Value[ReducedNbrInteractionFactors];
      TmpC[2] = this->M3Value[ReducedNbrInteractionFactors];
      TmpA[0] = this->M4Value[ReducedNbrInteractionFactors];
      TmpA[1] = this->M5Value[ReducedNbrInteractionFactors];
      TmpA[2] = this->M6Value[ReducedNbrInteractionFactors];

      TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
          Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	    if (Index <= i)
              {
	        vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                if (Index < i)
   	          vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
              }
           } 
	 vDestination[i] += this->EnergyShift * vSource[i];
	}

      if (this->OneBodyInteractionFactors != 0)
        for (int i = firstComponent; i < LastComponent; ++i)
	    for (int j = 0; j <= this->MaxMomentum; ++j) 
		vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];

      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
          ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (j = 0; j < TmpNbrInteraction; ++j)
               { 
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * vSource[i];
                  if (TmpIndexArray[j] < i)
                     vDestination[i] += vSource[TmpIndexArray[j]] * Conj(TmpCoefficientArray[j]);
               }
	      vDestination[i] += Shift * vSource[i];
	    }
        
          if (this->OneBodyInteractionFactors != 0)
           for (int i = firstComponent; i < LastComponent; ++i)
	     for (int j = 0; j <= this->MaxMomentum; ++j) 
		vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];
          delete TmpParticles;
	}
      else
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      for (j = 0; j < TmpNbrInteraction; ++j)
               {
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * vSource[i];
                  if (TmpIndexArray[j] < i)
                     vDestination[i] += vSource[TmpIndexArray[j]] * Conj(TmpCoefficientArray[j]);
               }
	      vDestination[i] += Shift * vSource[i];
	      ++Pos;
	    }
           if (this->OneBodyInteractionFactors != 0)
            for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	     for (int j = 0; j <= this->MaxMomentum; ++j) 
	       vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];

	  int Index;
	  int m1;
	  int m2;
	  int m3;
	  int m4;
          int m5;
          int m6;
          int TmpA[3];
          int TmpC[3];
	  Complex TmpInteraction;
	  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;

	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
		  {
	            TmpC[0] = this->M1Value[j];
	            TmpC[1] = this->M2Value[j];
	            TmpC[2] = this->M3Value[j];
	            TmpA[0] = this->M4Value[j];
	            TmpA[1] = this->M5Value[j];
	            TmpA[2] = this->M6Value[j];

		    TmpInteraction = this->InteractionFactors[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	                   if (Index <= i)
                             {
		                vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                                if (Index < i)
  		                  vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];

                             }
                         } 
		      }
		  }
                TmpC[0] = this->M1Value[ReducedNbrInteractionFactors];
                TmpC[1] = this->M2Value[ReducedNbrInteractionFactors];
                TmpC[2] = this->M3Value[ReducedNbrInteractionFactors];
                TmpA[0] = this->M4Value[ReducedNbrInteractionFactors];
                TmpA[1] = this->M5Value[ReducedNbrInteractionFactors];
                TmpA[2] = this->M6Value[ReducedNbrInteractionFactors];

		TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
	            Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
                    if (Coefficient2 != 0.0)
                     {
                      Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	              if (Index <= i)
                        {
		           vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                           if (Index < i)
                             vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                        }
                      } 
		    vDestination[i] += this->EnergyShift * vSource[i];
		  }
                if (this->OneBodyInteractionFactors != 0)
                 for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	          for (int j = 0; j <= this->MaxMomentum; ++j) 
		     vDestination[i] += this->OneBodyInteractionFactors[j] * TmpParticles->AdA(i, j) * vSource[i];
	      }
	  delete TmpParticles;
	}
   }
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnCylinderThreeBodyHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnCylinderThreeBodyHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEOnCylinderThreeBodyHamiltonian::FastMultiplicationMemory(long allowedMemory)
{

  this->NbrInteractionPerComponent = new int [this->Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension();
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(Complex)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(Complex)));
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * this->Particles->GetHilbertSpaceDimension()) + (Memory * (sizeof (int) + sizeof(Complex)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// WORKING VERSION -- FULL
// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element
/*
long AbstractQHEOnCylinderThreeBodyHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  int m5;
  int m6;
  int TmpA[3];
  int TmpC[3];
  int LastComponent = lastComponent + firstComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();

  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{

	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	    if (Index < Dim)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 
	}    
    }
  Memory = ((sizeof (int*) + sizeof (int) + 2 * sizeof(double*)) * TmpParticles->GetHilbertSpaceDimension() + 
	    Memory *  (sizeof (int) + sizeof(Complex)));
  delete TmpParticles;
  return Memory;
}

*/

long AbstractQHEOnCylinderThreeBodyHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  int m5;
  int m6;
  int TmpA[3];
  int TmpC[3];
  int LastComponent = lastComponent + firstComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();

  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{

	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 
	}    
    }
  Memory = ((sizeof (int*) + sizeof (int) + 2 * sizeof(double*)) * TmpParticles->GetHilbertSpaceDimension() + 
	    Memory *  (sizeof (int) + sizeof(Complex)));
  delete TmpParticles;
  return Memory;
}



// enable fast multiplication algorithm
//
void AbstractQHEOnCylinderThreeBodyHamiltonian::EnableFastMultiplication()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  
  gettimeofday (&(TotalStartingTime2), 0);
   cout<<endl<<"EnableFast..."<<endl;
  cout << "start enable" << endl;
  int ReducedSpaceDimension = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new Complex* [ReducedSpaceDimension];
  
  QHEParticlePrecalculationOperation Operation(this, false);
  //cout<<"After Precalculation"<<endl;
  Operation.ApplyOperation(this->Architecture);

  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout<<"End enableFast..."<<endl;
}


// WORKING VERSION -- FULL
// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
/*
void AbstractQHEOnCylinderThreeBodyHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  int NbrTranslation;
  int* TmpIndexArray;
  Complex* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();  
  int LastComponent = firstComponent + lastComponent;
  int m1, m2, m3, m4, m5, m6;
  int TmpC[3];
  int TmpA[3];      
  int SumIndices;
  int Pos;
  int ReducedNbrInteractionFactors;
  int PosIndex = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++PosIndex;
      PosMod = this->FastMultiplicationStep - PosMod;
    }

  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)  
    {

      this->InteractionPerComponentIndex[PosIndex] = new int [this->NbrInteractionPerComponent[PosIndex]];
      this->InteractionPerComponentCoefficient[PosIndex] = new Complex [this->NbrInteractionPerComponent[PosIndex]];      
      TmpIndexArray = this->InteractionPerComponentIndex[PosIndex];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[PosIndex];
      Pos = 0;

      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];

          Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
            if (Index < Dim)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors[j];
	       ++Pos;
              }
           } 
	}
     ++PosIndex;
   }

  delete TmpParticles;
}
*/

void AbstractQHEOnCylinderThreeBodyHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  int NbrTranslation;
  int* TmpIndexArray;
  Complex* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();  
  int LastComponent = firstComponent + lastComponent;
  int m1, m2, m3, m4, m5, m6;
  int TmpC[3];
  int TmpA[3];      
  int SumIndices;
  int Pos;
  int ReducedNbrInteractionFactors;
  int PosIndex = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++PosIndex;
      PosMod = this->FastMultiplicationStep - PosMod;
    }

  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)  
    {

      this->InteractionPerComponentIndex[PosIndex] = new int [this->NbrInteractionPerComponent[PosIndex]];
      this->InteractionPerComponentCoefficient[PosIndex] = new Complex [this->NbrInteractionPerComponent[PosIndex]];      
      TmpIndexArray = this->InteractionPerComponentIndex[PosIndex];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[PosIndex];
      Pos = 0;

      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];

          Coefficient2 = TmpParticles->ProdA(i, TmpA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, 3, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors[j];
	       ++Pos;
              }
           } 
	}
     ++PosIndex;
   }

  delete TmpParticles;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnCylinderThreeBodyHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Particles->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(Complex) * this->NbrInteractionPerComponent[i]);	  
	}
      File.close();
      return true;
    }
  else
    {
      return false;
    }
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractQHEOnCylinderThreeBodyHamiltonian::LoadPrecalculation (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Particles->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrInteractionPerComponent = new int [Tmp];
  File.read((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficient = new Complex* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new Complex [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(Complex) * this->NbrInteractionPerComponent[i]);	  
    }
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

