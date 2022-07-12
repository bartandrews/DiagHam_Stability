////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 29/06/2010                      //
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


#include "Hamiltonian/ParticleOnCylinderPermanentState.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderPermanentState::ParticleOnCylinderPermanentState(ParticleOnSphereWithSpin* particles, int nbrParticles, int maxMomentum,
										   double ratio,AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;

  if (precalculationFileName == 0)
    {
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
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnCylinderPermanentState::~ParticleOnCylinderPermanentState() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;
  delete[] this->M5Value;
  delete[] this->M6Value;


  if (this->FastMultiplicationFlag == true)
    {
      int ReducedDim = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      this->FastMultiplicationFlag = false;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnCylinderPermanentState::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnCylinderPermanentState::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderPermanentState::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 < m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, m5, m6);
		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->M5Value = new int [Pos];
      this->M6Value = new int [Pos];

      this->InteractionFactors = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
            for (int m6 = 0; m6 <= this->MaxMomentum; ++m6)
	      for (int m5 = 0; m5 < m6; ++m5)
 	       {
	         m4 = m1 + m2 + m3 - m6 - m5;
	         if ((m4 >= 0) && (m4 <= this->MaxMomentum))
		     {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		           this->M1Value[this->NbrInteractionFactors] = m1;
		           this->M2Value[this->NbrInteractionFactors] = m2;
		           this->M3Value[this->NbrInteractionFactors] = m3;
		           this->M4Value[this->NbrInteractionFactors] = m4;
		           this->M5Value[this->NbrInteractionFactors] = m5;
		           this->M6Value[this->NbrInteractionFactors] = m6;
                           cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactors;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient;
    }
  else //bosons
    {
    }
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// m5 = fifth index
// m6 = sixth index
// return value = numerical coefficient

Complex ParticleOnCylinderPermanentState::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int m5, int m6)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m6 - m4 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double GaussianExp;

  Complex Coefficient(0,0);

     GaussianExp = Xr * Xr + Xs * Xs + Xr * Xs;
     Coefficient.Re = exp(-GaussianExp) * (Xr - Xs);
     Coefficient.Im = 0.0;

     GaussianExp = Xrp * Xrp + Xsp * Xsp + Xrp * Xsp;
     Coefficient.Re *= (exp(-GaussianExp) * (Xrp - Xsp));
     Coefficient.Im = 0.0;
     return (-Coefficient * 16.0 * sqrt(M_PI) * sqrt(3.0 * M_PI)/(2.0 * M_PI * this->Ratio * this->NbrLzValue));
 
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnCylinderPermanentState::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Shift = this->EnergyShift;
  double Coefficient, Coefficient2;

  if (this->FastMultiplicationFlag == false)
    {
      //cout<<"Flagfalse"<<endl;

      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      int m5;
      int m6;
      int TmpA[3];
      int TmpC[3];
      int TmpSpinA[3];
      int TmpSpinC[3];

      Complex TmpInteraction;
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();

      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{

          //UPUPDOWN 

	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 0;
	  TmpSpinA[2] = 0;
          
	  TmpInteraction = this->InteractionFactors[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }


          //DOWNDOWNUP 

	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
	  TmpSpinA[2] = 1;
          
	  TmpInteraction = this->InteractionFactors[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<"DDU "<<TmpC[0]<<" "<<TmpC[1]<<" "<<TmpC[2]<<" "<<TmpA[0]<<" "<<TmpA[1]<<" "<<TmpA[2]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                } 
	    }

	} //j

      for (int i = firstComponent; i < LastComponent; ++i)
	{
	 vDestination[i] += this->EnergyShift * vSource[i];
	}

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
          ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
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
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
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

	  int Index;
	  int m1;
	  int m2;
	  int m3;
	  int m4;
          int m5;
          int m6;
          int TmpA[3];
          int TmpC[3];
          int TmpSpinA[3];
          int TmpSpinC[3];

	  Complex TmpInteraction;

	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		  {
                    //UPUPDOWN 

	            TmpC[0] = this->M1Value[j];
	            TmpC[1] = this->M2Value[j];
	            TmpC[2] = this->M3Value[j];
	            TmpA[0] = this->M4Value[j];
	            TmpA[1] = this->M5Value[j];
	            TmpA[2] = this->M6Value[j];

	            TmpSpinC[0] = 0;
	            TmpSpinC[1] = 0;
	            TmpSpinC[2] = 1;
	            TmpSpinA[0] = 1;
	            TmpSpinA[1] = 0;
             	    TmpSpinA[2] = 0;

		    TmpInteraction = this->InteractionFactors[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	                   if (Index <= i)
                             {
		                vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                                if (Index < i)
  		                  vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];

                             }
                         } 
		      }

                    //UPUPDOWN 

	            TmpC[0] = this->M1Value[j];
	            TmpC[1] = this->M2Value[j];
	            TmpC[2] = this->M3Value[j];
	            TmpA[0] = this->M4Value[j];
	            TmpA[1] = this->M5Value[j];
	            TmpA[2] = this->M6Value[j];

	            TmpSpinC[0] = 1;
	            TmpSpinC[1] = 1;
	            TmpSpinC[2] = 0;
	            TmpSpinA[0] = 0;
	            TmpSpinA[1] = 1;
             	    TmpSpinA[2] = 1;

		    TmpInteraction = this->InteractionFactors[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	                   if (Index <= i)
                             {
		                vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                                if (Index < i)
  		                  vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];

                             }
                         } 
		      }


		  }

		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    vDestination[i] += this->EnergyShift * vSource[i];
		  }
                
	      }
	  delete TmpParticles;
	}
   }


  return vDestination;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnCylinderPermanentState::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
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
  int TmpSpinA[3];
  int TmpSpinC[3];

  int LastComponent = lastComponent + firstComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();

  cout<<"Partial fast memory permanent "<<endl;

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

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 0;
          TmpSpinA[2] = 0;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
          TmpSpinA[2] = 1;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
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

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void ParticleOnCylinderPermanentState::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  int NbrTranslation;
  int* TmpIndexArray;
  Complex* TmpCoefficientArray;
  int* TmpNbrTranslationArray;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();  
  int LastComponent = firstComponent + lastComponent;
  int m1, m2, m3, m4, m5, m6;
  int TmpC[3];
  int TmpA[3]; 
  int TmpSpinC[3];
  int TmpSpinA[3]; 
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

  cout<<"Partial enable fast permanent "<<endl;

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

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinC[2] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 0;
          TmpSpinA[2] = 0;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactors[j];
	       ++Pos;
              }
           }


	  TmpC[0] = this->M1Value[j];
	  TmpC[1] = this->M2Value[j];
	  TmpC[2] = this->M3Value[j];
	  TmpA[0] = this->M4Value[j];
	  TmpA[1] = this->M5Value[j];
	  TmpA[2] = this->M6Value[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinC[2] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
          TmpSpinA[2] = 1;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 3);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 3, Coefficient);
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
