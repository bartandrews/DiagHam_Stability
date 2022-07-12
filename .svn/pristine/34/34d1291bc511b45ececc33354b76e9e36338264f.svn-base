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


#include "Hamiltonian/ParticleOnCylinderHaldaneRezayi.h"
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

ParticleOnCylinderHaldaneRezayi::ParticleOnCylinderHaldaneRezayi(ParticleOnSphereWithSpin* particles, int nbrParticles, int maxMomentum,
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
  this->HermitianSymmetryFlag=true;

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

ParticleOnCylinderHaldaneRezayi::~ParticleOnCylinderHaldaneRezayi() 
{
  delete[] this->InteractionFactorsIntra;
  delete[] this->M1ValueIntra;
  delete[] this->M2ValueIntra;
  delete[] this->M3ValueIntra;
  delete[] this->M4ValueIntra;

  delete[] this->InteractionFactorsInter;
  delete[] this->M1ValueInter;
  delete[] this->M2ValueInter;
  delete[] this->M3ValueInter;
  delete[] this->M4ValueInter;


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

void ParticleOnCylinderHaldaneRezayi::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactorsIntra;
  delete[] this->InteractionFactorsInter;

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

void ParticleOnCylinderHaldaneRezayi::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderHaldaneRezayi::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m3;
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];


      //********************** INTRA ************************************

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
 	       {
	         m3 = m1 + m2 - m4;
	         if ((m3 >= 0) && (m3 <= this->MaxMomentum))
                   if (m3 > m4)
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientIntra(m1, m2, m3, m4);
                                             //-this->EvaluateInteractionCoefficientIntra(m2, m1, m3, m4)
                                             //-this->EvaluateInteractionCoefficientIntra(m1, m2, m4, m3)
                                             //+this->EvaluateInteractionCoefficientIntra(m2, m1, m4, m3);

		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactorsIntra = 0;
      this->M1ValueIntra = new int [Pos];
      this->M2ValueIntra = new int [Pos];
      this->M3ValueIntra = new int [Pos];
      this->M4ValueIntra = new int [Pos];

      this->InteractionFactorsIntra = new Complex [Pos];
      cout << "nbr interaction intra = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

       for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
 	       {
	         m3 = m1 + m2 - m4;
	         if ((m3 >= 0) && (m3 <= this->MaxMomentum))
                  if (m3 > m4)
  	             {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactorsIntra[this->NbrInteractionFactorsIntra] = TmpCoefficient[Pos];
		           this->M1ValueIntra[this->NbrInteractionFactorsIntra] = m1;
		           this->M2ValueIntra[this->NbrInteractionFactorsIntra] = m2;
		           this->M3ValueIntra[this->NbrInteractionFactorsIntra] = m3;
		           this->M4ValueIntra[this->NbrInteractionFactorsIntra] = m4;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactorsIntra;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction intra = " << this->NbrInteractionFactorsIntra << endl;

      //********************** INTER ************************************

      Pos = 0;
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= this->MaxMomentum; ++m2)
	  for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
 	       {
	         m3 = m1 + m2 - m4;
	         if ((m3 >= 0) && (m3 <= this->MaxMomentum))
  	             {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientInter(m1, m2, m3, m4);
                                             //+this->EvaluateInteractionCoefficientInter(m2, m1, m4, m3);

		       if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		         MaxCoefficient = Norm(TmpCoefficient[Pos]);
		       ++Pos;
		     }
	        }

      this->NbrInteractionFactorsInter = 0;
      this->M1ValueInter = new int [Pos];
      this->M2ValueInter = new int [Pos];
      this->M3ValueInter = new int [Pos];
      this->M4ValueInter = new int [Pos];

      this->InteractionFactorsInter = new Complex [Pos];
      cout << "nbr interaction inter = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

       for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= this->MaxMomentum; ++m2)
	  for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
 	       {
	         m3 = m1 + m2 - m4;
	         if ((m3 >= 0) && (m3 <= this->MaxMomentum))
  	             {
		       if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		         {
		           this->InteractionFactorsInter[this->NbrInteractionFactorsInter] = TmpCoefficient[Pos];
		           this->M1ValueInter[this->NbrInteractionFactorsInter] = m1;
		           this->M2ValueInter[this->NbrInteractionFactorsInter] = m2;
		           this->M3ValueInter[this->NbrInteractionFactorsInter] = m3;
		           this->M4ValueInter[this->NbrInteractionFactorsInter] = m4;
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactorsInter;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction inter = " << this->NbrInteractionFactorsInter << endl;

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
// return value = numerical coefficient

Complex ParticleOnCylinderHaldaneRezayi::EvaluateInteractionCoefficientIntra(int m1, int m2, int m3, int m4)
{
  double Length = sqrt(2.0 * M_PI * this->NbrLzValue * this->Ratio);
  double kappa = 2.0 * M_PI/Length;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;	
  double Xr = 0.5 * (Xm1 - Xm2);
  double Xs = 0.5 * (Xm4 - Xm3);

  Complex Coefficient(0,0);

  //int rtimes2 = m1 - m2;
  //int rptimes2 = m4 -m3;
  //if ((rtimes2 * rtimes2 + rptimes2 * rptimes2) <= 18)
  //  {
      Coefficient.Re = exp(-Xr*Xr - Xs*Xs) * 4.0 * Xr * Xs;
      Coefficient.Im = 0.0;
  //  }
  return (Coefficient/sqrt(this->Ratio * this->NbrLzValue));
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnCylinderHaldaneRezayi::EvaluateInteractionCoefficientInter(int m1, int m2, int m3, int m4)
{
  double Length = sqrt(2.0 * M_PI * this->NbrLzValue * this->Ratio);
  double kappa = 2.0 * M_PI/Length;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;	
  double Xr = 0.5 * (Xm1 - Xm2);
  double Xs = 0.5 * (Xm4 - Xm3);

  Complex Coefficient(0,0);

  Coefficient.Re = exp(-Xr*Xr - Xs*Xs) * 8.0 * Xr * Xs;
  Coefficient.Im = 0.0;
  return (Coefficient/sqrt(this->Ratio * this->NbrLzValue));
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnCylinderHaldaneRezayi::HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
      int TmpA[2];
      int TmpC[2];
      int TmpSpinA[2];
      int TmpSpinC[2];

      Complex TmpInteraction;
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();

      for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
	{

          //Intra

	  TmpC[0] = this->M1ValueIntra[j];
	  TmpC[1] = this->M2ValueIntra[j];
	  TmpA[0] = this->M3ValueIntra[j];
	  TmpA[1] = this->M4ValueIntra[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 0;
          
	  TmpInteraction = this->InteractionFactorsIntra[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
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


	  TmpC[0] = this->M1ValueIntra[j];
	  TmpC[1] = this->M2ValueIntra[j];
	  TmpA[0] = this->M3ValueIntra[j];
	  TmpA[1] = this->M4ValueIntra[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 1;
          
	  TmpInteraction = this->InteractionFactorsIntra[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<" "<<TmpA[0]<<" "<<TmpA[1]<<endl;
                      if (Index < i)
                        vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];
                    }
                }
	    }

	} //j


      for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
	{

          //Inter

	  TmpC[0] = this->M1ValueInter[j];
	  TmpC[1] = this->M2ValueInter[j];
	  TmpA[0] = this->M3ValueInter[j];
	  TmpA[1] = this->M4ValueInter[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
          
	  TmpInteraction = this->InteractionFactorsInter[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
              if (Coefficient2 != 0.0)
                {
                   Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	           if (Index <= i)
                    {
		      vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                      //TmpParticles->PrintState(cout,i);cout<<" ";TmpParticles->PrintState(cout,Index);
                      //cout<<" "<<TmpC[0]<<" "<<TmpC[1]<<" "<<" "<<TmpA[0]<<" "<<TmpA[1]<<endl;
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
          int TmpA[2];
          int TmpC[2];
          int TmpSpinA[2];
          int TmpSpinC[2];

	  Complex TmpInteraction;

	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
		  {
                    //Intra

	            TmpC[0] = this->M1ValueIntra[j];
	            TmpC[1] = this->M2ValueIntra[j];
	            TmpA[1] = this->M3ValueIntra[j];
	            TmpA[0] = this->M4ValueIntra[j];

	            TmpSpinC[0] = 0;
	            TmpSpinC[1] = 0;
	            TmpSpinA[0] = 0;
	            TmpSpinA[1] = 0;

		    TmpInteraction = this->InteractionFactorsIntra[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	                   if (Index <= i)
                             {
		                vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                                if (Index < i)
  		                  vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];

                             }
                         } 
		      }


	            TmpC[0] = this->M1ValueIntra[j];
	            TmpC[1] = this->M2ValueIntra[j];
	            TmpA[1] = this->M3ValueIntra[j];
	            TmpA[0] = this->M4ValueIntra[j];

	            TmpSpinC[0] = 1;
	            TmpSpinC[1] = 1;
	            TmpSpinA[0] = 1;
	            TmpSpinA[1] = 1;

		    TmpInteraction = this->InteractionFactorsIntra[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	                   if (Index <= i)
                             {
		                vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                                if (Index < i)
  		                  vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];

                             }
                         } 
		      }


		  } //j


		for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
		  {
                    //Inter

	            TmpC[0] = this->M1ValueInter[j];
	            TmpC[1] = this->M2ValueInter[j];
	            TmpA[1] = this->M3ValueInter[j];
	            TmpA[0] = this->M4ValueInter[j];

	            TmpSpinC[0] = 1;
	            TmpSpinC[1] = 0;
	            TmpSpinA[0] = 0;
	            TmpSpinA[1] = 1;

		    TmpInteraction = this->InteractionFactorsInter[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
	                Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
                        if (Coefficient2 != 0.0)
                         {
                           Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	                   if (Index <= i)
                             {
		                vDestination[Index] += Coefficient * Coefficient2 * TmpInteraction * vSource[i];
                                if (Index < i)
  		                  vDestination[i] += Coefficient * Coefficient2 * Conj(TmpInteraction) * vSource[Index];

                             }
                         } 
		      }

		  } //j

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

long ParticleOnCylinderHaldaneRezayi::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient, Coefficient2;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  int TmpA[2];
  int TmpC[2];
  int TmpSpinA[2];
  int TmpSpinC[2];

  int LastComponent = lastComponent + firstComponent;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int Dim = TmpParticles->GetHilbertSpaceDimension();

  cout<<"Partial fast memory HR "<<endl;

  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
	{

	  TmpC[0] = this->M1ValueIntra[j];
	  TmpC[1] = this->M2ValueIntra[j];
	  TmpA[0] = this->M3ValueIntra[j];
	  TmpA[1] = this->M4ValueIntra[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 0;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	  TmpC[0] = this->M1ValueIntra[j];
	  TmpC[1] = this->M2ValueIntra[j];
	  TmpA[0] = this->M3ValueIntra[j];
	  TmpA[1] = this->M4ValueIntra[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 1;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	} //j    

      for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
	{

	  TmpC[0] = this->M1ValueInter[j];
	  TmpC[1] = this->M2ValueInter[j];
	  TmpA[0] = this->M3ValueInter[j];
	  TmpA[1] = this->M4ValueInter[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;
          
          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
  	    if (Index <= i)
             {
               ++Memory;
	       ++this->NbrInteractionPerComponent[i];
             }
           } 

	} //j    

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

void ParticleOnCylinderHaldaneRezayi::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
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
  int m1, m2, m3, m4;
  int TmpC[2];
  int TmpA[2]; 
  int TmpSpinC[2];
  int TmpSpinA[2]; 
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

      for (int j = 0; j < this->NbrInteractionFactorsIntra; ++j) 
	{
	  TmpC[0] = this->M1ValueIntra[j];
	  TmpC[1] = this->M2ValueIntra[j];
	  TmpA[0] = this->M3ValueIntra[j];
	  TmpA[1] = this->M4ValueIntra[j];

	  TmpSpinC[0] = 0;
	  TmpSpinC[1] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 0;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactorsIntra[j];
	       ++Pos;
              }
           }


	  TmpC[0] = this->M1ValueIntra[j];
	  TmpC[1] = this->M2ValueIntra[j];
	  TmpA[0] = this->M3ValueIntra[j];
	  TmpA[1] = this->M4ValueIntra[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 1;
	  TmpSpinA[0] = 1;
	  TmpSpinA[1] = 1;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactorsIntra[j];
	       ++Pos;
              }
           }

 
	}  //j


      for (int j = 0; j < this->NbrInteractionFactorsInter; ++j) 
	{
	  TmpC[0] = this->M1ValueInter[j];
	  TmpC[1] = this->M2ValueInter[j];
	  TmpA[0] = this->M3ValueInter[j];
	  TmpA[1] = this->M4ValueInter[j];

	  TmpSpinC[0] = 1;
	  TmpSpinC[1] = 0;
	  TmpSpinA[0] = 0;
	  TmpSpinA[1] = 1;

          Coefficient2 = TmpParticles->ProdA(i, TmpA, TmpSpinA, 2);
          if (Coefficient2 != 0.0)
           {
            Index = TmpParticles->ProdAd(TmpC, TmpSpinC, 2, Coefficient);
            if (Index <= i)
              {
	       TmpIndexArray[Pos] = Index;
	       TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * this->InteractionFactorsInter[j];
	       ++Pos;
              }
           }
 
	}  //j

     ++PosIndex;
   }

  delete TmpParticles;
}
