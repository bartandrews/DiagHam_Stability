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


#include "Hamiltonian/ParticleOnCylinderOneQuarterPfaffian.h"
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
  // gaffnianFlag = consider spinful Gaffnian instead of permanent times 221
  // nassFlag = consider NASS state at nu=4/3
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderOneQuarterPfaffian::ParticleOnCylinderOneQuarterPfaffian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum, double ratio, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
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
	
  this->EnergyShift = 0;

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

ParticleOnCylinderOneQuarterPfaffian::~ParticleOnCylinderOneQuarterPfaffian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;
  delete[] this->M5Value;
  delete[] this->M6Value;
  delete[] this->InteractionFactors2b;
  delete[] this->M1Value2b;
  delete[] this->M2Value2b;
  delete[] this->M3Value2b;
  delete[] this->M4Value2b;

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

void ParticleOnCylinderOneQuarterPfaffian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnCylinderOneQuarterPfaffian::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderOneQuarterPfaffian::EvaluateInteractionFactors()
{
  int Pos;
  int m4, m6;
  double MaxCoefficient;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {

      //2-body
     
      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient2b = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	        if (m3 > m4)
		  {
		    TmpCoefficient2b[Pos] = (this->EvaluateInteractionCoefficient2b(m1, m2, m3, m4)
			  		 + this->EvaluateInteractionCoefficient2b(m2, m1, m4, m3)
					 - this->EvaluateInteractionCoefficient2b(m1, m2, m4, m3)
					 - this->EvaluateInteractionCoefficient2b(m2, m1, m3, m4));
		    if (MaxCoefficient < Norm(TmpCoefficient2b[Pos]))
		      MaxCoefficient = Norm(TmpCoefficient2b[Pos]);
		    ++Pos;
		  }
	    }
      this->NbrInteractionFactors2b = 0;
      this->M1Value2b = new int [Pos];
      this->M2Value2b = new int [Pos];
      this->M3Value2b = new int [Pos];
      this->M4Value2b = new int [Pos];

      this->InteractionFactors2b = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
              if ((m4 >= 0) && (m4 <= this->MaxMomentum))
	        if (m3 > m4)
		  {
		    if  (Norm(TmpCoefficient2b[Pos]) > MaxCoefficient)
		      {
		        this->InteractionFactors2b[this->NbrInteractionFactors2b] = TmpCoefficient2b[Pos];
		        this->M1Value2b[this->NbrInteractionFactors2b] = m1;
		        this->M2Value2b[this->NbrInteractionFactors2b] = m2;
		        this->M3Value2b[this->NbrInteractionFactors2b] = m3;
		        this->M4Value2b[this->NbrInteractionFactors2b] = m4;
                        //cout<<TmpCoefficient[Pos].Re<<" "<<(m1+1)<<" "<<(m2+1)<<" "<<(m3+1)<<" "<<(m4+1)<<endl;
		        ++this->NbrInteractionFactors2b;
		      }
		    ++Pos;
		  }
	    }

     cout << "nbr interaction 2-body = " << this->NbrInteractionFactors2b << endl;
     cout << "====================================" << endl;

     delete[] TmpCoefficient2b;	

//----------------------------------------------------------------------------
      Pos = 0;    
      MaxCoefficient = 0.0;
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue * this->NbrLzValue)];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < m2; ++m3)
            for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
	      for (int m5 = 0; m5 < m4; ++m5)
 	       {
	         m6 = m1 + m2 + m3 - m4 - m5;
	         if ((m6 >= 0) && (m6 < m5))
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
	  for (int m3 = 0; m3 < m2; ++m3)
            for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
	      for (int m5 = 0; m5 < m4; ++m5)
 	       {
	         m6 = m1 + m2 + m3 - m4 - m5;
	         if ((m6 >= 0) && (m6 < m5))
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
                           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<TmpCoefficient[Pos]<<endl;
		           ++this->NbrInteractionFactors;
		         }
		       ++Pos;
		    }
	       }

     cout << "nbr interaction 3b = " << this->NbrInteractionFactors << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient;

//----------------------------------------------------------------------------------------------------------
    }
  else
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

Complex ParticleOnCylinderOneQuarterPfaffian::EvaluateInteractionCoefficient2b(int m1, int m2, int m3, int m4)
{

  double Length = sqrt(2.0 * M_PI * this->NbrLzValue * this->Ratio);
  double kappa = 2.0 * M_PI/Length;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;	

  Complex Coefficient(0,0);

  Coefficient.Re = exp(-0.5*pow(Xm1-Xm3,2.0)-0.5*pow(Xm1-Xm4,2.0)) * (pow(Xm1-Xm3,2.0)-pow(Xm1-Xm4,2.0)-1.0);
  Coefficient.Im = 0.0;
  return (Coefficient/sqrt(this->Ratio * this->NbrLzValue));
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

Complex ParticleOnCylinderOneQuarterPfaffian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int m5, int m6)
{

  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m4 - m6 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double A, Ap, Y, Yp, W, Wp, Sigma, Sigmap;

  Complex Coefficient(0, 0);

  Sigma = Xr * Xr + Xs * Xs + Xr * Xs;
  Sigmap = Xrp * Xrp + Xsp * Xsp + Xrp * Xsp;

  A = (Xr - Xs) * (Xr + 2.0 * Xs) * ( Xs + 2.0 * Xr);  
  Ap = (Xrp - Xsp) * (Xrp + 2.0 * Xsp) * ( Xsp + 2.0 * Xrp);  

  Y = -27.0 * Xr * Xs * (Xr + Xs);
  Yp = -27.0 * Xrp * Xsp * (Xrp + Xsp);  

  Coefficient.Re = A * Ap;

  Coefficient.Re += A * (2.0 - Sigma) * Ap * (2.0 - Sigmap); 

  Coefficient.Re += (9.0/5.0) * A * Y * Ap * Yp/729.0;

  Coefficient.Re += 10.0 * A * (1.0 - Sigma + 0.2 * Sigma * Sigma) * Ap * (1.0 - Sigmap + 0.2 * Sigmap * Sigmap);

  Coefficient.Re += (7.0/5.0) * A * Y * (1.0 - (2.0/7.0) * Sigma) * Ap * Yp * (1.0 - (2.0/7.0) * Sigmap)/81.0;

  Coefficient.Re += 0.2 * A * (-10.0 + 15.0 * Sigma - 6.0 * Sigma * Sigma + (2.0/3.0) * Sigma * Sigma * Sigma) * Ap * (-10.0 + 15.0 * Sigmap - 6.0 * Sigmap * Sigmap + (2.0/3.0) * Sigmap * Sigmap * Sigmap);
  
  Coefficient.Re += (1.0/(9.0 * 105.0)) * A * (Sigma * Sigma * Sigma - 27.0 * Y * Y) * Ap * (Sigmap * Sigmap * Sigmap - 27.0 * Yp * Yp)/1061425.0;

  Coefficient.Re *= exp(-Sigma -Sigmap);

  return (Coefficient * 16.0 * sqrt(M_PI) * sqrt(3.0 * M_PI)/(2.0 * M_PI * this->Ratio * this->NbrLzValue));
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnCylinderOneQuarterPfaffian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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
      int ReducedNbrInteractionFactors2b = this->NbrInteractionFactors2b - 1;

      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();

      //---------------- 3b --------------------------------
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

      //---------------- 2b --------------------------------
      for (int j = 0; j < ReducedNbrInteractionFactors2b; ++j) 
	{
	  m1 = this->M1Value2b[j];
	  m2 = this->M2Value2b[j];
	  m3 = this->M3Value2b[j];
	  m4 = this->M4Value2b[j];
	  TmpInteraction = this->InteractionFactors2b[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index <= i)
               { 
		 vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
                 if (Index < i)
                    vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
               }
	    }
	}


      //---------------- 3b --------------------------------
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
	}

      //---------------- 2b --------------------------------
      m1 = this->M1Value2b[ReducedNbrInteractionFactors2b];
      m2 = this->M2Value2b[ReducedNbrInteractionFactors2b];
      m3 = this->M3Value2b[ReducedNbrInteractionFactors2b];
      m4 = this->M4Value2b[ReducedNbrInteractionFactors2b];
      TmpInteraction = this->InteractionFactors2b[ReducedNbrInteractionFactors2b];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index <= i)
           {
	    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
            if (Index < i)
              vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
           }
	}

      for (int i = firstComponent; i < LastComponent; ++i)
	 vDestination[i] += this->EnergyShift * vSource[i];

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
	  int ReducedNbrInteractionFactors2b = this->NbrInteractionFactors2b - 1;

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
	      }


          //----------------------- 2b -----------------------------
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < ReducedNbrInteractionFactors2b; ++j) 
		  {
		    m1 = this->M1Value2b[j];
		    m2 = this->M2Value2b[j];
		    m3 = this->M3Value2b[j];
		    m4 = this->M4Value2b[j];
		    TmpInteraction = this->InteractionFactors2b[j];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
			Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
			if (Index <= i)
                         {
			  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
                          if (Index < i)
                            vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
                         }
		      }
		  }
		m1 = this->M1Value[ReducedNbrInteractionFactors2b];
		m2 = this->M2Value[ReducedNbrInteractionFactors2b];
		m3 = this->M3Value[ReducedNbrInteractionFactors2b];
		m4 = this->M4Value[ReducedNbrInteractionFactors2b];
		TmpInteraction = this->InteractionFactors2b[ReducedNbrInteractionFactors2b];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    if (Index <= i)
                     {
		       vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
                       if (Index < i)
                        vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
                     }
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

long ParticleOnCylinderOneQuarterPfaffian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
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

      //-------------------2b--------------------------
      for (int j = 0; j < this->NbrInteractionFactors2b; ++j) 
	{
	  m1 = this->M1Value2b[j];
	  m2 = this->M2Value2b[j];
	  m3 = this->M3Value2b[j];
	  m4 = this->M4Value2b[j];
	  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index <= i)
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i];
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

void ParticleOnCylinderOneQuarterPfaffian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
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

      //-------------2b----------------------- 
      for (int j = 0; j < this->NbrInteractionFactors2b; ++j) 
	{
	  m1 = this->M1Value2b[j];
	  m2 = this->M2Value2b[j];
	  m3 = this->M3Value2b[j];
	  m4 = this->M4Value2b[j];
	  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index <= i)
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors2b[j];
	      ++Pos;
	    }
	}

     ++PosIndex;
   }

  delete TmpParticles;
}
