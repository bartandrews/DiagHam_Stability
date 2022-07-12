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


#include "Hamiltonian/ParticleOnCylinderPairAmplitude.h"
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
// pairCenter = center of the pair
// pairAnisotropy = anisotropy of the pair
// pairMomentum = angular momentum of the pair
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderPairAmplitude::ParticleOnCylinderPairAmplitude(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   double ratio, int pairCenter, double pairAnisotropy, int pairMomentum, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->PairCenter = pairCenter;
  this->PairAnisotropy = pairAnisotropy;
  this->PairMomentum = pairMomentum;
 
  this->HermiteCoefficients = new int[this->PairMomentum + 1];
  this->ComputeHermiteCoefficients(this->PairMomentum);

  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;

  double Length = sqrt(2.0 * M_PI * this->NbrLzValue * this->Ratio);
  this->Kappa = 2.0 * M_PI/Length;

  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;

  this->OneBodyInteractionFactors = 0;

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

ParticleOnCylinderPairAmplitude::~ParticleOnCylinderPairAmplitude() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;

  if (this->OneBodyInteractionFactors != 0)
    delete[] this->OneBodyInteractionFactors;

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

void ParticleOnCylinderPairAmplitude::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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

void ParticleOnCylinderPairAmplitude::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderPairAmplitude::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m2, m4;
  Complex* TmpCoefficient = new Complex [this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m2 = this->PairCenter - m1; 	
	      m4 = this->PairCenter - m3;
	      if ((m2 >= 0) && (m2 <= this->MaxMomentum) && (m4 >= 0) && (m4 <= this->MaxMomentum))
		  {
		    TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4);

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
      this->InteractionFactors = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m2 = this->PairCenter - m1; 	
	      m4 = this->PairCenter - m3;
	      if ((m2 >= 0) && (m2 <= this->MaxMomentum) && (m4 >= 0) && (m4 <= this->MaxMomentum))
		  {
		    if  (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		      {
		        this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		        this->M1Value[this->NbrInteractionFactors] = m1;
		        this->M2Value[this->NbrInteractionFactors] = m2;
		        this->M3Value[this->NbrInteractionFactors] = m3;
		        this->M4Value[this->NbrInteractionFactors] = m4;
                        //cout<<"m1= "<<m1<<" m2="<<m2<<" m3="<<m3<<" m4="<<m4<<" : "<<TmpCoefficient[Pos]<<endl; 
		        ++this->NbrInteractionFactors;
		      }
		    ++Pos;
		  }
	    }
    }
  else //bosons
    {
/*
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m4 <= this->MaxMomentum))
 	       if (m3 > m4)
		 {
		  if (m1 != m2)
		    {
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		    }
		  else
		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3));
		  if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		    MaxCoefficient = Norm(TmpCoefficient[Pos]);
		  ++Pos;
		}
	      else
		if (m3 == m4)
		  {
		    if (m1 != m2)
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		    else
		      TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4);
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
      this->InteractionFactors = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m4 <= this->MaxMomentum))
	       if (m3 >= m4)
		 {
		  if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
		}
	     }
*/
    }
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}

/*
double  ParticleOnCylinderPairAmplitude::Hermite(int m, double x)
{
   double Coeff = 0.0;
   if (m==0)
     Coeff = 1.0;
   else if (m==1)
     Coeff = 2.0 * x;
   else if (m==2)
     Coeff = 4.0 * pow(x, 2.0) - 2.0;
   else if (m==3)
     Coeff = 8.0 * pow(x, 3.0) - 12.0 * x;
   else if (m==4)
     Coeff = 16.0 * pow(x, 4.0) - 48.0 * pow(x, 2.0) + 12.0;
   else if (m==5)
     Coeff = 32.0 * pow(x, 5.0) - 160.0 * pow(x, 3.0) + 120.0 * x;
   else if (m==6)
     Coeff = 64.0 * pow(x, 6.0) - 480.0 * pow(x, 4.0) + 720.0 * pow(x, 2.0)- 120.0;
   else if (m==7)
     Coeff = 128.0 * pow(x, 7.0) - 1344.0 * pow(x, 5.0) + 3360.0 * pow(x, 3.0) - 1680.0 * x;
   else if (m==8)
     Coeff = 256.0 * pow(x, 8.0) - 3584.0 * pow(x, 6.0) + 13440.0 * pow(x, 4.0) - 13440.0 * pow(x, 2.0) + 1680;
   else if (m==9)
     Coeff = 512.0 * pow(x, 9.0) - 9216.0 * pow(x, 7.0) + 4838.0 * pow(x, 5.0) - 80640.0 * pow(x, 3.0) + 30240.0 * x;
   else
     {
       cout << "Not implemented." << endl;
       exit(2);
     }

   return Coeff;
}
*/

void ParticleOnCylinderPairAmplitude::ComputeHermiteCoefficients(int MaxIndex)  
{
  int i,j;
 
  if (MaxIndex == 0)
    {
      this->HermiteCoefficients[0] = 1;
    }
  else if (MaxIndex == 1)
   {
     this->HermiteCoefficients[0] = 0;
     this->HermiteCoefficients[1] = 2;
   }
  else
   {
      int* Step0 = new int[MaxIndex + 1];
      for (int i = 0; i <= MaxIndex; i++)
        Step0[i] = 0;
      Step0[0] = 1;	

      int* Step1 = new int[MaxIndex + 1];
      for (int i = 0; i <= MaxIndex; i++)
        Step1[i] = 0;
      Step1[1] = 2;

      int* StepN = new int[MaxIndex+1]; 
      for (int i = 0; i <= MaxIndex; i++)
        StepN[i] = 0; 

      for (int M = 2; M <= MaxIndex; M++)
        {
           for (int i = 0; i <= MaxIndex; i++)
             StepN[i] = -2 * (M-1) * Step0[i];

           for (int i = 1; i <= MaxIndex; i++)
             StepN[i] += 2 * Step1[i-1];
 
           for (int i = 0; i <= MaxIndex; i++)
             Step0[i] = Step1[i];            
           for (int i = 0; i <= MaxIndex; i++)
             Step1[i] = StepN[i];            
        }

      for (int i = 0; i <= MaxIndex; i++)
         this->HermiteCoefficients[i] = StepN[i];
   }
}


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnCylinderPairAmplitude::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Xm1 = this->Kappa * m1;
  double Xm2 = this->Kappa * m2;
  double Xm3 = this->Kappa * m3;
  double Xm4 = this->Kappa * m4;

  double Xr = 0.5 * (Xm1 - Xm2);
  double Xrp = 0.5 * (Xm4 - Xm3);

  long double Fact = 1.0;
  for(int i = 1; i <= this->PairMomentum; i++)
    Fact *= i;
  
  Complex Coefficient(0.0, 0.0); 

  double HermitePol1 = 0;
  for (int i = 0; i <= this->PairMomentum; ++i)
    HermitePol1 += this->HermiteCoefficients[i] * pow(Xr * sqrt(2.0/this->PairAnisotropy), i);
  double HermitePol2 = 0;
  for (int i = 0; i <= this->PairMomentum; ++i)
    HermitePol2 += this->HermiteCoefficients[i] * pow(Xrp * sqrt(2.0/this->PairAnisotropy), i);

  Coefficient.Re = exp(-Xr*Xr/this->PairAnisotropy-Xrp*Xrp/this->PairAnisotropy) * HermitePol1 * HermitePol2 /(pow(2.0, this->PairMomentum) * Fact * sqrt(this->PairAnisotropy));

  return Coefficient;
}
