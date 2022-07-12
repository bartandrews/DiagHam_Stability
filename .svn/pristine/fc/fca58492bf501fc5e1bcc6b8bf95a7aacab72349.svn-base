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


#include "Hamiltonian/ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian.h"
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
// theta1 = hypermetric angle theta1
// theta2 = hypermetric angle theta2
// confinement = amplitude of the quadratic confinement potential
// electricFieldParameter = amplitude of the electric field along the cylinder
// bFieldfParameter = amplitude of the magnetic field (to set the energy scale)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   double ratio, double theta1, double theta2, double confinement, double electricFieldParameter, double bFieldParameter, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->HypermetricTheta1 = theta1;
  this->HypermetricTheta2 = theta2;
  this->Architecture = architecture;
  this->Confinement = confinement;
  this->ElectricField = electricFieldParameter;
  this->MagneticField = bFieldParameter;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;
  this->HermitianSymmetryFlag=true;

  this->OneBodyInteractionFactors = 0;
  if ((this->ElectricField != 0) || (this->Confinement != 0))
    {
      this->OneBodyInteractionFactors = new Complex [this->NbrLzValue];
      Complex Factor;
      double kappa = sqrt(2.0 * M_PI /(this->NbrLzValue * this->Ratio));
      for (int i = 0; i < this->NbrLzValue; ++i)
        { 
           Factor.Re = this->Confinement * pow(kappa * (i - 0.5 * this->MaxMomentum), 2.0);
           Factor.Im = 0.0;
           //add contribution from electric field
           Factor.Re += 0.194 * sqrt(this->MagneticField) * ((this->ElectricField/(1.0 + this->ElectricField)) * kappa * kappa * ((double)i - 0.5 * this->MaxMomentum) * ((double)i - 0.5 * this->MaxMomentum)); 
           Factor.Im += 0.0;
	   this->OneBodyInteractionFactors[i] = Factor;
        }
    }

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

ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::~ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;
  delete[] this->M5Value;
  delete[] this->M6Value;


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

void ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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

void ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m6;
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * (this->NbrLzValue - 1) * (this->NbrLzValue - 2)/6) * this->NbrLzValue * (this->NbrLzValue - 1)/2];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < m2; ++m3)
            for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
	      for (int m5 = 0; m5 < m4; ++m5)
 	       {
	         m6 = m1 + m2 + m3 - m4 - m5;
	         if ((m6 >= 0) && (m6 <= this->MaxMomentum))
  	           if (m6 < m5)
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
	         if ((m6 >= 0) && (m6 <= this->MaxMomentum))
  	           if (m6 < m5)
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

     cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient;
    }
  else //bosons
    {
      Complex* TmpCoefficient = new Complex [(this->NbrLzValue * (this->NbrLzValue + 1) * (this->NbrLzValue + 2)/6) * this->NbrLzValue * (this->NbrLzValue + 1)/2];

      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= m2; ++m3)
            for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
	      for (int m5 = 0; m5 <= m4; ++m5)
 	       {
	         m6 = m1 + m2 + m3 - m4 - m5;
	         if ((m6 >= 0) && (m6 <= this->MaxMomentum))
  	           if (m6 <= m5)
		     {
 		       TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientBosons(m1, m2, m3, m4, m5, m6) * this->NumberOfPermutations(m1, m2, m3) * this->NumberOfPermutations(m4, m5, m6);
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
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= m2; ++m3)
            for (int m4 = 0; m4 <= this->MaxMomentum; ++m4)
	      for (int m5 = 0; m5 <= m4; ++m5)
 	       {
	         m6 = m1 + m2 + m3 - m4 - m5;
	         if ((m6 >= 0) && (m6 <= this->MaxMomentum))
  	           if (m6 <= m5)
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

     cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
     cout << "====================================" << endl;
     delete[] TmpCoefficient;
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

Complex ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int m5, int m6)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m4 - m6 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double GaussianExp;

  Complex Coefficient(0,0);

  if ((this->HypermetricTheta1 != 0.0) || (this->HypermetricTheta2 != 0.0))
   {
     cout << "Hypermetric for fermions not implemented. " << endl;
     exit(2);
   }

  if (this->ElectricField == 0)
   {
     GaussianExp = Xr * Xr + Xs * Xs + Xr * Xs;
     Coefficient.Re = exp(-GaussianExp) * (Xr - Xs) * (2.0 * Xr + Xs) * (2.0 * Xs + Xr);
     Coefficient.Im = 0.0;

     GaussianExp = Xrp * Xrp + Xsp * Xsp + Xrp * Xsp;
     Coefficient.Re *= (exp(-GaussianExp) * (Xrp - Xsp) * (2.0 * Xrp + Xsp) * (2.0 * Xsp + Xrp));
     Coefficient.Im = 0.0;
     return (Coefficient * 16.0 * sqrt(M_PI) * sqrt(3.0 * M_PI)/(2.0 * M_PI * this->Ratio * this->NbrLzValue));
   }
  else
   {
     /*
     double alpha = sqrt(1.0 + this->ElectricField);
     Coefficient.Re = exp(-pow(Xm1-Xm3,2.0)/(2.0*pow(alpha,3.0))-pow(Xm1-Xm4,2.0)/(2.0 * pow(alpha,3.0))) * (pow(Xm1-Xm3,2.0)-alpha*alpha*pow(Xm1-Xm4,2.0)+alpha*alpha-alpha*alpha*alpha);
     Coefficient.Im = 0.0;
     return (Coefficient/sqrt(this->Ratio * this->NbrLzValue * alpha * alpha * alpha));
     */
     cout<<"Not implemented for electric fields!" << endl;
     exit(1);
   }
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a_m4 a_m5 a_m6 coupling term for bosons
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// m5 = fifth index
// m6 = sixth index
// return value = numerical coefficient

Complex ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian::EvaluateInteractionCoefficientBosons(int m1, int m2, int m3, int m4, int m5, int m6)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (2.0 * m1 - m2 - m3)/3.0;
  double Xs = kappa * (2.0 * m2 - m1 - m3)/3.0;
  double Xrp = kappa * (2.0 * m4 - m6 - m5)/3.0;
  double Xsp = kappa * (2.0 * m5 - m4 - m6)/3.0;
  double GaussianExp;

  Complex Coefficient(0,0);

  if (this->ElectricField == 0)
   {
     if ((this->HypermetricTheta1 == 0.0) && (this->HypermetricTheta2 == 0.0))
       {
         GaussianExp = Xr * Xr + Xs * Xs + Xr * Xs + Xrp * Xrp + Xsp * Xsp + Xrp * Xsp;
         Coefficient.Re = exp(-GaussianExp);
         Coefficient.Im = 0.0;
       }
     else //Hypermetric case
       {
         double GammaPlus = exp(-2.0*this->HypermetricTheta1) + 3.0 * exp(-2.0*this->HypermetricTheta2);
         double GammaMinus = exp(-2.0*this->HypermetricTheta1) - 3.0 * exp(-2.0*this->HypermetricTheta2);

	 GaussianExp = -0.25 * (GammaPlus * Xr * Xr + GammaPlus * Xs * Xs - 2.0 * GammaMinus * Xr * Xs)
		       -0.25 * (GammaPlus * Xrp * Xrp + GammaPlus * Xsp * Xsp - 2.0 * GammaMinus * Xrp * Xsp);
         Coefficient.Re = exp(GaussianExp) * exp(-this->HypermetricTheta1-this->HypermetricTheta2);
         Coefficient.Im = 0.0;
       } 

     return (Coefficient * (2.0/3.0) * sqrt(M_PI) * sqrt(3.0 * M_PI)/(2.0 * M_PI * this->Ratio * this->NbrLzValue));
   }
  else
   {
     /*
     double alpha = sqrt(1.0 + this->ElectricField);
     Coefficient.Re = exp(-pow(Xm1-Xm3,2.0)/(2.0*pow(alpha,3.0))-pow(Xm1-Xm4,2.0)/(2.0 * pow(alpha,3.0))) * (pow(Xm1-Xm3,2.0)-alpha*alpha*pow(Xm1-Xm4,2.0)+alpha*alpha-alpha*alpha*alpha);
     Coefficient.Im = 0.0;
     return (Coefficient/sqrt(this->Ratio * this->NbrLzValue * alpha * alpha * alpha));
     */
     cout<<"Not implemented for electric fields!" << endl;
     exit(1);
   }
}

// Get the number of permutations of annihilation/creation indices c_n1 c_n2 c_n3 for bosons

int ParticleOnCylinderThreeBodyLaplacianDeltaHamiltonian:: NumberOfPermutations(int n1, int n2, int n3)
{
  if ((n1 != n2) && (n2 != n3) && (n3 != n1) ) 
    return 6;
  else if ( ((n1==n2) && (n1 != n3)) || ((n1==n3) && (n1 != n2)) || ((n2==n3) && (n2 != n1)) ) 
    return 3;
  else
    return 1;
  return 0;
}
