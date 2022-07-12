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


#include "Hamiltonian/ParticleOnCylinderFourBodyDeltaHamiltonian.h"
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
// confinement = amplitude of the quadratic confinement potential
// electricFieldParameter = amplitude of the electric field along the cylinder
// bFieldfParameter = amplitude of the magnetic field (to set the energy scale)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderFourBodyDeltaHamiltonian::ParticleOnCylinderFourBodyDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   double ratio, double confinement, double electricFieldParameter, double bFieldParameter, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
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

ParticleOnCylinderFourBodyDeltaHamiltonian::~ParticleOnCylinderFourBodyDeltaHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->CreationIndices;
  delete[] this->AnnihilationIndices;


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

void ParticleOnCylinderFourBodyDeltaHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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

void ParticleOnCylinderFourBodyDeltaHamiltonian::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderFourBodyDeltaHamiltonian::EvaluateInteractionFactors()
{
  unsigned L8Mask = (1u<<8)-1;
  unsigned L16Mask = ((1u<<16)-1) - ((1u<<8)-1);
  unsigned L24Mask = ((1u<<24)-1) - ((1u<<16)-1);
  unsigned L32Mask = (~((1u<<24)-1));
  if (this->MaxMomentum >= 63) 
    {
       cout<<"Overflow"<<endl;
       exit(1);
    }
  int Norb = this->MaxMomentum + 1;
  int Pos = 0;
  int m8;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      cout<<"Fermions unsupported."<<endl;
      exit(1);
    }
  else //bosons
    {
  int NbrIndices = Norb * (Norb+1) * (Norb+2) * (Norb+3)/24;
  NbrIndices *= (Norb * (Norb+1) * (Norb+2)/6);
  cout<<"Total nbr of indices: "<<NbrIndices<<endl;
 
  long TmpMemory = NbrIndices * 2 * (2 * sizeof(unsigned) + sizeof(Complex)); 
  if (TmpMemory < 1024)
    cout  << "Need to allocate temporarily for matrix elements = " <<  TmpMemory << "b ";
  else
    if (TmpMemory < (1 << 20))
      cout  << "Need to allocate temporarily for matrix elements = " << (TmpMemory >> 10) << "kb ";
    else
      if (TmpMemory < (1 << 30))
	cout  << "Need to allocate temporarily for matrix elements = " << (TmpMemory >> 20) << "Mb ";
      else
	cout  << "Need to allocate temporarily for matrix elements = " << (TmpMemory >> 30) << "Gb ";
  cout << endl;
  //if (TmpMemory >= memory)
  //  {
  //    cout << "Insufficient memory...exiting."<<endl;
  //    exit(1);   
  //  }

  Complex* TmpCoefficient = new Complex [NbrIndices];
  unsigned* TmpIndicesCreation = new unsigned [NbrIndices];
  unsigned* TmpIndicesAnnihilation = new unsigned [NbrIndices];
  double MaxCoefficient=0.0;
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= m2; ++m3)
	    for (int m4 = 0; m4 <= m3; ++m4)
	      for (int m5 = 0; m5 <= this->MaxMomentum; ++m5)
	        for (int m6 = 0; m6 <= m5; ++m6)
	          for (int m7 = 0; m7 <= m6; ++m7)
		    {
		      m8 = m1 + m2 + m3 + m4 - m5 - m6 - m7;
		      if ((m8 >= 0) && (m8 <= this->MaxMomentum))
 	                if (m8 <= m7)
			  {
			       TmpIndicesCreation[Pos] = (((m1&L8Mask) | ((m2&L8Mask)<<8)) | ((m3&L8Mask)<<16)) | ((m4&L8Mask)<<24);
			       TmpIndicesAnnihilation[Pos] = (((m5&L8Mask)|((m6&L8Mask)<<8))|((m7&L8Mask)<<16))|((m8&L8Mask)<<24);
  		               TmpCoefficient[Pos] = this->EvaluateInteractionCoefficientBosons(m1, m2, m3, m4, m5, m6, m7, m8) * this->NumberOfPermutations(m1, m2, m3, m4) * this->NumberOfPermutations(m5, m6, m7, m8);
		               if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		                 MaxCoefficient = Norm(TmpCoefficient[Pos]);
   		  	       ++Pos;
			  }
	            }

   cout << "Total Nbr of Indices: " << Pos << endl;


   MaxCoefficient *= MACHINE_PRECISION;

  this->CreationIndices = new unsigned [Pos];
  this->AnnihilationIndices = new unsigned [Pos];
  this->InteractionFactors = new Complex [Pos]; 

  this->NbrInteractionFactors = 0;
    for (int i = 0; i < Pos; ++i)
     {
       if  (Norm(TmpCoefficient[i]) > MaxCoefficient)
	{
	 this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[i];
	 this->AnnihilationIndices[this->NbrInteractionFactors] = TmpIndicesAnnihilation[i];
	 this->CreationIndices[this->NbrInteractionFactors] = TmpIndicesCreation[i];

           //****** Decode m1,...,m8 *************
           //int m1 = this->CreationIndices[this->NbrInteractionFactors] & L8Mask;
           //int m2 = (this->CreationIndices[this->NbrInteractionFactors] & L16Mask)>>8;
           //int m3 = (this->CreationIndices[this->NbrInteractionFactors] & L24Mask)>>16;
           //int m4 = (this->CreationIndices[this->NbrInteractionFactors] & L32Mask)>>24;
      
           //int m5 = this->AnnihilationIndices[this->NbrInteractionFactors] & L8Mask;
           //int m6 = (this->AnnihilationIndices[this->NbrInteractionFactors] & L16Mask)>>8;
           //int m7 = (this->AnnihilationIndices[this->NbrInteractionFactors] & L24Mask)>>16;
           //m8 = (this->AnnihilationIndices[this->NbrInteractionFactors] & L32Mask)>>24;
           //***************************************
           //cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<" "<<m7<<" "<<m8<<" "<<this->InteractionFactors[this->NbrInteractionFactors]<<endl;
	 ++this->NbrInteractionFactors;
	}
    }

  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
  delete[] TmpIndicesCreation;
  delete[] TmpIndicesAnnihilation;
 }
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a^+_m4 a_m5 a_m6 a_m7 a_m8 coupling term
//
// m1,..,m8 = indices from left to right
// return value = numerical coefficient

Complex ParticleOnCylinderFourBodyDeltaHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int m5, int m6, int m7, int m8)
{
  Complex Result(0,0);
  return Result;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a^+_m3 a^+_m4 a_m5 a_m6 a_m7 a_m8 coupling term for bosons
//
// m1,..,m8 = indices from left to right
// return value = numerical coefficient

Complex ParticleOnCylinderFourBodyDeltaHamiltonian::EvaluateInteractionCoefficientBosons(int m1, int m2, int m3, int m4, int m5, int m6, int m7, int m8)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xr = kappa * (3.0 * m1 - m2 - m3 - m4)/4.0;
  double Xs = kappa * (3.0 * m2 - m1 - m3 - m4)/4.0;
  double Xt = kappa * (3.0 * m3 - m1 - m2 - m4)/4.0;

  double Xrp = kappa * (3.0 * m5 - m6 - m7 - m8)/4.0;
  double Xsp = kappa * (3.0 * m6 - m5 - m7 - m8)/4.0;
  double Xtp = kappa * (3.0 * m7 - m5 - m6 - m8)/4.0;

  double GaussianExp;

  Complex Coefficient(0.0 , 0.0);

  if (this->ElectricField == 0)
   {
     //Momentum 0, delta interaction
     GaussianExp = Xr * Xr + Xs * Xs + Xt * Xt + Xr * Xs + Xr * Xt + Xs * Xt;
     GaussianExp += Xrp * Xrp + Xsp * Xsp + Xtp * Xtp + Xrp * Xsp + Xrp * Xtp + Xsp * Xtp;
     Coefficient += exp(-GaussianExp);

     //Interaction for "Gaffnian" 3/7	

     //Momentum 2
     //double W2 = (8.0/3.0) * (Xr * Xr + Xs * Xs + Xt * Xt + Xr * Xs + Xr * Xt + Xs * Xt);
     //double Wp2 = (8.0/3.0) * (Xrp * Xrp + Xsp * Xsp + Xtp * Xtp + Xrp * Xsp + Xrp * Xtp + Xsp * Xtp);
     //Coefficient += (1.5 * (-1.0 + 0.5 * W2) * (-1.0 + 0.5 * Wp2)) * exp(-GaussianExp);

     //Momentum 3; had to divide by extra factor sqrt(3) wrt to PRX paper
     //double S3 = - Xs * Xt * (Xr + Xs + Xt) - Xr * Xt * (Xr + Xs + Xt) - Xr * Xs * (Xr + Xs + Xt) - Xr * Xs * Xt; 
     //double Sp3 = - Xsp * Xtp * (Xrp + Xsp + Xtp) - Xrp * Xtp * (Xrp + Xsp + Xtp) - Xrp * Xsp * (Xrp + Xsp + Xtp) - Xrp * Xsp * Xtp; 
     //Coefficient += (8.0 * S3 * Sp3)/sqrt(3.0) * exp(-GaussianExp);

     //Momentum 4, first sector
     //Coefficient += (1.0/64.0) * 0.3 * (20.0 - 20.0 * W2 + 3.0 * W2 * W2) * (20.0 - 20.0 * Wp2 + 3.0 * Wp2 * Wp2) * exp(-GaussianExp);

     //Momentum 4, second sector
     //double S4 = - Xr * Xs * Xt * (Xr + Xs + Xt);
     //double Sp4 = - Xrp * Xsp * Xtp * (Xrp + Xsp + Xtp);

     //Coefficient += 0.2 * ((-3.0/32.0) * W2 * W2 + (40.0/3.0) * S4) * ((-3.0/32.0) * Wp2 * Wp2 + (40.0/3.0) * Sp4) * exp(-GaussianExp);
     

     return (Coefficient * (M_PI/3.0) * sqrt(4.0 * M_PI)/pow(Length, 3.0));
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

// Get the number of permutations of annihilation/creation indices c_n1 c_n2 c_n3 c_n4 for bosons

// Get the number of permutations

int ParticleOnCylinderFourBodyDeltaHamiltonian:: NumberOfPermutations(int n1, int n2, int n3, int n4)
{
  if ( (n1!=n2) && (n1!=n3) && (n1!=n4) && (n2!=n3) && (n2!=n4) && (n3!=n4) )
    return 24; //all different
  else if ( ((n1==n2) && (n1!=n3) && (n1!=n4) && (n3!=n4)) || ((n1==n3) && (n1!=n2) && (n1!=n4) && (n2!=n4)) || ((n1==n4) && (n1!=n2) && (n1!=n3) && (n2!=n3)) || ((n2==n3) && (n2!=n1) && (n2!=n4) && (n1!=n4)) || ((n2==n4) && (n2!=n1) && (n2!=n3) && (n1!=n3)) || ((n3==n4) && (n3!=n1) && (n3!=n2) && (n1!=n2))) 
    return 12; //only two are equal
  else if ( ((n1==n2) && (n3==n4) && (n1!=n3)) || ((n1==n3) && (n2==n4) && (n1!=n2)) || ((n1==n4) && (n2==n3) && (n1!=n2)) )
    return 6; //two groups of two
  else if ( ((n1==n2) && (n1==n3) && (n1!=n4)) || ((n1==n2) && (n1==n4) && (n1!=n3)) || ((n1==n3) && (n3==n4) && (n1!=n2)) || ((n2==n3) && (n2==n4) && (n2!=n1)) )
    return 4; //three are equal
  else if ((n1==n2) && (n2==n3) && (n3==n4))
    return 1; //all four are equal
  return 0;
}
