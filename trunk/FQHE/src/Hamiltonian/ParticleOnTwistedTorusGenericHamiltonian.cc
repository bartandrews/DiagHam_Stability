////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 16/08/2011                      //
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
#include "Hamiltonian/ParticleOnTwistedTorusGenericHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;

#define M1_12 0.08333333333333333

static double MySqrArg;
#define GETSQR(a) ((MySqrArg=(a)) == 1.0 ? 1.0 : MySqrArg*MySqrArg)


// flag for switching extra output 
//#define VERBOSE_ONEBODY

// default constructor
//

ParticleOnTwistedTorusGenericHamiltonian::ParticleOnTwistedTorusGenericHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// angle = obliquity angle of torus (pi/2 = rectangular)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnTwistedTorusGenericHamiltonian::ParticleOnTwistedTorusGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   
										   double ratio, double angle, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->MaxMomentum = maxMomentum;
  this->Angle = M_PI*angle;
  this->Ratio = ratio;
  // calculate some convenient lattice geometry parameters
  this->CosTheta = cos(this->Angle);
  this->SinTheta = sqrt(1.0 - this->CosTheta * this->CosTheta);
  this->Lx = sqrt(2.0 * M_PI * (double)this->MaxMomentum * this->Ratio/this->SinTheta);
  this->Ly = sqrt(2.0 * M_PI * (double)this->MaxMomentum / (this->Ratio * this->SinTheta));
  this->Gx = 2.0 * M_PI / this->Lx;
  this->Gy = 2.0 * M_PI / this->Ly;
 
  this->LandauLevel = landauLevel;
  this->NbrPseudopotentials = nbrPseudopotentials;
  if (this->NbrPseudopotentials>0)
    {
      this->Pseudopotentials = pseudopotentials;
      this->LaguerreM=new Polynomial[NbrPseudopotentials];
      for (int i=0; i<NbrPseudopotentials; ++i)
	    this->LaguerreM[i]=LaguerrePolynomial(i);
    }
  else
    {
      this->Pseudopotentials = NULL;
      this->LaguerreM=NULL;
    }
  this->HaveCoulomb=haveCoulomb;
  if (HaveCoulomb)
    {
      if (this->LandauLevel>=0)
	{
	  // simple coulomb interactions
	  this->FormFactor=LaguerrePolynomial(this->LandauLevel);
	}
      else
	{
	  // coulomb interactions in graphene
	  this->FormFactor=0.5*(LaguerrePolynomial(abs(this->LandauLevel))+LaguerrePolynomial(abs(this->LandauLevel)-1));
	}
    }
  cout << "FormFactor=" << this->FormFactor << endl;
  if ((particles->GetHilbertSpaceDimension() > 0) && (noWignerEnergy == false))
    this->WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
  else 
    this->WignerEnergy = 0.0;
  cout << "Wigner Energy = " << WignerEnergy << endl;  
  
  this->HamiltonianShift = ((double) this->NbrParticles)*WignerEnergy;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex; 
  this->EvaluateInteractionFactors();

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

ParticleOnTwistedTorusGenericHamiltonian::~ParticleOnTwistedTorusGenericHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnTwistedTorusGenericHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrSectorSums = this->MaxMomentum;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	    	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
		for (int m2 = 0; m2 < m1; ++m2)
	  		++this->NbrSectorIndicesPerSum[(m1 + m2) % this->MaxMomentum];
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	   {
	     if (this->NbrSectorIndicesPerSum[i]  > 0)
	      {
	        this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	        this->NbrSectorIndicesPerSum[i] = 0;
	      }
	   }
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	    for (int m2 = 0; m2 < m1; ++m2)
	     {
	       int TmpSum = (m1 + m2) % this->MaxMomentum;
	       this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
	       this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
	       ++this->NbrSectorIndicesPerSum[TmpSum];    
	     }

      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	   {
	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	        {
	          int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	          int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	          for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		       {
		         int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		         int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		         int TotalMomentum=Index1+Index2-Index3-Index4;
		         if (TotalMomentum<0)
		            TotalMomentum+=this->MaxMomentum;
		         if( (TotalMomentum % this->MaxMomentum) == 0)
		          {
                     this->InteractionFactors[i][Index] = (this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4)
					      + this->EvaluateInteractionCoefficient(Index2, Index1, Index4, Index3)
					      - this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3)
					      - this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4));

		            TotalNbrInteractionFactors++;
		            ++Index;
		          }
		        }
	        }
	    }
    }
  // Bosonic case
  else
    {
      // Two-Body
      this->NbrSectorSums = this->MaxMomentum;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int ky1 = 0; ky1 < this->MaxMomentum; ++ky1)
	for (int ky2 = 0; ky2 <= ky1; ++ky2) 
	  ++this->NbrSectorIndicesPerSum[((ky1 + ky2) % this->MaxMomentum)];    
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int ky1 = 0; ky1 < this->MaxMomentum; ++ky1)
	for (int ky2 = 0; ky2 <= ky1; ++ky2) 
	  {
	    int TmpSum =  ((ky1 + ky2) % this->MaxMomentum);
	    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = ky1;
	    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = ky2;
	    ++this->NbrSectorIndicesPerSum[TmpSum];
	  }
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int TotalMomentum=Index1+Index2-Index3-Index4;
		  if (TotalMomentum<0)
		    TotalMomentum+=this->MaxMomentum;
		  if( (TotalMomentum % this->MaxMomentum) == 0)
		    {
			  Complex sumUFQHE = 0.0;
			  if (Index2 < Index1)
			    {
			      if (Index3 != Index4)
				{
				  sumUFQHE = (this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4)
					      + this->EvaluateInteractionCoefficient(Index2, Index1, Index4, Index3)
					      + this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3)
					      + this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4));
				}
			      else
				sumUFQHE = (this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4)
					    + this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4));			      
			    }
			  else
			    if (Index1 == Index2)
			      {
				if (Index3 != Index4)
				  sumUFQHE = (this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4)
					      + this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3));
				else
				  sumUFQHE = this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4);
			      }
			  this->InteractionFactors[i][Index] = sumUFQHE ;
			  
			
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}
	    }
	}
    }
  
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnTwistedTorusGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  Complex Sum(0.0, 0.0);
  double N1;
  double N2 = (double)(m1 - m4);
  double Q2, Qx, Qy;
  double Xj13 = this->Gy * (double)(m1 - m3);
  Complex Coefficient(1.0,0.0);
  double PrecisionPos, PrecisionNeg, Precision;
  Complex Phase;

  while (((fabs(Sum.Re) + fabs(Coefficient.Re)) != fabs(Sum.Re)) || ((fabs(Sum.Im) + fabs(Coefficient.Im)) != fabs(Sum.Im)))
    {
      Qx = 0.0;
      Qy = this->Gy * N2;
      Qy /= this->SinTheta;
      Q2 = Qx * Qx + Qy * Qy;

	  Coefficient.Re = this->GetVofQ(0.5 * Q2);
	  Coefficient.Im = 0.0;  	
	  if (Q2 == 0.0)
	  	Precision = 1.0;
	  else
  	   Precision = Coefficient.Re;

      N1 = 1.0;
      while ((fabs(Coefficient.Re) + fabs(Precision)) != fabs(Coefficient.Re))
	{
      //Sum over positive N1
       Qx = this->Gx * N1;
       Qy = this->Gy * N2 - this->Gx * N1 * this->CosTheta;
       Qy /= this->SinTheta;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionPos = this->GetVofQ(0.5 * Q2);          
       Phase.Re = cos(Qx * Xj13/this->SinTheta);
       Phase.Im =  -sin(Qx * Xj13/this->SinTheta);
       Coefficient += (PrecisionPos * Phase);

       //Sum over negative N1
       Qx = -this->Gx * N1;
       Qy = this->Gy * N2 + this->Gx * N1 * this->CosTheta;
       Qy /= this->SinTheta;
       Q2 = Qx * Qx + Qy * Qy;
	
       PrecisionNeg = this->GetVofQ(0.5 * Q2);   
       Phase.Re = cos(Qx * Xj13/this->SinTheta);
       Phase.Im = -sin(Qx * Xj13/this->SinTheta);
       Coefficient += (PrecisionNeg * Phase);
      //Increment N1
       N1 += 1.0;
       Precision = PrecisionPos + PrecisionNeg;
	}
      Sum += Coefficient;
      N2 += (double)this->MaxMomentum;
    }

  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = Sum;	    
  while (((fabs(Sum.Re) + fabs(Coefficient.Re)) != fabs(Sum.Re)) || ((fabs(Sum.Im) + fabs(Coefficient.Im)) != fabs(Sum.Im)))
    {
      Qx = 0.0;
      Qy = this->Gy * N2;
      Qy /= this->SinTheta;
      Q2 = Qx * Qx + Qy * Qy;

	  Coefficient.Re = this->GetVofQ(0.5 * Q2);
	  Coefficient.Im = 0.0;
	  if (Q2 == 0.0)
	  	Precision = 1.0;
	  else
  	   Precision = Coefficient.Re;

      N1 = 1.0;
      while ((fabs(Coefficient.Re) + fabs(Precision)) != fabs(Coefficient.Re))
	{
       //Sum over positive N1
       Qx = this->Gx * N1;
       Qy = this->Gy * N2 - this->Gx * N1 * this->CosTheta;
       Qy /= this->SinTheta;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionPos = this->GetVofQ(0.5 * Q2);
       Phase.Re = cos(Qx * Xj13/this->SinTheta);
       Phase.Im = -sin(Qx * Xj13/this->SinTheta);
       Coefficient += (PrecisionPos * Phase);

       //Sum over negative N1
       Qx = -this->Gx * N1;
       Qy = this->Gy * N2 + this->Gx * N1 * this->CosTheta;
       Qy /= this->SinTheta;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionNeg = this->GetVofQ(0.5 * Q2); 
       Phase.Re = cos(Qx * Xj13/this->SinTheta);
       Phase.Im = -sin(Qx * Xj13/this->SinTheta);
       Coefficient += (PrecisionNeg * Phase);
       //Increment N1
       N1 += 1.0;
       Precision = PrecisionPos + PrecisionNeg;
	}
      Sum += Coefficient;
      N2 -= (double)this->MaxMomentum;
    }
 return (Sum / (4.0 * M_PI * (double)this->MaxMomentum));
}


/*
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double  ParticleOnTwistedTorusGenericHamiltonian::RectangularEvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
  double InvRatio = 1.0/Ratio;

 
  //  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " : ";
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 1.0;
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = InvRatio * N1 * N1 + Ratio * N2 * N2;
	  Precision = 2.0 * exp(- PIOnM * Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += MaxMomentum;
    }
  N2 = (double) (m1 - m4 - MaxMomentum);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 1.0;
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = InvRatio * N1 * N1 + Ratio * N2 * N2;
	  Precision = 2.0 *  exp(- PIOnM * Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= MaxMomentum;
    }
  //  cout << Sum << endl;
  return (Sum / (4.0 * M_PI * MaxMomentum));
}



// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnTwistedTorusGenericHamiltonian::TwistedEvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  Complex Coefficient;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double PIOnMS = PIOnM / sin(this->Angle);
  double cosine = cos(this->Angle);
  double Factor =  ((double) (m1-m3)) * PIOnM * 2.0;
  Complex Sum = 0.0;
  double N2;
  double N1;
  double Q2;
  double Precision1;
  double Precision2;
  double InvRatio = 1.0/this->Ratio;


  N2 = (double) (m1 - m4);
  Coefficient = 1.0;
  while ((Norm(Sum) + Norm(Coefficient)) != Norm(Sum))
    {
      Q2 = Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = this->GetVofQ(PIOnMS*Q2);
	  Precision1 = (Precision2 = Coefficient.Re);
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnMS*Q2); // yields non-zero terms only for non-singular interactions
	  Precision1 = (Precision2 = 1.0);
	}
      N1 = 1.0;
      while ((Norm(Coefficient) + (fabs(Precision1) + fabs(Precision2))) != Norm(Coefficient))
	{
	  Q2 = InvRatio * N1 * N1 - 2 * N1 * N2 * cosine + Ratio * N2 * N2;
	  Precision1 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision1 * Phase(N1 * Factor);
	  
	  Q2 = InvRatio * N1 * N1 + 2 * N1 * N2 * cosine + Ratio * N2 * N2;
	  Precision2 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision2 * Phase(- N1 * Factor);

	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
    }

  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = 1.0;
  while ((Norm(Sum) + Norm(Coefficient)) != Norm(Sum))
    {
      Q2 = Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = this->GetVofQ(PIOnMS*Q2);
	  Precision1 = (Precision2 = Norm(Coefficient.Re));
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnMS*Q2); // yields non-zero terms only for non-singular interactions
	  Precision1 = (Precision2 = 1.0);
	}
      N1 = 1.0;
      while ((Norm(Coefficient) + fabs(Precision1) + fabs(Precision2)) != Norm(Coefficient))
	{
	  Q2 = InvRatio * N1 * N1 - 2 * N1 * N2 * cosine + Ratio * N2 * N2;
	  Precision1 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision1 * Phase(N1 * Factor);

	  Q2 = InvRatio * N1 * N1 + 2 * N1 * N2 * cosine + Ratio * N2 * N2;
	  Precision2 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision2 * Phase(- N1 * Factor);

	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (8.0 * M_PI * this->MaxMomentum));
}
*/

// get fourier transform of interaction
// Q2_half = one half of q² value
double ParticleOnTwistedTorusGenericHamiltonian::GetVofQ(double Q2_half)
{ 
  double Result;
  double Q2=2.0*Q2_half;
  if ((this->HaveCoulomb)&&(Q2_half!=0.0))
    {
      Result=GETSQR(this->FormFactor(Q2_half)) * (2.0 * M_PI)/ sqrt(Q2);
    }
  else
    Result=0.0;
    for (int i=0; i<NbrPseudopotentials; ++i)
      if (this->Pseudopotentials[i]!=0.0)
        Result += 2.0 * (2.0 * M_PI) * this->Pseudopotentials[i]*this->LaguerreM[i].PolynomialEvaluate(Q2);
  return Result * exp(-Q2_half);
}



// evaluate Wigner crystal energy per particle
//
// return value = Wigner crystal energy per particle

double ParticleOnTwistedTorusGenericHamiltonian::EvaluateWignerCrystalEnergy ()
{
  double TmpRatio = M_PI * this->Ratio;
  double TmpInvRatio = M_PI / this->Ratio;
  double Energy = this->MisraFunction(-0.5, TmpRatio);
  double Precision = Energy;
  int L1 = 2;
  while ((Energy + Precision) > Energy)
    {
      Precision = this->MisraFunction(-0.5, TmpRatio * L1 * L1);
      Energy += Precision;
      ++L1;
    }
  Energy *= 2.0;
  int L2 = 1;
  double PartialEnergy = Energy;
  while ((PartialEnergy + Energy) > Energy)
    {
      PartialEnergy = 2.0 * this->MisraFunction(-0.5, TmpInvRatio * L2 * L2);
      Precision = PartialEnergy;
      L1 = 1;
      while (((PartialEnergy + Precision) > PartialEnergy))// && ((fabs(PartialEnergy - Precision) + Energy) > Energy))
	{
	  Precision = 4.0 * this->MisraFunction(-0.5, TmpRatio * L1 * L1 + TmpInvRatio * L2 * L2);
	  PartialEnergy += Precision;
	  ++L1;	  
	}
      Energy += PartialEnergy;
      ++L2;
    }
  return 2.0 * (Energy - 2.0) / sqrt (2.0 * M_PI * this->MaxMomentum);
}

// evaluate Misra function (integral of t^n exp (-xt) between 1 and +inf)
//
// n = index of the Misra function
// x = point where the function has to be evaluated (> 0)
// return value = value of the n-Misra function at x

double ParticleOnTwistedTorusGenericHamiltonian::MisraFunction (double n, double x)
{
  int Count=0;
  int NbrSubdivision = 100000;
  double PreviousSum = this->PartialMisraFunction(n, x, 0.0, 1.0, NbrSubdivision);
  double NewSum = PreviousSum;
  PreviousSum *= 2.0;
  while (((fabs(PreviousSum - NewSum) / PreviousSum) > MACHINE_PRECISION) && (Count<5))
    {
      if ((fabs(PreviousSum - NewSum) / PreviousSum) < 1e-11)
	++Count;
      PreviousSum = NewSum;
      NbrSubdivision += 10000;
      NewSum = this->PartialMisraFunction(n, x, 0.0, 1.0, NbrSubdivision);
      //cout << " PreviousSum = " << PreviousSum << "   NewSum = " << NewSum << "  diff="<<PreviousSum-NewSum<<endl;
    }
  return 2.0 * (sqrt(M_PI * 0.25 / x) - NewSum);
}

// evaluate part of the integral needed in the Misra function (integral of t^n exp (-xt) between min and max)
//
// n = index of the Misra function
// x = point where the function has to be evaluated (> 0)
// min = lower bound of the integral
// max = upper bound of the integral
// nbrSubdivision = number of subdivision used for the integral
// return value = value of the integral

double ParticleOnTwistedTorusGenericHamiltonian::PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision)
{
  double Sum = 0.0;
  x *= -1.0;
  --nbrSubdivision;
  max  = (max - min) / ((double) nbrSubdivision);
  Sum += (0.5 + M1_12 * 2.0 * x * min * max )* exp(min * min * x);
  min += max;
  --nbrSubdivision;
  while (nbrSubdivision > 0)
    {
      Sum += exp(min * min * x);
      min += max;
      --nbrSubdivision;
    }
  Sum += (0.5 - M1_12 * 2.0 * x * min * max) * exp(min * min * x);
  Sum *= max;
  return Sum;
}
