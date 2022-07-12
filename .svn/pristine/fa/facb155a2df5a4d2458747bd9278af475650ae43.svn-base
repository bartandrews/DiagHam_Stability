////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//       class of hamiltonian associated to particles on a twisted torus      //
//             with coulomb interaction and magnetic translations             //
//                                                                            //
//                        last modification : 11/04/2012                      //
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

#include "Hamiltonian/ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Polynomial/SpecialPolynomial.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333

static double MySqrArg;
#define GETSQR(a) ((MySqrArg=(a)) == 1.0 ? 1.0 : MySqrArg*MySqrArg)

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the lengths of the two fundamental cycles of the torus L1 / L2
// angle = angle between the two fundamental cycles of the torus, along (L1 sin, L1 cos) and (0, L2)
// haveCoulomb = flag indicating whether a coulomb term is present
// landauLevel = landauLevel to be simulated (GaAs (>=0) or graphene (<0))
// nbrPseudopotentials = number of pseudopotentials indicated
// pseudopotentials = pseudopotential coefficients
// noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
        int nbrParticles, int maxMomentum, int xMomentum, double ratio, double angle, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy,
        AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Angle = angle;
  // calculate some convenient lattice geometry parameters
  this->CosTheta = cos(this->Angle);
  this->SinTheta = sqrt(1.0 - this->CosTheta * this->CosTheta);
  this->Lx = sqrt(2.0 * M_PI * (double)this->MaxMomentum * this->Ratio/this->SinTheta);
  this->Ly = sqrt(2.0 * M_PI * (double)this->MaxMomentum / (this->Ratio * this->SinTheta));
  this->Gx = 2.0 * M_PI / this->Lx;
  this->Gy = 2.0 * M_PI / this->Ly;
  cout << "Twisted torus cos(angle) = " << this->CosTheta << endl;

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

  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  cout << "Wigner Energy = " << WignerEnergy << endl;  
  this->EvaluateInteractionFactors();
  this->EnergyShift = ((double) this->NbrParticles)*WignerEnergy;
  this->CosinusTable = new double [this->MaxMomentum];
  this->SinusTable = new double [this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->CosinusTable[i] = cos(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
      this->SinusTable[i] = sin(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
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

ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::~ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
  if (this->FastMultiplicationFlag == true)
    {
      int ReducedDim = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	  delete[] this->InteractionPerComponentNbrTranslation[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      delete[] this->InteractionPerComponentNbrTranslation;
    }
  if (this->NbrPseudopotentials>0)
    delete [] this->LaguerreM;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->EnergyShift=shift+this->WignerEnergy;
}
  
// evaluate all interaction factors
//   

void ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
  Complex* TmpCoefficient = new Complex [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnTorusWithMagneticTranslations::FermionicStatistic)
    {
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
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
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  if  (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
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
    }
  else
    {
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					 + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					 + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (m1 == m2)
		    TmpCoefficient[Pos] *= 0.5;
		  if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
		    MaxCoefficient = Norm(TmpCoefficient[Pos]);
		  ++Pos;
		}
	      else
		{
		  if (m3 == m4)
		    {
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		      if (m1 == m2)
			TmpCoefficient[Pos] *= 0.5;
		      if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
			MaxCoefficient = Norm(TmpCoefficient[Pos]);
		      ++Pos;
		    }
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
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  if  (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      if (m1 == m2)
			TmpCoefficient[Pos] *= 0.5;
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
		}
	      else
		{
		  if (m3 == m4)
		    {
		      if  (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
			{
			  this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
			  this->M1Value[this->NbrInteractionFactors] = m1;
			  this->M2Value[this->NbrInteractionFactors] = m2;
			  this->M3Value[this->NbrInteractionFactors] = m3;
			  this->M4Value[this->NbrInteractionFactors] = m4;
			  if (m1 == m2)
			    TmpCoefficient[Pos] *= 0.5;
			  ++this->NbrInteractionFactors;
			}
		      ++Pos;
		    }
		}
	    }
    }
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
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

//ZP: The code below does not seem to work for correctly when the torus angle is non-zero

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

Complex ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
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

  N2 = (double) (m1 - m4);
  Coefficient = 1.0;
  while ((Norm(Sum) + Norm(Coefficient)) != Norm(Sum))
    {
      Q2 = this->Ratio * N2 * N2;
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
	  Q2 = this->InvRatio * N1 * N1 - 2.0 * N1 * N2 * cosine + this->Ratio * N2 * N2;
	  Precision1 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision1 * Phase(N1 * Factor);

	  Q2 = this->InvRatio * N1 * N1 + 2.0 * N1 * N2 * cosine + this->Ratio * N2 * N2;
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
      Q2 = this->Ratio * N2 * N2;
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
	  Q2 = this->InvRatio * N1 * N1 - 2.0 * N1 * N2 * cosine + this->Ratio * N2 * N2;
	  Precision1 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision1 * Phase(N1 * Factor);

	  Q2 = this->InvRatio * N1 * N1 + 2.0 * N1 * N2 * cosine + this->Ratio * N2 * N2;
	  Precision2 = this->GetVofQ(PIOnMS*Q2);
	  Coefficient += Precision2 * Phase(- N1 * Factor);

	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * this->MaxMomentum));
}


// get fourier transform of interaction
// Q2_half = one half of q² value
double ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half)
{
  double Result;
  double Q2=2.0*Q2_half;
  if ((this->HaveCoulomb)&&(Q2_half!=0.0))
    {
      //cout << "branch 1 : Ln="<<this->FormFactor.GetValue(Q2_half)<<" Ln2="<<GETSQR(this->FormFactor(Q2_half))<<", exp="<<exp(-Q2_half)<<" 1/Q="<<1.0/sqrt(Q2)<<" ";
      //this->FormFactor.PrintValue(cout, Q2_half)<<" ";
      Result=GETSQR(this->FormFactor(Q2_half)) / sqrt(Q2);
    }
  else
    Result=0.0;
  for (int i=0; i<NbrPseudopotentials; ++i)
    if (this->Pseudopotentials[i]!=0.0)
      Result += 2.0*this->Pseudopotentials[i]*this->LaguerreM[i].PolynomialEvaluate(Q2);
  //cout <<"V("<<2*Q2_half<<")="<<Result<<" LL="<<this->LandauLevel<<endl;
  return Result * exp(-Q2_half);
}


*/

// get fourier transform of interaction
// Q2_half = one half of q² value
double ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half)
{ 
  double Result;
  double Q2 = 2.0 * Q2_half;
  if ((this->HaveCoulomb) && (Q2 != 0.0))
    {
      Result = pow(this->FormFactor.PolynomialEvaluate(Q2_half), 2) * (2.0 * M_PI)/sqrt(Q2);
    }
  else
    Result = 0.0;

  //add pseudopotentials
  for (int i=0; i<NbrPseudopotentials; ++i)
     if (this->Pseudopotentials[i]!=0.0)
        Result += 2.0 * (2.0 * M_PI) * this->Pseudopotentials[i]*this->LaguerreM[i].PolynomialEvaluate(Q2);
   
  return Result * exp(-Q2_half);
}



// evaluate Wigner crystal energy per particle
//
// return value = Wigner crystal energy per particle

double ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateWignerCrystalEnergy ()
{
  double TmpRatio = M_PI * this->Ratio;
  double TmpInvRatio = M_PI * this->InvRatio;
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

double ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::MisraFunction (double n, double x)
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

double ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian::PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision)
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

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Particles->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i] << "    ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnTwistedTorusCoulombWithMagneticTranslationsHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Particles->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Particles->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][H.Particles->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

