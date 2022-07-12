////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Gunnar Möller                         //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                  coulomb interaction and magnetic translations             //
//                                                                            //
//                        last modification : 02/10/2003                      //
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


#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

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
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// layerSeparation = layer separation in units of magnetic length
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian
(ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
 double layerSeparation, AbstractArchitecture* architecture, int memory, char* precalculationFileName)
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
  this->LayerSeparation=layerSeparation;
//  double WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
  double WignerEnergy = 0.0;
  this->Architecture = architecture;
  cout << "Wigner Energy = " << WignerEnergy << endl;  
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;
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
	  int TmpMemory = this->FastMultiplicationMemory(memory);
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
	  cout << endl;
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

ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::~ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian() 
{
  delete[] M12IntraValue;
  for (int i=0; i<NbrM12IntraIndices; ++i)
    delete[] M34IntraValues[i];
  delete[] M34IntraValues;
  delete[] NbrM34IntraValues;

  delete[] M12InterValue;
  for (int i=0; i<NbrM12InterIndices; ++i)
    delete[] M34InterValues[i];
  delete[] M34InterValues;
  delete[] NbrM34InterValues;

  delete[] InteractionFactorsUpUp;
  delete[] InteractionFactorsDownDown;
  delete[] InteractionFactorsUpDown;

  delete[] OneBodyInteractionFactorsUpUp;
  delete[] OneBodyInteractionFactorsDownDown;  
  
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
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] InteractionFactorsUpUp;
  delete[] InteractionFactorsDownDown;
  delete[] InteractionFactorsUpDown;

  delete[] OneBodyInteractionFactorsUpUp;
  delete[] OneBodyInteractionFactorsDownDown;  

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
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  unsigned L16Mask = (1u<<16)-1;
  int Pos = 0;
  int M12Index = 0;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;  
  if (this->Particles->GetParticleStatistic() == ParticleOnTorusWithSpinAndMagneticTranslations::FermionicStatistic)
    {
      this->NbrM12IntraIndices = ((this->NbrLzValue-1) * (this->NbrLzValue - 2)) / 2;
      this->M12IntraValue = new unsigned [this->NbrM12IntraIndices];
      this->NbrM34IntraValues = new int [this->NbrM12IntraIndices];
      this->M34IntraValues = new unsigned* [this->NbrM12IntraIndices];
      for (int i=0; i<this->NbrM12IntraIndices; ++i) this->M34IntraValues[i] = new unsigned[this->NbrLzValue];
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
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
		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->LayerSeparation)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->LayerSeparation)
					   - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->LayerSeparation)
					   - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->LayerSeparation));
		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
	      }
	  }
      cout << "Max Nbr InteractionUpUp = " << Pos << endl;            
      MaxCoefficient *= MACHINE_PRECISION;
      InteractionFactorsUpUp = new double[Pos];
      InteractionFactorsDownDown = new double[Pos];
      M12Index = 0;
      Pos = 0;
      int TmpNbrInteractionFactors = 0;
      double Factor = - 4.0;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
	    this->M12IntraValue[M12Index] = (m1&L16Mask)|((m2&L16Mask)<<16);
	    this->NbrM34IntraValues[M12Index]=0;
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
		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactorsUpUp[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
		      this->InteractionFactorsDownDown[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M34IntraValues[M12Index][this->NbrM34IntraValues[M12Index]]
			= (m3&L16Mask)|((m4&L16Mask)<<16);
		      ++TmpNbrInteractionFactors;
		      ++this->NbrM34IntraValues[M12Index];
		    }
		  ++Pos;
		}
	    }
	    ++M12Index;
	  }
      cout << "Actual Nbr InteractionUpUp = " << TmpNbrInteractionFactors << endl;
      // matrix elements for different spin
      this->NbrM12InterIndices = (this->NbrLzValue-1) * (this->NbrLzValue-1);
      this->M12InterValue = new unsigned [this->NbrM12InterIndices];
      this->NbrM34InterValues = new int [this->NbrM12InterIndices];
      this->M34InterValues = new unsigned*[this->NbrM12InterIndices];      
      for (int i=0; i<this->NbrM12InterIndices; ++i) this->M34InterValues[i] = new unsigned[this->NbrLzValue];
      Pos = 0;
      M12Index = 0;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < this->MaxMomentum; ++m2)
	  {
	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	      {
		m4 = m1 + m2 - m3;
		if (m4 < 0)
		  m4 += this->MaxMomentum;
		else
		  if (m4 >= this->MaxMomentum)
		    m4 -= this->MaxMomentum;		
		TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->LayerSeparation);
		if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		  MaxCoefficient = fabs(TmpCoefficient[Pos]);
		++Pos;
	      }
	  }
      cout << "Max Nbr InteractionUpDown = " << Pos << endl;            
      MaxCoefficient *= MACHINE_PRECISION;
      InteractionFactorsUpDown = new double[Pos];
      M12Index = 0;
      Pos = 0;
      TmpNbrInteractionFactors = 0;
      Factor = - 4.0;  // check factor here... xxx 
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < this->MaxMomentum; ++m2)
	  {
	    this->M12InterValue[M12Index] = (m1&L16Mask)|((m2&L16Mask)<<16);
	    this->NbrM34InterValues[M12Index]=0;
	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	      {		
		m4 = m1 + m2 - m3;
		if (m4 < 0)
		  m4 += this->MaxMomentum;
		else
		  if (m4 >= this->MaxMomentum)
		    m4 -= this->MaxMomentum;
		if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  {
		    this->InteractionFactorsUpDown[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
		    this->M34InterValues[M12Index][this->NbrM34InterValues[M12Index]]
		      = (m3&L16Mask)|((m4&L16Mask)<<16);
		    ++TmpNbrInteractionFactors;
		    ++this->NbrM34InterValues[M12Index];
		  }
		++Pos;
	      }
	    ++M12Index;
	  }
      cout << "Actual Nbr InteractionUpDown = " << TmpNbrInteractionFactors << endl;
      // no one-body interactions:
      this->OneBodyInteractionFactorsUpUp = 0;
      this->OneBodyInteractionFactorsDownDown = 0;
    }
  else
    {
    }
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// layerSeparation = separation of layers
// return value = numerical coefficient

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, double layerSeparation)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q,Q2;
  double D2 = layerSeparation * layerSeparation;
  double Precision;
//  cout << "new coef====================================" << m1 << " "  << m2 << " "  << m3 << " "  << m4 << endl;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2 ) / sqrt(Q2+D2);
	  Precision = Coefficient;
	}
      else
	{
	  if (layerSeparation != 0.0)
	    {
	      Coefficient = 1.0/layerSeparation;
	      Precision = Coefficient;
	    }
	  else
	    {
	      Coefficient = 0.0; // could be infinite, here!
	      Precision = 1.0;
	    }
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Q=sqrt(Q2);
	  Precision = 2.0 * exp(- PIOnM * Q2) / sqrt(Q2 + D2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
    }
  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2 + D2);
	  Precision = Coefficient;
	}
      else
	{
	  if (layerSeparation != 0.0)
	    {
	      Coefficient = 1.0/layerSeparation;
	      Precision = Coefficient;
	    }
	  else
	    {
	      Coefficient = 0.0; // could be infinite, here!
	      Precision = 1.0;
	    }
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 *  exp(- PIOnM * Q2) / sqrt(Q2 + D2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * sqrt(2.0 * M_PI * this->MaxMomentum)));
}

// evaluate Wigner crystal energy per particle
//
// return value = Wigner crystal energy per particle

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::EvaluateWignerCrystalEnergy ()
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

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::MisraFunction (double n, double x)
{
  int NbrSubdivision = 100000;
  double PreviousSum = this->PartialMisraFunction(n, x, 0.0, 1.0, NbrSubdivision);
  double NewSum = PreviousSum;
  PreviousSum *= 2.0;
  while ((fabs(PreviousSum - NewSum) / PreviousSum) > MACHINE_PRECISION)
    {
      PreviousSum = NewSum;
      NbrSubdivision += 10000;
      NewSum = this->PartialMisraFunction(n, x, 0.0, 1.0, NbrSubdivision);
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

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision)
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

ostream& operator << (ostream& Str, ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian& H) 
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

