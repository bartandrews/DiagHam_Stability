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

#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// default constructor
//

ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian()
{
    
  this->Particles = 0;
  this->LzMax = 0;
  this->NbrLzValue = 0;
  this->NbrParticles = 0;
  this->MomentumModulo = 0;
  this->Ratio = 0.0;
  this->InvRatio = 0.0;
  this->MaxMomentum = 0;
  this->XMomentum = 0;
  
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;//true;
  this->Ratio = 0;  
  this->InvRatio = 0;
  this->SpinFluxUp = 0;
  this->SpinFluxDown = 0;
  this->HamiltonianShift = 0.0;
  this->Architecture = 0;
  this->PrecalculationShift = 0;  

  this->NbrPseudopotentialsUpUp = 0;
  this->PseudopotentialsUpUp = 0;
  this->NbrPseudopotentialsDownDown = 0;
  this->PseudopotentialsDownDown = 0;
  this->NbrPseudopotentialsUpDown = 0;
  this->PseudopotentialsUpDown = 0;
  
  this->MaxNbrPseudopotentials = 0;
  this->LaguerrePolynomials = 0;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;

}

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
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  this->LayerSeparation=layerSeparation;
//   double WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
  double WignerEnergy = 0.0;
  this->SpinFluxUp = 0.0;
  this->SpinFluxDown = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  cout << "Wigner Energy = " << WignerEnergy << endl;  
  this->PseudopotentialsUpUp = 0;
  this->PseudopotentialsDownDown = 0;
  this->PseudopotentialsUpDown = 0;
  this->LaguerrePolynomials = 0;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;

  this->EvaluateExponentialFactors();
  this->EvaluateInteractionFactors();

  this->HamiltonianShift = ((double) this->NbrParticles) * WignerEnergy;

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
  if (this->PseudopotentialsUpUp != 0)
    delete[] this->PseudopotentialsUpUp;
  if (this->PseudopotentialsDownDown != 0)
    delete[] this->PseudopotentialsDownDown;
  if (this->PseudopotentialsUpDown != 0)
    delete[] this->PseudopotentialsUpDown;
  if (this->LaguerrePolynomials != 0)
    delete[] this->LaguerrePolynomials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift += shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();

//   unsigned L16Mask = (1u<<16)-1;
//   int Pos = 0;
//   int M12Index = 0;
//   int m4;
//   double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
//   double MaxCoefficient = 0.0;  

  if (this->Particles->GetParticleStatistic() == ParticleOnTorusWithSpinAndMagneticTranslations::FermionicStatistic)
    {
      // upup-upup 
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupup[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupup[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // downdown-downdown
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsdowndown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsdowndown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // updown-updown
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->LayerSeparation)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->LayerSeparation));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->LayerSeparation)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->LayerSeparation));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupdown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupdown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}




//       this->NbrM12IntraIndices = ((this->NbrLzValue-1) * (this->NbrLzValue - 2)) / 2;
//       this->M12IntraValue = new unsigned [this->NbrM12IntraIndices];
//       this->NbrM34IntraValues = new int [this->NbrM12IntraIndices];
//       this->M34IntraValues = new unsigned* [this->NbrM12IntraIndices];
//       for (int i=0; i<this->NbrM12IntraIndices; ++i) this->M34IntraValues[i] = new unsigned[this->NbrLzValue];
//       for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
// 	for (int m2 = 0; m2 < m1; ++m2)
// 	  {
// 	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
// 	      {
// 		m4 = m1 + m2 - m3;
// 		if (m4 < 0)
// 		  m4 += this->MaxMomentum;
// 		else
// 		  if (m4 >= this->MaxMomentum)
// 		    m4 -= this->MaxMomentum;
// 		if (m3 > m4)
// 		  {
// 		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
// 					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
// 					   - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
// 					   - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
// 		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
// 		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
// 		    ++Pos;
// 		  }
// 	      }
// 	  }
//       cout << "Max Nbr InteractionUpUp = " << Pos << endl;            
//       MaxCoefficient *= MACHINE_PRECISION;
//       InteractionFactorsUpUp = new double[Pos];
//       InteractionFactorsDownDown = new double[Pos];
//       M12Index = 0;
//       Pos = 0;
//       int TmpNbrInteractionFactors = 0;
//       for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
// 	for (int m2 = 0; m2 < m1; ++m2)
// 	  {
// 	    this->M12IntraValue[M12Index] = (m1&L16Mask)|((m2&L16Mask)<<16);
// 	    this->NbrM34IntraValues[M12Index]=0;
// 	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
// 	      {		
// 		m4 = m1 + m2 - m3;
// 		if (m4 < 0)
// 		  m4 += this->MaxMomentum;
// 		else
// 		  if (m4 >= this->MaxMomentum)
// 		    m4 -= this->MaxMomentum;
// 		if (m3 > m4)
// 		{
// 		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
// 		    {
// 		      this->InteractionFactorsUpUp[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
// 		      this->InteractionFactorsDownDown[TmpNbrInteractionFactors] = TmpCoefficient[Pos];
// 		      this->M34IntraValues[M12Index][this->NbrM34IntraValues[M12Index]]
// 			= (m3&L16Mask)|((m4&L16Mask)<<16);
// 		      ++TmpNbrInteractionFactors;
// 		      ++this->NbrM34IntraValues[M12Index];
// 		    }
// 		  ++Pos;
// 		}
// 	    }
// 	    ++M12Index;
// 	  }
//       cout << "Actual Nbr InteractionUpUp = " << TmpNbrInteractionFactors << endl;
//       // matrix elements for different spin
//       this->NbrM12InterIndices = (this->NbrLzValue-1) * (this->NbrLzValue-1);
//       this->M12InterValue = new unsigned [this->NbrM12InterIndices];
//       this->NbrM34InterValues = new int [this->NbrM12InterIndices];
//       this->M34InterValues = new unsigned*[this->NbrM12InterIndices];      
//       for (int i=0; i<this->NbrM12InterIndices; ++i) this->M34InterValues[i] = new unsigned[this->NbrLzValue];
//       Pos = 0;
//       M12Index = 0;
//       for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
// 	for (int m2 = 0; m2 < this->MaxMomentum; ++m2)
// 	  {
// 	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
// 	      {
// 		m4 = m1 + m2 - m3;
// 		if (m4 < 0)
// 		  m4 += this->MaxMomentum;
// 		else
// 		  if (m4 >= this->MaxMomentum)
// 		    m4 -= this->MaxMomentum;
// 		if ((m1!=m2)||(m3!=m4))
// 		  {
// 		    TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->LayerSeparation)
// 		      + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->LayerSeparation);
// 		  }
// 		else
// 		  {
// 		    TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->LayerSeparation);
// 		  }
// 		if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
// 		  MaxCoefficient = fabs(TmpCoefficient[Pos]);
// 		++Pos;
// 	      }
// 	  }
//       cout << "Max Nbr InteractionUpDown = " << Pos << endl;            
//       MaxCoefficient *= MACHINE_PRECISION;
//       InteractionFactorsUpDown = new double[Pos];
//       M12Index = 0;
//       Pos = 0;
//       TmpNbrInteractionFactors = 0;
//       for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
// 	for (int m2 = 0; m2 < this->MaxMomentum; ++m2)
// 	  {
// 	    this->M12InterValue[M12Index] = (m1&L16Mask)|((m2&L16Mask)<<16);
// 	    this->NbrM34InterValues[M12Index]=0;
// 	    for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
// 	      {		
// 		m4 = m1 + m2 - m3;
// 		if (m4 < 0)
// 		  m4 += this->MaxMomentum;
// 		else
// 		  if (m4 >= this->MaxMomentum)
// 		    m4 -= this->MaxMomentum;
// 		if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
// 		  {
// 		    // swap 3,4, introduce additional minus sign
// 		    this->InteractionFactorsUpDown[TmpNbrInteractionFactors] = -1.0*TmpCoefficient[Pos];
// 		    this->M34InterValues[M12Index][this->NbrM34InterValues[M12Index]]
// 		      = (m4&L16Mask)|((m3&L16Mask)<<16); // rather than (m3&L16Mask)|((m4&L16Mask)<<16);
// 		    ++TmpNbrInteractionFactors;
// 		    ++this->NbrM34InterValues[M12Index];
// 		  }
// 		++Pos;
// 	      }
// 	    ++M12Index;
// 	  }
//       cout << "Actual Nbr InteractionUpDown = " << TmpNbrInteractionFactors << endl;
//       // no one-body interactions:
//       this->OneBodyInteractionFactorsUpUp = 0;
//       this->OneBodyInteractionFactorsDownDown = 0;
    }
  else
    {
      cout << "Bosonic statistics not defined yet for ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian"<<endl;
      exit(1);
    }
  cout << "====================================" << endl;
//  delete[] TmpCoefficient;
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
//   if ((m1==m2)&&(m1==m3)&&(m1==m4)) return 0.5;
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;// this->GetVofQ(PIOnM*Q2, layerSeparation); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
    }
  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = 1.0;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient =  this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0; // this->GetVofQ(PIOnM*Q2, layerSeparation); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * this->MaxMomentum));
}


// get fourier transform of interaction
// Q2_half = one half of q² value
// layerSeparation = layer separation

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsHamiltonian::GetVofQ(double Q2_half, double layerSeparation)
{
  double Q=sqrt(2.0*Q2_half);
  return exp(-Q2_half-Q*layerSeparation)/Q;
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

