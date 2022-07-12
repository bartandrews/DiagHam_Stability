////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                   coulombian interaction implemented for DMRG              //
//                                                                            //
//                        last modification : 07/11/2002                      //
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
#include "Hamiltonian/ParticleOnTorusCoulombDMRGReadyHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Matrix/RealMatrix.h"

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

ParticleOnTorusCoulombDMRGReadyHamiltonian::ParticleOnTorusCoulombDMRGReadyHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum,
										       double ratio)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
  cout << "Wigner Energy = " << WignerEnergy << endl;
  this->EvaluateInteractionFactors();
  cout << "fast = " << this->FastMultiplicationMemory() << endl;
  this->EnableFastMultiplication();
}

// destructor
//

ParticleOnTorusCoulombDMRGReadyHamiltonian::~ParticleOnTorusCoulombDMRGReadyHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
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
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusCoulombDMRGReadyHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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
  this->Particles = (ParticleOnTorus*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusCoulombDMRGReadyHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusCoulombDMRGReadyHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusCoulombDMRGReadyHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnTorusCoulombDMRGReadyHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnTorusCoulombDMRGReadyHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent) 
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    vDestination[i] = 0.0;
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Shift = ((double) this->NbrParticles) * this->WignerEnergy;
  if (this->FastMultiplicationFlag == false)
    {
      double Coefficient;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = this->M4Value[j];
	  TmpInteraction = this->InteractionFactors[j];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	    }
	}
      m1 = this->M1Value[ReducedNbrInteractionFactors];
      m2 = this->M2Value[ReducedNbrInteractionFactors];
      m3 = this->M3Value[ReducedNbrInteractionFactors];
      m4 = this->M4Value[ReducedNbrInteractionFactors];
      TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index < Dim)
	    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
//	  vDestination[i] += Shift * vSource[i];
	}
    }
  else
    {
      double Coefficient;
      int* TmpIndexArray;
      double* TmpCoefficientArray; 
      int j;
      int TmpNbrInteraction;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	  TmpIndexArray = this->InteractionPerComponentIndex[i];
	  TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	  Coefficient = vSource[i];
	  for (j = 0; j < TmpNbrInteraction; ++j)
	    vDestination[TmpIndexArray[j]] += TmpCoefficientArray[j] * Coefficient;
//	  vDestination[i] += Shift * Coefficient;
	}
   }
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = 0.0;
      vDestination.Im(i) = 0.0;
    }
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusCoulombDMRGReadyHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnTorusCoulombDMRGReadyHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      RealMatrix* TmpMatrix = new RealMatrix(this->Particles->GetHilbertSpaceDimension(), this->Particles->GetHilbertSpaceDimension(), true);
      this->Particles->Ad(i, *TmpMatrix);
//      cout << "i = " << i << endl;
//       cout << *TmpMatrix << endl;
     TmpList += TmpMatrix;
    }
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnTorusCoulombDMRGReadyHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      RealMatrix* TmpMatrix = new RealMatrix(this->Particles->GetHilbertSpaceDimension(), this->Particles->GetHilbertSpaceDimension(), true);
      this->Particles->Ad(i, *TmpMatrix);
      TmpList += TmpMatrix;
    }
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnTorusCoulombDMRGReadyHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
/*  this->MaxMomentum = 50;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->Ratio = 25.0 / 4.0;
  this->InvRatio = 1.0 / this->Ratio;*/

  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

/*  cout << this->EvaluateInteractionCoefficient(2, 0, 2, 0) << endl;
//  this->InvRatio = this->Ratio;
//  this->Ratio = 1.0 / this->InvRatio;
  cout << this->EvaluateInteractionCoefficient(0, 2, 2, 0) << endl;
  exit(0);*/

  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
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
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
		  ++Pos;
		}
	    }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= 1e-6;//MACHINE_PRECISION;
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
		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      /*		cout << this->M1Value[this->NbrInteractionFactors] << " " << this->M2Value[this->NbrInteractionFactors] 
					<< " " << this->M3Value[this->NbrInteractionFactors] 
					<< " " << this->InteractionFactors[this->NbrInteractionFactors] << endl;*/
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
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
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
		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
	    }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
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
	      if (m3 >= m4)
		{
		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      /*		cout << this->M1Value[this->NbrInteractionFactors] << " " << this->M2Value[this->NbrInteractionFactors] 
					<< " " << this->M3Value[this->NbrInteractionFactors] 
					<< " " << this->InteractionFactors[this->NbrInteractionFactors] << endl;*/
/*		      if ((m1 == 1) && (m2 == 0) && (m3 == 3) && (m4 == 2))
			{
			  cout << "found " << TmpCoefficient[Pos] << endl;
			}*/
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
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

double ParticleOnTorusCoulombDMRGReadyHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
//  cout << "new coef====================================" << m1 << " "  << m2 << " "  << m3 << " "  << m4 << endl;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;
	  Precision = 1.0;
	}
//      cout << Coefficient << endl;
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * exp(- PIOnM * Q2) / sqrt(Q2);
/*	  if ((Coefficient * (Precision * cos (N1 * Factor)) < 0) && (fabs((Coefficient + (Precision * cos (N1 * Factor)))) < 1e-3 * (fabs(Coefficient))))
	    cout << " possible lost of precision at " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " " << Coefficient << " " << (Precision * cos (N1 * Factor))
		 << " " << (Coefficient + (Precision * cos (N1 * Factor))) << endl;*/
	  Coefficient += Precision * cos (N1 * Factor);
//	  cout << Coefficient << " " << (Precision * cos (N1 * Factor)) << endl;;
//	  Precision = fabs(Precision);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
//      cout << "partial sum = " << Sum << " " << Coefficient << endl;
    }
  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;
	  Precision = 1.0;
	}
//      cout << Coefficient << endl;
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 *  exp(- PIOnM * Q2) / sqrt(Q2);
	  Coefficient += Precision * cos (N1 * Factor);
//	  cout << Coefficient << " " << (Precision * cos (N1 * Factor)) << endl;;
//	  Precision = fabs(Precision);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
//      cout << "partial sum = " << Sum << " " << Coefficient << endl;
    }
  return (Sum / (2.0 * sqrt(2.0 * M_PI * this->MaxMomentum)));
/*  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
      else
	Coefficient = 0.0;
      Precision = Coefficient;
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * exp(- PIOnM * Q2) / sqrt(Q2);
	  Coefficient += Precision * cos (N1 * Factor);
//	  Precision = fabs(Precision);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
    }
  N2 = (double) (m3 - m4 - this->MaxMomentum);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
      else
	Coefficient = 0.0;
      Precision = Coefficient;
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 *  exp(- PIOnM * Q2) / sqrt(Q2);
	  Coefficient += Precision * cos (N1 * Factor);
//	  Precision = fabs(Precision);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
//  Sum += (drand48() - 0.5) * Sum * 1e-5;
  return -(Sum / (2.0 * sqrt(2.0 * M_PI * this->MaxMomentum)));*/
}

// evaluate Wigner crystal energy per particle
//
// return value = Wigner crystal energy per particle

double ParticleOnTorusCoulombDMRGReadyHamiltonian::EvaluateWignerCrystalEnergy ()
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

double ParticleOnTorusCoulombDMRGReadyHamiltonian::MisraFunction (double n, double x)
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
//      cout << " PreviousSum = " << PreviousSum << "   NewSum = " << NewSum << endl;
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

double ParticleOnTorusCoulombDMRGReadyHamiltonian::PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision)
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

// test the amount of memory needed for fast multiplication algorithm
//
// return value = amount of memory needed

int ParticleOnTorusCoulombDMRGReadyHamiltonian::FastMultiplicationMemory()
{
  int Index;
  double Coefficient;
  int memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  this->NbrInteractionPerComponent = new int [this->Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    this->NbrInteractionPerComponent[i] = 0;
  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
    {
      m1 = this->M1Value[j];
      m2 = this->M2Value[j];
      m3 = this->M3Value[j];
      m4 = this->M4Value[j];
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      ++memory;
	      ++this->NbrInteractionPerComponent[i];
	    }
	}   
    }
  memory = ((sizeof (int*) + sizeof (int) + 2 * sizeof(double*)) * this->Particles->GetHilbertSpaceDimension() + 
	    memory *  (sizeof (int) + 2 * sizeof(double)));
  return memory;
}

// enable fast multiplication algorithm
//

void ParticleOnTorusCoulombDMRGReadyHamiltonian::EnableFastMultiplication()
{
  int Index;
  double Coefficient;
  int m1;
  int m2;
  int m3;
  int m4;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  this->InteractionPerComponentIndex = new int* [this->Particles->GetHilbertSpaceDimension()];
  this->InteractionPerComponentCoefficient = new double* [this->Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];      
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
      Pos = 0;
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = this->M4Value[j];
	  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
//	  cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Index << endl;
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
	      ++Pos;
	    }
	}
    }
  this->FastMultiplicationFlag = true;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnTorusCoulombDMRGReadyHamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnTorusCoulombDMRGReadyHamiltonian& H) 
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

