////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of approximated hamiltonian associated to fermions          //
//                    on a torus with coulombian interaction                  //
//                                                                            //
//                        last modification : 09/09/2002                      //
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


#include "Hamiltonian/ApproximateFermionOnTorusCoulombHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// constructor from default datas
//
// fermions = Hilbert space associated to the system
// nbrFermions = number of fermions
// maxMomentum = maximum Lz value reached by a fermion in the state
// ratio = ratio between the width in the x direction and the width in the y direction

ApproximateFermionOnTorusCoulombHamiltonian::ApproximateFermionOnTorusCoulombHamiltonian(FermionOnTorus* fermions, int nbrFermions, int maxMomentum,
								     double ratio)
{
  this->Fermions = fermions;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrFermions = nbrFermions;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
  cout << "Wigner Energy = " << WignerEnergy << endl;
/*  this->InteractionCoefficient = (this->EvaluateInteractionCoefficient(1, 2, 0, 3)
				  + this->EvaluateInteractionCoefficient(2, 1, 3, 0)
				  - this->EvaluateInteractionCoefficient(1, 2, 3, 0)
				  - this->EvaluateInteractionCoefficient(2, 1, 0, 3));
  this->EvaluateDiagonalMatrixElements();*/
  this->EvaluateInteractionFactors();
  int Fast = this->FastMultiplicationMemory();
  cout << "fast = " << Fast << endl;
  if ((Fast > 0) && (Fast < (1 << 28)))
    {
      this->EnableFastMultiplication();
    }
  else
    {
      delete[] this->NbrInteractionPerComponent;
    }
}

// destructor
//

ApproximateFermionOnTorusCoulombHamiltonian::~ApproximateFermionOnTorusCoulombHamiltonian() 
{
//  delete[] this->DiagonalCoefficients;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
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

void ApproximateFermionOnTorusCoulombHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
  this->Fermions = (FermionOnTorus*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ApproximateFermionOnTorusCoulombHamiltonian::GetHilbertSpace ()
{
  return this->Fermions;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ApproximateFermionOnTorusCoulombHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Fermions->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ApproximateFermionOnTorusCoulombHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ApproximateFermionOnTorusCoulombHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Fermions->GetHilbertSpaceDimension();
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

Complex ApproximateFermionOnTorusCoulombHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Fermions->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
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

RealVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Fermions->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Fermions->GetHilbertSpaceDimension();
  double Shift = ((double) this->NbrFermions) * this->WignerEnergy;
  if (this->FastMultiplicationFlag == false)
    {
      double Coefficient;
      int Index;
      double TmpInteractionCoefficient;
      int ReducedNbrInteractionFactors3 = this->MaxMomentum - 3;
      int ReducedNbrInteractionFactors2 = this->MaxMomentum - 2;
      int ReducedNbrInteractionFactors1 = this->MaxMomentum - 1;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  TmpInteractionCoefficient = this->InteractionCoefficient * vSource[i];
	  for (int j = 0; j < ReducedNbrInteractionFactors3; ++j) 
	    {
	      Index = this->Fermions->AdAdAA(i, j + 2, j + 1, j + 3, j, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteractionCoefficient;
	    }
	  Index = this->Fermions->AdAdAA(i, ReducedNbrInteractionFactors1, ReducedNbrInteractionFactors2, 0, ReducedNbrInteractionFactors3, Coefficient);
	  if (Index < Dim)
	    vDestination[Index] += Coefficient * TmpInteractionCoefficient;
	  Index = this->Fermions->AdAdAA(i, 0, ReducedNbrInteractionFactors1, 1, ReducedNbrInteractionFactors2, Coefficient);
	  if (Index < Dim)
	    vDestination[Index] += Coefficient * TmpInteractionCoefficient;
	  Index = this->Fermions->AdAdAA(i, 1, 0, 2, ReducedNbrInteractionFactors1, Coefficient);
	  if (Index < Dim)
	    vDestination[Index] += Coefficient * TmpInteractionCoefficient;
	  vDestination[i] += this->DiagonalCoefficients[i] * vSource[i];
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
	  vDestination[i] += Shift * Coefficient;
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

ComplexVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Fermions->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Fermions->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ApproximateFermionOnTorusCoulombHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ApproximateFermionOnTorusCoulombHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ApproximateFermionOnTorusCoulombHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ApproximateFermionOnTorusCoulombHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

/*  this->MaxMomentum = 15;
  double Max = 0.0;
  int Val = 0;
  cout << (this->MaxMomentum * this->MaxMomentum) << " " << this->MaxMomentum << " " << this->MaxMomentum << endl;
  for (; Val < this->MaxMomentum; ++Val)
    for (int i = 0; i < Val; ++i)
      for (int j = 0; j < this->MaxMomentum; ++j)
	{
	  int k = Val + i - j;
	  if (k >= this->MaxMomentum)
	    k -= this->MaxMomentum;
	  else
	    if (k < 0)
	      k += this->MaxMomentum;
	  if (k < j)
	    {	
	      double Tmp = fabs(this->EvaluateInteractionCoefficient(Val, i, j, k) + this->EvaluateInteractionCoefficient(i, Val, k, j)
				- this->EvaluateInteractionCoefficient(i, Val, j, k) - this->EvaluateInteractionCoefficient(Val, i, k, j));
	      if (Tmp > Max)
		Max = Tmp;
//	  cout << i << " " << j << " " << fabs(this->EvaluateInteractionCoefficient(Val, i, j, k) + this->EvaluateInteractionCoefficient(i, Val, k, j)
//					       - this->EvaluateInteractionCoefficient(i, Val, j, k) - this->EvaluateInteractionCoefficient(Val, i, k, j)) << endl;
	    }
	}
  cout << " Max = " << Max << endl;
  cout << "----------------------------------------------------" << endl;
  Max *= 1e-2;

  for (Val = 0; Val < this->MaxMomentum; ++Val)
    {
      int Count = 0;
      for (int i = 0; i < Val; ++i)
	for (int j = 0; j < this->MaxMomentum; ++j)
	  {
	    int k = Val + i - j;
	    if (k >= this->MaxMomentum)
	    k -= this->MaxMomentum;
	    else
	      if (k < 0)
		k += this->MaxMomentum;	
	    double Tmp = fabs(this->EvaluateInteractionCoefficient(Val, i, j, k) + this->EvaluateInteractionCoefficient(i, Val, k, j)
			      - this->EvaluateInteractionCoefficient(i, Val, j, k) - this->EvaluateInteractionCoefficient(Val, i, k, j));
	  if (k < j)
	    {	
	      if (Tmp > Max)
		{
		  cout << Val << " " << i << " " << j << " " << k << " = " << Tmp << endl;
		  ++Count;
		}
	    }
	  }
      cout << "total = " << Count << endl;
      cout << "----------------------------------------------------" << endl;
    }
  exit(0);*/

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

double ApproximateFermionOnTorusCoulombHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
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
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * exp(- PIOnM * Q2) / sqrt(Q2);
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
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 *  exp(- PIOnM * Q2) / sqrt(Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * sqrt(2.0 * M_PI * this->MaxMomentum)));
}

// evaluate Wigner crystal energy per fermion
//
// return value = Wigner crystal energy per fermion

double ApproximateFermionOnTorusCoulombHamiltonian::EvaluateWignerCrystalEnergy ()
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

double ApproximateFermionOnTorusCoulombHamiltonian::MisraFunction (double n, double x)
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

double ApproximateFermionOnTorusCoulombHamiltonian::PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision)
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

int ApproximateFermionOnTorusCoulombHamiltonian::FastMultiplicationMemory()
{
  int Index;
  int ReducedMaxMomentum = this->MaxMomentum - 1;
  int Lim = 1;
  double Coefficient;
  int memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  this->NbrInteractionPerComponent = new int [this->Fermions->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
    this->NbrInteractionPerComponent[i] = 0;
  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
    {
      m1 = this->M1Value[j];
      m2 = this->M2Value[j];
      m3 = this->M3Value[j];
      m4 = this->M4Value[j];
      if ((((abs(m1 - m3) <= Lim) || (abs(m1 - m3) >= ReducedMaxMomentum)) && ((abs(m2 - m4) <= Lim) || (abs(m2 - m4) >= ReducedMaxMomentum))) ||
	  (((abs(m1 - m4) <= Lim) || (abs(m1 - m4) >= ReducedMaxMomentum)) && ((abs(m2 - m3) <= Lim) || (abs(m2 - m3) >= ReducedMaxMomentum))))
	for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
	  {
	    Index = this->Fermions->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	    if (Index < this->Fermions->GetHilbertSpaceDimension())
	      {
		++memory;
		++this->NbrInteractionPerComponent[i];
	      }
	  }    
    }
  memory = ((sizeof (int*) + sizeof (int) + 2 * sizeof(double*)) * this->Fermions->GetHilbertSpaceDimension() + 
	    memory *  (sizeof (int) + 2 * sizeof(double)));
  return memory;
}

// enable fast multiplication algorithm
//

void ApproximateFermionOnTorusCoulombHamiltonian::EnableFastMultiplication()
{
  int Index;
  double Coefficient;
  int m1;
  int m2;
  int m3;
  int m4;
  int* TmpIndexArray;
  int ReducedMaxMomentum = this->MaxMomentum - 1;
  int Lim = 1;
  double* TmpCoefficientArray;
  int Pos;
  this->InteractionPerComponentIndex = new int* [this->Fermions->GetHilbertSpaceDimension()];
  this->InteractionPerComponentCoefficient = new double* [this->Fermions->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
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
	  if ((((abs(m1 - m3) <= Lim) || (abs(m1 - m3) == ReducedMaxMomentum)) && ((abs(m2 - m4) <= Lim) || (abs(m2 - m4) >= ReducedMaxMomentum))) ||
	      (((abs(m1 - m4) <= Lim) || (abs(m1 - m4) == ReducedMaxMomentum)) && ((abs(m2 - m3) <= Lim) || (abs(m2 - m3) >= ReducedMaxMomentum))))
	    {
	      Index = this->Fermions->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < this->Fermions->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Pos] = Index;
		  TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
		  ++Pos;
		}
	    }
	}
    }
  this->FastMultiplicationFlag = true;
}

// evaluate all matrix diagonal elements
//   

void ApproximateFermionOnTorusCoulombHamiltonian::EvaluateDiagonalMatrixElements()
{
/*  double Shift = ((double) this->NbrFermions) * this->WignerEnergy;
  int ReducedMaxMomentum = this->MaxMomentum - 1;  
  double* TmpCoef = new double [ReducedMaxMomentum];
  for (int i = 1; i < this->MaxMomentum; ++i)
    {
      TmpCoef[i - 1] = (this->EvaluateInteractionCoefficient(0, i, i, 0)
		    + this->EvaluateInteractionCoefficient(i, 0, 0, i)
		    - this->EvaluateInteractionCoefficient(0, i, 0, i)
		    - this->EvaluateInteractionCoefficient(i, 0, i, 0));
      cout << TmpCoef[i - 1] << endl;
    }
  double Coefficient;
  double Coefficient2;
  this->DiagonalCoefficients = new double [this->Fermions->GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
    {
      Coefficient = Shift;
      for (int k = 0; k < ReducedMaxMomentum; ++k)
	{
	  Coefficient2 = 0;
	  int j = 0;
	  for (; j < (ReducedMaxMomentum - k); ++j)
	    Coefficient2 += TmpCoef[j] * this->Fermions->AdA(i, k + j + 1);
	  for (; j < ReducedMaxMomentum; ++j)
	    Coefficient2 += TmpCoef[j] * this->Fermions->AdA(i, j - ReducedMaxMomentum + k);	  
	  Coefficient += Coefficient2 * this->Fermions->AdA(i, k);	  
	}
      this->DiagonalCoefficients[i] = Coefficient;
    }*/
/*  double Coef1 = (this->EvaluateInteractionCoefficient(0, 1, 1, 0)
		  + this->EvaluateInteractionCoefficient(1, 0, 0, 1)
		  - this->EvaluateInteractionCoefficient(0, 1, 0, 1)
		  - this->EvaluateInteractionCoefficient(1, 0, 1, 0));
  double Coef2 = (this->EvaluateInteractionCoefficient(0, 2, 2, 0)
		  + this->EvaluateInteractionCoefficient(2, 0, 0, 2)
		  - this->EvaluateInteractionCoefficient(0, 2, 0, 2)
		  - this->EvaluateInteractionCoefficient(2, 0, 2, 0));
  double Coef3 = (this->EvaluateInteractionCoefficient(0, 3, 3, 0)
		  + this->EvaluateInteractionCoefficient(3, 0, 0, 3)
		  - this->EvaluateInteractionCoefficient(0, 3, 0, 3)
		  - this->EvaluateInteractionCoefficient(3, 0, 3, 0));
  cout << Coef1 << " " << Coef2 << " " << Coef3 << endl;
  this->DiagonalCoefficients = new double [this->Fermions->GetHilbertSpaceDimension()];
  double Coefficient;
  int ReducedMaxMomentum = this->MaxMomentum - 3;
  for (int i = 0; i < this->Fermions->GetHilbertSpaceDimension(); ++i)
    {
      Coefficient = Shift;
      for (int k = 0; k < ReducedMaxMomentum; ++k)
	{
	  Coefficient += (Coef1 * this->Fermions->AdA(i, k + 1) + Coef2 *  this->Fermions->AdA(i, k + 2) +
			  Coef3 *  this->Fermions->AdA(i, k + 3)) * this->Fermions->AdA(i, k);
	}
      Coefficient += (Coef1 * this->Fermions->AdA(i, ReducedMaxMomentum + 1) + Coef2 *  this->Fermions->AdA(i, ReducedMaxMomentum + 2) +
		      Coef3 *  this->Fermions->AdA(i, 0)) * this->Fermions->AdA(i, ReducedMaxMomentum);
      Coefficient += (Coef1 * this->Fermions->AdA(i, ReducedMaxMomentum + 2) + Coef2 *  this->Fermions->AdA(i, 0) +
		      Coef3 *  this->Fermions->AdA(i, 1)) * this->Fermions->AdA(i, ReducedMaxMomentum + 1);
      Coefficient += (Coef1 * this->Fermions->AdA(i, 0) + Coef2 *  this->Fermions->AdA(i, 1) +
		      Coef3 *  this->Fermions->AdA(i, 2)) * this->Fermions->AdA(i, ReducedMaxMomentum + 2);
      cout << Coefficient << " ";
      this->DiagonalCoefficients[i] = Coefficient;
    }*/
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ApproximateFermionOnTorusCoulombHamiltonian& H) 
{
  RealVector TmpV2 (H.Fermions->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Fermions->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Fermions->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Fermions->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Fermions->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Fermions->GetHilbertSpaceDimension(); j++)
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

MathematicaOutput& operator << (MathematicaOutput& Str, ApproximateFermionOnTorusCoulombHamiltonian& H) 
{
  RealVector TmpV2 (H.Fermions->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Fermions->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Fermions->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Fermions->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Fermions->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Fermions->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Fermions->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Fermions->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Fermions->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Fermions->GetHilbertSpaceDimension() - 1][H.Fermions->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

