////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of internal interaction for                    //
//                 particles on torus with coulombian interaction             //
//                                                                            //
//                        last modification : 08/11/2002                      //
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
#include "Interaction/InternalInteraction/ParticleOnTorusCoulombInternalInteraction.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::endl;


// constructor
//
// fermion = flag to indicate if particles are fermions
// nbrSite = total number of sites
// ratio = ratio between the width in the x direction and the width in the y direction
// precision = precision on interaction coefficient 

ParticleOnTorusCoulombInternalInteraction::ParticleOnTorusCoulombInternalInteraction(bool fermion, int nbrSite, double ratio, 
										     double precision)
{
  this->FermionFlag = fermion;
  this->SystemSize = nbrSite;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Precision = precision;
}

// destructor
//

ParticleOnTorusCoulombInternalInteraction::~ParticleOnTorusCoulombInternalInteraction()
{
}

// evaluate and add to hamitonian internal interaction terms
//
// hamiltonian = reference on hamiltonian matrix representation
// operators = list of operators used to construct interaction terms
// return value = reference on hamiltonian matrix representation

Matrix& ParticleOnTorusCoulombInternalInteraction::AddInteraction (Matrix& hamiltonian, List<Matrix*>& operators)
{
  int MaxMomentum = operators.GetNbrElement();
  double* TmpCoefficient = new double [MaxMomentum * (MaxMomentum + 1) * MaxMomentum];
  double MaxCoefficient = 0.0;
  int m4;
  int Pos = 0;
  if (this->FermionFlag == true)
    {
      for (int m1 = 0; m1 < MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->SystemSize;
	      else
		if (m4 >= this->SystemSize)
		  m4 -= this->SystemSize;
	      if ((m3 > m4))// && ((m1 * MaxMomentum + m2) >= (m3 * MaxMomentum + m4)))
		{
//		  cout << m1 << " " << m2 << " " << m3 << " " << m4 << endl;
		  TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
		  ++Pos;
		}
	    }
      int NbrInteractionFactors = 0;
      Pos = 0;
      MaxCoefficient *= this->Precision;
      for (int m1 = 0; m1 < MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->SystemSize;
	      else
		if (m4 >= this->SystemSize)
		  m4 -= this->SystemSize;
	      if ((m3 > m4))
		{
		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      ++NbrInteractionFactors;
//		      cout << m1 << " " << m2 << " " << m3 << " " << m4 << endl;
		      ((RealSymmetricMatrix&) hamiltonian).AddAAAtAt((RealMatrix&) (*(operators[m1])), (RealMatrix&) (*(operators[m2])),
								     (RealMatrix&) (*(operators[m3])), (RealMatrix&) (*(operators[m4])), TmpCoefficient[Pos]);
		    }
		  ++Pos;
		}
	    }
    }
  else
    {
    }
  delete[] TmpCoefficient;
  return hamiltonian;
}


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusCoulombInternalInteraction::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->SystemSize);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double CoefPrecision;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
	  CoefPrecision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;
	  CoefPrecision = 1.0;
	}
      while ((fabs(Coefficient) + CoefPrecision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  CoefPrecision = 2.0 * exp(- PIOnM * Q2) / sqrt(Q2);
	  Coefficient += CoefPrecision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->SystemSize;
    }
  N2 = (double) (m1 - m4 - this->SystemSize);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = exp(- PIOnM * Q2) / sqrt(Q2);
	  CoefPrecision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;
	  CoefPrecision = 1.0;
	}
      while ((fabs(Coefficient) + CoefPrecision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  CoefPrecision = 2.0 *  exp(- PIOnM * Q2) / sqrt(Q2);
	  Coefficient += CoefPrecision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->SystemSize;
    }
  return (Sum / (2.0 * sqrt(2.0 * M_PI * this->SystemSize)));
}
