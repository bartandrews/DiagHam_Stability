////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                               delta interaction                            //
//                                                                            //
//                        last modification : 23/06/2003                      //
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


#include "Hamiltonian/ParticleOnTorusDeltaHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Polynomial/Polynomial.h"
#include "Polynomial/SpecialPolynomial.h"

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
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusDeltaHamiltonian::ParticleOnTorusDeltaHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum,
        double ratio, double angle, int level, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
    this->Particles = particles;
    this->LzMax = maxMomentum - 1;
    this->NbrLzValue = this->LzMax + 1;
    this->NbrParticles = nbrParticles;
    this->FastMultiplicationFlag = false;
    this->Ratio = ratio;
    this->InvRatio = 1.0 / ratio;
    this->Theta = angle;
    this->CosTheta = cos(angle);
    this->SinTheta = sin(angle);
    this->Lx = maxMomentum;
    this->Ly = maxMomentum / ratio;
    this->Architecture = architecture;
    this->MixingLevel = level;
    this->divisor = 1;
    for(int j = 1; j <= level; j++)
    {
        this->divisor *= (1/sqrt(2*j));
    }

    this->LPolynomial = LaguerrePolynomial(MixingLevel);
    long MinIndex;
    long MaxIndex;
    this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
    this->PrecalculationShift = (int) MinIndex;

    this->EvaluateInteractionFactors();
    this->HamiltonianShift = 0.0;



    if (precalculationFileName == 0)
    {
        if (memory > 0)
        {
            long TmpMemory = this->FastMultiplicationMemory(memory);
            if (TmpMemory < 1024)
                cout  << "fast = " <<  TmpMemory << "b ";
            else if (TmpMemory < (1 << 20))
                cout  << "fast = " << (TmpMemory >> 10) << "kb ";
            else if (TmpMemory < (1 << 30))
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

ParticleOnTorusDeltaHamiltonian::~ParticleOnTorusDeltaHamiltonian()
{
}

// This is the original DiagHam function.
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

// double ParticleOnTorusDeltaHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
// {
//     double Coefficient = 1.0;
//     double PIOnM = M_PI / ((double) this->NbrLzValue);
//     double Factor =  - ((double) (m1-m4)) * PIOnM * 2.0;
//     double Sum = 0.0;
//     double N2 = (double) (m1 - m3);
//     double N1;
//     double Q2;
//     double Precision;
// //  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " : ";
//     while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
//     {
//         N1 = 1.0;
//         Q2 = this->Ratio * N2 * N2;
//         if (N2 != 0.0)
//         {
//             Coefficient = exp(- PIOnM * Q2);
//             Precision = Coefficient;
//         }
//         else
//         {
//             Coefficient = 1.0;
//             Precision = 1.0;
//         }
//         while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
//         {
//             Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
//             Precision = 2.0 * exp(- PIOnM * Q2);
//             Coefficient += Precision * cos (N1 * Factor);
//             N1 += 1.0;
//         }
//         Sum += Coefficient;
//         N2 += this->NbrLzValue;
//     }
//     N2 = (double) (m1 - m3 - this->NbrLzValue);
//     Coefficient = Sum;
//     while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
//     {
//         N1 = 1.0;
//         Q2 = this->Ratio * N2 * N2;
//         if (N2 != 0.0)
//         {
//             Coefficient = exp(- PIOnM * Q2);
//             Precision = Coefficient;
//         }
//         else
//         {
//             Coefficient = 1.0;
//             Precision = 1.0;
//         }
//         while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
//         {
//             Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
//             Precision = 2.0 *  exp(- PIOnM * Q2);
//             Coefficient += Precision * cos (N1 * Factor);
//             N1 += 1.0;
//         }
//         Sum += Coefficient;
//         N2 -= this->NbrLzValue;
//     }
// //  cout << Sum << endl;
//     return (Sum / (4.0 * M_PI * this->NbrLzValue));
// }


double ParticleOnTorusDeltaHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
    return EvaluateInteractionCoefficientWithMixing(m1,m2,m3,m4);

}

double ParticleOnTorusDeltaHamiltonian::EvaluateInteractionCoefficientWithMixing(int m1, int m2, int m3, int m4)
{
    double Coefficient = 1.0;
    double PIOnM = M_PI / ((double) this->NbrLzValue);
    double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
    double c0 = 1.0;
    double c4 = 0.0;
    //double Cos4 = this->CosTheta*this->CosTheta*this->CosTheta*this->CosTheta;
    //double Sin4 = this->SinTheta*this->SinTheta*this->SinTheta*this->SinTheta;
    double Multiplier = 0.0;
    double Ly2 = this->Ly*this->Ly;
    double Lx2 = this->Lx*this->Lx;
    double Ell2 = this->Lx*this->Ly / (2 * M_PI * this->NbrLzValue);
    double Sin2Cos2 = this->SinTheta*this->SinTheta*this->CosTheta*this->CosTheta;
    double NormalizationFactor = 1;// / (1 - 2*Sin2Cos2);
    Complex ZQ,ZQbar = 1;

    double Sum = 0.0;
    double N2 = (double) (m1 - m4);
    double N1;
    double Q2;
    double q2l2by2;
    double Precision;
//  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " : ";
    while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
        N1 = 1.0;
        Q2 = this->Ratio * N2 * N2;
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

            Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
            q2l2by2 = (N2*N2/Ly2 + N1*N1/Lx2)*Ell2*M_PI*M_PI*2;
            ZQ = 1;
            ZQbar = 1;
            for(int i = 1; i <= this->MixingLevel; i++)
            {
                ZQ *= 2*sqrt(Ell2)*M_PI*(I()*N1/this->Lx + N2/this->Ly);
                ZQbar *= 2*sqrt(Ell2)*M_PI*(-I()*N1/this->Lx + N2/this->Ly);
            }

            Multiplier = 0.0;
            Multiplier += c0*c0;
            Multiplier += c4*c4*this->LPolynomial(q2l2by2);
            Multiplier += c0*c4*this->divisor*Real(ZQ + ZQbar);
            Multiplier *= NormalizationFactor;
            Multiplier *= Multiplier;

            Precision = 2.0 * exp(- PIOnM * Q2);
            Coefficient += Multiplier*Precision * cos (N1 * Factor);
            N1 += 1.0;
        }

        Sum += Coefficient;
        N2 += this->NbrLzValue;
    }
    N2 = (double) (m1 - m4 - this->NbrLzValue);
    Coefficient = Sum;
    while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
        N1 = 1.0;
        Q2 = this->Ratio * N2 * N2;
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
            Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
            q2l2by2 = (N1*N1/Lx2 + N2*N2/Ly2)*Ell2*M_PI*M_PI*2;
            ZQ = 1;
            ZQbar =1;
            for(int i = 1; i <= this->MixingLevel; i++)
            {
                ZQ *= 2*sqrt(Ell2)*M_PI*(I()*N1/this->Lx + N2/this->Ly);
                ZQbar *= 2*sqrt(Ell2)*M_PI*(-I()*N1/this->Lx + N2/this->Ly);
            }

            Multiplier = 0.0;
            Multiplier += c0*c0;
            Multiplier += c4*c4*this->LPolynomial(q2l2by2);
            Multiplier += c0*c4*this->divisor*Real(ZQ + ZQbar);
            Multiplier *= NormalizationFactor;
            Multiplier *= Multiplier;

            Precision = 2.0 *  exp(- PIOnM * Q2);
            Coefficient += Multiplier * Precision * cos (N1 * Factor);
            N1 += 1.0;
        }
        Sum += Coefficient;
        N2 -= this->NbrLzValue;
    }
    //std::cout << Sum << endl;
    return (Sum / (4.0 * M_PI * this->NbrLzValue));


}

// double ParticleOnTorusDeltaHamiltonian::EvaluateInteractionCoefficientWithMixing(int m1, int m2, int m3, int m4)
// {
//     double Coefficient = 1.0;
//     double PIOnM = M_PI / ((double) this->NbrLzValue);
//     double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
//     double Cos4 = this->CosTheta*this->CosTheta*this->CosTheta*this->CosTheta;
//     double Sin4 = this->SinTheta*this->SinTheta*this->SinTheta*this->SinTheta;
//     double Multiplier = 0.0;
//     double Ly2 = this->Ly*this->Ly;
//     double Lx2 = this->Lx*this->Lx;
//     double Ell2 = this->Lx*this->Ly / (2 * M_PI * this->NbrLzValue);
//     double Sin2Cos2 = this->SinTheta*this->SinTheta*this->CosTheta*this->CosTheta;
//     double NormalizationFactor = 1 / (1 - 2*Sin2Cos2);
//     Complex ZQ,ZQbar = 1;

//     double Sum = 0.0;
//     double N2 = (double) (m1 - m4);
//     double N1;
//     double Q2;
//     double q2l2by2;
//     double Precision;
// //  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " : ";
//     while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
//     {
//         N1 = 1.0;
//         Q2 = this->Ratio * N2 * N2;
//         if (N2 != 0.0)
//         {
//             Coefficient = exp(- PIOnM * Q2);
//             Precision = Coefficient;
//         }
//         else
//         {
//             Coefficient = 1.0;
//             Precision = 1.0;
//         }
//         while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
//         {

//             Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
//             q2l2by2 = (N2*N2/Ly2 + N1*N1/Lx2)*Ell2*M_PI*M_PI*2;
//             ZQ = 1;
//             ZQbar = 1;
//             for(int i = 1; i <= this->MixingLevel; i++)
//             {
//                 ZQ *= 2*sqrt(Ell2)*M_PI*(I()*N1/this->Lx + N2/this->Ly);
//                 ZQbar *= 2*sqrt(Ell2)*M_PI*(-I()*N1/this->Lx + N2/this->Ly);
//             }

//             Multiplier = 0.0;
//             Multiplier += Cos4;
//             Multiplier += Sin4*this->LPolynomial(q2l2by2);
//             Multiplier += Sin2Cos2*this->divisor*Real(ZQ + ZQbar);
//             Multiplier *= NormalizationFactor;
//             Multiplier *= Multiplier;

//             Precision = 2.0 * exp(- PIOnM * Q2);
//             Coefficient += Multiplier*Precision * cos (N1 * Factor);
//             N1 += 1.0;
//         }

//         Sum += Coefficient;
//         N2 += this->NbrLzValue;
//     }
//     N2 = (double) (m1 - m4 - this->NbrLzValue);
//     Coefficient = Sum;
//     while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
//     {
//         N1 = 1.0;
//         Q2 = this->Ratio * N2 * N2;
//         if (N2 != 0.0)
//         {
//             Coefficient = exp(- PIOnM * Q2);
//             Precision = Coefficient;
//         }
//         else
//         {
//             Coefficient = 1.0;
//             Precision = 1.0;
//         }
//         while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
//         {
//             Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
//             q2l2by2 = (N1*N1/Lx2 + N2*N2/Ly2)*Ell2*M_PI*M_PI*2;
//             ZQ = 1;
//             ZQbar =1;
//             for(int i = 1; i <= this->MixingLevel; i++)
//             {
//                 ZQ *= 2*sqrt(Ell2)*M_PI*(I()*N1/this->Lx + N2/this->Ly);
//                 ZQbar *= 2*sqrt(Ell2)*M_PI*(-I()*N1/this->Lx + N2/this->Ly);
//             }

//             Multiplier = 0.0;
//             Multiplier += Cos4;
//             Multiplier += Sin4*this->LPolynomial(q2l2by2);
//             Multiplier += Sin2Cos2*this->divisor*Real(ZQ + ZQbar);
//             Multiplier *= NormalizationFactor;
//             Multiplier *= Multiplier;

//             Precision = 2.0 *  exp(- PIOnM * Q2);
//             Coefficient += Multiplier * Precision * cos (N1 * Factor);
//             N1 += 1.0;
//         }
//         Sum += Coefficient;
//         N2 -= this->NbrLzValue;
//     }
//     //std::cout << Sum << endl;
//     return (Sum / (4.0 * M_PI * this->NbrLzValue));


// }
// Evaluate Interaction Coefficients for cos^2 |0> + sin^2 |4>
/*
double ParticleOnTorusDeltaHamiltonian::EvaluateLL4MixingCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Cos4 = this->CosTheta*this->CosTheta*this->CosTheta*this->CosTheta;
  double Sin4 = this->SinTheta*this->SinTheta*this->SinTheta*this->SinTheta;
  double Multiplier = 0.0;
  double Ly2 = this->Ly*this->Ly;
  double Lx2 = this->Lx*this->Lx;
  double MomentumFactor = 2 * M_PI / this->Ly;
  double Ell2 = this->Lx*this->Ly / (2 * M_PI * this->NbrLzValue);
  double Pi4 = M_PI*M_PI*M_PI*M_PI;
  double Sin2Theta2 = 4*this->SinTheta*this->SinTheta*this->CosTheta*this->CosTheta;
  double NormalizationFactor = 1 / (1 - (1/2)*Sin2Theta2);

  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double q2l2by2;
  double Precision;
//  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " : ";
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
  {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
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
      Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
      q2l2by2 = (N2*N2/Lx2 + N1*N1/Ly2)*(N2*N2/Lx2 + N1*N1/Ly2)*Ell2*M_PI*M_PI*2;
      Multiplier = NormalizationFactor*(Cos4 + (1/24)*Sin4*(q2l2by2*q2l2by2*q2l2by2*q2l2by2
                                         - 16*q2l2by2*q2l2by2*q2l2by2
                                         + 72*q2l2by2*q2l2by2
                                         - 96*q2l2by2 +24)

                          + (1/2)*sqrt(2)/sqrt(3)*(Sin2Theta2)*(N2*N2*N2*N2/(Ly2*Ly2) + N1*N1*N1*N1/(Lx2*Lx2) - 6*N1*N1*N2*N2/(Lx2*Ly2)));
      Multiplier = Multiplier * Multiplier;
      Precision = 2.0 * exp(- PIOnM * Q2);
      Coefficient += Multiplier*Precision * cos (N1 * Factor);
      N1 += 1.0;
    }

    Sum += Coefficient;
    N2 += this->NbrLzValue;
  }
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = Sum;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
  {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
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
      q2l2by2 = (N2*N2/Lx2 + N1*N1/Ly2)*(N2*N2/Lx2 + N1*N1/Ly2)*Ell2*M_PI*M_PI*2;
      Multiplier = NormalizationFactor*(Cos4 + (1/24)*Sin4*(q2l2by2*q2l2by2*q2l2by2*q2l2by2
                                         - 16*q2l2by2*q2l2by2*q2l2by2
                                         + 72*q2l2by2*q2l2by2
                                         - 96*q2l2by2 +24)

                          + (1/2)*sqrt(2)/sqrt(3)*(Sin2Theta2)*(N2*N2*N2*N2/(Ly2*Ly2) + N1*N1*N1*N1/(Lx2*Lx2) - 6*N1*N1*N2*N2/(Lx2*Ly2)));
      Multiplier = Multiplier * Multiplier;
      Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
      Precision = 2.0 *  exp(- PIOnM * Q2);
      Coefficient += Multiplier * Precision * cos (N1 * Factor);
      N1 += 1.0;
    }
    Sum += Coefficient;
    N2 -= this->NbrLzValue;
  }
  //std::cout << Sum << endl;
  return (Sum / (4.0 * M_PI * this->NbrLzValue));

}
*/


// Evaluate Interaction Coefficients for cos^2 * |0> + sin^2 * |2>
/*
double ParticleOnTorusDeltaHamiltonian::EvaluateInteractionCoefficientWithMixing(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Cos4 = this->CosTheta*this->CosTheta*this->CosTheta*this->CosTheta;
  double Sin4 = this->SinTheta*this->SinTheta*this->SinTheta*this->SinTheta;
  double Multiplier = 0.0;
  double Ly2 = this->Ly*this->Ly;
  double Lx2 = this->Lx*this->Lx;
  double MomentumFactor = 2 * M_PI / this->Ly;
  double Ell2 = this->Lx*this->Ly / (2 * M_PI * this->NbrLzValue);
  double Pi4 = M_PI*M_PI*M_PI*M_PI;
  double Sin2Theta2 = 4*this->SinTheta*this->SinTheta*this->CosTheta*this->CosTheta;
  double NormalizationFactor = 1 / (1 - (1/2)*Sin2Theta2);


  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
//  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " : ";
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
  {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
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
      Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
      Multiplier = NormalizationFactor*(Cos4 + Sin4*(1
                                 + 2*Ell2*Ell2*Pi4*(N2*N2/Lx2 + N1*N1/Ly2)*(N2*N2/Lx2 + N1*N1/Ly2)
                                 - 4*Ell2*M_PI*M_PI*(N2*N2/Lx2 + N1*N1/Ly2)
                                 )
                          +sqrt(2)/2*Ell2*M_PI*M_PI*(Sin2Theta2)*(N2*N2/Ly2 - N1*N1/Lx2));
      Multiplier = Multiplier * Multiplier;
      std::cout << Multiplier << "\n";
      Precision = 2.0 * exp(- PIOnM * Q2);
      Coefficient += Multiplier*Precision * cos (N1 * Factor);
      N1 += 1.0;
    }

    Sum += Coefficient;
    N2 += this->NbrLzValue;
  }
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = Sum;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
  {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
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
      Multiplier = NormalizationFactor*(Cos4 + Sin4*(1
                                 + 2*Ell2*Ell2*Pi4*(N2*N2/Lx2 + N1*N1/Ly2)*(N2*N2/Lx2 + N1*N1/Ly2)
                                 - 4*Ell2*M_PI*M_PI*(N2*N2/Lx2 + N1*N1/Ly2)
                                 )
                          +sqrt(2)/2*Ell2*M_PI*M_PI*(Sin2Theta2)*(N2*N2/Ly2 - N1*N1/Lx2));
      Multiplier = Multiplier * Multiplier;
      Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
      Precision = 2.0 *  exp(- PIOnM * Q2);
      Coefficient += Multiplier * Precision * cos (N1 * Factor);
      N1 += 1.0;
    }
    Sum += Coefficient;
    N2 -= this->NbrLzValue;
  }
  //std::cout << Sum << endl;
  return (Sum / (4.0 * M_PI * this->NbrLzValue));

}
*/
