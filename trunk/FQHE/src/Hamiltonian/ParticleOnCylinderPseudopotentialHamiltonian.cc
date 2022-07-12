////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a cylinder with      //
//                        pseudopotential interaction                         //
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


#include "Hamiltonian/ParticleOnCylinderPseudopotentialHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
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
// confinement = amplitude of the quadratic confinement potential
// lineCharge = use line charge instead of parabolic confinement
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = array containing pseudopotential values
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderPseudopotentialHamiltonian::ParticleOnCylinderPseudopotentialHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   double ratio, double confinement, bool lineCharge, int nbrPseudopotentials, double* pseudopotentials, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;

  this->NbrPseudopotentials = nbrPseudopotentials;
  this->Pseudopotentials = new double[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->Pseudopotentials[i] = pseudopotentials[i];
  this->LaguerrePolynomials =new Polynomial[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->Architecture = architecture;
  this->Confinement = confinement;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;


  this->OneBodyInteractionFactors = new Complex [this->NbrLzValue];
  Complex Factor;
  double kappa = sqrt(2.0 * M_PI /(this->NbrLzValue * this->Ratio));
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double Height = sqrt(2.0 * M_PI * this->NbrLzValue / this->Ratio);

  for (int i = 0; i < this->NbrLzValue; ++i)
   { 
       if (lineCharge == false)
         {
           //Parabolic confinement     
           Factor.Re = this->Confinement * pow(kappa * (i - 0.5 * this->MaxMomentum), 2.0);
           Factor.Im = 0.0;
         }
       else
        {
           double error;

           //Realistic confinement
           Factor.Re = this->Confinement * this->LineChargeMatrixElement(i, this->MaxMomentum, Length, Height, error);
           Factor.Im = 0.0;
           if (fabs(error) > 1e-6)
             cout << "Insufficient accuracy in line charge " << endl;
         } 
       this->OneBodyInteractionFactors[i] += Factor;

       if (Norm(Factor) != 0.0)
         cout << "One body: i= " << i << " " << Factor << endl;

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

ParticleOnCylinderPseudopotentialHamiltonian::~ParticleOnCylinderPseudopotentialHamiltonian() 
{
  delete[] this->Pseudopotentials;
  delete[] this->LaguerrePolynomials;
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;

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

void ParticleOnCylinderPseudopotentialHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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

void ParticleOnCylinderPseudopotentialHamiltonian::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderPseudopotentialHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
  Complex* TmpCoefficient = new Complex [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m4 <= this->MaxMomentum))
  	        if (m3 > m4)
		  {
 		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
                                           + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
                                           - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
                                           - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));

		     //cout << m1 << " " << m2 << " " << m3 << " " << m4 << " : " << this->EvaluateInteractionCoefficient(m1, m2, m3, m4) << endl;
		     //cout << m1 << " " << m2 << " " << m4 << " " << m3 << " : " << this->EvaluateInteractionCoefficient(m1,m2,m4,m3) << endl;
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
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
              if ((m4 >= 0) && (m4 <= this->MaxMomentum))
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
  else //bosons
    {
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m4 <= this->MaxMomentum))
		{
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
		      if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
			MaxCoefficient = Norm(TmpCoefficient[Pos]);
		      ++Pos;
		    }
		  else
		    {
		      if (m3 == m4)
			{
			  if (m1 != m2)
			    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
						   + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
			  else
			    TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4);
			  if (MaxCoefficient < Norm(TmpCoefficient[Pos]))
			    MaxCoefficient = Norm(TmpCoefficient[Pos]);
			  ++Pos;
			}
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
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m4 <= this->MaxMomentum))
	       if (m3 >= m4)
		{
		  if (Norm(TmpCoefficient[Pos]) > MaxCoefficient)
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

Complex ParticleOnCylinderPseudopotentialHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;	
  double error;

  Complex Coefficient(0, 0);

  Coefficient.Re = this->PseudopotentialMatrixElement(Xm1-Xm4, Xm1-Xm3, this->NbrPseudopotentials, this->Pseudopotentials, this->LaguerrePolynomials, error);
  Coefficient.Im = 0.0;

  if (fabs(error) > 1e-6)
    {
      cout << "Warning: large error in matrix elements! " ;
    }

  return (Coefficient/Length);
}


#ifdef HAVE_GSL  

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace CoulombMatEl
{

struct f_params {
  double Xj14;
  double Xj13;
  int NbrPseudopotentials;
  double* Pseudopotentials;
  Polynomial* LaguerrePolynomials;
};

double Integrand(double qx, void *p)
{
  f_params &params= *reinterpret_cast<f_params *>(p);

  double q2 = qx * qx + params.Xj14 * params.Xj14;

  double Vq = 0.0;
  for (int i = 0; i < params.NbrPseudopotentials; ++i)
     if (params.Pseudopotentials[i] != 0.0)
        Vq += params.Pseudopotentials[i] * params.LaguerrePolynomials[i].PolynomialEvaluate(q2);

  return (exp(-0.5*q2) * 2.0 * cos(qx * params.Xj13) * Vq);  
}


}

namespace LineChargeMatEl
{

struct f_params {
  int index;
  int Nphi;
  double L;
  double H;
};

double LineChargeIntegrand(double x, void *p)
{
  f_params &params= *reinterpret_cast<f_params *>(p);

  double Xm = (params.index - 0.5 * params.Nphi) * 2.0 * M_PI/params.L;

  double LogNum, LogDen, IntegrandPos, IntegrandNeg;

  LogNum = x - 0.5 * params.H + sqrt(pow(x - 0.5 * params.H, 2.0) + pow(params.L/(2.0 * M_PI),2.0));
  LogDen = x + 0.5 * params.H + sqrt(pow(x + 0.5 * params.H, 2.0) + pow(params.L/(2.0 * M_PI),2.0));

  IntegrandPos = exp(-pow(x - Xm,2.0)) * log(fabs(LogNum/LogDen));

  LogNum = -x - 0.5 * params.H + sqrt(pow(-x - 0.5 * params.H, 2.0) + pow(params.L/(2.0 * M_PI),2.0));
  LogDen = -x + 0.5 * params.H + sqrt(pow(-x + 0.5 * params.H, 2.0) + pow(params.L/(2.0 * M_PI),2.0));

  IntegrandNeg = exp(-pow(-x - Xm,2.0)) * log(fabs(LogNum/LogDen));

  return ( ((IntegrandPos + IntegrandNeg) * params.Nphi ) / (sqrt(M_PI) ) );
}


}
#endif

double ParticleOnCylinderPseudopotentialHamiltonian::PseudopotentialMatrixElement(double xj14, double xj13, int nbrPseudopotentials, double* pseudopotentials, Polynomial* laguerrePolynomials, double &error)
{
#ifdef HAVE_GSL

  gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (1000000);

  double lower_limit = 0.0;
  double abs_error = 1.0e-10;
  double rel_error = 1.0e-10;
  double result, finalresult;

  CoulombMatEl::f_params params;
  params.Xj14 = xj14;
  params.Xj13 = xj13;
  params.NbrPseudopotentials = nbrPseudopotentials;
  params.Pseudopotentials = pseudopotentials;
  params.LaguerrePolynomials = laguerrePolynomials;
  
  gsl_function F;

  F.function = &CoulombMatEl::Integrand;
  F.params = reinterpret_cast<void *>(&params);
    
  gsl_integration_qagiu (&F, lower_limit,
			 abs_error, rel_error, 10000, work_ptr, &result,
			 &error);

  finalresult = result;

  gsl_integration_workspace_free (work_ptr);
  
  //cout << "result          = " << result << endl;
  //cout << "estimated error = " << error << endl;
  //cout << "intervals =  " << work_ptr->size << endl;

  return finalresult;
#endif
}

double ParticleOnCylinderPseudopotentialHamiltonian::LineChargeMatrixElement(int index, int MaxMomentum, double Length, double Height, double &error)
{
#ifdef HAVE_GSL

  gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (1000000);

  double lower_limit = 0.0;
  double abs_error = 1.0e-10;
  double rel_error = 1.0e-10;
  double result, finalresult;

  LineChargeMatEl::f_params params;
  params.index = index;
  params.Nphi = MaxMomentum;
  params.L = Length;
  params.H = Height;

  gsl_function F;

  F.function = &LineChargeMatEl::LineChargeIntegrand;
  F.params = reinterpret_cast<void *>(&params);
    
  gsl_integration_qagiu (&F, lower_limit,
			 abs_error, rel_error, 10000, work_ptr, &result,
			 &error);

  finalresult = result;

  gsl_integration_workspace_free (work_ptr);
  
  //cout << "result          = " << result << endl;
  //cout << "estimated error = " << error << endl;
  //cout << "intervals =  " << work_ptr->size << endl;

  return finalresult;
#endif
}
