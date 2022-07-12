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


#include "Hamiltonian/ParticleOnCylinderCoulombHamiltonian.h"
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
// fillingFactor = filling factor of the FQHE state
// landauLevel = LL index
// confinement = amplitude of the quadratic confinement potential
// lineCharge = use line charge instead of parabolic confinement
// electricFieldParameter = amplitude of the electric field along the cylinder
// bFieldfParameter = amplitude of the magnetic field (to set the energy scale)
// deltaV1 = tweak of V1 pseudopotential
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderCoulombHamiltonian::ParticleOnCylinderCoulombHamiltonian(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   double ratio, double fillingFactor, int landauLevel, double confinement, bool lineCharge, double electricFieldParameter, double bFieldParameter, double deltaV1, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->FillingFactor = fillingFactor;
  this->LLIndex = landauLevel;
  this->Architecture = architecture;
  this->Confinement = confinement;
  this->ElectricField = electricFieldParameter;
  this->MagneticField = bFieldParameter;
  this->DeltaV1 = deltaV1;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;


  //this->OneBodyInteractionFactors = 0;

  //add the Hartree terms
  this->OneBodyInteractionFactors = new Complex [this->NbrLzValue];
  Complex Factor;
  double kappa = sqrt(2.0 * M_PI /(this->NbrLzValue * this->Ratio));
  for (int i = 0; i < this->NbrLzValue; ++i)
   { 
     Factor.Re = 0.0;
     Factor.Im = 0.0;
     for (int j = 0; j < this->NbrLzValue; ++j)
      {
        Factor += this->EvaluateInteractionCoefficient(i, j, j, i);
      }
     Factor *= (-this->FillingFactor);
     this->OneBodyInteractionFactors[i] = Factor;
   }

  if ((this->ElectricField != 0.0) || (this->Confinement != 0.0))
    {
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
             Factor.Re = this->Confinement * this->LineChargeMatrixElement(i, this->NbrParticles, this->MaxMomentum, Length, Height, error);
             Factor.Im = 0.0;
             if (fabs(error) > 1e-6)
               cout << "Insufficient accuracy in line charge " << endl;
           } 

           //add the contribution from electric field
           Factor.Re += 0.194 * sqrt(this->MagneticField) * ((this->ElectricField/(1.0 + this->ElectricField)) * kappa * kappa * ((double)i - 0.5 * this->MaxMomentum) * ((double)i - 0.5 * this->MaxMomentum)); 
           Factor.Im += 0.0;
	   this->OneBodyInteractionFactors[i] += Factor;

           if (Norm(Factor) != 0.0)
             cout << "One body: i= " << i << " " << this->OneBodyInteractionFactors[i] << endl;
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

ParticleOnCylinderCoulombHamiltonian::~ParticleOnCylinderCoulombHamiltonian() 
{
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

void ParticleOnCylinderCoulombHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
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

void ParticleOnCylinderCoulombHamiltonian::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnCylinderCoulombHamiltonian::EvaluateInteractionFactors()
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

Complex ParticleOnCylinderCoulombHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double Xm3 = kappa * m3;
  double Xm4 = kappa * m4;	
  double error;

  Complex Coefficient(0,0);

  if (this->ElectricField == 0.0)
   {
     Coefficient.Re = exp(-0.5*(Xm1-Xm4)*(Xm1-Xm4)) * this->CoulombMatrixElement(Xm1-Xm4,Xm1-Xm3,error);
     Coefficient.Im = 0.0;
     return (Coefficient/(Length * sqrt(2.0 * M_PI)));
   }
  else
   {
     double alpha = sqrt(1.0 + this->ElectricField);
     Coefficient.Re = exp(-pow(Xm1-Xm4,2.0)/(2.0 * pow(alpha,3.0))) * this->CoulombMatrixElement(Xm1-Xm4,Xm1-Xm3,error);;
     Coefficient.Im = 0.0;
     return (Coefficient/(Length * sqrt(2.0 * M_PI * alpha * alpha * alpha)));
   }
}


#ifdef HAVE_GSL  

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace CoulombMatEl
{

struct f_params {
  double Xj14;
  double Xj13;
  double ElectricField;
  int LLIndex;
  double DeltaV1;
};

double Integrand(double qx, void *p)
{
  f_params &params= *reinterpret_cast<f_params *>(p);

  if (params.ElectricField == 0.0)
   {
    double CoulombFormFactor = 1.0/sqrt(qx*qx+params.Xj14*params.Xj14);
    if (params.LLIndex == 1)
      CoulombFormFactor *= pow(1.0-0.5*(qx*qx+params.Xj14*params.Xj14), 2.0);
    else if (params.LLIndex >= 2)
      {
        cout << "Only considering up to LL=1" << endl;
        exit(1);
      }
    if (params.DeltaV1 != 0.0)
      CoulombFormFactor -= 2.0 * params.DeltaV1 * (qx * qx + params.Xj14*params.Xj14) * 2.0 * M_PI;

    return (exp(-0.5*qx*qx) * 2.0 * cos(qx*params.Xj13) * CoulombFormFactor);   
   }
  else //applied electric field
   {
    double alpha = sqrt(1.0 + params.ElectricField);

    double CoulombFormFactor = 1.0/sqrt(qx*qx+params.Xj14*params.Xj14);
    if (params.LLIndex == 1)
      CoulombFormFactor *= pow(1.0-0.5*(qx*qx/alpha+params.Xj14*params.Xj14/pow(alpha,3.0)), 2.0);
    else if (params.LLIndex >= 2)
      {
        cout << "Only considering up to LL=1" << endl;
        exit(1);
      }
    if (params.DeltaV1 != 0.0)
      CoulombFormFactor -= 2.0 * params.DeltaV1 * (qx * qx + params.Xj14*params.Xj14) * 2.0 * M_PI;

    return (exp(-0.5*qx*qx/alpha) * 2.0 * cos(qx*params.Xj13/(alpha*alpha)) * CoulombFormFactor);
   }
}

double SingularIntegrand(double qx, void *p)
{
  f_params &params= *reinterpret_cast<f_params *>(p);

  if (params.ElectricField == 0.0)
   {
    double CoulombFormFactor = 1.0/qx;
    if (params.LLIndex == 1)
      CoulombFormFactor *= pow(1.0-0.5*(qx*qx), 2.0);
    else if (params.LLIndex >= 2)
      {
        cout << "Only considering up to LL=1" << endl;
        exit(1);
      }
    if (params.DeltaV1 != 0.0)
      CoulombFormFactor -= 2.0 * params.DeltaV1 * (qx * qx) * 2.0 * M_PI;

    return (2.0 * (exp(-0.5*qx*qx) * cos(qx*params.Xj13)-1.0) * CoulombFormFactor);   
   }
  else //applied electric field
   {
    double alpha = sqrt(1.0 + params.ElectricField);

    double CoulombFormFactor = 1.0/qx;
    if (params.LLIndex == 1)
      CoulombFormFactor *= pow(1.0-0.5*(qx*qx/alpha), 2.0);
    else if (params.LLIndex >= 2)
      {
        cout << "Only considering up to LL=1" << endl;
        exit(1);
      }
    if (params.DeltaV1 != 0.0)
      CoulombFormFactor -= 2.0 * params.DeltaV1 * (qx * qx) * 2.0 * M_PI;

    return (2.0 * (exp(-0.5*qx*qx/alpha) * cos(qx*params.Xj13/(alpha*alpha)) - 1.0) * CoulombFormFactor);
   }
}

}

namespace LineChargeMatEl
{

struct f_params {
  int index;
  int Ne;
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

double ParticleOnCylinderCoulombHamiltonian::CoulombMatrixElement(double xj14, double xj13, double &error)
{
#ifdef HAVE_GSL

  gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (1000000);

  double lower_limit = 0.0;
  double abs_error = 1.0e-10;
  double rel_error = 1.0e-10;
  double result, finalresult;

  CoulombMatEl::f_params params;
  params.Xj14=xj14;
  params.Xj13=xj13;
  params.ElectricField = this->ElectricField;
  params.LLIndex = this->LLIndex;
  params.DeltaV1 = this->DeltaV1;
  
  //cout<<"Xj14 "<<params.Xj14<<" Xj13 "<<params.Xj13<<" Efield "<<params.ElectricField<<endl;

  gsl_function F;

  if (params.Xj14 != 0.0)
    {
      F.function = &CoulombMatEl::Integrand;
      F.params = reinterpret_cast<void *>(&params);
    
      gsl_integration_qagiu (&F, lower_limit,
			 abs_error, rel_error, 10000, work_ptr, &result,
			 &error);

      finalresult = result;
     }
   else
    {
      lower_limit = 1.0;
      F.function = &CoulombMatEl::Integrand;
      F.params = reinterpret_cast<void *>(&params);
    
      gsl_integration_qagiu (&F, lower_limit,
			 abs_error, rel_error, 10000, work_ptr, &result,
			 &error);

      finalresult = result;

      F.function = &CoulombMatEl::SingularIntegrand;
      F.params = reinterpret_cast<void *>(&params);

      gsl_integration_qags (&F, 0.0, 1.0, 0, 1e-10, 10000,
                             work_ptr, &result, &error); 
 
      finalresult += result;
    }

  //cout<<"result "<<result<<endl;

  gsl_integration_workspace_free (work_ptr);
  
  //cout << "result          = " << result << endl;
  //cout << "estimated error = " << error << endl;
  //cout << "intervals =  " << work_ptr->size << endl;

  return finalresult;
#endif
}

double ParticleOnCylinderCoulombHamiltonian::LineChargeMatrixElement(int index, int NbrParticles, int MaxMomentum, double Length, double Height, double &error)
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
  params.Ne = NbrParticles;
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
