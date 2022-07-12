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


#include "Hamiltonian/ParticleOnCylinderOrbitalProjection.h"
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
// orbitalIndex = index of the orbital to be projected out
// anisotropy = shape (anisotropy) parameter
// x0,y0 = position in real space
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCylinderOrbitalProjection::ParticleOnCylinderOrbitalProjection(ParticleOnSphere* particles, int nbrParticles, int maxMomentum,
										   double ratio, int orbitalIndex, double anisotropy, double x0, double y0, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->OrbitalIndex = orbitalIndex;
  this->Anisotropy = anisotropy;
  this->X0 = x0;
  this->Y0 = y0;
  this->Architecture = architecture;
  this->EnergyShift = 0.0;

  this->LaguerreM = new Polynomial[this->OrbitalIndex + 1];
  for (int i = 0; i < (this->OrbitalIndex + 1); ++i)
    this->LaguerreM[i] = LaguerrePolynomial(i);
 
  this->NbrOneBodyInteractionFactors = this->NbrLzValue * this->NbrLzValue;
  this->OneBodyInteractionFactors = new Complex [this->NbrOneBodyInteractionFactors];
  this->OneBodyM1Values = new int [this->NbrOneBodyInteractionFactors];
  this->OneBodyM2Values = new int [this->NbrOneBodyInteractionFactors];


  double kappa = sqrt(2.0 * M_PI /(this->NbrLzValue * this->Ratio));
  this->NbrOneBodyInteractionFactors = 0;
  for (int i = 0; i < this->NbrLzValue; ++i)
    for (int j = 0; j < this->NbrLzValue; ++j)
      //if (i==j)
      { 
        this->OneBodyM1Values[this->NbrOneBodyInteractionFactors] = i;
        this->OneBodyM2Values[this->NbrOneBodyInteractionFactors] = j;
        this->OneBodyInteractionFactors[this->NbrOneBodyInteractionFactors] = this->EvaluateInteractionCoefficient(i, j);
        //if (Norm(this->OneBodyInteractionFactors[this->NbrOneBodyInteractionFactors])>1e-8)
        //cout<<"i= "<<i<<" j= "<<j<<" "<<this->OneBodyInteractionFactors[this->NbrOneBodyInteractionFactors]<<endl;
        this->NbrOneBodyInteractionFactors++;
      }
  cout<<"Done calculating one body terms, total nbr = " << this->NbrOneBodyInteractionFactors << endl;


  if (precalculationFileName != 0)
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnCylinderOrbitalProjection::~ParticleOnCylinderOrbitalProjection() 
{
  delete[] this->OneBodyInteractionFactors;
  delete[] this->OneBodyM1Values;
  delete[] this->OneBodyM2Values;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnCylinderOrbitalProjection::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnSphere*) hilbertSpace;
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnCylinderOrbitalProjection::ShiftHamiltonian (double shift)
{
  this->EnergyShift = shift;
}
  
// evaluate all interaction factors
//   
void ParticleOnCylinderOrbitalProjection::EvaluateInteractionFactors()
{
}

// evaluate the numerical coefficient  in front of the a+_m1 a_m2 coupling term
//
// m1 = first index
// m2 = second index
// return value = numerical coefficient

Complex ParticleOnCylinderOrbitalProjection::EvaluateInteractionCoefficient(int m1, int m2)
{
  double Length = sqrt(2.0 * M_PI * this->Ratio * this->NbrLzValue);
  double kappa = 2.0 * M_PI/Length;
  double Xm1 = kappa * m1;
  double Xm2 = kappa * m2;
  double error;

  Complex Coefficient, Phase;

  Phase.Re = cos((Xm1-Xm2) * this->Y0);
  Phase.Im = -sin((Xm1-Xm2) * this->Y0);

  Coefficient = exp(-0.25*(Xm1-Xm2)*(Xm1-Xm2)/this->Anisotropy) * Phase * this->OrbitalProjectionMatrixElement(Xm1, Xm2, this->X0, this->Y0, this->OrbitalIndex, this->Anisotropy, this->LaguerreM, kappa, this->MaxMomentum, error);
 
  return (-Coefficient/Length);
}


#ifdef HAVE_GSL  

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

namespace OrbitalProjectionMatEl
{

struct f_params {
  double Xj1;
  double Xj2;
  double Anisotropy;
  double X0;
  double Y0;
  double Kappa;
  double MaxMomentum;
  int OrbitalIndex;
  Polynomial* LaguerreM;
};

double Integrand(double qx, void *p)
{
  f_params &params= *reinterpret_cast<f_params *>(p);

  double Q2 = params.Anisotropy * qx * qx + (params.Xj1 - params.Xj2) * (params.Xj1 - params.Xj2)/params.Anisotropy;

  return (exp(-0.25*qx*qx*params.Anisotropy) * params.LaguerreM[params.OrbitalIndex].PolynomialEvaluate(0.5 * Q2) * 2.0 * cos(qx * (params.X0 + 0.5 * (params.Xj1 + params.Xj2 - params.Kappa * params.MaxMomentum))));   
}

}


#endif

double ParticleOnCylinderOrbitalProjection::OrbitalProjectionMatrixElement(double xj1, double xj2, double x0, double y0, int orbitalIndex, double anisotropy, Polynomial* laguerreM, double kappa, int maxMomentum, double &error)
{
#ifdef HAVE_GSL

  gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (1000000);

  double lower_limit = 0.0;
  double abs_error = 1.0e-8;
  double rel_error = 1.0e-8;
  double result;

  OrbitalProjectionMatEl::f_params params;
  params.Xj1=xj1;
  params.Xj2=xj2;
  params.Anisotropy = anisotropy;
  params.X0 = x0;
  params.Y0 = y0; 
  params.Kappa = kappa;
  params.MaxMomentum = maxMomentum;
  params.OrbitalIndex = orbitalIndex;
  params.LaguerreM = laguerreM;

  gsl_function F;

  F.function = &OrbitalProjectionMatEl::Integrand;
  F.params = reinterpret_cast<void *>(&params);
    
  gsl_integration_qagiu (&F, lower_limit,
			 abs_error, rel_error, 10000, work_ptr, &result,
			 &error);

  gsl_integration_workspace_free (work_ptr);

  //cout << "result          = " << result << endl;
  //cout << "estimated error = " << error << endl;
  //cout << "intervals =  " << work_ptr->size << endl;

  return result;
#endif
}


ComplexVector& ParticleOnCylinderOrbitalProjection::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  int Index;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();

  if (this->OneBodyInteractionFactors != 0)
    for (int i = firstComponent; i < LastComponent; ++i)
      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j) 
        {
           Index = TmpParticles->AdA(i, this->OneBodyM1Values[j], this->OneBodyM2Values[j], Coefficient);
 
           if (Index < Dim)
            {
              vDestination[Index] += Coefficient * this->OneBodyInteractionFactors[j] * vSource[i];
              //if ((Coefficient * this->OneBodyInteractionFactors[j] * vSource[i]) != 0)
              //{
              //cout<<"Attempt "; TmpParticles->PrintState(cout,i); cout<<" i= "<<this->OneBodyM1Values[j]<<" j= "<<this->OneBodyM2Values[j]<<" Index= "<<Index<< " Coeff= "<<Coefficient<<" ";
              //TmpParticles->PrintState(cout,Index);
              //cout<<endl;
              //}
            }
        }           

  delete TmpParticles;

  //cout<<vDestination;
  //exit(1);
  return vDestination;
}

