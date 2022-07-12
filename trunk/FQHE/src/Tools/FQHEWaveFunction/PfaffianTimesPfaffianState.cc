////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class implementing a paired CF wave function on the sphere          //
//                                                                            //
//                        last modification : 18/05/2007                      //
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
#include "PfaffianTimesPfaffianState.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

PfaffianTimesPfaffianState::PfaffianTimesPfaffianState()
{
  this->NbrParticlesPerLayer = 0;
  this->NbrParameters = 0;
}

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1 (same coefficients are used for up and down)
  // MRCoeff = coefficient of MR term
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised


PfaffianTimesPfaffianState::PfaffianTimesPfaffianState(int nbrParticles, int nbrLandauLevels,
						       double * givenCFCoefficients, double MRcoeff,
						       bool correctPrefactors)
{
  int jastrowPower=2;
  if (nbrParticles%4!=0)
    {
      cout << "Error: PfaffianTimesPfaffianState requires an even number of particles in each layer!"<<endl;
      exit(1);
    }
  this->NbrParticlesPerLayer = nbrParticles/2;
  this->NbrLandauLevels = nbrLandauLevels;
  this->NbrParameters = this->NbrLandauLevels+1; // inherited field
  this->AbsEffectiveFlux = 1;
  this->MRCoeff = MRcoeff;
  int NbrEffectiveFlux = -1;
  this->Orbitals1 = new JainCFOnSphereOrbitals(NbrParticlesPerLayer, nbrLandauLevels, NbrEffectiveFlux, jastrowPower);
  this->Orbitals2 = new JainCFOnSphereOrbitals(NbrParticlesPerLayer, nbrLandauLevels, NbrEffectiveFlux, jastrowPower);
  this->TrialParameters= new double [NbrParameters];
  for (int i=0; i<NbrLandauLevels; ++i) this->TrialParameters[i] = givenCFCoefficients[i];
  this->TrialParameters[this->NbrLandauLevels] = this->MRCoeff;
  this->ElementNorm=1.0;
  this->Flag.Initialize();
  this->Slater1 = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  this->Slater2 = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  
  this->J11 = new Complex[this->NbrParticlesPerLayer];
  this->J22 = new Complex[this->NbrParticlesPerLayer];
  
  this->gAlpha1 = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha1[i]= new Complex[NbrParticlesPerLayer*NbrParticlesPerLayer];

  this->gAlpha2 = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha2[i]= new Complex[NbrParticlesPerLayer*NbrParticlesPerLayer];

  if (correctPrefactors)
    {
      int p=Orbitals1->GetJastrowPower();
      FactorialCoefficient Coef;
      for (int n=0; n<NbrLandauLevels; ++n)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialDivide(NbrEffectiveFlux+p*(nbrParticles-1)+2,NbrEffectiveFlux+2*p*(nbrParticles-1)+1);
	  Coef.PartialFactorialMultiply(p*(nbrParticles-1)+n+2,2*p*(nbrParticles-1)+n+1);
	  this->TrialParameters[n]*=Coef.GetNumericalValue()*Coef.GetNumericalValue();
	  //cout << "Correction["<<n<<"]="<<Coef.GetNumericalValue()<<"from: " << nbrEffectiveFlux+p*(nbrParticles-1)+2<<","<<nbrEffectiveFlux+2*p*(nbrParticles-1)+1<<","<< p*(nbrParticles-1)+n+2<<","<<2*p*(nbrParticles-1)+n+1<<endl;
	}
    }
}

// copy constructor
//
// function = reference on the wave function to copy

PfaffianTimesPfaffianState::PfaffianTimesPfaffianState(const PfaffianTimesPfaffianState& function)
{
  this->NbrParticlesPerLayer = function.NbrParticlesPerLayer;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->NbrParameters = function.NbrParameters;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->Flag = function.Flag;  
  this->Orbitals1 = function.Orbitals1;
  this->Orbitals2 = function.Orbitals2;
  this->TrialParameters=function.TrialParameters;
  this->ElementNorm=function.ElementNorm;

  this->Slater1=function.Slater1;
  this->Slater2=function.Slater2;
  
  this->J11 = new Complex[this->NbrParticlesPerLayer];
  this->J22 = new Complex[this->NbrParticlesPerLayer];
  
  this->gAlpha1 = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha1[i]= new Complex[NbrParticlesPerLayer*NbrParticlesPerLayer];
  
  this->gAlpha2 = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha2[i]= new Complex[NbrParticlesPerLayer*NbrParticlesPerLayer];

  
}

// destructor
//

PfaffianTimesPfaffianState::~PfaffianTimesPfaffianState()
{
  if ( (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete Orbitals1;
      delete Orbitals2;
      delete [] TrialParameters;
    }
  for (int i=0; i< this->NbrLandauLevels; ++i)
    {
      delete [] gAlpha1[i];
      delete [] gAlpha2[i];
    }
  delete [] gAlpha1;
  delete [] gAlpha2;

  delete [] J11;
  delete [] J22;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianTimesPfaffianState::Clone ()
{
  return new PfaffianTimesPfaffianState(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianTimesPfaffianState::operator ()(RealVector& x)
{
  RealVector part=x.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();
  Complex tmp, tmp2;

  // initialize Slater determinant 
  for (int i=0;i<this->NbrParticlesPerLayer;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  tmp2=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    {
	      tmp+=this->TrialParameters[n]*this->gAlpha1[n][i*this->NbrParticlesPerLayer+j];
	      tmp2+=this->TrialParameters[n]*this->gAlpha2[n][i*this->NbrParticlesPerLayer+j];
	    }
	  tmp+=MRCoeff/Orbitals1->JastrowFactorElement(i,j);
	  tmp2+=MRCoeff/Orbitals2->JastrowFactorElement(i,j);
	  tmp *=  - this->ElementNorm * J11[i]*J11[j];
	  tmp2 *=  - this->ElementNorm * J22[i]*J22[j];
	  Slater1->SetMatrixElement(i,j,tmp);
	  Slater2->SetMatrixElement(i,j,tmp2);
	}
    }
  //cout << *Matrix << endl;
  return Slater1->Pfaffian()*Slater2->Pfaffian() * Interpolation * JInter / J1 / J2 ;
  // return JInter * J1 * J2;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex PfaffianTimesPfaffianState::CalculateFromSpinorVariables(ComplexVector& uv)
{
  ComplexVector part=uv.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->OrbitalValues1 = Orbitals1->CalculateFromSpinorVariables(part);
  part = uv.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->OrbitalValues2 = Orbitals2->CalculateFromSpinorVariables(part);
  
  this->EvaluateTables();
  Complex tmp, tmp2;

  // initialize Slater determinant 
  for (int i=0;i<this->NbrParticlesPerLayer;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  tmp2=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    {
	      tmp+=this->TrialParameters[n]*this->gAlpha1[n][i*this->NbrParticlesPerLayer+j];
	      tmp2+=this->TrialParameters[n]*this->gAlpha2[n][i*this->NbrParticlesPerLayer+j];
	    }
	  tmp+=MRCoeff/Orbitals1->JastrowFactorElement(i,j);
	  tmp2+=MRCoeff/Orbitals2->JastrowFactorElement(i,j);
	  tmp *=  - this->ElementNorm * J11[i]*J11[j];
	  tmp2 *=  - this->ElementNorm * J22[i]*J22[j];
	  Slater1->SetMatrixElement(i,j,tmp);
	  Slater2->SetMatrixElement(i,j,tmp2);
	}
    }
  //cout << *Matrix << endl;
  return Slater1->Pfaffian()*Slater2->Pfaffian()*Interpolation * JInter / J1 / J2;
}

// get a value of the wavefunction for the last set of coordinates, but with different variational coefficients
// coefficients: array of variational coefficients f_0, f_1, ...
//                    length has to be identical to the initial set of parameters!
// singular: new prefactor of the Moore-Read part
//
Complex PfaffianTimesPfaffianState::GetForOtherParameters( double *coefficients)
{
  Complex tmp, tmp2;

  // initialize Slater determinant 
  for (int i=0;i<this->NbrParticlesPerLayer;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  tmp2=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    {
	      tmp+=coefficients[n]*this->gAlpha1[n][i*this->NbrParticlesPerLayer+j];
	      tmp2+=coefficients[n]*this->gAlpha2[n][i*this->NbrParticlesPerLayer+j];
	    }
	  tmp+=coefficients[this->NbrLandauLevels]/Orbitals1->JastrowFactorElement(i,j);
	  tmp2+=coefficients[this->NbrLandauLevels]/Orbitals2->JastrowFactorElement(i,j);
	  tmp *=  - this->ElementNorm * J11[i]*J11[j];
	  tmp2 *=  - this->ElementNorm * J22[i]*J22[j];
	  Slater1->SetMatrixElement(i,j,tmp);
	  Slater2->SetMatrixElement(i,j,tmp2);
	}
    }
  //cout << *Matrix << endl;
  return Slater1->Pfaffian()*Slater2->Pfaffian()*Interpolation * JInter / J1 / J2;
}


// do many evaluations, storing the result in the vector results given in the call
// x: positions to evaluate the wavefuntion in
// format for passing parameters in the matrix coefficients: coefficients[nbrSet][LandauLevel],
// the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
void PfaffianTimesPfaffianState::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
{
  RealVector part=x.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();

  Complex tmp, tmp2;
  int numParamSets=results.GetVectorDimension();
  double *tmpCoefficients;
  for (int s=0; s<numParamSets; ++s)
    {
      tmpCoefficients = coefficients[s];
      /*      
      cout<<"tmpCoefficients="<<tmpCoefficients[0];
      for (int k=1; k<this->NbrLandauLevels; ++k)
	cout<<" "<<tmpCoefficients[k];
      cout<<"; "<<tmpCoefficients[NbrLandauLevels]<<endl;
      */
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<i;++j)
	    {
	      tmp=0.0;
	      tmp2=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		{
		  tmp+=tmpCoefficients[n]*this->gAlpha1[n][i*this->NbrParticlesPerLayer+j];
		  tmp2+=tmpCoefficients[n]*this->gAlpha2[n][i*this->NbrParticlesPerLayer+j];
		}
 	      tmp+=tmpCoefficients[this->NbrLandauLevels]/Orbitals1->JastrowFactorElement(i,j);
 	      tmp2+=tmpCoefficients[this->NbrLandauLevels]/Orbitals2->JastrowFactorElement(i,j);
	      tmp *=  - this->ElementNorm * J11[i]*J11[j];
	      tmp2 *=  - this->ElementNorm * J22[i]*J22[j];
	      Slater1->SetMatrixElement(i,j,tmp);
	      Slater2->SetMatrixElement(i,j,tmp2);
	    }
	}
      //cout << *Matrix << endl;
       tmp = Slater1->Pfaffian()*Slater2->Pfaffian()*Interpolation * JInter / J1 / J2;
       results.Re(s)=Real(tmp);
       results.Im(s)=Imag(tmp);
    }
}  


void PfaffianTimesPfaffianState::SetTrialParameters(double * coefficients)
{
  for (int n=0; n<this->NbrParameters; ++n)
    this->TrialParameters[n]=coefficients[n];
  this->MRCoeff =  this->TrialParameters[this->NbrLandauLevels];
}

// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void PfaffianTimesPfaffianState::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      //cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)1.0/this->NbrParticlesPerLayer);
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)1.0/this->NbrParticlesPerLayer);
      else 
	this->ElementNorm*= pow(det,(double)-1.0/this->NbrParticlesPerLayer);
      det=Norm((*this)(x));
      //cout <<"N'="<< this->ElementNorm << endl;
    }
}


// normalize the wave-function over an average number of MC positions

void PfaffianTimesPfaffianState::AdaptAverageMCNorm(int thermalize, int average)
{
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(2*this->NbrParticlesPerLayer);
  this->AdaptNorm(Particles->GetPositions());
  Complex TmpMetropolis, TrialValue = (*this)(Particles->GetPositions());  
  double PreviousSamplingAmplitude = SqrNorm(TrialValue);
  double CurrentSamplingAmplitude = PreviousSamplingAmplitude;
  int NextCoordinates=0;
  // do some MC moves: accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
  for (int i = 0; i < thermalize; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;      
    }
  this->AdaptNorm(Particles->GetPositions());
  double SumTrialValues=0.0;
  for (int i = 0; i < average; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
    }  
  this->ElementNorm*= pow(SumTrialValues/average,(double)-2.0/this->NbrParticlesPerLayer);
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
// assumes that Orbitals1 and Orbitals2 have already been initialized prior to call
void PfaffianTimesPfaffianState::EvaluateTables()
{
  int offset, alpha;  
  Complex tmp, tmp2;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  this->J1=1.0;
  this->J2=1.0;
  if (Orbitals1->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J11[i]=1.0;
	  for(int j=0;j<i;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticlesPerLayer;j++)
	    {
	      tmp= Orbitals1->JastrowFactorElement(i,j);
	      J11[i] *=tmp;
	      J1*=tmp;
	    }
	}
    }
  else // if some interpolation occurred, the true values of the J11's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J11[i]=1.0;
	  for(int j=0;j<i;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	  for(int j=i+1;j<NbrParticlesPerLayer;j++)
	    {
	      tmp = ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	      J11[i] *=tmp;
	      J1*=tmp;
	    }

	}
    }

  if (Orbitals2->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J22[i]=1.0;
	  for(int j=0;j<i;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticlesPerLayer;j++)
	    {
	      tmp = Orbitals2->JastrowFactorElement(i,j);
	      J22[i] *= tmp;
	      J2 *= tmp;
	    }
	}
    }
  else // if some interpolation occurred, the true values of the J22's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J22[i]=1.0;
	  for(int j=0;j<i;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	  for(int j=i+1;j<NbrParticlesPerLayer;j++)
	    {
	      tmp = ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	      J22[i] *= tmp;
	      J2 *= tmp;
	    }
	}
    }

  //inter-layer terms
  this->JInter=1.0;
  for (int i=0; i<this->NbrParticlesPerLayer; ++i)
    for (int j=0; j<this->NbrParticlesPerLayer; ++j)
      {
	this->JInter *= ((Orbitals1->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals1->SpinorV(i)));
      }

      
  // evaluate sums over orbitals m for each LL:
  for (int i=0;i<this->NbrParticlesPerLayer;i++)
    for(int j=0;j<i;j++)
      {
	alpha=0;
	for (int n=0;n<this->NbrLandauLevels;n++)
	  {
	    tmp=0.0;
	    tmp2=0.0;
	    offset=2*n*(n+this->AbsEffectiveFlux+1)+this->AbsEffectiveFlux;    
	    for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	      {
		tmp+=this->fsgn((m2+AbsEffectiveFlux)/2)*OrbitalValues1[alpha][i]*OrbitalValues1[offset-alpha][j];
		tmp2+=this->fsgn((m2+AbsEffectiveFlux)/2)*OrbitalValues2[alpha][i]*OrbitalValues2[offset-alpha][j];
		alpha++;
	      }
	    this->gAlpha1[n][i*this->NbrParticlesPerLayer+j] = tmp;
	    this->gAlpha2[n][i*this->NbrParticlesPerLayer+j] = tmp2;
	  }
      }
}
