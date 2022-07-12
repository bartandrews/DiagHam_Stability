////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 23/01/2008                      //
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


#include "HalperinSamplingFunction.h"
#include "ParticleOnSphereCollection.h"
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;


// constructor
// nbrParticles = number of particles in system
// exponentK = power raised for up-spins
// exponentL = power raised for down-spins
// exponentM = power raised for inter-spin factors
// nbrUp = number of particles with up spin (negative number -> choose half the particles)
HalperinSamplingFunction::HalperinSamplingFunction(int nbrParticles, int exponentK, int exponentL, int exponentM, int nbrUp)
{
  if ((nbrUp<0)&&(nbrParticles&1))
     {
       cout << "HalperinSamplingFunction sampling function requires an even number of particles, or else give number of up particles" << endl;
       exit(1);
     }
  this->NbrParticles = nbrParticles;
  if (nbrUp<0)
    this->NbrUp=nbrParticles/2;
  else
    this->NbrUp=nbrUp;
  this->Exponent_K=exponentK;
  this->Exponent_L=exponentL;
  this->Exponent_M=exponentM;
  this->System=NULL;
  this->ElementNorm=1.0;
}
  

// virtual destructor
HalperinSamplingFunction::~HalperinSamplingFunction()
{
}


// register basic system of particles
// this function needs to be called before any of the other routines are functional
void HalperinSamplingFunction::RegisterSystem(AbstractParticleCollection *system)
{
  this->System=system;
  if (((ParticleOnSphereCollection*)System)->GetNbrParticles() != this->NbrParticles)
    {
      cout << "Number of particles in system not compatible in sampling function";
      exit(1);
    }
  // pointers to spinor coordinates (external)
  ((ParticleOnSphereCollection*)System)->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}



// method for ratio of probabilities with respect to the last configuration
// allows for more rapid calculation due to cancellation of factors
double HalperinSamplingFunction::GetTransitionRatio()
{
  double ratio1=1.0, ratio2=1.0;
  int tomove = System->GetMovedNbr();
  ((ParticleOnSphereCollection*)System)->GetPreviousPos(LastU,LastV);
  if (tomove<this->NbrUp)
    {
      // intra-spin terms:
      for (int i=0;i<tomove;i++)
	{
	  ratio1 *= Norm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	    Norm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
	}
      for (int i=tomove+1;i<this->NbrUp;i++)
	{
	  ratio1 *= Norm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	    Norm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
	}
      ratio1 = pow( ratio1, 2.0*this->Exponent_K);

      // inter-spin terms:
      for (int i=this->NbrUp;i<this->NbrParticles;i++)
	{
	  ratio2 *= Norm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	    Norm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
	}
      ratio2 = pow( ratio2, 2.0*this->Exponent_M);
      
    }
  else
    {
      // intra-spin terms:
      for (int i=this->NbrUp;i<tomove;i++)
	{
	  ratio1 *= Norm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	    Norm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
	}
      for (int i=tomove+1;i<this->NbrParticles;i++)
	{
	  ratio1 *= Norm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	    Norm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
	}
      ratio1 = pow( ratio1, 2.0*this->Exponent_L);
      
      // inter-spin terms:
      for (int i=0;i<this->NbrUp;i++)
	{
	  ratio2 *= Norm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	    Norm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
	}
      ratio2 = pow( ratio2, 2.0*this->Exponent_M);
    }  
  return ratio1*ratio2;
}


// get the full function value for a system of particles
Complex HalperinSamplingFunction::GetFunctionValue()
{
  Complex Term1=1.0, Term2=1.0, Term3=1.0;
  for (int j=1; j<this->NbrUp; ++j)
    for (int i=0;i<j;i++)
      Term1 *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  Complex Result=Term1;
  for (int i=1; i<this->Exponent_K; ++i)
    Result *= Term1;
  for (int j=this->NbrUp+1; j<this->NbrParticles; ++j)
    for (int i=this->NbrUp;i<j;i++)
      Term2 *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  for (int i=0; i<this->Exponent_L; ++i)
    Result *= Term2;
  for (int j=0; j<this->NbrUp; ++j)
    for (int i=this->NbrUp; i<this->NbrParticles; i++)
      Term3 *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  for (int i=0; i<this->Exponent_M; ++i)
    Result *= Term3;
  return Result;
}


// call this method to scale the sampling function (needed to normalize the function)
// scale = total scaling factor
void HalperinSamplingFunction::ScaleByFactor(double scale)
{
  double NbrDown = NbrParticles-NbrUp;
  double Factors = (double)NbrUp*(NbrUp-1)*this->Exponent_K/2.0;
  Factors += NbrDown*(NbrDown-1)*this->Exponent_L/2.0;
  Factors += (double)NbrUp*NbrDown*this->Exponent_M;
  this->ElementNorm *= std::pow(scale,1.0/Factors);
}
