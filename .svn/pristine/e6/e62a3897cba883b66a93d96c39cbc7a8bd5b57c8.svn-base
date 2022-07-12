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


#include "LaughlinSamplingFunction.h"
#include "ParticleOnSphereCollection.h"
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;


// constructor
LaughlinSamplingFunction::LaughlinSamplingFunction(int nbrParticles, int exponent)
{
  this->NbrParticles=nbrParticles;
  this->Exponent=exponent;
  this->System=NULL;
  this->ElementNorm=1.0;
}
  

// virtual destructor
LaughlinSamplingFunction::~LaughlinSamplingFunction()
{
}


// register basic system of particles
// this function needs to be called before any of the other routines are functional
void LaughlinSamplingFunction::RegisterSystem(AbstractParticleCollection *system)
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
double LaughlinSamplingFunction::GetTransitionRatio()
{
  double ratio=1.0;
  int tomove = System->GetMovedNbr();
  ((ParticleOnSphereCollection*)System)->GetPreviousPos(LastU,LastV);
  for (int i=0;i<tomove;i++)
    {
      ratio *= SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	SqrNorm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
    }
  for (int i=tomove+1;i<this->NbrParticles;i++)
    {
      ratio *= SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	SqrNorm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
    }
  double Base=ratio;
  for (int i=1; i<this->Exponent; ++i)
    ratio *= Base;
  return ratio;
}


// get the full function value for a system of particles
Complex LaughlinSamplingFunction::GetFunctionValue()
{
  Complex Result=1.0;
  for (int j=1; j<this->NbrParticles; ++j)
    for (int i=0;i<j;i++)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  Complex Base=Result;
  for (int i=1; i<this->Exponent; ++i)
    Result *= Base;
  return Result;
}


// call this method to scale the sampling function (needed to normalize the function)
// scale = total scaling factor
void LaughlinSamplingFunction::ScaleByFactor(double scale)
{
  double factors = (double)NbrParticles*(NbrParticles-1)*this->Exponent/2.0;
  this->ElementNorm *= std::pow(scale,1.0/factors);
}
