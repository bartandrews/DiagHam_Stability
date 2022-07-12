////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2016 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 28/07/2016                      //
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


#include "LaughlinSamplingFunctionOnDisk.h"
#include "AbstractParticleCollectionOnDisk.h"
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

// constructor
LaughlinSamplingFunctionOnDisk::LaughlinSamplingFunctionOnDisk(int nbrParticles, int exponent, double defectAngle, double spin)
{
  this->NbrParticles=nbrParticles;
  this->Exponent=exponent;
  if (defectAngle==1.0)
    {
      std::cout << "Defect angle has to be between 0 and 1"<<std::endl;
      exit(1);
    }
  this->DefectAngle = defectAngle;
  this->Gamma = 1.0 - this->DefectAngle;
  this->InvGammaSqr = 1.0/(this->Gamma * this->Gamma);
  this->Spin = spin;
  this->System=NULL;
  this->ElementNorm=1.0;
  this->LogScale = 0.0;
}
  

// virtual destructor
LaughlinSamplingFunctionOnDisk::~LaughlinSamplingFunctionOnDisk()
{
}


// register basic system of particles
// this function needs to be called before any of the other routines are functional
void LaughlinSamplingFunctionOnDisk::RegisterSystem(AbstractParticleCollectionOnDisk *system)
{
  this->System=system;
  if (System->GetNbrParticles() != this->NbrParticles)
    {
      cout << "Number of particles in system not compatible in sampling function";
      exit(1);
    }
  // pointers to spinor coordinates (external)
  System->GetCoordinates(CoordinatesZ);
}



// method for ratio of probabilities with respect to the last configuration
// allows for more rapid calculation due to cancellation of factors
double LaughlinSamplingFunctionOnDisk::GetTransitionRatio()
{
  double ratio=1.0;
  int tomove = System->GetMovedNbr();
  System->GetPreviousPos(LastZ);
  for (int i=0;i<tomove;i++)
    {
      ratio *= SqrNorm(CoordinatesZ[i]-CoordinatesZ[tomove])/SqrNorm(CoordinatesZ[i]-LastZ);
    }
  for (int i=tomove+1;i<this->NbrParticles;i++)
    {
      ratio *= SqrNorm(CoordinatesZ[i]-CoordinatesZ[tomove])/SqrNorm(CoordinatesZ[i]-LastZ);
    }
  double Base=ratio;
  for (int i=1; i<this->Exponent; ++i)
    ratio *= Base;
  if (this->DefectAngle!=0.0)
    {
      if (this->Spin!=0.0)
	ratio *= std::pow(SqrNorm(CoordinatesZ[tomove])/SqrNorm(LastZ), this->DefectAngle*this->Spin);
      ratio *= std::exp(this->InvGammaSqr*( -std::pow(SqrNorm(CoordinatesZ[tomove]), this->Gamma) + std::pow(SqrNorm(LastZ), this->Gamma)));
    }
  else
    ratio *= std::exp(-SqrNorm(CoordinatesZ[tomove])+SqrNorm(LastZ));
  return ratio;
}


// get the full function value for a system of particles
Complex LaughlinSamplingFunctionOnDisk::GetFunctionValue()
{
  Complex Result=1.0;
  for (int j=1; j<this->NbrParticles; ++j)
    for (int i=0;i<j;i++)
      Result *= this->ElementNorm*(CoordinatesZ[i]-CoordinatesZ[j]);
  Complex Base=Result;
  for (int i=1; i<this->Exponent; ++i)
    Result *= Base;
  if (this->DefectAngle!=0.0)
    {
      if (this->Spin!=0.0)
	for (int i = 0; i < this->NbrParticles; ++i)
	  Result *= std::pow(SqrNorm(CoordinatesZ[i]), 0.5*this->DefectAngle*this->Spin);
      double SumExp = 0.0;
      for (int i = 0; i < this->NbrParticles; ++i)
	SumExp += std::pow(SqrNorm(CoordinatesZ[i]), 0.5*this->Gamma);
      Result *= std::exp(-0.5*this->InvGammaSqr*SumExp + this->LogScale);
    }
  else
    {
      double SumSqr=0.0;
      for (int i = 0; i < this->NbrParticles; ++i)
	SumSqr += SqrNorm(CoordinatesZ[i]);
      Result *= std::exp(-0.5*SumSqr+this->LogScale);
    }
  return Result;
}


// call this method to scale the sampling function (needed to normalize the function)
// scale = total scaling factor
void LaughlinSamplingFunctionOnDisk::ScaleByFactor(double scale)
{
  double factors = (double)NbrParticles*(NbrParticles-1)*this->Exponent/2.0;
  this->ElementNorm *= pow(scale,1.0/factors);
}



// set function value to one for present particle positions
void LaughlinSamplingFunctionOnDisk::AdaptNorm()
{
  int countdown=100;
  double SumSqr=0.0;
  for (int i = 0; i < this->NbrParticles; ++i)
    SumSqr += SqrNorm(CoordinatesZ[i]);
  this->LogScale = 0.5*SumSqr; // assume this is a typical configuration, so set the exponential factor to one
  double norm=Norm(this->GetFunctionValue());
  while ((((norm<.1)||(norm>10.0))||(std::isnan((double)norm)))&&(countdown-- > 0))
    {
      if ((norm>1e300)||std::isnan((double)norm))
	this->ScaleByFactor(1e-300);
      else if (norm==0.0)
	this->ScaleByFactor(1e300);
      else
	this->ScaleByFactor(1.0/norm);
      norm=Norm(this->GetFunctionValue());
    }
  if (countdown <= 0)
    cout << "Problem with scaling of sampling function"<<endl; 
}
