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


#include "TrivialSamplingFunction.h"
#include "ParticleOnSphereCollection.h"
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;


// constructor
// waveFunction = wavefunction to be sampled from
TrivialSamplingFunction::TrivialSamplingFunction(Abstract1DComplexFunction *waveFunction)
{
  this->WaveFunction=waveFunction;
  this->LastAmplitude=1e-300;
  this->System=NULL;
}
  

// virtual destructor
TrivialSamplingFunction::~TrivialSamplingFunction()
{
}


// register basic system of particles
// this function needs to be called before any of the other routines are functional
void TrivialSamplingFunction::RegisterSystem(AbstractParticleCollection *system)
{  
  this->System=system;
  if (system->GetCollectionType()==AbstractParticleCollection::OnSphereCollection)
    {
      this->Positions=((ParticleOnSphereCollection*)System)->GetPositions();
    }
  else
    {
      cout << "Unknown type of particle collection in TrivialSamplingFunction";
    }
      
  // pointers to spinor coordinates (external)
  // ((ParticleOnSphereCollection*)System)->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}


// signal that the last move was accepted
void TrivialSamplingFunction::AcceptedMove()
{
  this->LastAmplitude=this->TentativeNewValue;
}


// method for ratio of probabilities with respect to the last configuration
// allows for more rapid calculation due to cancellation of factors
double TrivialSamplingFunction::GetTransitionRatio()
{
  this->TentativeNewValue = SqrNorm((*WaveFunction)(Positions));
  //  cout << "Positions at present time:"<<endl<<Positions<<endl;
  //  cout << "Evaluated wavefunction: "<<TentativeNewValue<<endl;
  double Result=TentativeNewValue/LastAmplitude;  
  return Result;
}


// get the full function value for a system of particles
Complex TrivialSamplingFunction::GetFunctionValue()
{
  return (*WaveFunction)(Positions);
}


// call this method to scale the sampling function (needed to normalize the function)
// scale = total scaling factor
void TrivialSamplingFunction::ScaleByFactor(double scale)
{
}
