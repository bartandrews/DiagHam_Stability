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


#include "AbstractMCSamplingFunction.h"

#include <iostream>
#include <cmath>
using std::cout;
using std::endl;


// virtual destructor
AbstractMCSamplingFunction::~AbstractMCSamplingFunction()
{
}

// register basic system of particles
void AbstractMCSamplingFunction::RegisterSystem(AbstractParticleCollection *system)
{
  cout << "Attention: AbstractMCSamplingFunction::RegisterSystem - should be overridden in derived classes!"<<endl;
}



// set function value to one for present particle positions
void AbstractMCSamplingFunction::AdaptNorm()
{
  int countdown=100;
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

// signal that the last move was accepted
void AbstractMCSamplingFunction::AcceptedMove()
{
}

// set function value for a typical average of particle positions in MC sampling
void AbstractMCSamplingFunction::AdaptAverageMCNorm(int thermalize, int average)
{
  this->AdaptNorm();
  double acceptanceProbability;
  // do some MC moves: accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
  for (int i = 0; i < thermalize; ++i)
    {
      this->GetSystem()->Move();
      acceptanceProbability = this->GetTransitionRatio();
      if ((acceptanceProbability > 1.0) ||  (this->GetSystem()->GetRandomNumber() < acceptanceProbability))
	{
	  // nothing to do for now...
	  this->AcceptedMove();
	}
      else
	{
	  this->GetSystem()->RestoreMove();
	}
    }
  this->AdaptNorm();
  double SumSamplingValues=0.0;
  for (int i = 0; i < average; ++i)
    {
      this->GetSystem()->Move();
      acceptanceProbability = this->GetTransitionRatio();
      if ((acceptanceProbability > 1.0) ||  (this->GetSystem()->GetRandomNumber() < acceptanceProbability))
	{
	  // nothing to do for now...
	  this->AcceptedMove();
	}
      else
	{
	  this->GetSystem()->RestoreMove();
	}            
      SumSamplingValues+=Norm(this->GetFunctionValue());
    }
  SumSamplingValues/=average;
  this->ScaleByFactor(1.0/SumSamplingValues);
}
