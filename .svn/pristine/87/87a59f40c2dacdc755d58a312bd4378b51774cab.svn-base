////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of QHE particle wave function evaluation operation         //
//                                                                            //
//                        last modification : 29/07/2004                      //
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
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


// constructor 
//
// space = pointer to the Hilbert space to use
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = indicate which coordinates will be changed during next time step (-1 if no time coherence has to be used)

QHEParticleWaveFunctionOperation::QHEParticleWaveFunctionOperation (AbstractQHEParticle* space, RealVector* state, RealVector* position, 
								    AbstractFunctionBasis* basis, int nextCoordinates)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->NextCoordinates = nextCoordinates;
  if (this->NextCoordinates != -1)
    space->InitializeWaveFunctionEvaluation(true);
  else
    space->InitializeWaveFunctionEvaluation(false);
  this->HilbertSpace = (AbstractQHEParticle*) space->Clone();
  this->State = state;
  this->Position = position;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Basis = basis;    
  this->Scalars = 0;
  this->States = 0;
  this->NbrScalars = 1;
}

// constructor for multiple wave function evaluations
//
// space = pointer to the Hilbert space to use
// states = array of vectors corresponding to the states in the Fock basis
// nbrStates = number of states in the states array
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = indicate which coordinates will be change during next time step (-1 if no time coherence has to be used)

QHEParticleWaveFunctionOperation::QHEParticleWaveFunctionOperation(AbstractQHEParticle* space, RealVector* states, int nbrStates, RealVector* position, 
								   AbstractFunctionBasis* basis, int nextCoordinates)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->NextCoordinates = nextCoordinates;
  if (this->NextCoordinates != -1)
    space->InitializeWaveFunctionEvaluation(true);
  else
    space->InitializeWaveFunctionEvaluation(false);
  this->HilbertSpace = (AbstractQHEParticle*) space->Clone();
  this->State = 0;
  this->Position = position;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Basis = basis;    
  this->NbrScalars = nbrStates;
  this->Scalars = new Complex[this->NbrScalars];
  this->States = new RealVector[this->NbrScalars];
  for (int i = 0; i < this->NbrScalars; ++i)
    {
      this->Scalars[i] = 0.0;
      this->States[i] = states[i];
    }
}

// copy constructor 
//
// operation = reference on operation to copy

QHEParticleWaveFunctionOperation::QHEParticleWaveFunctionOperation(const QHEParticleWaveFunctionOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->State = operation.State;
  this->HilbertSpace = (AbstractQHEParticle*) operation.HilbertSpace->Clone();
  this->Position = operation.Position;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Basis = operation.Basis;
  this->Scalar = operation.Scalar;
  this->NextCoordinates = operation.NextCoordinates;
  this->NbrScalars = operation.NbrScalars;
  if (this->NbrScalars > 1)
    {
      this->Scalars = new Complex[this->NbrScalars];
      this->States = new RealVector[this->NbrScalars];
      for (int i = 0; i < this->NbrScalars; ++i)
	{
	  this->Scalars[i] = operation.Scalars[i];
	  this->States[i] = operation.States[i];
	}
    }
  else
    {
      this->Scalars = 0;
      this->States = 0;
    }
}
  
// destructor
//

QHEParticleWaveFunctionOperation::~QHEParticleWaveFunctionOperation()
{
  if (this->NbrScalars > 1)
    {
      delete[] this->Scalars;
      delete[] this->States;
    }
  delete this->HilbertSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void QHEParticleWaveFunctionOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* QHEParticleWaveFunctionOperation::Clone()
{
  return new QHEParticleWaveFunctionOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool QHEParticleWaveFunctionOperation::RawApplyOperation()
{
  if (this->NextCoordinates == -1)
    if (this->NbrScalars <= 1)
      {
	this->Scalar = this->HilbertSpace->EvaluateWaveFunction(*(this->State), *(this->Position), *(this->Basis),
								this->FirstComponent, this->NbrComponent);
      }
    else
      {
	this->HilbertSpace->EvaluateWaveFunctions(this->States, this->NbrScalars, *(this->Position), *(this->Basis),
						  this->Scalars, this->FirstComponent, this->NbrComponent);
      }
  else
    this->Scalar = this->HilbertSpace->EvaluateWaveFunctionWithTimeCoherence(*(this->State), *(this->Position), 
									     *(this->Basis), this->NextCoordinates, 
									     this->FirstComponent, this->NbrComponent);
  return true;
}

