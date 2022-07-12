////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(2) spin on a sphere           //
//                including all possible two-body coupling terms              //
//                                                                            //
//                        last modification : 12/06/2009                      //
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
#include "Hamiltonian/AbstractQHEOnSphereWithSpinFullHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <fstream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

AbstractQHEOnSphereWithSpinFullHamiltonian::~AbstractQHEOnSphereWithSpinFullHamiltonian()
{
  if (this->InteractionFactorsUpUpUpUp != 0)
    {
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	delete[] this->InteractionFactorsUpUpUpUp[i];
      delete[] this->InteractionFactorsUpUpUpUp;
    }
  if (this->InteractionFactorsUpDownUpUp != 0)
    {
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	delete[] this->InteractionFactorsUpDownUpUp[i];
      delete[] this->InteractionFactorsUpDownUpUp;
    }
  if (this->InteractionFactorsDownDownUpUp != 0)
    {
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	delete[] this->InteractionFactorsDownDownUpUp[i];
      delete[] this->InteractionFactorsDownDownUpUp;
    }
  if (this->InteractionFactorsUpUpUpDown != 0)
    {
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	delete[] this->InteractionFactorsUpUpUpDown[i];
      delete[] this->InteractionFactorsUpUpUpDown;
    }
  if (this->InteractionFactorsUpDownUpDown != 0)
    {
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	delete[] this->InteractionFactorsUpDownUpDown[i];
      delete[] this->InteractionFactorsUpDownUpDown;
    }
  if (this->InteractionFactorsDownDownUpDown != 0)
    {
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	delete[] this->InteractionFactorsDownDownUpDown[i];
      delete[] this->InteractionFactorsDownDownUpDown;
    }
  if (this->InteractionFactorsUpUpDownDown != 0)
    {
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	delete[] this->InteractionFactorsUpUpDownDown[i];
      delete[] this->InteractionFactorsUpUpDownDown;
    }
  if (this->InteractionFactorsUpDownDownDown != 0)
    {
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	delete[] this->InteractionFactorsUpDownDownDown[i];
      delete[] this->InteractionFactorsUpDownDownDown;
    }
  if (this->InteractionFactorsDownDownDownDown != 0)
    {
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	delete[] this->InteractionFactorsDownDownDownDown[i];
      delete[] this->InteractionFactorsDownDownDownDown;
    }
  
  if (this->NbrUpUpSectorSums > 0)
    {
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	if ( this->UpUpSectorIndicesPerSum[i] != 0 ) delete[] this->UpUpSectorIndicesPerSum[i];
      delete[] this->UpUpSectorIndicesPerSum;
      delete[] this->NbrUpUpSectorIndicesPerSum;
    }
  if (this->NbrUpDownSectorSums > 0)
    {
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	if ( this->UpDownSectorIndicesPerSum[i] != 0 ) delete[] this->UpDownSectorIndicesPerSum[i];
      delete[] this->UpDownSectorIndicesPerSum;
      delete[] this->NbrUpDownSectorIndicesPerSum;
    }
  if (this->NbrDownDownSectorSums > 0)
    {
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	if ( this->DownDownSectorIndicesPerSum[i] != 0 ) delete[] this->DownDownSectorIndicesPerSum[i];
      delete[] this->DownDownSectorIndicesPerSum;
      delete[] this->NbrDownDownSectorIndicesPerSum;
    }

  if (this->OneBodyInteractionFactorsUpUp != 0)
    {
      delete[] this->OneBodyInteractionFactorsUpUp;
      delete[] this->OneBodyMValuesUpUp;
    }
  if (this->OneBodyInteractionFactorsUpDown != 0)
    {
      delete[] this->OneBodyInteractionFactorsUpDown;
      delete[] this->OneBodyMValuesUpDown;
    }
  if (this->OneBodyInteractionFactorsDownUp != 0)
    {
      delete[] this->OneBodyInteractionFactorsDownUp;
      delete[] this->OneBodyMValuesDownUp;
    }
  if (this->OneBodyInteractionFactorsDownDown != 0)
    {
      delete[] this->OneBodyInteractionFactorsDownDown;
      delete[] this->OneBodyMValuesDownDown;
    }
}


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereWithSpinFullHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  delete TmpParticles;
  return Memory;
}

