////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of hamiltonian associated to particles on a cylinder with        //
//        SU(2) spin with opposite magnetic field for each species            //
//      a generic interaction defined by its pseudopotential, pairing         //
//               and translation breaking one body potentials                 //
//                                                                            //
//                        last modification : 19/09/2016                      //
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


#include "Hamiltonian/ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// lzmax = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// maxMomentumTransfer = maxixum monentum transfer that can appear in a one body operator
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyOffDiagonalPotentialUpUp = off-diagonal contribution of the one-body potential for particles with spin up, the first entry is the annihilation index, 
//                                   the second entry is the momentum tranfer
// onebodyOffDiagonalPotentialDownDown = off-diagonal contribution of the one-body potential for particles with spin down, the first entry is the annihilation index, 
//                                       the second entry is the momentum tranfer
// onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// onebodyOffDiagonalPotentialPairing = off diagonal contribution of the one-body pairing term, the first entry is the index of the rightmost creation operator, 
//	    				  the second entry is the momentum tranfer
// chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
// averageNumberParticles = average number of particles in the system
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta(ParticleOnSphereWithSpin* particles, int lzmax, double ratio, int maxMomentumTransfer, 
																					     double** pseudoPotential, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
																					     Complex** onebodyOffDiagonalPotentialUpUp, Complex** onebodyOffDiagonalPotentialDownDown,
																					     Complex* onebodyPotentialPairing, Complex** onebodyOffDiagonalPotentialPairing,
																					     double chargingEnergy, double averageNumberParticles, 
																					     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrParticles = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / this->Ratio;
  this->MaximumMomentumTransfer = maxMomentumTransfer;

  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;

  this->NbrPseudopotentialsUpUp = this->LzMax + 1;
  while ((this->NbrPseudopotentialsUpUp > 0) && (pseudoPotential[0][this->NbrPseudopotentialsUpUp - 1] == 0.0))
    {
      --this->NbrPseudopotentialsUpUp;
    }
  if (this->NbrPseudopotentialsUpUp > 0)
    {
      this->PseudopotentialsUpUp = new double[this->NbrPseudopotentialsUpUp];
      for (int i = 0; i < this->NbrPseudopotentialsUpUp; ++i)
	{
	  this->PseudopotentialsUpUp[i] = pseudoPotential[0][i];
	}
    } 
  else
    {
      this->PseudopotentialsUpUp = 0;
    }
  this->NbrPseudopotentialsDownDown = this->LzMax + 1;
  while ((this->NbrPseudopotentialsDownDown > 0) && (pseudoPotential[1][this->NbrPseudopotentialsDownDown - 1] == 0.0))
    {
      --this->NbrPseudopotentialsDownDown;
    }
  if (this->NbrPseudopotentialsDownDown > 0)
    {
      this->PseudopotentialsDownDown = new double[this->NbrPseudopotentialsDownDown];
      for (int i = 0; i < this->NbrPseudopotentialsDownDown; ++i)
	{
	  this->PseudopotentialsDownDown[i] = pseudoPotential[1][i];
	}
    } 
  else
    {
      this->PseudopotentialsDownDown = 0;
    }
  this->NbrPseudopotentialsUpDown = this->LzMax + 1;
  while ((this->NbrPseudopotentialsUpDown > 0) && (pseudoPotential[2][this->NbrPseudopotentialsUpDown - 1] == 0.0))
    {
      --this->NbrPseudopotentialsUpDown;
    }
  if (this->NbrPseudopotentialsUpDown > 0)
    {
      this->PseudopotentialsUpDown = new double[this->NbrPseudopotentialsUpDown];
      for (int i = 0; i < this->NbrPseudopotentialsUpDown; ++i)
	{
	  this->PseudopotentialsUpDown[i] = pseudoPotential[2][i];
	}
    } 
  else
    {
      this->PseudopotentialsUpDown = 0;
    }

  this->ChargingEnergy = chargingEnergy;
  this->EvaluateInteractionFactors();
  this->AverageNumberParticles = averageNumberParticles;
  this->HamiltonianShift =  this->ChargingEnergy * this->AverageNumberParticles * this->AverageNumberParticles;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = new double [this->LzMax + 1];
  for (int i = 0; i <= this->LzMax; ++i)
    this->OneBodyInteractionFactorsupup[i] = this->ChargingEnergy * (1.0 - (2.0 * this->AverageNumberParticles));
  if (onebodyPotentialUpUp != 0)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupup[i] += onebodyPotentialUpUp[i];
    }
  this->OneBodyInteractionFactorsdowndown = new double [this->LzMax + 1];
  for (int i = 0; i <= this->LzMax; ++i)
    this->OneBodyInteractionFactorsdowndown[i] = this->ChargingEnergy * (1.0 - (2.0 * this->AverageNumberParticles));
  if (onebodyPotentialDownDown != 0)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsdowndown[i] += onebodyPotentialDownDown[i];
    }
  this->OneBodyInteractionFactorsupdown = 0;
  double TmpSign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    TmpSign = -1.0;
  this->OneBodyOffDiagonalInteractionFactorsupup = 0;
  if (onebodyOffDiagonalPotentialUpUp != 0)
    {
      this->OneBodyOffDiagonalInteractionFactorsupup = new Complex*[this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyOffDiagonalInteractionFactorsupup[i] = new Complex[this->MaximumMomentumTransfer];
	  for (int j = 0; j < this->MaximumMomentumTransfer; ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsupup[i][j] = TmpSign * onebodyOffDiagonalPotentialUpUp[i][j];
	    }
	}
    }
  this->OneBodyOffDiagonalInteractionFactorsdowndown = 0;
  if (onebodyOffDiagonalPotentialDownDown != 0)
    {
      this->OneBodyOffDiagonalInteractionFactorsdowndown = new Complex*[this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyOffDiagonalInteractionFactorsdowndown[i] = new Complex[this->MaximumMomentumTransfer];
	  for (int j = 0; j < this->MaximumMomentumTransfer; ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsdowndown[i][j] = TmpSign * onebodyOffDiagonalPotentialDownDown[i][j];
	    }
	}
    }
  this->OneBodyInteractionFactorsPairing = 0;
  if (onebodyPotentialPairing != 0)
    {
      this->OneBodyInteractionFactorsPairing = new Complex [this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyInteractionFactorsPairing[i] = onebodyPotentialPairing[i];
	}
    }
  this->OneBodyOffDiagonalInteractionFactorsPairing = 0;
  if (onebodyOffDiagonalPotentialPairing != 0)
    {
      this->OneBodyOffDiagonalInteractionFactorsPairing = new Complex* [this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyOffDiagonalInteractionFactorsPairing[i] = new Complex[2 * this->MaximumMomentumTransfer + 1];
	  for (int j = 0; j < this->MaximumMomentumTransfer; ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsPairing[i][j] = onebodyOffDiagonalPotentialPairing[i][j];
	    }
	  for (int j = this->MaximumMomentumTransfer + 1; j <= (2 * this->MaximumMomentumTransfer); ++j)
	    {
	      this->OneBodyOffDiagonalInteractionFactorsPairing[i][j] = onebodyOffDiagonalPotentialPairing[i][j];
	    }
	}
    }
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
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::~ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta() 
{
  if (this->PseudopotentialsUpUp != 0)
    delete[] this->PseudopotentialsUpUp;
  if (this->PseudopotentialsDownDown != 0)
    delete[] this->PseudopotentialsDownDown;
  if (this->PseudopotentialsUpDown != 0)
    delete[] this->PseudopotentialsUpDown;
}

// evaluate all interaction factors
//   

void ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;

  this->NbrInterSectorSums = 2 * this->LzMax + 1;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrInterSectorIndicesPerSum[m1 + m2];      
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
      this->NbrInterSectorIndicesPerSum[i] = 0;
    }
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      {
	this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = 2 * this->LzMax - 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
	  }

      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp));
		  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
								+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
								- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
								- this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown));
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      if (this->NbrPseudopotentialsUpDown > 0)
	{
	  this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  double Factor = 1.0;
		  int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		      this->InteractionFactorsupdown[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown)
								  - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown));
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	}
      else
	{
	  this->InteractionFactorsupdown = 0;
	  this->NbrInterSectorSums = 0;
	}
    }
  else
    {
      this->NbrIntraSectorSums = 2 * this->LzMax + 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    this->IntraSectorIndicesPerSum[m1 + m2][this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1] = m1;
	    this->IntraSectorIndicesPerSum[m1 + m2][1 + (this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
	  }

      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  if (m1 != m2)
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
									+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp)
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp);
			  this->InteractionFactorsdowndown[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown);
			}
		    }
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      if (this->NbrPseudopotentialsUpDown > 0)
	{
	  this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  double Factor = 2.0;
		  int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		      this->InteractionFactorsupdown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown)
								  + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown));
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	}
      else
	{
	  this->InteractionFactorsupdown = 0;
	  this->NbrInterSectorSums = 0;
	}
    }


 cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
 cout << "====================================" << endl;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// return value = numerical coefficient

double ParticleOnCylinderWithSpinTimeReversalSymmetricGenericHamiltonianAndPairingAllMomenta::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials)
{
  double Kappa2Factor = (2.0 * M_PI * this->Ratio / ((double) (this->LzMax + 1)));
  double Coefficient = pseudopotentials[0] * exp (-0.25 * Kappa2Factor * (((double) ((m1 - m2) * (m1 - m2))) + ((double) ((m3 - m4) * (m3 - m4))))) / sqrt(2.0 * M_PI);
  if (nbrPseudopotentials > 1)
    {
      Coefficient += pseudopotentials[1] * (((double) ((m1 - m2) * (m4 - m3)))
					    * exp (-0.25 * Kappa2Factor * (((double) ((m1 - m2) * (m1 - m2))) 
									  + ((double) ((m3 - m4) * (m3 - m4)))))) * Kappa2Factor / sqrt(this->Ratio * ((double) (this->LzMax + 1)));
    }
  return Coefficient;
}

