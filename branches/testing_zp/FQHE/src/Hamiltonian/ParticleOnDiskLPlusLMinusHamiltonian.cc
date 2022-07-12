////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a disk where         //
//             the hamiltonian is reduced to the L^+ L^- operator             //
//                                                                            //
//                        last modification : 04/07/2008                      //
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


#include "Hamiltonian/ParticleOnDiskLPlusLMinusHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// nbrFluxQuanta = number of flux quanta
// architecture = architecture to use for precalculation
// lFactor = multiplicative factor in front of the L^+L^- operator in the Hamiltonian
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// fixedLz = true if the contribution of the of the Lz^2 has to be computed from the total Lz, false if it has to be computed using the two body operators
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnDiskLPlusLMinusHamiltonian::ParticleOnDiskLPlusLMinusHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrFluxQuanta,
									   AbstractArchitecture* architecture, double lFactor,  long memory, bool onDiskCacheFlag,
									   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = true;
  this->LFactor = lFactor;
  this->Architecture = architecture;
  this->HamiltonianShift = 0;
  this->EvaluateInteractionFactors();  
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  if (precalculationFileName == 0)
    {
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
  this->L2Operator = 0;
}

// destructor
//

ParticleOnDiskLPlusLMinusHamiltonian::~ParticleOnDiskLPlusLMinusHamiltonian() 
{
}

// evaluate all interaction factors
//   

void ParticleOnDiskLPlusLMinusHamiltonian::EvaluateInteractionFactors()
{
  RealMatrix CoefficientsLMinus (this->LzMax + 1, this->LzMax + 1);
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt((double) (i + 1)) * ((double) (this->NbrFluxQuanta - i));
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsLMinus(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt((double) i) * ((double) (this->NbrFluxQuanta + 1 - i));
      //      double TmpCoefficient = sqrt((double) i);
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsLMinus(j, i) *= 2.0 * TmpCoefficient;
    }

  RealMatrix CoefficientsLPlus (this->LzMax + 1, this->LzMax + 1);
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt((double) (i + 1));
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsLPlus(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt((double) i);
      //      double TmpCoefficient = sqrt((double) i);
      for (int j = 0; j <= this->LzMax; ++j)
	CoefficientsLPlus(j, i) *= 0.5 * TmpCoefficient;
    }

  double Factor = this->LFactor;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      Factor *= -1.0;
      this->NbrInteractionFactors = this->LzMax * (this->LzMax - 1) + 1;
    }
  else
    this->NbrInteractionFactors = this->LzMax * (this->LzMax + 1) + 1;
  
  this->M1Value = new int [this->NbrInteractionFactors];
  this->M2Value = new int [this->NbrInteractionFactors];
  this->M3Value = new int [this->NbrInteractionFactors];
  this->InteractionFactors = new double [this->NbrInteractionFactors];
  this->NbrInteractionFactors = 0;


  for (int m3 = 1; m3 <= this->LzMax; ++m3)
    {
      if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	  (m3 != 2))
	{
	  this->InteractionFactors[this->NbrInteractionFactors] = Factor * (CoefficientsLMinus(0, m3) + CoefficientsLPlus(0, m3));
	  this->M1Value[this->NbrInteractionFactors] = m3 - 1;
	  this->M2Value[this->NbrInteractionFactors] = 1;
	  this->M3Value[this->NbrInteractionFactors] = m3;
	  ++this->NbrInteractionFactors;
	}
    }
  for (int m4 = 1; m4 < this->LzMax; ++m4)
    {
      int m3= 1;
      for (; m3 < m4; ++m3)
	{
	  if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	      (m3 != (m4 + 2)))
	    {
	      this->InteractionFactors[this->NbrInteractionFactors] = Factor * (CoefficientsLMinus(m4, m3) + CoefficientsLPlus(m4, m3));
	      this->M1Value[this->NbrInteractionFactors] = m3 - 1;
	      this->M2Value[this->NbrInteractionFactors] = m4 + 1;
	      this->M3Value[this->NbrInteractionFactors] = m3;
	      ++this->NbrInteractionFactors;
	    }
	}
      if (this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic)
	{
	  this->InteractionFactors[this->NbrInteractionFactors] = Factor * (CoefficientsLMinus(m4, m3) + CoefficientsLPlus(m4, m3));
	  this->M1Value[this->NbrInteractionFactors] = m3 - 1;
	  this->M2Value[this->NbrInteractionFactors] = m4 + 1;
	  this->M3Value[this->NbrInteractionFactors] = m3;
	  ++this->NbrInteractionFactors;	  
	}
      ++m3;
      for (; m3 <= this->LzMax; ++m3)
	{
	  if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	      (m3 != (m4 + 2)))
	    {
	      this->InteractionFactors[this->NbrInteractionFactors] = Factor * (CoefficientsLMinus(m4, m3) + CoefficientsLPlus(m4, m3));
	      this->M1Value[this->NbrInteractionFactors] = m3 - 1;
	      this->M2Value[this->NbrInteractionFactors] = m4 + 1;
	      this->M3Value[this->NbrInteractionFactors] = m3;
	      ++this->NbrInteractionFactors;
	    }
	}
    }  
  

//  Factor = 0.0 * this->LFactor;
  Factor = this->LFactor;
   if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
     Factor *= -1.0;
  this->NbrOneBodyInteractionFactors = this->LzMax + 1;
  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
  this->OneBodyInteractionFactors = new double[this->NbrOneBodyInteractionFactors];
  this->OneBodyMValues[0] = 0;
  this->OneBodyNValues[0] = 0;
  //  this->OneBodyInteractionFactors[0] = Factor * ((double) (this->NbrFluxQuanta * this->NbrFluxQuanta));
  this->OneBodyInteractionFactors[0] = Factor * (CoefficientsLMinus(0, 1) + CoefficientsLPlus(1, 0));// + Coefficients(1, 0));
  for (int i = 1; i < this->LzMax; ++i)
    {
      this->OneBodyMValues[i] = i;
      this->OneBodyNValues[i] = i;
//      this->OneBodyInteractionFactors[i] = Factor * ((double) ((i + 1) * (this->NbrFluxQuanta - i) * (this->NbrFluxQuanta - i)));
      this->OneBodyInteractionFactors[i] = Factor * (CoefficientsLMinus(i, i + 1) + CoefficientsLPlus(i - 1, i));
    }
  this->OneBodyMValues[this->LzMax] = this->LzMax;
  this->OneBodyNValues[this->LzMax] = this->LzMax;
//  this->OneBodyInteractionFactors[this->LzMax] = Factor * ((double) ((this->LzMax + 1) * (this->NbrFluxQuanta - this->LzMax) * (this->NbrFluxQuanta - this->LzMax)));
  this->OneBodyInteractionFactors[this->LzMax] = Factor * (CoefficientsLMinus(this->LzMax, this->LzMax - 1) + CoefficientsLPlus(this->LzMax - 1, this->LzMax));
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

