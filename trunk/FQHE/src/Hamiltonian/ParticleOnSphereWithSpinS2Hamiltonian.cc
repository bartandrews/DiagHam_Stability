////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with spin   //
//  where the hamiltonian is reduced to a simple total square spin momentum   //
//                                                                            //
//                        last modification : 27/07/2007                      //
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


#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "GeneralTools/StringTools.h"

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
// totalLz = twice the projected momentum total value
// totalSz = twice the projected spin total value
// architecture = architecture to use for precalculation
// s2Factor = multiplicative factor in front of the S^2 operator in the Hamiltonian
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// fixedLz = flag indicating whether Sz needs to be evaluated

ParticleOnSphereWithSpinS2Hamiltonian::ParticleOnSphereWithSpinS2Hamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, int totalLz, int totalSz,
									     AbstractArchitecture* architecture, double s2Factor, long memory, bool onDiskCacheFlag,
									     char* precalculationFileName, bool fixedSz)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->TotalLz = totalLz;
  this->TotalSz = totalSz;
  this->NbrLzValue = this->LzMax + 1;
  this->FixedSz = fixedSz;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = true;
  this->S2Factor = s2Factor;
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  if (this->FixedSz)
    this->HamiltonianShift = this->S2Factor * (0.25 * ((double) (this->TotalSz * this->TotalSz)) + 0.5 * ((double) this->NbrParticles));
  else
    this->HamiltonianShift = this->S2Factor * 0.5 * ((double) this->NbrParticles);
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout  << "fast = ";
	  PrintMemorySize(cout, TmpMemory)<<endl;
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
}

// destructor
//

ParticleOnSphereWithSpinS2Hamiltonian::~ParticleOnSphereWithSpinS2Hamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinS2Hamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinS2Hamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinS2Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereWithSpinS2Hamiltonian::ShiftHamiltonian (double shift)
{
  if (this->FixedSz == true)
    this->HamiltonianShift = shift + this->S2Factor * (0.25 * ((double) (this->TotalSz * this->TotalSz)) + 0.5 * ((double) this->NbrParticles));
  else
    this->HamiltonianShift = shift + this->S2Factor * 0.5 * ((double) this->NbrParticles);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinS2Hamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinS2Hamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereWithSpinS2Hamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereWithSpinS2Hamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereWithSpinS2Hamiltonian::EvaluateInteractionFactors()
{
  double Factor = this->S2Factor;
  double Factor2 = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphereWithSpin::FermionicStatistic)
    Factor2 *= -1.0;
  

  if (this->FixedSz==false)
    {
      // have intra - indices
      this->NbrM12IntraIndices = this->LzMax * (this->LzMax + 1) / 2;
      if (this->Particles->GetParticleStatistic() != ParticleOnSphereWithSpin::FermionicStatistic)
	this->NbrM12IntraIndices += this->LzMax + 1;
      this->M1IntraValue = new int [this->NbrM12IntraIndices];
      this->M2IntraValue = new int [this->NbrM12IntraIndices];
      this->M3IntraValues = new int* [this->NbrM12IntraIndices];
      this->NbrM3IntraValues  = new int [this->NbrM12IntraIndices];
      this->M12InteractionFactorsupup  = new double [this->NbrM12IntraIndices];
      this->M12InteractionFactorsdowndown  = new double [this->NbrM12IntraIndices];
      // number of inter-indices - twice as many for this case!
      this->NbrM12InterIndices = 2 * (this->LzMax + 1) * (this->LzMax + 1);
      
    }
  else
    {
      // no intra - indices
      this->NbrM12IntraIndices = 0;
      this->M1IntraValue = 0;
      this->M2IntraValue = 0;
      this->M3IntraValues = 0;
      this->NbrM3IntraValues  = 0;
      this->M12InteractionFactorsupup = 0;
      this->M12InteractionFactorsdowndown = 0;
      // number of inter-indices
      this->NbrM12InterIndices = (this->LzMax + 1) * (this->LzMax + 1);
    }
  
  this->M1InterValue = new int [this->NbrM12InterIndices];
  this->M2InterValue = new int [this->NbrM12InterIndices];
  this->M3InterValues = new int* [this->NbrM12InterIndices];
  this->NbrM3InterValues  = new int [this->NbrM12InterIndices];
  this->M12InteractionFactorsupdown = new double [this->NbrM12InterIndices];

  // reset counters (will be incremented during initialization)
  this->NbrM12InterIndices = 0;
  this->NbrM12IntraIndices = 0;
    
  for (int m3 = 0; m3 <= this->LzMax; ++m3)
    for (int m4 = 0; m4 <= this->LzMax; ++m4)
      {
	this->M12InteractionFactorsupdown[this->NbrM12InterIndices] = Factor;
	this->M1InterValue[this->NbrM12InterIndices] = m3;
	this->M2InterValue[this->NbrM12InterIndices] = m4;
	this->M3InterValues[this->NbrM12InterIndices] = new int [1];
	this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
	this->M3InterValues[this->NbrM12InterIndices][0] =  m4;
	++this->NbrM12InterIndices;
      }

  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  
  if (this->FixedSz==false)
    {
      for (int m3 = 1; m3 <= this->LzMax; ++m3)
	for (int m4 = 0; m4 < m3; ++m4)
	  {
	    this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = 0.5*Factor*Factor2;
	    this->M12InteractionFactorsdowndown[this->NbrM12IntraIndices] = 0.5*Factor*Factor2;
	    this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	    this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	    this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	    this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	    this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3;
	    ++this->NbrM12IntraIndices;
	  }
      // diagonal terms: one and two-body
      if (this->Particles->GetParticleStatistic() != ParticleOnSphereWithSpin::FermionicStatistic)
	{
	  for (int m3 = 0; m3 <= this->LzMax; ++m3)
	    {
	      this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = 0.25*Factor*Factor2;
	      this->M12InteractionFactorsdowndown[this->NbrM12IntraIndices] = 0.25*Factor*Factor2;
	      this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	      this->M2IntraValue[this->NbrM12IntraIndices] = m3;
	      this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	      this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	      this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3;
	      ++this->NbrM12IntraIndices;
	    }
	  this->OneBodyInteractionFactorsupup = new double[this->LzMax+1];
	  this->OneBodyInteractionFactorsdowndown = new double[this->LzMax+1];
	  for (int m=0; m<=this->LzMax; ++m)
	    {
	      this->OneBodyInteractionFactorsupup[m]=0.25*Factor;
	      this->OneBodyInteractionFactorsdowndown[m]=0.25*Factor;
	    }
	}
      for (int m3 = 0; m3 <= this->LzMax; ++m3)
	for (int m4 = 0; m4 <= this->LzMax; ++m4)
	  {
	    this->M12InteractionFactorsupdown[this->NbrM12InterIndices] = -0.5*Factor*Factor2;
	    this->M1InterValue[this->NbrM12InterIndices] = m3;
	    this->M2InterValue[this->NbrM12InterIndices] = m4;
	    this->M3InterValues[this->NbrM12InterIndices] = new int [1];
	    this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
	    this->M3InterValues[this->NbrM12InterIndices][0] =  m3;
	    ++this->NbrM12InterIndices;
	  }
    }

  cout << "nbr interaction = " << this->NbrM12InterIndices << endl;
  cout << "====================================" << endl;
}

