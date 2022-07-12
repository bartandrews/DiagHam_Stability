////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with spin   //
// where the hamiltonian is reduced to a simple total square angular momentum //
//                                                                            //
//                        last modification : 06/07/2007                      //
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


#include "Hamiltonian/ParticleOnSphereWithSU4SpinL2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "Hamiltonian/ParticleOnSphereWithSU4SpinS2Hamiltonian.h"
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
// architecture = architecture to use for precalculation
// l2Factor = multiplicative factor in front of the L^2 operator in the Hamiltonian
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSU4SpinL2Hamiltonian::ParticleOnSphereWithSU4SpinL2Hamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int lzmax, int totalLz,
									     AbstractArchitecture* architecture, double l2Factor, long memory, bool onDiskCacheFlag,
									     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->TotalLz = totalLz;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = true;
  this->L2Factor = l2Factor;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.25 * this->L2Factor * ((double) (this->TotalLz * this->TotalLz));
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

	  PrintMemorySize(cout,TmpMemory);

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

ParticleOnSphereWithSU4SpinL2Hamiltonian::~ParticleOnSphereWithSU4SpinL2Hamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSU4SpinL2Hamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphereWithSU4Spin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSU4SpinL2Hamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSU4SpinL2Hamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereWithSU4SpinL2Hamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift + 0.25 * this->L2Factor * ((double) (this->TotalLz * this->TotalLz));
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSU4SpinL2Hamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex ParticleOnSphereWithSU4SpinL2Hamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereWithSU4SpinL2Hamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereWithSU4SpinL2Hamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereWithSU4SpinL2Hamiltonian::EvaluateInteractionFactors()
{
  // no entangled terms mixing spin and pseudospin
  this->M12InteractionFactorsupdmEnt = NULL;
  this->M12InteractionFactorsumdpEnt = NULL;
  
  RealMatrix Coefficients (this->LzMax + 1, this->LzMax + 1);
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax + 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	Coefficients(i, j) = TmpCoefficient;
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      double TmpCoefficient = sqrt(0.25 * ((double) ((((this->LzMax + 2) * this->LzMax) - (((2 * i) - this->LzMax) * ((2 * i) - this->LzMax - 2))))));
      for (int j = 0; j <= this->LzMax; ++j)
	Coefficients(j, i) *= 0.5 * TmpCoefficient;
    }

  // cout << "Coefficients="<<endl<<Coefficients;

  //  this->L2Factor  = 1.0;

  double Factor = 2.0*this->L2Factor;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      Factor *= -1.0;
      this->NbrM12IntraIndices = this->LzMax * (this->LzMax - 1) + 1;
    }
  else this->NbrM12IntraIndices = this->LzMax * this->LzMax + 1;
  this->M1IntraValue = new int [this->NbrM12IntraIndices];
  this->M2IntraValue = new int [this->NbrM12IntraIndices];
  this->M3IntraValues = new int* [this->NbrM12IntraIndices];
  this->NbrM3IntraValues  = new int [this->NbrM12IntraIndices];
  this->M12InteractionFactorsupup = new double [this->NbrM12IntraIndices];
  this->M12InteractionFactorsumum = new double [this->NbrM12IntraIndices];
  this->M12InteractionFactorsdpdp = new double [this->NbrM12IntraIndices];
  this->M12InteractionFactorsdmdm = new double [this->NbrM12IntraIndices];
  
  this->NbrM12IntraIndices = 0;

  for (int m3 = 1; m3 <= this->LzMax; ++m3)
    {
      if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	  (m3 != 2))
	{
	  this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * Coefficients(0, m3);
	  this->M12InteractionFactorsumum[this->NbrM12IntraIndices] = Factor * Coefficients(0, m3);
	  this->M12InteractionFactorsdpdp[this->NbrM12IntraIndices] = Factor * Coefficients(0, m3);
	  this->M12InteractionFactorsdmdm[this->NbrM12IntraIndices] = Factor * Coefficients(0, m3);
	  this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	  this->M2IntraValue[this->NbrM12IntraIndices] = 0;
	  this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	  this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	  this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3 - 1;
	  ++this->NbrM12IntraIndices;
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
	      this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M12InteractionFactorsumum[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M12InteractionFactorsdpdp[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M12InteractionFactorsdmdm[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	      this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	      this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	      this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	      this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3 - 1;
	      ++this->NbrM12IntraIndices;
	    }
	}
      if (this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic)
	{
	  this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsumum[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsdpdp[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsdmdm[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	  this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	  this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	  this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	  this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	  this->M3IntraValues[this->NbrM12IntraIndices][0] = m3 - 1;
	  ++this->NbrM12IntraIndices;	  
	}
      ++m3;
      for (; m3 <= this->LzMax; ++m3)
	{
 	  if ((this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic) ||
	      (m3 != (m4 + 2)))
            {
	      this->M12InteractionFactorsupup[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M12InteractionFactorsumum[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M12InteractionFactorsdpdp[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M12InteractionFactorsdmdm[this->NbrM12IntraIndices] = Factor * Coefficients(m4, m3);
	      this->M1IntraValue[this->NbrM12IntraIndices] = m3;
	      this->M2IntraValue[this->NbrM12IntraIndices] = m4;
	      this->M3IntraValues[this->NbrM12IntraIndices] = new int [1];
	      this->NbrM3IntraValues[this->NbrM12IntraIndices] = 1;
	      this->M3IntraValues[this->NbrM12IntraIndices][0] =  m3 - 1;
	      ++this->NbrM12IntraIndices;
	    }
	}
    }

  this->NbrM12InterIndices = (this->LzMax + 3) * this->LzMax;
  this->M1InterValue = new int [this->NbrM12InterIndices];
  this->M2InterValue = new int [this->NbrM12InterIndices];
  this->M3InterValues = new int* [this->NbrM12InterIndices];
  this->NbrM3InterValues  = new int [this->NbrM12InterIndices];

  this->M12InteractionFactorsupum = new double [4 * this->NbrM12InterIndices];
  this->M12InteractionFactorsupdp = new double [4 * this->NbrM12InterIndices];
  this->M12InteractionFactorsupdm = new double [4 * this->NbrM12InterIndices];
  this->M12InteractionFactorsumdp = new double [4 * this->NbrM12InterIndices];
  this->M12InteractionFactorsumdm = new double [4 * this->NbrM12InterIndices];
  this->M12InteractionFactorsdpdm = new double [4 * this->NbrM12InterIndices];
  
  this->NbrM12InterIndices = 0;
  int TmpNbrM12InterIndices = 0;
  
  //   this->L2Factor  = 0.0;
  //   Factor = 2.0 * this->L2Factor;
  //   if (this->Particles->GetParticleStatistic() == ParticleOnSphereWithSU4Spin::FermionicStatistic)
  //     Factor *= -1.0;
  
  for (int m3 = 1; m3 <= this->LzMax; ++m3)
    {
      this->M12InteractionFactorsupum[TmpNbrM12InterIndices] = Factor * Coefficients(0, m3);
      this->M12InteractionFactorsupdp[TmpNbrM12InterIndices] = Factor * Coefficients(0, m3);
      this->M12InteractionFactorsupdm[TmpNbrM12InterIndices] = Factor * Coefficients(0, m3);
      this->M12InteractionFactorsumdp[TmpNbrM12InterIndices] = Factor * Coefficients(0, m3);
      this->M12InteractionFactorsumdm[TmpNbrM12InterIndices] = Factor * Coefficients(0, m3);
      this->M12InteractionFactorsdpdm[TmpNbrM12InterIndices] = Factor * Coefficients(0, m3);
      ++TmpNbrM12InterIndices;
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = 0;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 - 1;
      ++this->NbrM12InterIndices;
    }
  for (int m4 = 1; m4 < this->LzMax; ++m4)
    {
      int m3= 0;
      this->M12InteractionFactorsupum[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
      this->M12InteractionFactorsupdp[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
      this->M12InteractionFactorsupdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
      this->M12InteractionFactorsumdp[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
      this->M12InteractionFactorsumdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
      this->M12InteractionFactorsdpdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
      ++TmpNbrM12InterIndices;
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = m4;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 + 1;
      ++this->NbrM12InterIndices;
      ++m3;
      for (; m3 < this->LzMax; ++m3)
	{
	  this->M12InteractionFactorsupum[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
	  this->M12InteractionFactorsupdp[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
	  this->M12InteractionFactorsupdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
	  this->M12InteractionFactorsumdp[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
	  this->M12InteractionFactorsumdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
	  this->M12InteractionFactorsdpdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, m4);
	  ++TmpNbrM12InterIndices;
	  this->M12InteractionFactorsupum[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsupdp[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsupdm[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsumdp[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsumdm[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
	  this->M12InteractionFactorsdpdm[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
	  ++TmpNbrM12InterIndices;
	  this->M1InterValue[this->NbrM12InterIndices] = m3;
	  this->M2InterValue[this->NbrM12InterIndices] = m4;
	  this->M3InterValues[this->NbrM12InterIndices] = new int [2];
	  this->NbrM3InterValues[this->NbrM12InterIndices] = 2;
	  this->M3InterValues[this->NbrM12InterIndices][0] =  m3 + 1;
	  this->M3InterValues[this->NbrM12InterIndices][1] =  m3 - 1;
	  ++this->NbrM12InterIndices;
	}
      this->M12InteractionFactorsupum[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
      this->M12InteractionFactorsupdp[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
      this->M12InteractionFactorsupdm[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
      this->M12InteractionFactorsumdp[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
      this->M12InteractionFactorsumdm[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
      this->M12InteractionFactorsdpdm[TmpNbrM12InterIndices] = Factor * Coefficients(m4, m3);
      ++TmpNbrM12InterIndices;
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = m4;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 - 1;
      ++this->NbrM12InterIndices;
    }
  for (int m3 = 0; m3 < this->LzMax; ++m3)
    {
      this->M12InteractionFactorsupum[TmpNbrM12InterIndices] = Factor * Coefficients(m3, this->LzMax);
      this->M12InteractionFactorsupdp[TmpNbrM12InterIndices] = Factor * Coefficients(m3, this->LzMax);
      this->M12InteractionFactorsupdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, this->LzMax);
      this->M12InteractionFactorsumdp[TmpNbrM12InterIndices] = Factor * Coefficients(m3, this->LzMax);
      this->M12InteractionFactorsumdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, this->LzMax);
      this->M12InteractionFactorsdpdm[TmpNbrM12InterIndices] = Factor * Coefficients(m3, this->LzMax);
      ++TmpNbrM12InterIndices;
      this->M1InterValue[this->NbrM12InterIndices] = m3;
      this->M2InterValue[this->NbrM12InterIndices] = this->LzMax;
      this->M3InterValues[this->NbrM12InterIndices] = new int [1];
      this->NbrM3InterValues[this->NbrM12InterIndices] = 1;
      this->M3InterValues[this->NbrM12InterIndices][0] =  m3 + 1;
      ++this->NbrM12InterIndices;
    }

  //   this->L2Factor  = 0.0;

  if (this->LzMax>0)
    {     
      this->OneBodyInteractionFactorsupup = new double[this->LzMax + 1];
      this->OneBodyInteractionFactorsumum = new double[this->LzMax + 1];
      this->OneBodyInteractionFactorsdpdp = new double[this->LzMax + 1];
      this->OneBodyInteractionFactorsdmdm = new double[this->LzMax + 1];
      this->OneBodyInteractionFactorsupup[0] = this->L2Factor * Coefficients(0, 1);
      this->OneBodyInteractionFactorsumum[0] = this->L2Factor * Coefficients(0, 1);
      this->OneBodyInteractionFactorsdpdp[0] = this->L2Factor * Coefficients(0, 1);
      this->OneBodyInteractionFactorsdmdm[0] = this->L2Factor * Coefficients(0, 1);
    
      for (int i = 1; i < this->LzMax; ++i)
	{
	  this->OneBodyInteractionFactorsupup[i] = this->L2Factor * (Coefficients(i, i + 1) + Coefficients(i - 1, i));
	  this->OneBodyInteractionFactorsumum[i] = this->L2Factor * (Coefficients(i, i + 1) + Coefficients(i - 1, i));
	  this->OneBodyInteractionFactorsdpdp[i] = this->L2Factor * (Coefficients(i, i + 1) + Coefficients(i - 1, i));
	  this->OneBodyInteractionFactorsdmdm[i] = this->L2Factor * (Coefficients(i, i + 1) + Coefficients(i - 1, i));
	  
	}

      this->OneBodyInteractionFactorsupup[this->LzMax] = this->L2Factor * Coefficients(this->LzMax - 1, this->LzMax);
      this->OneBodyInteractionFactorsumum[this->LzMax] = this->L2Factor * Coefficients(this->LzMax - 1, this->LzMax);
      this->OneBodyInteractionFactorsdpdp[this->LzMax] = this->L2Factor * Coefficients(this->LzMax - 1, this->LzMax);
      this->OneBodyInteractionFactorsdmdm[this->LzMax] = this->L2Factor * Coefficients(this->LzMax - 1, this->LzMax);
    }
  else
    {
      this->OneBodyInteractionFactorsupup = 0;
      this->OneBodyInteractionFactorsumum = 0;
      this->OneBodyInteractionFactorsdpdp = 0;
      this->OneBodyInteractionFactorsdmdm = 0;
    }

  cout << "nbr interaction = " << ((4 * (this->NbrM12IntraIndices + this->LzMax)) + 4 + 6 * this->NbrM12InterIndices) << endl;
  cout << "====================================" << endl;
  
}
