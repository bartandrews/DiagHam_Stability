////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2006 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//     class of hamiltonian associated to particles on a sphere with spin     //
//     and interacting with a delta interaction in the s and p-wave channels  //
//                                                                            //
//                        last modification : 27/01/2006                      //
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


#include "Hamiltonian/ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian (ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, double v0, double v1,
   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;
  this->V0=v0;
  this->V1=v1;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
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
	  long TmpMemoryNeeded = this->FastMultiplicationMemory(memory);
	  if (TmpMemoryNeeded < 1024)
	    cout  << "fast = " <<  TmpMemoryNeeded << "b ";
	  else
	    if (TmpMemoryNeeded < (1 << 20))
	      cout  << "fast = " << (TmpMemoryNeeded >> 10) << "kb ";
	    else
	      if (TmpMemoryNeeded < (1 << 30))
		cout  << "fast = " << (TmpMemoryNeeded >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemoryNeeded >> 30) << "Gb ";
	  
	  if (this->DiskStorageFlag == false)
	    {
	      cout << "Enabling fast calculation!" << endl;
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

ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::~ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian() 
{
  delete[] this->UUInteractionFactors;
  delete[] this->UUM1Value;
  delete[] this->UUM2Value;
  delete[] this->UUM3Value;
  delete[] this->UDInteractionFactors;
  delete[] this->UDM1Value;
  delete[] this->UDM2Value;
  delete[] this->UDM3Value;
  if (this->FastMultiplicationFlag == true)
    {
      if (this->DiskStorageFlag == false)
	{
	  long MinIndex;
	  long MaxIndex;
	  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
	  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
	  int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	    ++ReducedDim;
	  for (int i = 0; i < ReducedDim; ++i)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	  delete[] this->InteractionPerComponentIndex;
	  delete[] this->InteractionPerComponentCoefficient;
	}
      else
	{
	  remove (this->DiskStorageFileName);
	}
      delete[] this->NbrInteractionPerComponent;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->UUInteractionFactors;
  delete[] this->UUM1Value;
  delete[] this->UUM2Value;
  delete[] this->UUM3Value;
  delete[] this->UDInteractionFactors;
  delete[] this->UDM1Value;
  delete[] this->UDM2Value;
  delete[] this->UDM3Value;
  
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}


// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

  // reset interaction coefficients
  // v0 interaction in s-wave Delta-Interaction channel
  // v1 interaction in p-wave Laplacian-Delta-Interaction channel
void ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::SetInteraction(double v0, double v1)
{
  this->V0=v0;
  this->V1=v1;
  this->EvaluateInteractionFactors();

}

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
Complex ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  cout << "Fake function: MatrixElement" << endl;
  return Complex();
}
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
Complex ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  cout << "Fake function: MatrixElement" << endl;
  return Complex();
}



// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value
void ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift=shift;
}


// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::EvaluateInteractionFactors()
{
  int Lim;
  int Min;
  int Pos = 0;
//  int NbrNonZero = 0;
//  cout << "this->LzMax=" << this->LzMax << endl;
  double TmpV = (((double) this->LzMax) + 1.0);
  TmpV = 1.0; //(TmpV * TmpV) / ((0.5 * ((double) this->LzMax)) * ((2.0 * (double) this->LzMax) - 1.0));
  
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);

  int J0, J1; 
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double MaxCoefficient = 0.0;
  // make sure that we have Fermionic Statistics -> we're not dealing with Spin 1 bosons, yet!
  if (this->Particles->GetParticleStatistic() != ParticleOnSphere::FermionicStatistic)
    {
      cout << "No Bosons with Spin defined, so far!" << endl;
      exit(1);
    }
  else // we do have Fermions, then!
    {
      // treat V0-Delta-Interaction, and UD-part of V1-Interaction, first (different spin channels):
      J0 = 2 * this->LzMax;
      J1 = 2 * (this->LzMax - 1);
      for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	for (int m2 =  -this->LzMax; m2 <= this->LzMax; m2 += 2)
	  {
	    Lim = m1 + m2 + this->LzMax;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < -this->LzMax)
	      Min = -this->LzMax;
	    for (int m3 = Min; m3 <= Lim; m3 += 2)
	      {
		m4 = m1 + m2 - m3;
		TmpCoefficient[Pos] = 0.0;
		Clebsch.InitializeCoefficientIterator(m1, m2);
		TmpCoefficient[Pos] += V0 * TmpV * Clebsch.GetCoefficient(m1, m2, J0) * Clebsch.GetCoefficient(m3, m4, J0);
		if ( (m1 != m2) && ( m3 != m4) && (abs(m1 + m2) <= J1) && (abs(m3 + m4) <= J1))
		  TmpCoefficient[Pos] += V1 * TmpV * Clebsch.GetCoefficient(m1, m2, J1) * Clebsch.GetCoefficient(m3, m4, J1);
		if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  MaxCoefficient = fabs(TmpCoefficient[Pos]);
		++Pos;
	      }
	  }
      this->UDNbrInteractionFactors = 0;
      this->UDM1Value = new int [Pos];
      this->UDM2Value = new int [Pos];
      this->UDM3Value = new int [Pos];
      this->UDInteractionFactors = new double [Pos];
      //      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      double Factor = - 1.0; // (0.5 * ((double) this->LzMax));  // A revoir! xxx
      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
	  {
	    Lim = m1 + m2;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < 0)
	      Min = 0;
	    for (int m3 = Min; m3 <= Lim; ++m3)
	      {
		if (( fabs(TmpCoefficient[Pos]) > MaxCoefficient))//non-negligleable matrix element
		  {
		    this->UDInteractionFactors[this->UDNbrInteractionFactors] = Factor * TmpCoefficient[Pos];
		    this->UDM1Value[this->UDNbrInteractionFactors] = m1;
		    this->UDM2Value[this->UDNbrInteractionFactors] = m2;
		    this->UDM3Value[this->UDNbrInteractionFactors] = m3;
		    /*cout <<"UD: "<< this->UDM1Value[this->UDNbrInteractionFactors] << " " << this->UDM2Value[this->UDNbrInteractionFactors] 
		      << " " << this->UDM3Value[this->UDNbrInteractionFactors] 
		      << " " << this->UDInteractionFactors[this->UDNbrInteractionFactors] << endl;*/
		    ++this->UDNbrInteractionFactors;
		  }
		++Pos;
	      }
	  }


      // Now, look at V1-Delta-Interaction in same spin channels:
      Pos=0;
      MaxCoefficient=0.0;
      J1 = 2 * (this->LzMax - 1);
      for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	for (int m2 =  -this->LzMax; m2  < m1; m2 += 2)
	  {
	    Lim = m1 + m2 + this->LzMax;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < -this->LzMax)
	      Min = -this->LzMax;
	    for (int m3 = Min; m3 <= Lim; m3 += 2)
	      {
		Clebsch.InitializeCoefficientIterator(m1, m2);
		m4 = m1 + m2 - m3;
		if ( m4 < m3 )
		  {
		    TmpCoefficient[Pos] = 0.0;
		    TmpCoefficient[Pos] += TmpV * Clebsch.GetCoefficient(m1, m2, J1) * Clebsch.GetCoefficient(m3, m4, J1);
		    if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
	      }
	  }
	    
      // now store factors in required format:
      this->UUNbrInteractionFactors = 0;
      this->UUM1Value = new int [Pos];
      this->UUM2Value = new int [Pos];
      this->UUM3Value = new int [Pos];
      this->UUInteractionFactors = new double [Pos];
//      cout << "nbr UU interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      Factor = - this->V1  * 2.0; // factor 4.0 is for symmetry issues!
      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	for (int m2 = 0; m2 < m1; ++m2) // m1; ++m2)
	  {
	    Lim = m1 + m2;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < 0)
	      Min = 0;
	    for (int m3 = Min; m3 <= Lim; ++m3)
	      {
		if ( 2 * m3 > m1 + m2 )
		  {
		    if ( fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		      {
			this->UUInteractionFactors[this->UUNbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			this->UUM1Value[this->UUNbrInteractionFactors] = m1;
			this->UUM2Value[this->UUNbrInteractionFactors] = m2;
			this->UUM3Value[this->UUNbrInteractionFactors] = m3;
			/* cout << "UU: " << this->UUM1Value[this->UUNbrInteractionFactors] << " " << this->UUM2Value[this->UUNbrInteractionFactors] 
			     << " " << this->UUM3Value[this->UUNbrInteractionFactors] 
			     << " " << this->UUInteractionFactors[this->UUNbrInteractionFactors] << endl; */
			++this->UUNbrInteractionFactors;
		      }
		    ++Pos;
		  }
	      }
	  }
    }
  delete[] TmpCoefficient;
}



  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
RealVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
RealVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
			       int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}


  // multiply a vector by the current hamiltonian for all indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
RealVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
RealVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
											int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient=0.0;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      // calculate same spin channel:
      int ReducedNbrInteractionFactors = this->UUNbrInteractionFactors - 1;
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      for (int j = 0; j < ReducedNbrInteractionFactors; ++j)
	{
	  m1 = this->UUM1Value[j];
	  m2 = this->UUM2Value[j];
	  m3 = this->UUM3Value[j];
	  TmpInteraction = this->UUInteractionFactors[j];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      // Interactions of Fermions with Spin up:
	      Index = this->Particles->AduAduAuAu(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      
	      // Interactions of Fermions with Spin down:
	      Index = this->Particles->AddAddAdAd(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	    }
	}
      m1 = this->UUM1Value[ReducedNbrInteractionFactors];
      m2 = this->UUM2Value[ReducedNbrInteractionFactors];
      m3 = this->UUM3Value[ReducedNbrInteractionFactors];
      TmpInteraction = this->UUInteractionFactors[ReducedNbrInteractionFactors];
      m4 = m1 + m2 - m3;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	      // Interactions of Fermions with Spin up:
	      Index = this->Particles->AduAduAuAu(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      // Interactions of Fermions with Spin down:
	      Index = this->Particles->AddAddAdAd(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      // apply artificial shift, if so wanted
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      
      // now treat mixed spin cases:
      for (int j = 0; j < this->UDNbrInteractionFactors; ++j) 
	{
	  m1 = this->UDM1Value[j];
	  m2 = this->UDM2Value[j];
	  m3 = this->UDM3Value[j];
	  TmpInteraction = this->UDInteractionFactors[j];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      // Interactions of Fermions Interactions mixing both spin species:
	      Index = this->Particles->AddAduAdAu(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{ 
	  int* TmpIndexArray; 
	  double* TmpCoefficientArray;  
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      return this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      return this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;

}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int l = 0; l < this->BufferSize; ++l)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int i = 0; i < Lim; ++i)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  return vDestination;
}



// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSource = array of vectors to be multiplied
// vDestination = array of vectors at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelAddMultiplyPartialFastMultiply 
(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      Coefficient=vSource[l];
      vDestination[l] += this->HamiltonianShift * Coefficient;
      for (j = 0; j < TmpNbrInteraction; ++j)
	vDestination[ TmpIndexArray[j] ] += TmpCoefficientArray[j]  * Coefficient;
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	// calculate same spin channel:
	int ReducedNbrInteractionFactors = this->UUNbrInteractionFactors - 1;

	for (int j = 0; j < ReducedNbrInteractionFactors; ++j)
	  {
	    m1 = this->UUM1Value[j];
	    m2 = this->UUM2Value[j];
	    m3 = this->UUM3Value[j];
	    TmpInteraction = this->UUInteractionFactors[j];
	    m4 = m1 + m2 - m3;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		// Interactions of Fermions with Spin up:
		Index = this->Particles->AduAduAuAu(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      
		// Interactions of Fermions with Spin down:
		Index = this->Particles->AddAddAdAd(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      }
	  }
	m1 = this->UUM1Value[ReducedNbrInteractionFactors];
	m2 = this->UUM2Value[ReducedNbrInteractionFactors];
	m3 = this->UUM3Value[ReducedNbrInteractionFactors];
	TmpInteraction = this->UUInteractionFactors[ReducedNbrInteractionFactors];
	m4 = m1 + m2 - m3;
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    // Interactions of Fermions with Spin up:
	    Index = this->Particles->AduAduAuAu(i, m1, m2, m3, m4, Coefficient);
	    if (Index < Dim)
	      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	    // Interactions of Fermions with Spin down:
	    Index = this->Particles->AddAddAdAd(i, m1, m2, m3, m4, Coefficient);
	    if (Index < Dim)
	      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	    // apply artificial shift, if so wanted
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	  }
	
	// now treat mixed spin cases:
	for (int j = 0; j < this->UDNbrInteractionFactors; ++j) 
	  {
	    m1 = this->UDM1Value[j];
	    m2 = this->UDM2Value[j];
	    m3 = this->UDM3Value[j];
	    TmpInteraction = this->UDInteractionFactors[j];
	    m4 = m1 + m2 - m3;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		// Interactions of Fermions Interactions mixing both spin species:
		Index = this->Particles->AddAduAdAu(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      }
	  }

      }
  delete TmpParticles;
  return vDestination;
}



  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
ComplexVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
ComplexVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination.Re(i) = 0.0;
      vDestination.Im(i) = 0.0;
    }
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
ComplexVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
ComplexVector& ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				   int firstComponent, int nbrComponent)
{
  cout << "Attention, ComplexVector& LowLevelAddMultiply not fully implemented!" << endl;
  return vDestination;  // not fully implemented!
}





// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Particles->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i] << "    ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Particles->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Particles->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][H.Particles->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}


  // test the amount of memory needed for fast multiplication algorithm
  //
  // return value = amount of memory needed
long ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
   }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      this->FastMultiplicationStep = 1;
    }
  
  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;

  return Memory;
}



// enable fast multiplication algorithm
//
void ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int Index;
  double Coefficient;
  int m1;
  int m2;
  int m3;
  int m4;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];

  int TotalPos = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]]; 
      this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
      Pos = 0;
      for (int j = 0; j < this->UUNbrInteractionFactors; ++j) 
	{
	  m1 = this->UUM1Value[j];
	  m2 = this->UUM2Value[j];
	  m3 = this->UUM3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = this->Particles->AddAddAdAd(i + this->PrecalculationShift, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->UUInteractionFactors[j];
	      ++Pos;
	    }
	  Index = this->Particles->AduAduAuAu(i + this->PrecalculationShift, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->UUInteractionFactors[j];
	      ++Pos;
	    }
	}
	 
      for (int j = 0; j < this->UDNbrInteractionFactors; ++j) 
	{
	  m1 = this->UDM1Value[j];
	  m2 = this->UDM2Value[j];
	  m3 = this->UDM3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = this->Particles->AddAduAdAu(i + this->PrecalculationShift, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->UDInteractionFactors[j];
	      ++Pos;
	    }
	}
      ++TotalPos;
    }
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;

}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  if (this->FastMultiplicationStep == 1)
    {
      this->DiskStorageFlag = false;
      this->DiskStorageFileName = 0;
      cout << "Enabling fast calculation!" << endl;
      this->EnableFastMultiplication();
      return;
    }
  cout << "Enabling disk storage for fast calculation!" << endl;
  this->DiskStorageFlag = true;
  this->DiskStorageFileName = new char [strlen(fileName) + 8];
  sprintf (this->DiskStorageFileName, "%s.ham", fileName);

  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->DiskStorageStart = (int) MinIndex;
  int DiskStorageEnd = 1 + (int) MaxIndex;

  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  this->MaxNbrInteractionPerComponent = 0;

  int TotalPos = 0;
  ofstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::out);
 
  File.write((char*) &(EffectiveHilbertSpaceDimension), sizeof(int));
  File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
  File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * EffectiveHilbertSpaceDimension);

  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      FileJump += (long) this->NbrInteractionPerComponent[i];
      if (this->MaxNbrInteractionPerComponent < this->NbrInteractionPerComponent[i])
	this->MaxNbrInteractionPerComponent = this->NbrInteractionPerComponent[i];
    }
  FileJump *= sizeof(int);

  TmpIndexArray = new int [this->MaxNbrInteractionPerComponent];
  TmpCoefficientArray = new double [this->MaxNbrInteractionPerComponent];      

  for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
    {
      if (this->NbrInteractionPerComponent[TotalPos] > 0)
	{
	  Pos = 0;
	  for (int j = 0; j < UUNbrInteractionFactors; ++j)
	    {
	      m1 = this->UUM1Value[j];
	      m2 = this->UUM2Value[j];
	      m3 = this->UUM3Value[j];
	      m4 = m1 + m2 - m3;
	      // Interactions of Fermions with Spin up:
	      Index = this->Particles->AduAduAuAu(i, m1, m2, m3, m4, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Pos] = Index;
		  TmpCoefficientArray[Pos] = Coefficient * this->UUInteractionFactors[j];
		  ++Pos;
		}
	      // Interactions of Fermions with Spin down:
	      Index = this->Particles->AddAddAdAd(i, m1, m2, m3, m4, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Pos] = Index;
		  TmpCoefficientArray[Pos] = Coefficient * this->UUInteractionFactors[j];
		  ++Pos;
		}
	    }
      
	  // now treat mixed spin cases:
	  for (int j = 0; j < this->UDNbrInteractionFactors; ++j) 
	    {
	      m1 = this->UDM1Value[j];
	      m2 = this->UDM2Value[j];
	      m3 = this->UDM3Value[j];
	      m4 = m1 + m2 - m3;
	      // Interactions of Fermions Interactions mixing both spin species:
	      Index = this->Particles->AddAduAdAu(i, m1, m2, m3, m4, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  TmpIndexArray[Pos] = Index;
		  TmpCoefficientArray[Pos] = Coefficient * this->UDInteractionFactors[j];
		  ++Pos;
		}
	    }
	  File.write((char*) TmpIndexArray, sizeof(int) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump -= sizeof(int) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(FileJump, ios::cur);
	  File.write((char*) TmpCoefficientArray, sizeof(double) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump += sizeof(double) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(-FileJump, ios::cur);	  
	}
      ++TotalPos;
    }
  delete[] TmpIndexArray;
  delete[] TmpCoefficientArray;
  File.close();

  this->FastMultiplicationFlag = true;
  this->BufferSize = this->Memory / ((this->MaxNbrInteractionPerComponent * (sizeof(int) + sizeof(double))) + sizeof(int*) + sizeof(double*));

  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}




// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < this->UUNbrInteractionFactors; ++j) 
	{
	  m1 = this->UUM1Value[j];
	  m2 = this->UUM2Value[j];
	  m3 = this->UUM3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = TmpParticles->AddAddAdAd(i, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
	    }
	  Index = TmpParticles->AduAduAuAu(i, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
	    }
	}
      for (int j = 0; j < this->UDNbrInteractionFactors; ++j) 
	{
	  m1 = this->UDM1Value[j];
	  m2 = this->UDM2Value[j];
	  m3 = this->UDM3Value[j];
	  m4 = m1 + m2 - m3;
	  // Interactions of Fermions with mixed Spins:
	  Index = this->Particles->AddAduAdAu(i, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
	    }
	}
    }
  delete TmpParticles;

  return Memory;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Particles->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
	}
      File.close();
      return true;
    }
  else
    {
      return false;
    }
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool ParticleOnSphereWithSpinDeltaLaplacianDeltaHamiltonian::LoadPrecalculation (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Particles->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrInteractionPerComponent = new int [Tmp];
  File.read((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficient = new double* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
    }
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}
