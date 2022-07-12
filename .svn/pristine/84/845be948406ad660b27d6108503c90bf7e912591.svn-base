////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//     SU(2) spin and a generic interaction defined by its pseudopotential    //
//     the hamiltonian is not assumed to be SU(2) symmetric                   //
//                                                                            //
//                        last modification : 05/10/2012                      //
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


#include "Hamiltonian/ParticleOnSphereWithSpinFullGenericHamiltonian.h"
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


ParticleOnSphereWithSpinFullGenericHamiltonian::ParticleOnSphereWithSpinFullGenericHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as UpUpUpUp, UpDownUpUp, DownDownUpUp, UpUpUpDown, ... (9 terms) )
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, from down to up
// onebodyPotentialDownUp =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, from up to down
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinFullGenericHamiltonian::ParticleOnSphereWithSpinFullGenericHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
										       double** pseudoPotential, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown, double* onebodyPotentialUpDown, double* onebodyPotentialDownUp,
										       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										       char* precalculationFileName)
{
  // For AbstractQHEOnSphereWithSpinHamiltonian destructor
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;

  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->PseudoPotentials = new double* [10];
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  for (int j = 0; j < 10; ++j)
    {
      this->PseudoPotentials[j] = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
    }
//  test : print pseudopotentials
//  for (int j = 0; j < 10; ++j)
//    {
//      for (int i = 0; i < this->NbrLzValue; ++i)
//          cout << this->PseudoPotentials[j][i] << " ";
//      cout << endl;
//    }
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  
  this->OneBodyInteractionFactorsUpUp = 0;
  this->NbrOneBodyInteractionFactorsUpUp = 0;
  this->OneBodyMValuesUpUp = 0;
  if (onebodyPotentialUpUp != 0)
    {
      this->OneBodyInteractionFactorsUpUp = new double [this->NbrLzValue];
      this->NbrOneBodyInteractionFactorsUpUp = this->LzMax+1;
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsUpUp[i] = onebodyPotentialUpUp[i];
      this->OneBodyMValuesUpUp = new int [this->NbrLzValue];
      for (int i=0; i <= this->LzMax; i++)
        OneBodyMValuesUpUp[i]=i;
    }

  this->OneBodyInteractionFactorsUpDown = 0;
  this->NbrOneBodyInteractionFactorsUpDown = 0;
  this->OneBodyMValuesUpDown = 0;
  if (onebodyPotentialUpDown != 0)
    {
      this->NbrOneBodyInteractionFactorsUpDown = this->LzMax+1;
      this->OneBodyInteractionFactorsUpDown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsUpDown[i] = onebodyPotentialUpDown[i];
      this->OneBodyMValuesUpDown = new int [this->NbrLzValue];
      for (int i=0; i <= this->LzMax; i++)
        OneBodyMValuesUpDown[i]=i;
    }

  this->OneBodyInteractionFactorsDownUp = 0;
  this->NbrOneBodyInteractionFactorsDownUp = 0;
  this->OneBodyMValuesDownUp = 0;
  if (onebodyPotentialDownUp != 0)
    {
      this->NbrOneBodyInteractionFactorsDownUp = this->LzMax+1;
      this->OneBodyInteractionFactorsDownUp = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsDownUp[i] = onebodyPotentialDownUp[i];
      this->OneBodyMValuesDownUp = new int [this->NbrLzValue];
      for (int i=0; i <= this->LzMax; i++)
        OneBodyMValuesDownUp[i]=i;
    }

  this->OneBodyInteractionFactorsDownDown = 0;
  this->NbrOneBodyInteractionFactorsDownDown = 0;
  this->OneBodyMValuesDownDown = 0;
  if (onebodyPotentialDownDown != 0)
    {
      this->NbrOneBodyInteractionFactorsDownDown = this->LzMax+1;
      this->OneBodyInteractionFactorsDownDown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsDownDown[i] = onebodyPotentialDownDown[i];
      this->OneBodyMValuesDownDown = new int [this->NbrLzValue];
      for (int i=0; i <= this->LzMax; i++)
        OneBodyMValuesDownDown[i]=i;
    }
    
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
//  for( int i=0; i<5; ++i) 
//	cout << "NbrInteractionPerComponent["<<i<<"] = " << this->NbrInteractionPerComponent[i] << endl;
}

// destructor
//
ParticleOnSphereWithSpinFullGenericHamiltonian::~ParticleOnSphereWithSpinFullGenericHamiltonian() 
{
  for (int j = 0; j < 10; ++j)
    delete[] this->PseudoPotentials[j];
  delete[] this->PseudoPotentials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinFullGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinFullGenericHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinFullGenericHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereWithSpinFullGenericHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinFullGenericHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex ParticleOnSphereWithSpinFullGenericHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereWithSpinFullGenericHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereWithSpinFullGenericHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   
void ParticleOnSphereWithSpinFullGenericHamiltonian::EvaluateInteractionFactors()
{
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);

  const int NbrInteractionSectors = 9;
  int J = 2 * this->LzMax - 2;
  double ClebschCoef;
  long TotalNbrInteractionFactors = 0;
  
  // Factors which multiply the InteractionsFactors after Clebsch-Gordan have been computed
  double Factors[NbrInteractionSectors]={-4.,-2.,-4.,-2.,-2.,-2.,-4.,-2.,-4.};
  
  // Sign = 1 if LzMax is even and 0 otherwise
  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;
  double TmpCoefficient2 = 0.0;

 
  // Initialize UpDownSectorIndicesPerSum
  this->NbrUpDownSectorSums = 2 * this->LzMax + 1;
  this->NbrUpDownSectorIndicesPerSum = new int[this->NbrUpDownSectorSums];
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    this->NbrUpDownSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrUpDownSectorIndicesPerSum[m1 + m2];  //Initialize NbrUpDownSectorIndicesPerSum to {1,2,3,...,LzMax+1,...,3,2,1}    
  this->UpDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    {
      this->UpDownSectorIndicesPerSum[i] = new int[2 * this->NbrUpDownSectorIndicesPerSum[i]];      
      this->NbrUpDownSectorIndicesPerSum[i] = 0;
    }
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      {
	this->UpDownSectorIndicesPerSum[(m1 + m2)][this->NbrUpDownSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	this->UpDownSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrUpDownSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	++this->NbrUpDownSectorIndicesPerSum[(m1 + m2)];
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      // Initialize UpUpSectorIndicesPerSum
      this->NbrUpUpSectorSums = 2 * this->LzMax + 1;
      this->NbrUpUpSectorIndicesPerSum = new int[this->NbrUpUpSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	this->NbrUpUpSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrUpUpSectorIndicesPerSum[(m1 + m2)]; //Initialize NbrUpUpSectorIndicesPerSum to {0,1,1,2,2,3,3,...,(IntegerPart((LzMax+1)/2)),....,3,3,2,2,1,1,0}
      this->UpUpSectorIndicesPerSum = new int* [this->NbrUpUpSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	{
	  this->UpUpSectorIndicesPerSum[i] = new int[2 * this->NbrUpUpSectorIndicesPerSum[i]];      
	  this->NbrUpUpSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    this->UpUpSectorIndicesPerSum[(m1 + m2)][this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	    this->UpUpSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	    ++this->NbrUpUpSectorIndicesPerSum[(m1 + m2)];
	  }

      // Initialize DownDownSectorIndicesPerSum
      this->NbrDownDownSectorSums = 2 * this->LzMax + 1;
      this->NbrDownDownSectorIndicesPerSum = new int[this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	this->NbrDownDownSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrDownDownSectorIndicesPerSum[(m1 + m2)];
      this->DownDownSectorIndicesPerSum = new int* [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	{
	  this->DownDownSectorIndicesPerSum[i] = new int[2 * this->NbrDownDownSectorIndicesPerSum[i]];      
	  this->NbrDownDownSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    this->DownDownSectorIndicesPerSum[(m1 + m2)][this->NbrDownDownSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	    this->DownDownSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrDownDownSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	    ++this->NbrDownDownSectorIndicesPerSum[(m1 + m2)];
	  }

      // Allocate memory for InteractionFactors
      this->InteractionFactorsUpUpUpUp = new double* [this->NbrUpUpSectorSums];
      this->InteractionFactorsUpDownUpUp = new double* [this->NbrUpUpSectorSums];
      this->InteractionFactorsDownDownUpUp = new double* [this->NbrUpUpSectorSums];
      this->InteractionFactorsUpUpUpDown = new double* [this->NbrUpDownSectorSums];
      this->InteractionFactorsUpDownUpDown = new double* [this->NbrUpDownSectorSums];
      this->InteractionFactorsDownDownUpDown = new double* [this->NbrUpDownSectorSums];
      this->InteractionFactorsUpUpDownDown = new double* [this->NbrDownDownSectorSums];
      this->InteractionFactorsUpDownDownDown = new double* [this->NbrDownDownSectorSums];
      this->InteractionFactorsDownDownDownDown = new double* [this->NbrDownDownSectorSums];

      //Compute interaction factors
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	{
          int Index1=0,Index2=0,Index3=0;

          this->InteractionFactorsUpUpUpUp[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]];
          this->InteractionFactorsUpDownUpUp[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]];
          this->InteractionFactorsDownDownUpUp[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]];
          this->InteractionFactorsUpUpUpDown[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]];
          this->InteractionFactorsUpDownUpDown[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]];
          this->InteractionFactorsDownDownUpDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]];
          this->InteractionFactorsUpUpDownDown[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]];
          this->InteractionFactorsUpDownDownDown[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]];
          this->InteractionFactorsDownDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]];

          
          // UpUp input particles
	  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;

              // UpUpUpUp sector
	      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsUpUpUpUp[i][Index1] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) == Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsUpUpUpUp[i][Index1] += this->PseudoPotentials[0][J >> 1] * TmpCoefficient;
			}
		    }
		  this->InteractionFactorsUpUpUpUp[i][Index1] *= Factors[0];
		  ++TotalNbrInteractionFactors;
		  ++Index1; 
                }
              // UpDownUpUp sector
	      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;

		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsUpDownUpUp[i][Index2] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) == Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  TmpCoefficient2 = ClebschCoef * Clebsch.GetCoefficient(m4, m3, J);
			  this->InteractionFactorsUpDownUpUp[i][Index2] += this->PseudoPotentials[1][J >> 1] * (TmpCoefficient-TmpCoefficient2);
			}
                    }

		  this->InteractionFactorsUpDownUpUp[i][Index2] *= Factors[1];
                  ++TotalNbrInteractionFactors;
		  ++Index2;
	        }

              // DownDownUpUp sector
	      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsDownDownUpUp[i][Index3] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) == Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsDownDownUpUp[i][Index3] += this->PseudoPotentials[2][J >> 1] * TmpCoefficient;
			}
		    }
		  this->InteractionFactorsDownDownUpUp[i][Index3] *= Factors[2];
		  ++TotalNbrInteractionFactors;
		  ++Index3;
	       }	
	    }


       Index1=0;Index2=0;Index3=0;
       // UpDown input particles
       for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
	 {
	   int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	   int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
    
           //UpUpUpDown sector
	   for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
             {
	       int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
	       int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
	       Clebsch.InitializeCoefficientIterator(m1, m2);
	       this->InteractionFactorsUpUpUpDown[i][Index1] = 0.0;
	       while (Clebsch.Iterate(J, ClebschCoef))
		 {
		   if (((J >> 1) & 1) == Sign)
                     {
		       TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		       TmpCoefficient2 = ClebschCoef * Clebsch.GetCoefficient(m4, m3, J);
		       this->InteractionFactorsUpUpUpDown[i][Index1] += this->PseudoPotentials[3][J >> 1] * (TmpCoefficient-TmpCoefficient2);
		     }
		 }
	       this->InteractionFactorsUpUpUpDown[i][Index1] *= Factors[3];
	       ++TotalNbrInteractionFactors;              
               ++Index1;
	     }
          
	   // UpDownUpDown sector
	   for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
	     {
	       int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
	       int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
	       Clebsch.InitializeCoefficientIterator(m1, m2);
	       this->InteractionFactorsUpDownUpDown[i][Index2] = 0.0;
	       while (Clebsch.Iterate(J, ClebschCoef))
		 {
                   TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
                   TmpCoefficient2 = ClebschCoef * Clebsch.GetCoefficient(m4, m3, J);
                   this->InteractionFactorsUpDownUpDown[i][Index2] += (this->PseudoPotentials[4][J >> 1] * TmpCoefficient - this->PseudoPotentials[5][J >> 1] * TmpCoefficient2);
		 }
	       this->InteractionFactorsUpDownUpDown[i][Index2] *= Factors[4];
	       ++TotalNbrInteractionFactors;
	       ++Index2;
	    }	

	  // DownDownUpDown sector
	  for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
	     {
	       int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
	       int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
	       Clebsch.InitializeCoefficientIterator(m1, m2);
	       this->InteractionFactorsDownDownUpDown[i][Index3] = 0.0;
	       while (Clebsch.Iterate(J, ClebschCoef))
		 {
		   if (((J >> 1) & 1) == Sign)
		     {
		       TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		       TmpCoefficient2 = ClebschCoef * Clebsch.GetCoefficient(m4, m3, J);
		       this->InteractionFactorsDownDownUpDown[i][Index3] += this->PseudoPotentials[6][J >> 1] * (TmpCoefficient-TmpCoefficient2);
		     }
		 }
	       this->InteractionFactorsDownDownUpDown[i][Index3] *= Factors[5];
	       ++TotalNbrInteractionFactors;
	       ++Index3;
	     }	
	 }

      Index1=0;Index2=0;Index3=0;
      // DownDown input particles
      for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
	{
	  int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	  int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;

	   // UpUpDownDown sector
	   for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
	     {
	       int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
	       int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
	       Clebsch.InitializeCoefficientIterator(m1, m2);
	       this->InteractionFactorsUpUpDownDown[i][Index1] = 0.0;
	       while (Clebsch.Iterate(J, ClebschCoef))
		 {
		   if (((J >> 1) & 1) == Sign)
		     {
		       TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		       this->InteractionFactorsUpUpDownDown[i][Index1] += this->PseudoPotentials[7][J >> 1] * TmpCoefficient;
		     }
		 }
	       this->InteractionFactorsUpUpDownDown[i][Index1] *= Factors[6];
	       ++TotalNbrInteractionFactors;
	       ++Index1;
	     }

	   // UpDownDownDown sector
	   for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
	     {
	       int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
	       int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
	       Clebsch.InitializeCoefficientIterator(m1, m2);
	       this->InteractionFactorsUpDownDownDown[i][Index2] = 0.0;
	       while (Clebsch.Iterate(J, ClebschCoef))
		 {
		   if (((J >> 1) & 1) == Sign)
		     {
		       TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		       TmpCoefficient2 = ClebschCoef * Clebsch.GetCoefficient(m4, m3, J);
		       this->InteractionFactorsUpDownDownDown[i][Index2] += this->PseudoPotentials[8][J >> 1] * (TmpCoefficient-TmpCoefficient2);
		     }
		 }
	       this->InteractionFactorsUpDownDownDown[i][Index2] *= Factors[7];
	       ++TotalNbrInteractionFactors;
	       ++Index2;
	    }	

	   // DownDownDownDown sector
	   for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
	     {
	       int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
	       int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
	       Clebsch.InitializeCoefficientIterator(m1, m2);
	       this->InteractionFactorsDownDownDownDown[i][Index3] = 0.0;
	       while (Clebsch.Iterate(J, ClebschCoef))
		 {
		   if (((J >> 1) & 1) == Sign)
		     {
		       TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		       this->InteractionFactorsDownDownDownDown[i][Index3] += this->PseudoPotentials[9][J >> 1] * TmpCoefficient;
		     }
		 }
	       this->InteractionFactorsDownDownDownDown[i][Index3] *= Factors[8];
	       ++TotalNbrInteractionFactors;
	       ++Index3;
	    }	
	 } 
      }
  }

  // Not yet implemented for bosons
  else 
    {
      cout << "Error: ParticleOnSphereWithSpinFullHamiltonian cannot handle bosons yet." << endl;
      exit(1);
    }
 cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
 cout << "====================================" << endl;
}
