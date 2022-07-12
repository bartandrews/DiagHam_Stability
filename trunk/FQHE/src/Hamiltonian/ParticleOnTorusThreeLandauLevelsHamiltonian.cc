////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 29/06/2010                      //
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


#include "Hamiltonian/ParticleOnTorusThreeLandauLevelsHamiltonian.h"
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

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define OneOverSqrt2 0.707106781186547524400844362104849039284


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusThreeLandauLevelsHamiltonian::ParticleOnTorusThreeLandauLevelsHamiltonian(ParticleOnSphereWithSU3Spin* particles, int nbrParticles, int maxMomentum, double cyclotronenergy,
												   double ratio, bool haveDelta, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  cout<<"Torus with spin constructor"<<endl;
  
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Lx = sqrt(2.0 * M_PI * (double)this->NbrLzValue * this->Ratio);
  this->Ly = sqrt(2.0 * M_PI * (double)this->NbrLzValue * this->InvRatio);
  this->Gx = 2.0 * M_PI / this->Lx;
  this->Gy = 2.0 * M_PI / this->Ly;
  cout << "-------------------------------------------------"<<endl;
  cout << "-     Geometry: Lx = " << this->Lx << " , Ly = " << this->Ly<< "  ;" << endl;
  cout << "-------------------------------------------------"<<endl;

  this->CyclotronEnergy=cyclotronenergy;
  
  this->HaveDelta = haveDelta;
  
  if (this->HaveDelta)
    this->FiniteQZeroComponent = 1.0;
  else
    this->FiniteQZeroComponent = 0.0;
  
  cout << "Cyclotron energy: "<<this->CyclotronEnergy << endl;
  
  this->OneBodyInteractionFactors11 = 0;
  this->OneBodyInteractionFactors22 = 0;
  this->OneBodyInteractionFactors33 = 0;
  if (this->CyclotronEnergy != 0)
    {
      this->OneBodyInteractionFactors11 = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->OneBodyInteractionFactors11[i] = 0.0;
      
      this->OneBodyInteractionFactors22 = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->OneBodyInteractionFactors22[i] = this->CyclotronEnergy;
      
      this->OneBodyInteractionFactors33 = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->OneBodyInteractionFactors33[i] = 2.0*this->CyclotronEnergy;
    }
  
  this->OneBodyInteractionFactors12 = 0;
  this->OneBodyInteractionFactors23 = 0;
  this->OneBodyInteractionFactors13 = 0;
  
  this->HermitianSymmetryFlag = false;
  
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  
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

ParticleOnTorusThreeLandauLevelsHamiltonian::~ParticleOnTorusThreeLandauLevelsHamiltonian() 
{
  if (this->OneBodyInteractionFactors11 != 0)
    {
      delete[] this->OneBodyInteractionFactors11;
    }
  if (this->OneBodyInteractionFactors22 != 0)
    {
      delete[] this->OneBodyInteractionFactors22;
    }
    if (this->OneBodyInteractionFactors33 != 0)
    {
      delete[] this->OneBodyInteractionFactors33;
    }
}


// evaluate all interaction factors
//   

void ParticleOnTorusThreeLandauLevelsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  
  this->NbrInterSectorSums = this->NbrLzValue;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrInterSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];      
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
      this->NbrInterSectorIndicesPerSum[i] = 0;
    }
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      {
	int TmpIndex = (m1 + m2) % this->NbrLzValue;
	this->InterSectorIndicesPerSum[TmpIndex][this->NbrInterSectorIndicesPerSum[TmpIndex] << 1] = m1;
	this->InterSectorIndicesPerSum[TmpIndex][1 + (this->NbrInterSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	++this->NbrInterSectorIndicesPerSum[TmpIndex];
      }
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      cout<<"Currently only implemented for bosons"<<endl;
      exit(1);
      /*
	this->NbrIntraSectorSums = this->NbrLzValue;
	this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
	for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpIndex = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpIndex][this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpIndex][1 + (this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpIndex];
	  }

      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupupupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
      */
    }
  else //Bosons
    {
      this->NbrIntraSectorSums = this->NbrLzValue;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpIndex = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpIndex][this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpIndex][1 + (this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpIndex];
	  }
      
//*************************************************************************************************
//*************************************************************************************************
//~~~~~~~~~~~~~~~~~~~~~~START FILLING IN THE INTERACTION TERMS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*************************************************************************************************
//*************************************************************************************************

//------------------------------------- INTRA TERMS --------------------------------------
    
    //****************** 1 1 1 1 *************************
	
      this->InteractionFactors1111 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1111[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors1111[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 0));
			}
		      else
			{
			  this->InteractionFactors1111[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 0, 0));
		          
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors1111[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 0));
			}
		      else
			{
			  this->InteractionFactors1111[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      
      //****************** 2 2 2 2 *************************
      this->InteractionFactors2222 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors2222[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors2222[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 1, 1));
			}
		      else
			{
			  this->InteractionFactors2222[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 1));
                          
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors2222[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 1, 1));
			}
		      else
			{
			  this->InteractionFactors2222[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 3 3 3 3 *************************
      this->InteractionFactors3333 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors3333[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors3333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 2, 2));
			}
		      else
			{
			  this->InteractionFactors3333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 2, 2));
                          
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors3333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 2, 2));
			}
		      else
			{
			  this->InteractionFactors3333[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 2, 2);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


      //****************** 1 1 2 2 *************************
      this->InteractionFactors1122 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1122[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors1122[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 1, 1));
			}
		      else
			{
			  this->InteractionFactors1122[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 1));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors1122[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 1, 1));
			}
		      else
			{
			  this->InteractionFactors1122[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** 2 2 1 1 *************************
      this->InteractionFactors2211 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors2211[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors2211[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 0));
			}
		      else
			{
			  this->InteractionFactors2211[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 0, 0));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors2211[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 0));
			}
		      else
			{
			  this->InteractionFactors2211[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 1 3 3 *************************
      this->InteractionFactors1133 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1133[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors1133[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 2, 2));
			}
		      else
			{
			  this->InteractionFactors1133[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 2, 2));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors1133[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 2, 2));
			}
		      else
			{
			  this->InteractionFactors1133[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 2, 2);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** 3 3 1 1 *************************
      this->InteractionFactors3311 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors3311[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors3311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 0, 0));
			}
		      else
			{
			  this->InteractionFactors3311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 0, 0));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors3311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 0, 0));
			}
		      else
			{
			  this->InteractionFactors3311[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 0);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 2 3 3 *************************
      this->InteractionFactors2233 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors2233[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors2233[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 2, 2));
			}
		      else
			{
			  this->InteractionFactors2233[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 2, 2));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors2233[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 2, 2));
			}
		      else
			{
			  this->InteractionFactors2233[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 2, 2);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** 3 3 2 2 *************************
      this->InteractionFactors3322 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors3322[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors3322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 1, 1));
			}
		      else
			{
			  this->InteractionFactors3322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 1, 1));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors3322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 1, 1));
			}
		      else
			{
			  this->InteractionFactors3322[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 1, 1);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


//-------------------------------------------INTER--------------------------------------------------------

      //****************** 1 2 1 2 *************************
      this->InteractionFactors1212 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1212[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors1212[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 1)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 2 1 3 *************************
      this->InteractionFactors1213 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1213[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors1213[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 2)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 2, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 2, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 2)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 2 2 3 *************************
      this->InteractionFactors1223 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1223[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors1223[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 1, 2)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 2, 1)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 2, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 1, 2)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 3 1 2 *************************
      this->InteractionFactors1312 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1312[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors1312[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 0, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 0, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 2, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 0, 1)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 3 1 3 *************************
      this->InteractionFactors1313 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1313[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors1313[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 0, 2)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 0, 2, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 2, 2, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 0, 2)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 3 2 3 *************************
      this->InteractionFactors1323 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1323[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors1323[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 1, 2)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 0, 2, 1)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 2, 2, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 1, 2)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 3 1 2 *************************
      this->InteractionFactors2312 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors2312[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors2312[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 0, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 1, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 2, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 0, 1)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 3 1 3 *************************
      this->InteractionFactors2313 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors2313[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors2313[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 0, 2)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 1, 2, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 2, 2, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 0, 2)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 3 2 3 *************************
      this->InteractionFactors2323 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors2323[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactors2323[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 1, 2)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 1, 2, 1)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 2, 2, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 1, 2)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


//-------------------------------------------MIXED--------------------------------------------------------

      //****************** 1 1 1 2 *************************
      this->InteractionFactors1112 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1112[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors1112[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 1, 0));
		    }
		  else
		   {
		      this->InteractionFactors1112[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 1 1 3 *************************
      this->InteractionFactors1113 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1113[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors1113[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 2, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 2, 0));
		    }
		  else
		   {
		      this->InteractionFactors1113[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 2, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 1 2 3 *************************
      this->InteractionFactors1123 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1123[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors1123[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 2, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 1, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 2, 1));
		    }
		  else
		   {
		      this->InteractionFactors1123[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 2, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
    //****************** 2 2 1 2*************************
      this->InteractionFactors2212 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors2212[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors2212[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 1, 0));
		    }
		  else
		   {
		      this->InteractionFactors2212[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 2 1 3*************************
      this->InteractionFactors2213 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors2213[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors2213[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 2, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 2, 0));
		    }
		  else
		   {
		      this->InteractionFactors2213[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 2, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 2 2 3*************************
      this->InteractionFactors2223 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors2223[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors2223[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 2, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 1, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 2, 1));
		    }
		  else
		   {
		      this->InteractionFactors2223[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 2, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

//****************** 3 3 1 2*************************
      this->InteractionFactors3312 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors3312[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors3312[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 1, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 0, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 2, 1, 0));
		    }
		  else  
		   {
		      this->InteractionFactors3312[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 1, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 3 3 1 3*************************
      this->InteractionFactors3313 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors3313[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors3313[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 2, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 0, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 2, 2, 0));
		    }
		  else
		   {
		      this->InteractionFactors3313[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 0, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 2, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 3 3 2 3*************************
      this->InteractionFactors3323 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors3323[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m3 != m4)
		    {
		      this->InteractionFactors3323[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 1, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 2, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 2, 1, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 2, 2, 1));
		    }
		  else
		   {
		      this->InteractionFactors3323[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 2, 2, 1, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 2, 2, 2, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	

      //****************** 1 2 1 1 *************************
      this->InteractionFactors1211 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1211[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors1211[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 0, 0));
		    }
		  else
		   {
		      this->InteractionFactors1211[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 1 3 1 1 *************************
      this->InteractionFactors1311 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1311[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors1311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 0, 0, 0));
		    }
		  else
		   {
		      this->InteractionFactors1311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 0, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	//****************** 2 3 1 1 *************************
      this->InteractionFactors2311 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors2311[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors2311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 1, 0, 0));
		    }
		  else
		   {
		      this->InteractionFactors2311[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 0, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** 1 2 2 2 *************************
      this->InteractionFactors1222 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1222[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors1222[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 1, 1));
		    }
		  else
		   {
		      this->InteractionFactors1222[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 1, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

//****************** 1 3 2 2 *************************
      this->InteractionFactors1322 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1322[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors1322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 0, 1, 1));
		    }
		  else
		   {
		      this->InteractionFactors1322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 1, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
      
      
//****************** 2 3 2 2 *************************
      this->InteractionFactors2322 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors2322[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors2322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 1, 1, 1));
		    }
		  else
		   {
		      this->InteractionFactors2322[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 1, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
	
	//****************** 1 2 3 3 *************************
      this->InteractionFactors1233 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1233[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors1233[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 2, 2));
		    }
		  else
		   {
		      this->InteractionFactors1233[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 2, 2));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

//****************** 1 3 3 3 *************************
      this->InteractionFactors1333 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1333[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors1333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 0, 2, 2));
		    }
		  else
		   {
		      this->InteractionFactors1333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 0, 2, 2));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
      
      
//****************** 2 3 3 3 *************************
      this->InteractionFactors2333 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors2333[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          if (m1 != m2)
		    {
		      this->InteractionFactors2333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 2, 2)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 2, 1, 2, 2));
		    }
		  else
		   {
		      this->InteractionFactors2333[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 2, 2, 2)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 2, 1, 2, 2));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
	
//*************************************************************************************************
//*************************************************************************************************
//~~~~~~~~~~~~~~~~~~~~~~END FILLING IN THE INTERACTION TERMS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*************************************************************************************************
//*************************************************************************************************

    } //Bosons

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;

}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// l1,...,l4 = LL index
// return value = numerical coefficient

Complex ParticleOnTorusThreeLandauLevelsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int ll1, int ll2, int ll3, int ll4)
{
  Complex Sum(0.0, 0.0);
  double N1;
  double N2 = (double)(m1 - m4);
  double Q2, Qx, Qy;
  double Xj13 = Gy * (double)(m1 - m3);
  Complex Coefficient(1.0,0.0);
  Complex PrecisionPos, PrecisionNeg;
  Complex Phase;

  while ((Norm(Sum) + Norm(Coefficient)) != Norm(Sum))
    {
      Qx = 0.0;
      Qy = this->Gy * N2;
      Q2 = Qx * Qx + Qy * Qy;

      if (Q2 != 0.0)
	{
	  Coefficient = exp(- 0.5 * Q2) * this->Potential(Q2, Qx, Qy, ll1, ll2, ll3, ll4);
	  PrecisionPos = Coefficient;
          PrecisionNeg = Coefficient;
	}
      else
	{
	  Coefficient = this->FiniteQZeroComponent * this->FormFactor(0, 0, 0, ll1, ll4) * this->FormFactor(0, 0, 0, ll2, ll3);
	  //Coefficient.Im = 0.0;
	  PrecisionPos.Re = 1.0;
	  PrecisionPos.Im = 0.0;
	  PrecisionNeg.Re = 1.0;
	  PrecisionNeg.Im = 0.0;
	}

      N1 = 1.0;
      while ((Norm(Coefficient) + Norm(PrecisionPos) + Norm(PrecisionNeg)) != Norm(Coefficient))
	{
      //Sum over positive N1
       Qx = this->Gx * N1;
       Qy = this->Gy * N2;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionPos = exp(-0.5 * Q2)*  this->Potential(Q2, Qx, Qy, ll1, ll2, ll3, ll4);          	
       Phase.Re = cos(Qx * Xj13);
       Phase.Im =  sin(Qx * Xj13);
       Coefficient += (PrecisionPos * Phase);

       //Sum over negative N1
       Qx = -this->Gx * N1;
       Qy = this->Gy * N2;
       Q2 = Qx * Qx + Qy * Qy;
	
       PrecisionNeg = exp(-0.5 * Q2)*  this->Potential(Q2, Qx, Qy,  ll1, ll2, ll3, ll4);	  
       Phase.Re = cos(Qx * Xj13);
       Phase.Im = sin(Qx * Xj13);
       Coefficient += (PrecisionNeg * Phase);

      //Increment N1
       N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += (double)this->NbrLzValue;
    }

  //cout<<"Tmp Sum "<<Sum<<endl;
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = Sum;	    
  while ((Norm(Sum) + Norm(Coefficient)) != Norm(Sum))
    {
      Qx = 0.0;
      Qy = this->Gy * N2;
      Q2 = Qx * Qx + Qy * Qy;

      if (Q2 != 0.0)
	{
	  Coefficient = exp(-0.5 * Q2)*  this->Potential(Q2, Qx, Qy,  ll1, ll2, ll3, ll4);
	  PrecisionPos = Coefficient;
          PrecisionNeg = Coefficient;
	}
      else
	{
	  Coefficient = this->FiniteQZeroComponent * this->FormFactor(0, 0, 0, ll1, ll4) * this->FormFactor(0, 0, 0, ll2, ll3);
          //Coefficient.Im = 0.0;
	  PrecisionPos.Re = 1.0;
	  PrecisionPos.Im = 0.0;
	  PrecisionNeg.Re = 1.0;
	  PrecisionNeg.Im = 0.0;
	}

      N1 = 1.0;
      while ((Norm(Coefficient) + Norm(PrecisionPos) + Norm(PrecisionNeg)) != Norm(Coefficient))
	{
       //Sum over positive N1
       Qx = this->Gx * N1;
       Qy = this->Gy * N2;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionPos = exp(-0.5 * Q2) *  this->Potential(Q2, Qx, Qy, ll1, ll2, ll3, ll4);	  
       Phase.Re = cos(Qx * Xj13);
       Phase.Im = sin(Qx * Xj13);
       Coefficient += (PrecisionPos * Phase);

       //Sum over negative N1
       Qx = -this->Gx * N1;
       Qy = this->Gy * N2;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionNeg = exp(-0.5 * Q2) *  this->Potential(Q2, Qx, Qy,  ll1, ll2, ll3, ll4);	  
       Phase.Re = cos(Qx * Xj13);
       Phase.Im = sin(Qx * Xj13);
       Coefficient += (PrecisionNeg * Phase);
       //Increment N1
       N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= (double)this->NbrLzValue;
    }
 
  if (this->HaveDelta)
    Sum /= ((double)this->NbrLzValue);
  else
    Sum /= (4.0 * M_PI * (double)this->NbrLzValue);

  return Sum;
}

// Returns the Potential times the form factor
// Q2 = momentum squared (Qx^2+Qy^2)
// Qx = qx component
// Qy = qy component 
// l1,...,l4 = LL index
Complex ParticleOnTorusThreeLandauLevelsHamiltonian::Potential(double Q2, double Qx, double Qy, int ll1, int ll2, int ll3, int ll4)
{
    Complex Potential;

    if (this->HaveDelta)
     {
       Potential.Re = 1.0;
       Potential.Im = 0.0;
     }
    else
     {
       Potential.Re = 2.0 * M_PI/sqrt(Q2);
       Potential.Im = 0.0;
     }

    Potential *= (FormFactor(Q2, Qx, Qy, ll1, ll4) * FormFactor(Q2, -Qx, -Qy, ll2, ll3)) ;
   
    return Potential;
}

// Returns the formfactor
// Input : momentum + LL index
Complex ParticleOnTorusThreeLandauLevelsHamiltonian::FormFactor(double Q2, double Qx, double Qy, int ll1, int ll2)
{
  Complex ff(1.0, 0.0);
  
  if (ll1 == 0)
    {
      if(ll2 == 1)
	{
	  ff.Re = Qy * OneOverSqrt2;
	  ff.Im = Qx * OneOverSqrt2;
	}
      else
	{
	  if(ll2 == 2)
	    {
	      ff.Re =  0.5* OneOverSqrt2 * ( Qy*Qy- Qx*Qx) ;
	      ff.Im = Qx*Qy*OneOverSqrt2;
	    }  
	}
    }
  else
    {
      if (ll1 == 1)
	{
	  if(ll2 == 0)
	    {
	      ff.Re = -Qy * OneOverSqrt2;
	      ff.Im = Qx * OneOverSqrt2; 
	    }
	  else
	    {
	      if(ll2 == 1)
		{
		  ff.Re = 1.0 - 0.5 * Q2;
		  ff.Im = 0.0;
		}
	      else
		{
		  ff.Re =  Qy*0.5*(2.0 - 0.5 * Q2) ;
		  ff.Im = Qx*0.5*(2.0 - 0.5 * Q2);
		}
	    }
	}
      else
	{
	  if (ll2 == 0)
	    {
	      ff.Re =  0.5* OneOverSqrt2 * ( Qy*Qy- Qx*Qx) ;
	      ff.Im = -Qx*Qy*OneOverSqrt2;
	    }
	  else
	    {
	      if(ll2 == 1)
		{
		  ff.Re =  -Qy*0.5*(2.0 - 0.5 * Q2) ;
		  ff.Im = Qx*0.5*(2.0 - 0.5 * Q2);
		}
	      else
		{
		  ff.Re = 1.0 -  Q2 + 0.125*Q2*Q2;
		  ff.Im = 0.0;
		}
	    }
	}
    }
    
  if (((ll1+ll2) & 0x1) != 0)
    ff *= -1.0;
  
    
  return ff;   
}
