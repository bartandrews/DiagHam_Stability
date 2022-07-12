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


#include "Hamiltonian/ParticleOnTorusTwoLandauLevelsHamiltonian.h"
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

ParticleOnTorusTwoLandauLevelsHamiltonian::ParticleOnTorusTwoLandauLevelsHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int maxMomentum, double cyclotronenergy,
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

  this->OneBodyInteractionFactorsupup = 0;
  if (this->CyclotronEnergy < 0)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->OneBodyInteractionFactorsupup[i] = this->CyclotronEnergy;
    }

  this->OneBodyInteractionFactorsdowndown = 0;
  if (this->CyclotronEnergy > 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->OneBodyInteractionFactorsdowndown[i] = this->CyclotronEnergy;
    }

  this->OneBodyInteractionFactorsupdown = 0;

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

ParticleOnTorusTwoLandauLevelsHamiltonian::~ParticleOnTorusTwoLandauLevelsHamiltonian() 
{
  cout<<"Destructor"<<endl;
  if (this->InteractionFactorsupupupup != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupupupup[i];
	  delete[] this->InteractionFactorsupupdowndown[i];
	  delete[] this->InteractionFactorsdowndownupup[i];
	  delete[] this->InteractionFactorsdowndowndowndown[i];
	}
      delete[] this->InteractionFactorsupupupup;
      delete[] this->InteractionFactorsupupdowndown;
      delete[] this->InteractionFactorsdowndownupup;
      delete[] this->InteractionFactorsdowndowndowndown;


      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupupupdown[i];
	  delete[] this->InteractionFactorsupdownupup[i];
	  delete[] this->InteractionFactorsupdowndowndown[i];
	  delete[] this->InteractionFactorsdowndownupdown[i];
	  delete[] this->InteractionFactorsupdownupdown[i];
	}
	  delete[] this->InteractionFactorsupupupdown;
	  delete[] this->InteractionFactorsupdownupup;
	  delete[] this->InteractionFactorsupdowndowndown;
	  delete[] this->InteractionFactorsdowndownupdown;
	  delete[] this->InteractionFactorsupdownupdown;
    }

  if (this->OneBodyInteractionFactorsupup != 0)
    {
      delete[] this->OneBodyInteractionFactorsupup;
    }
  if (this->OneBodyInteractionFactorsdowndown != 0)
    {
      delete[] this->OneBodyInteractionFactorsdowndown;
    }
  if (this->OneBodyInteractionFactorsupdown != 0)
    {
      delete[] this->OneBodyInteractionFactorsupdown;
    }

  cout<<"Exit destructor"<<endl;
}


// evaluate all interaction factors
//   

void ParticleOnTorusTwoLandauLevelsHamiltonian::EvaluateInteractionFactors()
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

      //****************** UP UP UP UP *************************
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
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
		  if (m1 != m2)
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupupupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 0));
			}
		      else
			{
			  this->InteractionFactorsupupupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 0, 0));
		          
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupupupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 0));
 		          
			}
		      else
			{
			  this->InteractionFactorsupupupup[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 0);
		          
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


      //****************** DOWN DOWN DOWN DOWN *************************
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactorsdowndowndowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 1, 1));
			}
		      else
			{
			  this->InteractionFactorsdowndowndowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 1));
                          
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsdowndowndowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 1, 1));
			}
		      else
			{
			  this->InteractionFactorsdowndowndowndown[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 1, 1);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


      //****************** UP UP DOWN DOWN *************************
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactorsupupdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 1, 1));
			}
		      else
			{
			  this->InteractionFactorsupupdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 1));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupupdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 1, 1));
			}
		      else
			{
			  this->InteractionFactorsupupdowndown[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 1, 1);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** DOWN DOWN UP UP *************************
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactorsdowndownupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 0));
			}
		      else
			{
			  this->InteractionFactorsdowndownupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 0, 0));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsdowndownupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 0));
			}
		      else
			{
			  this->InteractionFactorsdowndownupup[i][Index] = this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 0);
			}
		    }	  
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


//-------------------------------------------INTER--------------------------------------------------------

      //****************** UP DOWN UP DOWN *************************
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
 	          this->InteractionFactorsupdownupdown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 1)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 1, 0)
                                                                   + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 1)
                                                                   );
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


//-------------------------------------------MIXED--------------------------------------------------------

      //****************** UP UP UP DOWN *************************
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupupupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		      this->InteractionFactorsupupupdown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 0, 0, 0, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 0, 0, 1, 0));
		    }
		  else
		   {
		      this->InteractionFactorsupupupdown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 0, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 0, 1, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}


      //****************** UP DOWN UP UP *************************
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		      this->InteractionFactorsupdownupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 0)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 0, 0));
		    }
		  else
		   {
		      this->InteractionFactorsupdownupup[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 0, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 0, 0));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** UP DOWN DOWN DOWN *************************
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		      this->InteractionFactorsupdowndowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 1, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 0, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 0, 1, 1));
		    }
		  else
		   {
		      this->InteractionFactorsupdowndowndown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 0, 1, 1, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 0, 1, 1));
	  	   }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}

      //****************** DOWN DOWN UP DOWN*************************
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		      this->InteractionFactorsdowndownupdown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 0)
								    + this->EvaluateInteractionCoefficient(m4, m3, m1, m2, 1, 1, 0, 1)
								    + this->EvaluateInteractionCoefficient(m4, m3, m2, m1, 1, 1, 1, 0));
		    }
		  else
		   {
		      this->InteractionFactorsdowndownupdown[i][Index] = (this->EvaluateInteractionCoefficient(m3, m4, m1, m2, 1, 1, 0, 1)
								    + this->EvaluateInteractionCoefficient(m3, m4, m2, m1, 1, 1, 1, 0));
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

Complex ParticleOnTorusTwoLandauLevelsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int ll1, int ll2, int ll3, int ll4)
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
       Phase.Im =  -sin(Qx * Xj13);
       Coefficient += (PrecisionPos * Phase);

       //Sum over negative N1
       Qx = -this->Gx * N1;
       Qy = this->Gy * N2;
       Q2 = Qx * Qx + Qy * Qy;
	
       PrecisionNeg = exp(-0.5 * Q2)*  this->Potential(Q2, Qx, Qy,  ll1, ll2, ll3, ll4);	  
       Phase.Re = cos(Qx * Xj13);
       Phase.Im = -sin(Qx * Xj13);
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
       Phase.Im = -sin(Qx * Xj13);
       Coefficient += (PrecisionPos * Phase);

       //Sum over negative N1
       Qx = -this->Gx * N1;
       Qy = this->Gy * N2;
       Q2 = Qx * Qx + Qy * Qy;

       PrecisionNeg = exp(-0.5 * Q2) *  this->Potential(Q2, Qx, Qy,  ll1, ll2, ll3, ll4);	  
       Phase.Re = cos(Qx * Xj13);
       Phase.Im = -sin(Qx * Xj13);
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
Complex ParticleOnTorusTwoLandauLevelsHamiltonian::Potential(double Q2, double Qx, double Qy, int ll1, int ll2, int ll3, int ll4)
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
Complex ParticleOnTorusTwoLandauLevelsHamiltonian::FormFactor(double Q2, double Qx, double Qy, int ll1, int ll2)
{

  Complex ff(1.0, 0.0);

  if ((ll1 == 1) && (ll2 == 1))
    {
      ff.Re = 1.0 - 0.5 * Q2;
      ff.Im = 0.0;
    } 
  else
   {
     if ((ll1 == 1) && (ll2 == 0))
      {
       ff.Re = Qy * OneOverSqrt2;
       ff.Im = Qx * OneOverSqrt2; 
      } 
     else 
      if ((ll1 == 0) && (ll2 == 1)) 
       {
        ff.Re = -Qy * OneOverSqrt2;
        ff.Im = Qx * OneOverSqrt2; 
       }
   }

 return ff;   
}
